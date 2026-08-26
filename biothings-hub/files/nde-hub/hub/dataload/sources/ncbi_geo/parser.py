import datetime
import json
import os
import sqlite3
import time
import traceback
from collections import Counter
from pathlib import Path

import dateutil.parser

from .sample_charateristics import insert_value, parse_sample_characteristics, parse_series_sample_characteristics

try:
    from config import logger
except ImportError:
    import logging

    logger = logging.getLogger(__name__)


GSM_SUMMARY_CACHE_VERSION = 1
GSM_SUMMARY_CACHE_FILENAME = ".gsm_summary_cache.sqlite3"
GSM_SUMMARY_CACHE_COMMIT_INTERVAL = 500
GSE_PARSE_LOG_INTERVAL = 1000

_GSM_SUMMARY_FIELDS = frozenset(
    {
        "!Sample_geo_accession",
        "!Sample_type",
        "!Sample_library_source",
    }
)
_GSM_CHARACTERISTICS_PREFIX = "!Sample_characteristics"


def _store_soft_value(result, key, value):
    """Store a SOFT field while preserving the parser's scalar/list behavior."""
    if key in result:
        if isinstance(result[key], list):
            result[key].append(value)
        else:
            result[key] = [result[key], value]
    else:
        result[key] = value


def parse_gsm_summary(filepath):
    """Read only the GSM fields that contribute to a parent GSE document."""
    result = {}
    with open(filepath, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("!Sample_") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            if key not in _GSM_SUMMARY_FIELDS and not key.startswith(_GSM_CHARACTERISTICS_PREFIX):
                continue
            _store_soft_value(result, key, value.strip())
    return result


def _freeze_value(value):
    """Return a hashable equality key for a parsed aggregate value."""
    if isinstance(value, dict):
        return "dict", tuple(sorted((key, _freeze_value(item)) for key, item in value.items()))
    if isinstance(value, list):
        return "list", tuple(_freeze_value(item) for item in value)
    if isinstance(value, tuple):
        return "tuple", tuple(_freeze_value(item) for item in value)
    try:
        hash(value)
    except TypeError:
        return type(value).__name__, repr(value)
    return type(value).__name__, value


class UniqueValueAccumulator:
    """Preserve first-seen order and scalar/list output with O(1) deduplication."""

    def __init__(self, output):
        self.output = output
        self._seen = {}

    def add(self, output, key, value):
        if output is not self.output:
            raise ValueError("Accumulator used with a different output dictionary")

        marker = _freeze_value(value)
        seen = self._seen.setdefault(key, set())
        if marker in seen:
            return
        seen.add(marker)

        if key not in output:
            output[key] = value
        elif isinstance(output[key], list):
            output[key].append(value)
        else:
            output[key] = [output[key], value]


class GSMSummaryCache:
    """Persistent summaries of the small GSM subset needed by the GSE parser."""

    _DDL = """
        CREATE TABLE IF NOT EXISTS gsm_summary (
            gsm_id TEXT PRIMARY KEY,
            source_size INTEGER NOT NULL,
            source_mtime_ns INTEGER NOT NULL,
            cache_version INTEGER NOT NULL,
            summary TEXT NOT NULL
        )
    """

    def __init__(self, data_folder):
        self.path = os.path.join(data_folder, GSM_SUMMARY_CACHE_FILENAME)
        self.conn = None
        self.hits = 0
        self.misses = 0
        self.stale = 0
        self.writes = 0
        self.parse_seconds = 0.0
        self._pending_writes = 0
        try:
            self.conn = sqlite3.connect(self.path, timeout=60)
            self.conn.execute("PRAGMA busy_timeout = 60000")
            self.conn.execute(self._DDL)
            self.conn.commit()
            logger.info(
                "GSM summary cache: path=%s version=%s",
                self.path,
                GSM_SUMMARY_CACHE_VERSION,
            )
        except (OSError, sqlite3.Error) as e:
            logger.warning("GSM summary cache disabled at %s: %s", self.path, e)
            if self.conn is not None:
                self.conn.close()
                self.conn = None

    def get(self, gsm_id, filepath):
        stat = os.stat(filepath)
        if self.conn is not None:
            row = self.conn.execute(
                "SELECT source_size, source_mtime_ns, cache_version, summary "
                "FROM gsm_summary WHERE gsm_id = ?",
                (gsm_id,),
            ).fetchone()
            if row is not None:
                source_size, source_mtime_ns, cache_version, summary = row
                if (
                    source_size == stat.st_size
                    and source_mtime_ns == stat.st_mtime_ns
                    and cache_version == GSM_SUMMARY_CACHE_VERSION
                ):
                    try:
                        value = json.loads(summary)
                    except json.JSONDecodeError:
                        self.stale += 1
                    else:
                        self.hits += 1
                        return value
                else:
                    self.stale += 1

        self.misses += 1
        started = time.monotonic()
        value = parse_gsm_summary(filepath)
        self.parse_seconds += time.monotonic() - started
        self._put(gsm_id, stat, value)
        return value

    def _put(self, gsm_id, stat, summary):
        if self.conn is None:
            return
        try:
            self.conn.execute(
                "INSERT OR REPLACE INTO gsm_summary "
                "(gsm_id, source_size, source_mtime_ns, cache_version, summary) VALUES (?, ?, ?, ?, ?)",
                (
                    gsm_id,
                    stat.st_size,
                    stat.st_mtime_ns,
                    GSM_SUMMARY_CACHE_VERSION,
                    json.dumps(summary, ensure_ascii=False, separators=(",", ":")),
                ),
            )
            self.writes += 1
            self._pending_writes += 1
            if self._pending_writes >= GSM_SUMMARY_CACHE_COMMIT_INTERVAL:
                self.conn.commit()
                self._pending_writes = 0
        except sqlite3.Error as e:
            logger.warning("Could not update GSM summary cache for %s: %s", gsm_id, e)
            self.conn.rollback()
            self._pending_writes = 0

    def close(self):
        if self.conn is None:
            return
        try:
            self.conn.commit()
        finally:
            self.conn.close()
            self.conn = None


def get_full_name(name):
    parts = name.split(",")
    while len(parts) < 3:
        parts.append("")
    if len(parts) > 3:
        parts = [parts[0], "".join(parts[1:-1]), parts[-1]]
    first, middle, last = parts
    full_name = " ".join([part.strip() for part in [first, middle, last] if part.strip()])
    return full_name


def build_species(names, taxids):
    """DefinedTerms for GEO organism names, carrying the NCBI taxid listed alongside them."""
    if not names:
        return []
    names = names if isinstance(names, list) else [names]
    taxids = taxids if isinstance(taxids, list) else ([] if taxids is None else [taxids])

    species = []
    for i, name in enumerate(names):
        term = {"@type": "DefinedTerm", "name": name}
        if i < len(taxids) and str(taxids[i]).strip().isdigit():
            term["identifier"] = str(taxids[i]).strip()
        if term not in species:
            species.append(term)
    return species


def parse_soft_series(filepath):
    """
    Parse a GEO SOFT series file into a dictionary.
    Each key is the SOFT field (e.g., '!Series_title'), value is a list if repeated, or a string.
    """
    result = {}
    with open(filepath, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("^"):
                continue
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()
                _store_soft_value(result, key, value)
    return result


def parse_gsm(data_folder):
    """
    Parse a GEO SOFT platform file into a dictionary.
    Each key is the SOFT field (e.g., '!Platform_title'), value is a list if repeated, or a string.
    """

    gse_dir = os.path.join(data_folder, "gsm")
    records = get_records(gse_dir)
    with open(Path(__file__).resolve().parent / "sex_mappings.json", "r") as f:
        sex_mapping = json.load(f)
    with open(Path(__file__).resolve().parent / "gsm_nde_mapping.json", "r") as f:
        nde_mapping = json.load(f)
    with open(Path(__file__).resolve().parent / "mapping_dict.json", "r") as f:
        sample_mapping = json.load(f)

    for item in records:
        if not item.get("!Sample_geo_accession"):
            continue
        _id = item.get("!Sample_geo_accession")
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + _id
        output = {
            "@context": "http://schema.org/",
            "@type": "Sample",
            "_id": _id.casefold(),
            "identifier": _id,
            "url": url,
            "distribution": [{"@type": "DataDownload", "contentUrl": url}],
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "NCBI GEO",
                "url": "https://www.ncbi.nlm.nih.gov/geo/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
            "additionalType": "ExperimentalRunSample",
            "conditionsOfAccess": "Open",
        }

        if name := item.get("!Sample_title"):
            output["name"] = name

        if date_published := item.get("!Sample_status"):
            date_str = date_published.replace("Public on ", "")
            try:
                dt = dateutil.parser.parse(date_str, ignoretz=True).date().isoformat()
                output["datePublished"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_str}': {e}")

        if date_created := item.get("!Sample_submission_date"):
            try:
                dt = dateutil.parser.parse(date_created, ignoretz=True).date().isoformat()
                output["dateCreated"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_created}': {e}")

        if date_modified := item.get("!Sample_last_update_date"):
            try:
                dt = dateutil.parser.parse(date_modified, ignoretz=True).date().isoformat()
                output["dateModified"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_modified}': {e}")

        if sample_type := item.get("!Sample_type"):
            if isinstance(sample_type, list):
                output["sampleType"] = [{"name": s, "@type": "DefinedTerm"} for s in sample_type]
            else:
                output["sampleType"] = [{"name": sample_type, "@type": "DefinedTerm"}]

        if description := item.get("!Sample_description"):
            if isinstance(description, list):
                output["description"] = " ".join(description)
            else:
                output["description"] = description

        if sample_process := item.get("!Sample_data_processing"):
            if isinstance(sample_process, list):
                output["sampleProcess"] = " ".join(sample_process)
            else:
                output["sampleProcess"] = sample_process

        if author := item.get("!Sample_contact_name"):
            output["author"] = {"@type": "Person", "name": get_full_name(author)}
            if affiliation := item.get("!Sample_contact_institute"):
                output["author"]["affiliation"] = {"@type": "Organization", "name": affiliation}

        species = []
        for key in sorted(k for k in item if k.startswith("!Sample_organism_ch")):
            channel = key[len("!Sample_organism_ch") :]
            for term in build_species(item[key], item.get(f"!Sample_taxid_ch{channel}")):
                if term not in species:
                    species.append(term)
        if species:
            output["species"] = species

        if instrument := item.get("!Sample_instrument_model"):
            if isinstance(instrument, list):
                output["instrument"] = [{"@type": "DefinedTerm", "name": i} for i in instrument]
            else:
                output["instrument"] = {"@type": "DefinedTerm", "name": instrument}

        mts = []

        if mt := item.get("!Sample_library_selection"):
            if isinstance(mt, list):
                for m in mt:
                    mts.append({"@type": "DefinedTerm", "name": m})
            else:
                mts.append({"@type": "DefinedTerm", "name": mt})

        if mt := item.get("!Sample_library_strategy"):
            if isinstance(mt, list):
                for m in mt:
                    mts.append({"@type": "DefinedTerm", "name": m})
            else:
                mts.append({"@type": "DefinedTerm", "name": mt})

        if mts:
            output["measurementTechnique"] = mts

        if sample_type := item.get("!Sample_library_source"):
            if "sampleType" not in output:
                output["sampleType"] = []
            if isinstance(sample_type, list):
                output["sampleType"].extend([{"name": s, "@type": "DefinedTerm"} for s in sample_type])
            else:
                output["sampleType"].append({"name": sample_type, "@type": "DefinedTerm"})

        if same_as := item.get("!Sample_relation"):
            output["sameAs"] = same_as

        if is_basis_for := item.get("!Sample_series_id"):
            if isinstance(is_basis_for, list):
                output["isBasisFor"] = [
                    {
                        "@type": "Dataset",
                        "identifier": sid,
                        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + sid,
                    }
                    for sid in is_basis_for
                ]
            else:
                output["isBasisFor"] = {
                    "@type": "Dataset",
                    "identifier": is_basis_for,
                    "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + is_basis_for,
                }

        for key, value in item.items():
            if key.startswith("!Sample_characteristics"):
                parse_sample_characteristics(output, value, sample_mapping, nde_mapping, sex_mapping)

            if "protocol" in key.lower() or key.startswith("!Sample_molecule"):
                if isinstance(value, list):
                    sample_process = " ".join(value)
                else:
                    sample_process = value
                if "sampleProcess" in output and output["sampleProcess"]:
                    if isinstance(output["sampleProcess"], list):
                        output["sampleProcess"] = " ".join(output["sampleProcess"])
                    output["sampleProcess"] = output["sampleProcess"] + " " + sample_process
                else:
                    output["sampleProcess"] = sample_process

            if key.startswith("!Sample_supplementary_file"):
                if isinstance(value, list):
                    for v in value:
                        insert_value(output, "distribution", {"@type": "DataDownload", "contentUrl": v})
                else:
                    insert_value(output, "distribution", {"@type": "DataDownload", "contentUrl": value})

        yield output


def get_records(data_folder):
    for root, dirs, _ in os.walk(data_folder):
        for dir in dirs:
            logger.info(dir)
            dirpath = os.path.join(root, dir)
            for file in os.listdir(dirpath):
                # logger.info(file)
                if file.endswith(".txt"):
                    fpath = os.path.join(dirpath, file)
                    try:
                        item = parse_soft_series(fpath)
                        yield item
                    except Exception as e:
                        logger.error(f"Error parsing {fpath}: {e}")


def find_gsm_file(data_folder, acc):
    # Extract prefix (GSE/GSM) and numeric part
    prefix = acc[:3]
    num = acc[3:]
    # Pad numeric part to at least 3 digits for nnn, or 6 for full
    padded = num.zfill(6)
    subdir = prefix + padded[:3] + "nnn"
    subdir_path = os.path.join(data_folder, subdir)
    file_path = os.path.join(subdir_path, f"{acc}.txt")
    if os.path.exists(file_path):
        return file_path
    else:
        logger.warning(f"GSM file not found: {file_path}")
        return None


def parse_series_sample(
    item,
    sample,
    sample_mapping,
    nde_mapping,
    sex_mapping,
    accumulator,
    unmapped_subproperties,
):
    """
    Given a GEO Series item and an output dictionary, parse for the sample field in the dataset
    """
    if not item.get("!Sample_geo_accession"):
        return
    sample_elements = sample["aggregateElement"]
    accumulator.add(
        sample_elements,
        "url",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + item.get("!Sample_geo_accession"),
    )
    if sample_type := item.get("!Sample_type"):
        if isinstance(sample_type, list):
            for s in sample_type:
                accumulator.add(sample_elements, "sampleType", {"name": s, "@type": "DefinedTerm"})
        else:
            accumulator.add(sample_elements, "sampleType", {"name": sample_type, "@type": "DefinedTerm"})

    if sample_type := item.get("!Sample_library_source"):
        if isinstance(sample_type, list):
            for s in sample_type:
                accumulator.add(sample_elements, "sampleType", {"name": s, "@type": "DefinedTerm"})
        else:
            accumulator.add(sample_elements, "sampleType", {"name": sample_type, "@type": "DefinedTerm"})

    for key, value in item.items():
        if key.startswith("!Sample_characteristics"):
            parse_series_sample_characteristics(
                sample_elements,
                value,
                sample_mapping,
                nde_mapping,
                sex_mapping,
                add_value=accumulator.add,
                unmapped_subproperties=unmapped_subproperties,
            )


def _log_gse_parser_stats(metrics, gsm_cache, *, finished=False):
    cache_requests = gsm_cache.hits + gsm_cache.misses
    hit_rate = 100.0 * gsm_cache.hits / cache_requests if cache_requests else 0.0
    top_unmapped = ", ".join(
        f"{name}={count}" for name, count in metrics["unmapped_subproperties"].most_common(10)
    )
    logger.info(
        "GSE parser%s: gse=%s gsm_references=%s gsm_cache_hits=%s gsm_cache_misses=%s "
        "gsm_cache_stale=%s gsm_cache_hit_rate=%.1f%% gsm_cache_writes=%s "
        "gsm_uncached_parse_seconds=%.1fs gse_transform_seconds=%.1fs "
        "unmapped_characteristics=%s top_unmapped=[%s]",
        " finished" if finished else "",
        metrics["gse"],
        metrics["gsm_references"],
        gsm_cache.hits,
        gsm_cache.misses,
        gsm_cache.stale,
        hit_rate,
        gsm_cache.writes,
        gsm_cache.parse_seconds,
        metrics["transform_seconds"],
        sum(metrics["unmapped_subproperties"].values()),
        top_unmapped,
    )


def parse_gse(data_folder):
    """Yield GSE documents while reusing persistent, source-validated GSM summaries."""
    gsm_cache = GSMSummaryCache(data_folder)
    metrics = {
        "gse": 0,
        "gsm_references": 0,
        "transform_seconds": 0.0,
        "unmapped_subproperties": Counter(),
    }
    try:
        yield from _parse_gse(data_folder, gsm_cache, metrics)
    finally:
        gsm_cache.close()
        _log_gse_parser_stats(metrics, gsm_cache, finished=True)


def _parse_gse(data_folder, gsm_cache, metrics):
    """
    Parse a GEO SOFT platform file into a dictionary.
    Each key is the SOFT field (e.g., '!Platform_title'), value is a list if repeated, or a string.
    """
    with open(Path(__file__).resolve().parent / "sex_mappings.json", "r") as f:
        sex_mapping = json.load(f)
    with open(Path(__file__).resolve().parent / "gse_nde_mapping.json", "r") as f:
        nde_mapping = json.load(f)
    with open(Path(__file__).resolve().parent / "mapping_dict.json", "r") as f:
        sample_mapping = json.load(f)

    gse_dir = os.path.join(data_folder, "gse")
    records = get_records(gse_dir)
    for item in records:
        record_started = time.monotonic()
        if not item.get("!Series_geo_accession"):
            continue
        _id = item.get("!Series_geo_accession")
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + _id
        output = {
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": _id.casefold(),
            "identifier": _id,
            "url": url,
            "distribution": {"@type": "DataDownload", "contentUrl": url},
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "NCBI GEO",
                "url": "https://www.ncbi.nlm.nih.gov/geo/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
            "conditionsOfAccess": "Open",
        }

        if gsm_ids := item.get("!Series_sample_id"):
            sample = {
                "@type": "SampleCollection",
                "itemListElement": [],
                "aggregateElement": {},
                "numberOfItems": {"@type": "QuantitativeValue", "value": 0, "unitText": "sample"},
            }
            gsm_dir = os.path.join(data_folder, "gsm")
            if not isinstance(gsm_ids, list):
                gsm_ids = [gsm_ids]
            accumulator = UniqueValueAccumulator(sample["aggregateElement"])
            metrics["gsm_references"] += len(gsm_ids)
            for gsm_id in gsm_ids:
                sample["itemListElement"].append(
                    {
                        "@type": "Sample",
                        "identifier": gsm_id,
                        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + gsm_id,
                        "_id": gsm_id.casefold(),
                    }
                )
                sample["numberOfItems"]["value"] += 1
                gsm_file = find_gsm_file(gsm_dir, gsm_id)
                try:
                    sample_item = gsm_cache.get(gsm_id, gsm_file)
                    parse_series_sample(
                        sample_item,
                        sample,
                        sample_mapping,
                        nde_mapping,
                        sex_mapping,
                        accumulator,
                        metrics["unmapped_subproperties"],
                    )
                except Exception as e:
                    logger.error(f"Error parsing GSM file {gsm_file}: {e}")
                    logger.error(traceback.format_exc())
                    continue

            if sample.get("itemListElement"):
                # logger.info(f"gse_id: {output['_id']}")
                # logger.info(f"sample: {sample}")
                output["sample"] = sample

        if name := item.get("!Series_title"):
            output["name"] = name

        if date_published := item.get("!Series_status"):
            date_str = date_published.replace("Public on ", "")
            try:
                dt = dateutil.parser.parse(date_str, ignoretz=True).date().isoformat()
                output["datePublished"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_str}': {e}")

        if date_created := item.get("!Series_submission_date"):
            try:
                dt = dateutil.parser.parse(date_created, ignoretz=True).date().isoformat()
                output["dateCreated"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_created}': {e}")

        if date_modified := item.get("!Series_last_update_date"):
            try:
                dt = dateutil.parser.parse(date_modified, ignoretz=True).date().isoformat()
                output["dateModified"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_modified}': {e}")

        if pmids := item.get("!Series_pubmed_id"):
            if isinstance(pmids, list):
                output["pmids"] = ",".join(pmids)
            else:
                output["pmids"] = pmids

        if descriptions := item.get("!Series_summary"):
            if isinstance(descriptions, list):
                output["description"] = " ".join(descriptions)
            else:
                output["description"] = descriptions

        if description := item.get("!Series_overall_design"):
            if isinstance(description, list):
                output["description"] = (output.get("description", "") + " " + " ".join(description)).strip()
            else:
                output["description"] = (output.get("description", "") + " " + description).strip()

        if mt := item.get("!Series_type"):
            if isinstance(mt, list):
                output["measurementTechnique"] = [{"@type": "DefinedTerm", "name": m} for m in mt]
            else:
                output["measurementTechnique"] = {"@type": "DefinedTerm", "name": mt}

        if authors := item.get("!Series_contributor"):
            if isinstance(authors, list):
                output["author"] = [{"@type": "Person", "name": get_full_name(a)} for a in authors]
            else:
                output["author"] = [{"@type": "Person", "name": get_full_name(authors)}]

        if institutes := item.get("!Series_contact_institute"):
            if isinstance(institutes, list):
                output["sourceOrganization"] = [{"@type": "Organization", "name": i} for i in institutes]
            else:
                output["sourceOrganization"] = [{"@type": "Organization", "name": institutes}]

        if species := build_species(item.get("!Series_platform_organism"), item.get("!Series_platform_taxid")):
            output["species"] = species

        metrics["gse"] += 1
        metrics["transform_seconds"] += time.monotonic() - record_started
        if metrics["gse"] % GSE_PARSE_LOG_INTERVAL == 0:
            _log_gse_parser_stats(metrics, gsm_cache)
        yield output
