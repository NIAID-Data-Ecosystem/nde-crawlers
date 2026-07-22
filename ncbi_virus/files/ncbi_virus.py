#!/usr/bin/env python3
"""
NCBI Virus DataCollection crawler for the NIAID Data Ecosystem.

The crawler reads NCBI Datasets Virus data reports and aggregates nucleotide
sequence records into one DataCollection per exact viral NCBI Taxonomy ID.

Data source
-----------
Virus reports are retrieved from:
    https://api.ncbi.nlm.nih.gov/datasets/v2/virus/taxon/{taxon}/dataset_report

The all-virus root taxon can be too large for a single report query, so the
default production path discovers the immediate child taxa under Viruses
(taxon 10239) and streams each child query separately.
"""

from __future__ import annotations

import copy
import datetime as dt
import hashlib
import json
import logging
import os
import re
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from typing import Any, Iterable, Iterator, Optional

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")

# ---------------------------------------------------------------------------
# API configuration
# ---------------------------------------------------------------------------

API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ROOT_VIRUS_TAXON = "10239"
USER_AGENT = "nde-ncbi-virus-crawler/0.1"
MAX_RETRIES = 4
RETRY_BACKOFF = 10

# Current top-level children under Viruses, used only if taxonomy discovery
# fails. These values are intentionally broad seeds; exact output records are
# still grouped by each report's virus.tax_id.
FALLBACK_SEED_TAXA = [
    "2840022",  # Adnaviria
    "2731341",  # Duplodnaviria
    "3700708",  # Efunaviria
    "2731342",  # Floreoviria
    "3700715",  # Pleomoviria
    "2559587",  # Riboviria
    "2842242",  # Ribozyviria
    "3412574",  # Singelaviria
    "2732004",  # Varidnaviria
    "3700717",  # Volvereviria
    "3700676",  # Basaltiviridae
    "3700688",  # Lomiviridae
    "3700703",  # Xigoviridae
    "3700704",  # Yamazakiviridae
    "3403061",  # Viruses incertae sedis
    "186616",  # environmental samples
    "12333",  # unclassified bacterial viruses
    "12429",  # unclassified viruses
]

PAGE_SIZE = int(os.environ.get("NCBI_VIRUS_PAGE_SIZE", "1000"))
MAX_PAGES_PER_SEED = int(os.environ.get("NCBI_VIRUS_MAX_PAGES_PER_SEED", "0"))
MAX_RECORDS_PER_SEED = int(os.environ.get("NCBI_VIRUS_MAX_RECORDS_PER_SEED", "0"))
MAX_SEED_TAXA = int(os.environ.get("NCBI_VIRUS_MAX_SEED_TAXA", "0"))
MAX_EXAMPLES = int(os.environ.get("NCBI_VIRUS_MAX_EXAMPLES", "5"))
MAX_AGGREGATE_VALUES = int(os.environ.get("NCBI_VIRUS_MAX_AGGREGATE_VALUES", "50"))
REQUEST_DELAY = float(os.environ.get("NCBI_VIRUS_REQUEST_DELAY", "0.34"))
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
DEDUP_ACCESSIONS = os.environ.get("NCBI_VIRUS_DEDUP_ACCESSIONS", "").lower() in {"1", "true", "yes"}
USE_SQL_CACHE = os.environ.get("NCBI_VIRUS_USE_SQL_CACHE", "true").lower() in {"1", "true", "yes"}
RESUME_INCOMPLETE_CACHE = os.environ.get("NCBI_VIRUS_RESUME_INCOMPLETE_CACHE", "true").lower() in {
    "1",
    "true",
    "yes",
}
CACHE_EXPIRE_DAYS = int(os.environ.get("NCBI_VIRUS_CACHE_EXPIRE_DAYS", "30"))
CACHE_DIR = os.environ.get("NCBI_VIRUS_CACHE_DIR", "/cache/ncbi_virus")
CACHE_DB = os.environ.get("NCBI_VIRUS_CACHE_DB", "ncbi_virus.db")
CACHE_SCHEMA_VERSION = "1"

# ---------------------------------------------------------------------------
# Static metadata from the approved mapping sheet
# ---------------------------------------------------------------------------

SOURCE_NAME = "NCBI Virus"
SOURCE_URL = "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/"
SOURCE_DESCRIPTION = (
    "Integrative resource for retrieval, display, and analysis of curated "
    "viral sequences and large sequence datasets."
)
SOURCE_CATALOG_RECORD = (
    "https://data.niaid.nih.gov/resources?id=dde_b3bb013bc893ca80"
)
LICENSE = "http://opendefinition.org/licenses/odc-odbl/"
USAGE_INFO = {"@type": "CreativeWork", "url": "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/"}

ABOUT = [
    {"@type": "DefinedTerm", "name": "Nucleotide sequence"},
    {"@type": "DefinedTerm", "name": "Genome sequence"},
    {"@type": "DefinedTerm", "name": "Protein sequence metadata"},
]

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "name": "Virology",
        "identifier": "topic_0781",
        "url": "http://edamontology.org/topic_0781",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Genomics",
        "identifier": "topic_0622",
        "url": "http://edamontology.org/topic_0622",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Infectious disease",
        "identifier": "topic_0804",
        "url": "http://edamontology.org/topic_0804",
        "inDefinedTermSet": "EDAM",
    },
]

MEASUREMENT_TECHNIQUE = [
    {
        "@type": "DefinedTerm",
        "name": "Genome sequencing",
    },
    {
        "@type": "DefinedTerm",
        "name": "Nucleotide sequencing",
    },
]

VARIABLE_MEASURED = [
    {"@type": "DefinedTerm", "name": "Viral genome sequence"},
    {"@type": "DefinedTerm", "name": "Viral nucleotide sequence"},
    {"@type": "DefinedTerm", "name": "Protein annotation"},
]

ENCODING_FORMAT = [
    {"@type": "DefinedTerm", "name": "FASTA"},
    {"@type": "DefinedTerm", "name": "GenBank"},
    {"@type": "DefinedTerm", "name": "JSON Lines"},
    {"@type": "DefinedTerm", "name": "TSV"},
]

HOWTO_STEPS = [
    (
        "Step 1. Discover broad viral seed taxa from the NCBI Datasets "
        "taxonomy filtered_subtree endpoint under Viruses (taxon 10239)."
    ),
    (
        "Step 2. Use the NCBI Datasets Virus API to query records by each "
        "seed taxon, paginating dataset_report responses."
    ),
    (
        "Step 3. Group each source record by its exact virus.tax_id and "
        "aggregate counts, release/update dates, hosts, locations, samples, "
        "submitters, linked BioProjects/SRA accessions, and representative "
        "sequence records."
    ),
    (
        "Step 4. Fill in static source-level fields from the NCBI Virus "
        "resource catalog mapping, including access conditions, license, "
        "topicCategory, measurementTechnique, variableMeasured, usageInfo, "
        "and crawler provenance."
    ),
]

INVALID_VALUES = {
    "",
    "-",
    "?",
    "na",
    "n/a",
    "nan",
    "none",
    "not applicable",
    "not collected",
    "not determined",
    "not provided",
    "not recorded",
    "null",
    "unknown",
    "unspecified",
}

DATE_RE = re.compile(r"(?P<year>\d{4})(?:-(?P<month>\d{1,2})(?:-(?P<day>\d{1,2}))?)?")


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def _clean_string(value: Any) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, bool):
        return str(value).lower()
    text = str(value).strip()
    if not text or text.lower() in INVALID_VALUES:
        return None
    return re.sub(r"\s+", " ", text)


def _clean_enum(value: Any, prefix: str = "") -> Optional[str]:
    text = _clean_string(value)
    if not text:
        return None
    if prefix and text.startswith(prefix):
        text = text[len(prefix) :]
    return text.replace("_", " ").strip().title()


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _add_counter(counter: Counter[str], value: Any, limit: int = MAX_AGGREGATE_VALUES) -> None:
    text = _clean_string(value)
    if not text:
        return
    counter[text] += 1
    if len(counter) > limit * 2:
        for key, _count in counter.most_common()[limit:]:
            counter.pop(key, None)


def _counter_values(counter: Counter[str], limit: int = MAX_AGGREGATE_VALUES) -> list[str]:
    return [value for value, _count in counter.most_common(limit)]


def _add_limited_dict(
    target: dict[str, dict[str, Any]],
    key: Any,
    value: Optional[dict[str, Any]],
    limit: int = MAX_AGGREGATE_VALUES,
) -> None:
    text_key = _clean_string(key)
    if not text_key or not value or text_key in target or len(target) >= limit:
        return
    target[text_key] = value


def _dedupe(items: Iterable[Any]) -> list[Any]:
    seen: set[str] = set()
    result = []
    for item in items:
        key = json.dumps(item, sort_keys=True) if isinstance(item, (dict, list)) else str(item)
        if key in seen:
            continue
        seen.add(key)
        result.append(item)
    return result


def _make_property(name: str, value: Any, property_id: Optional[str] = None) -> Optional[dict[str, Any]]:
    if isinstance(value, list):
        cleaned = [_clean_string(v) for v in value]
        value = [v for v in cleaned if v]
        if not value:
            return None
    elif isinstance(value, bool):
        value = value
    elif isinstance(value, (int, float)):
        value = value
    else:
        value = _clean_string(value)
        if not value:
            return None

    prop = {"@type": "PropertyValue", "name": name, "value": value}
    if property_id:
        prop["propertyID"] = property_id
    return prop


def _taxonomy_identifier(tax_id: Any) -> Optional[str]:
    text = _clean_string(tax_id)
    return f"taxonomy:{text}" if text else None


def _taxonomy_url(tax_id: Any) -> Optional[str]:
    text = _clean_string(tax_id)
    return f"https://www.ncbi.nlm.nih.gov/taxonomy/{text}" if text else None


def _taxonomy_term(tax_id: Any, name: Any, alternate_names: Optional[Iterable[Any]] = None) -> Optional[dict[str, Any]]:
    clean_name = _clean_string(name)
    clean_id = _clean_string(tax_id)
    if not clean_name and not clean_id:
        return None

    term: dict[str, Any] = {"@type": "DefinedTerm"}
    if clean_name:
        term["name"] = clean_name
    if clean_id:
        term["identifier"] = f"taxonomy:{clean_id}"
        term["url"] = f"https://www.ncbi.nlm.nih.gov/taxonomy/{clean_id}"
        term["inDefinedTermSet"] = "NCBI Taxonomy"

    alts = [_clean_string(v) for v in alternate_names or []]
    alts = [v for v in alts if v and v != clean_name]
    if alts:
        term["alternateName"] = _dedupe(alts)
    return term


def _place(name: Any, alternate_name: Any = None) -> Optional[dict[str, Any]]:
    clean_name = _clean_string(name)
    alt = _clean_string(alternate_name)
    if not clean_name and not alt:
        return None
    place: dict[str, Any] = {"@type": "Place"}
    if clean_name:
        place["name"] = clean_name
    if alt:
        place["alternateName"] = alt
    return place


def _ncbi_url(kind: str, identifier: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/{kind}/{urllib.parse.quote(identifier)}"


def _report_api_url(taxon: str, page_size: int = 1) -> str:
    params = urllib.parse.urlencode({"page_size": page_size})
    return f"{API_BASE}/virus/taxon/{urllib.parse.quote(str(taxon))}/dataset_report?{params}"


def _virus_search_url(taxon: str, name: str) -> str:
    query_value = f"{name}, taxid:{taxon}"
    return (
        "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?"
        f"VirusLineage_ss={urllib.parse.quote(query_value)}"
    )


def _parse_iso_date(value: Any) -> Optional[str]:
    text = _clean_string(value)
    if not text:
        return None
    try:
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        return dt.datetime.fromisoformat(text).date().isoformat()
    except Exception:
        match = DATE_RE.search(text)
        if match:
            return _date_from_match(match)
    return None


def _date_from_match(match: re.Match[str]) -> str:
    year = int(match.group("year"))
    month = match.group("month")
    day = match.group("day")
    if month and day:
        return f"{year:04d}-{int(month):02d}-{int(day):02d}"
    if month:
        return f"{year:04d}-{int(month):02d}"
    return f"{year:04d}"


def _date_key(value: str, is_end: bool = False) -> tuple[int, int, int]:
    parts = value.split("-")
    year = int(parts[0])
    month = int(parts[1]) if len(parts) > 1 else (12 if is_end else 1)
    day = int(parts[2]) if len(parts) > 2 else (31 if is_end else 1)
    return year, month, day


def _collection_date_bounds(value: Any) -> tuple[Optional[str], Optional[str]]:
    text = _clean_string(value)
    if not text:
        return None, None

    pieces = [piece.strip() for piece in re.split(r"\s*/\s*|\s+to\s+", text, flags=re.I) if piece.strip()]
    if not pieces:
        pieces = [text]

    starts = []
    ends = []
    for piece in pieces:
        match = DATE_RE.search(piece)
        if not match:
            continue
        parsed = _date_from_match(match)
        starts.append(parsed)
        ends.append(parsed)

    if not starts:
        return None, None
    return starts[0], ends[-1]


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------


def _fetch_json(path: str, params: Optional[dict[str, Any]] = None) -> Any:
    params = dict(params or {})
    if NCBI_API_KEY:
        params.setdefault("api_key", NCBI_API_KEY)

    if path.startswith("http"):
        url = path
    else:
        url = f"{API_BASE}{path}"

    if params:
        separator = "&" if "?" in url else "?"
        url = f"{url}{separator}{urllib.parse.urlencode(params, doseq=True)}"

    last_error: Optional[Exception] = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": USER_AGENT,
                },
            )
            with urllib.request.urlopen(req, timeout=180) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", "replace")
            last_error = RuntimeError(f"HTTP {exc.code} for {url}: {body[:500]}")
            if exc.code not in {429, 500, 502, 503, 504} or attempt == MAX_RETRIES:
                raise last_error
        except Exception as exc:
            last_error = exc
            if attempt == MAX_RETRIES:
                raise

        sleep_for = RETRY_BACKOFF * attempt
        logger.warning("Request failed on attempt %s/%s; sleeping %ss", attempt, MAX_RETRIES, sleep_for)
        time.sleep(sleep_for)

    raise RuntimeError(f"Failed to fetch {url}") from last_error


def _sleep_between_requests() -> None:
    if REQUEST_DELAY > 0:
        time.sleep(REQUEST_DELAY)


def _configured_taxa() -> list[str]:
    taxa: list[str] = []

    taxa_file = os.environ.get("NCBI_VIRUS_TAXA_FILE")
    if taxa_file:
        try:
            with open(taxa_file, encoding="utf-8") as fh:
                for line in fh:
                    taxa.extend(part.strip() for part in line.split(","))
        except FileNotFoundError:
            logger.warning("Configured NCBI_VIRUS_TAXA_FILE does not exist: %s", taxa_file)

    taxa_env = os.environ.get("NCBI_VIRUS_TAXA")
    if taxa_env:
        taxa.extend(part.strip() for part in taxa_env.split(","))

    cleaned = []
    for taxon in taxa:
        clean = _clean_string(taxon)
        if clean:
            cleaned.append(clean)
    return _dedupe(cleaned)


def _discover_seed_taxa() -> list[str]:
    configured = _configured_taxa()
    if configured:
        logger.info("Using %d configured NCBI Virus seed taxa", len(configured))
        return configured[:MAX_SEED_TAXA] if MAX_SEED_TAXA else configured

    root_taxon = os.environ.get("NCBI_VIRUS_ROOT_TAXON", ROOT_VIRUS_TAXON)
    try:
        data = _fetch_json(
            f"/taxonomy/taxon/{urllib.parse.quote(str(root_taxon))}/filtered_subtree",
            {"levels": 2},
        )
        edges = data.get("edges") or {}
        root = edges.get(str(root_taxon)) or {}
        taxa = [str(tax_id) for tax_id in root.get("visible_children") or []]
        taxa = [taxon for taxon in taxa if taxon != str(root_taxon)]
        if taxa:
            logger.info("Discovered %d NCBI Virus seed taxa under %s", len(taxa), root_taxon)
            return taxa[:MAX_SEED_TAXA] if MAX_SEED_TAXA else taxa
    except Exception as exc:
        logger.warning("Failed to discover NCBI Virus seed taxa: %s", exc)

    logger.warning("Using fallback NCBI Virus seed taxa")
    return FALLBACK_SEED_TAXA[:MAX_SEED_TAXA] if MAX_SEED_TAXA else FALLBACK_SEED_TAXA


def _fetch_report_page(taxon: str, page_token: Optional[str] = None) -> dict[str, Any]:
    params: dict[str, Any] = {"page_size": min(PAGE_SIZE, 1000)}
    if page_token:
        params["page_token"] = page_token
    return _fetch_json(f"/virus/taxon/{urllib.parse.quote(str(taxon))}/dataset_report", params)


def _fetch_reports_for_seed(taxon: str) -> Iterator[dict[str, Any]]:
    page_token: Optional[str] = None
    page_count = 0
    record_count = 0

    while True:
        logger.info("Fetching NCBI Virus reports for taxon %s page %s", taxon, page_count + 1)
        data = _fetch_report_page(taxon, page_token)
        reports = data.get("reports") or []
        total = data.get("total_count")
        if page_count == 0:
            logger.info("Taxon %s report total_count: %s", taxon, total)

        if not reports:
            break

        for report in reports:
            yield report
            record_count += 1
            if MAX_RECORDS_PER_SEED and record_count >= MAX_RECORDS_PER_SEED:
                logger.info("Reached NCBI_VIRUS_MAX_RECORDS_PER_SEED=%s for %s", MAX_RECORDS_PER_SEED, taxon)
                return

        page_count += 1
        if MAX_PAGES_PER_SEED and page_count >= MAX_PAGES_PER_SEED:
            logger.info("Reached NCBI_VIRUS_MAX_PAGES_PER_SEED=%s for %s", MAX_PAGES_PER_SEED, taxon)
            return

        page_token = data.get("next_page_token") or data.get("nextPageToken")
        if not page_token:
            break
        _sleep_between_requests()


# ---------------------------------------------------------------------------
# Report-to-schema helpers
# ---------------------------------------------------------------------------


def _biosample_accession(report: dict[str, Any]) -> Optional[str]:
    value = report.get("biosample") or report.get("bio_sample")
    if isinstance(value, dict):
        for key in ("accession", "accession_version", "identifier"):
            text = _clean_string(value.get(key))
            if text:
                return text
        return None
    return _clean_string(value)


def _source_database(report: dict[str, Any]) -> Optional[str]:
    return _clean_string(report.get("source_database"))


def _report_accession(report: dict[str, Any]) -> Optional[str]:
    return _clean_string(report.get("accession") or (report.get("nucleotide") or {}).get("accession_version"))


def _lineage_text(lineage: Any) -> Optional[str]:
    pieces = []
    for node in _as_list(lineage):
        if not isinstance(node, dict):
            continue
        name = _clean_string(node.get("name"))
        tax_id = _clean_string(node.get("tax_id"))
        if name and tax_id:
            pieces.append(f"{name} (taxonomy:{tax_id})")
        elif name:
            pieces.append(name)
    return " > ".join(pieces) if pieces else None


def _sample_location(report: dict[str, Any]) -> list[dict[str, Any]]:
    location = report.get("location") or {}
    places = []
    geo = _place(location.get("geographic_location"))
    region = _place(location.get("geographic_region"))
    state = _place(location.get("geographic_location"), location.get("usa_state"))
    for place in (state if state and location.get("usa_state") else geo, region):
        if place:
            places.append(place)
    return _dedupe(places)


def _sample_associated_genotype(report: dict[str, Any]) -> list[str]:
    values = []
    host = report.get("host") or {}
    names = host.get("infraspecific_names") or {}
    for label, key in (
        ("Breed", "breed"),
        ("Cultivar", "cultivar"),
        ("Ecotype", "ecotype"),
        ("Strain", "strain"),
    ):
        text = _clean_string(names.get(key))
        if text:
            values.append(f"{label}: {text}")
    return values


def _sample_from_report(report: dict[str, Any]) -> Optional[dict[str, Any]]:
    sample: dict[str, Any] = {"@type": "Sample"}

    isolate = report.get("isolate") or {}
    isolate_name = _clean_string(isolate.get("name"))
    if isolate_name:
        sample["name"] = isolate_name

    biosample = _biosample_accession(report)
    if biosample:
        sample["identifier"] = biosample
        sample["url"] = _ncbi_url("biosample", biosample)

    sample_types = []
    isolate_source = _clean_string(isolate.get("source"))
    mol_type = _clean_string(report.get("mol_type"))
    if isolate_source:
        sample_types.append({"@type": "DefinedTerm", "name": isolate_source})
    if mol_type:
        sample_types.append({"@type": "DefinedTerm", "name": mol_type})
    if sample_types:
        sample["sampleType"] = _dedupe(sample_types)

    collection_date = _clean_string(isolate.get("collection_date"))
    if collection_date:
        sample["dateCollected"] = collection_date
        start, end = _collection_date_bounds(collection_date)
        if start or end:
            sample["temporalCoverage"] = {
                "@type": "TemporalInterval",
                "temporalType": "sample collection",
            }
            if start:
                sample["temporalCoverage"]["startDate"] = start
            if end:
                sample["temporalCoverage"]["endDate"] = end

    places = _sample_location(report)
    if places:
        sample["locationOfOrigin"] = places

    host_names = (report.get("host") or {}).get("infraspecific_names") or {}
    sex = _clean_string(host_names.get("sex"))
    if sex:
        sample["sex"] = sex

    associated = _sample_associated_genotype(report)
    if associated:
        sample["associatedGenotype"] = associated

    host_isolate = _clean_string(host_names.get("isolate"))
    if host_isolate:
        sample["alternateName"] = host_isolate

    purpose = _clean_enum(report.get("purpose_of_sampling"), "PURPOSE_OF_SAMPLING_")
    if purpose and purpose != "Unknown":
        sample["experimentalPurpose"] = {"@type": "DefinedTerm", "name": purpose}

    lab_host = _clean_string(report.get("lab_host"))
    if lab_host:
        sample["cellType"] = {"@type": "DefinedTerm", "name": lab_host}
        sample["sampleProcess"] = f"Passaged in {lab_host}"

    return sample if len(sample) > 1 else None


def _example_work_from_report(report: dict[str, Any]) -> Optional[dict[str, Any]]:
    accession = _report_accession(report)
    if not accession:
        return None

    nucleotide = report.get("nucleotide") or {}
    isolate = report.get("isolate") or {}
    virus = report.get("virus") or {}
    host = report.get("host") or {}
    submitter = report.get("submitter") or {}

    identifiers = [
        accession,
        _clean_string(nucleotide.get("accession_version")),
        _clean_string(nucleotide.get("sequence_hash")),
    ]
    identifiers = _dedupe([value for value in identifiers if value])

    example: dict[str, Any] = {
        "@type": "CreativeWork",
        "identifier": identifiers[0] if len(identifiers) == 1 else identifiers,
        "name": _clean_string(isolate.get("name")) or _clean_string(nucleotide.get("title")) or accession,
        "url": _ncbi_url("nuccore", accession),
        "encodingFormat": copy.deepcopy(ENCODING_FORMAT),
    }

    source_db = _source_database(report)
    if source_db:
        example["isPartOf"] = {"@type": "DataCatalog", "name": source_db}

    additional = []
    for name, value in (
        ("Virus lineage", _lineage_text(virus.get("lineage"))),
        ("Host lineage", _lineage_text(host.get("lineage"))),
        ("Completeness", report.get("completeness")),
        ("Is annotated", report.get("is_annotated")),
        ("Sequence length", report.get("length")),
        ("Gene count", report.get("gene_count")),
        ("Protein count", report.get("protein_count")),
        ("Mature peptide count", report.get("mature_peptide_count")),
        ("Molecule type", report.get("mol_type")),
        ("Purpose of sampling", _clean_enum(report.get("purpose_of_sampling"), "PURPOSE_OF_SAMPLING_")),
        ("Lab host", report.get("lab_host")),
        ("Is lab host", report.get("is_lab_host")),
        ("Is vaccine strain", report.get("is_vaccine_strain")),
        ("Segment", report.get("segment")),
        ("Submitter country", submitter.get("country")),
        ("BioSample", _biosample_accession(report)),
    ):
        prop = _make_property(name, value)
        if prop:
            additional.append(prop)

    bioprojects = [_clean_string(v) for v in _as_list(report.get("bioprojects"))]
    bioprojects = [v for v in bioprojects if v]
    if bioprojects:
        prop = _make_property("BioProject accessions", bioprojects)
        if prop:
            additional.append(prop)

    sra_accessions = [_clean_string(v) for v in _as_list(report.get("sra_accessions"))]
    sra_accessions = [v for v in sra_accessions if v]
    if sra_accessions:
        prop = _make_property("SRA accessions", sra_accessions)
        if prop:
            additional.append(prop)

    if additional:
        example["additionalProperty"] = additional
    return example


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


class TaxonAccumulator:
    COUNTER_FIELDS = (
        "virus_common_names",
        "virus_strains",
        "pangolin",
        "keywords",
        "source_databases",
        "completeness",
        "sample_sources",
        "mol_types",
        "purposes",
        "lab_hosts",
        "sample_sexes",
        "sample_genotypes",
        "sample_alternate_names",
        "bioprojects",
        "biosamples",
        "sra_accessions",
        "accessions",
    )
    DICT_FIELDS = (
        "hosts",
        "spatial_coverage",
        "sample_locations",
        "authors",
    )
    LIST_FIELDS = (
        "example_works",
        "sample_examples",
    )
    SCALAR_FIELDS = (
        "tax_id",
        "record_count",
        "virus_name",
        "earliest_release",
        "latest_update",
        "earliest_collection",
        "latest_collection",
    )

    def __init__(self, tax_id: str) -> None:
        self.tax_id = tax_id
        self.record_count = 0
        self.virus_name: Optional[str] = None
        self.virus_common_names: Counter[str] = Counter()
        self.virus_strains: Counter[str] = Counter()
        self.pangolin: Counter[str] = Counter()
        self.keywords: Counter[str] = Counter()
        self.source_databases: Counter[str] = Counter()
        self.completeness: Counter[str] = Counter()
        self.hosts: dict[str, dict[str, Any]] = {}
        self.spatial_coverage: dict[str, dict[str, Any]] = {}
        self.sample_locations: dict[str, dict[str, Any]] = {}
        self.sample_sources: Counter[str] = Counter()
        self.mol_types: Counter[str] = Counter()
        self.purposes: Counter[str] = Counter()
        self.lab_hosts: Counter[str] = Counter()
        self.sample_sexes: Counter[str] = Counter()
        self.sample_genotypes: Counter[str] = Counter()
        self.sample_alternate_names: Counter[str] = Counter()
        self.bioprojects: Counter[str] = Counter()
        self.biosamples: Counter[str] = Counter()
        self.sra_accessions: Counter[str] = Counter()
        self.accessions: Counter[str] = Counter()
        self.authors: dict[str, dict[str, Any]] = {}
        self.example_works: list[dict[str, Any]] = []
        self.sample_examples: list[dict[str, Any]] = []
        self.earliest_release: Optional[str] = None
        self.latest_update: Optional[str] = None
        self.earliest_collection: Optional[str] = None
        self.latest_collection: Optional[str] = None

    def to_cache_dict(self) -> dict[str, Any]:
        data = {field: getattr(self, field) for field in self.SCALAR_FIELDS}
        for field in self.COUNTER_FIELDS:
            data[field] = dict(getattr(self, field))
        for field in self.DICT_FIELDS:
            data[field] = getattr(self, field)
        for field in self.LIST_FIELDS:
            data[field] = getattr(self, field)
        return data

    @classmethod
    def from_cache_dict(cls, data: dict[str, Any]) -> "TaxonAccumulator":
        accumulator = cls(str(data["tax_id"]))
        accumulator.record_count = int(data.get("record_count") or 0)
        accumulator.virus_name = data.get("virus_name")
        accumulator.earliest_release = data.get("earliest_release")
        accumulator.latest_update = data.get("latest_update")
        accumulator.earliest_collection = data.get("earliest_collection")
        accumulator.latest_collection = data.get("latest_collection")

        for field in cls.COUNTER_FIELDS:
            setattr(accumulator, field, Counter(data.get(field) or {}))
        for field in cls.DICT_FIELDS:
            setattr(accumulator, field, data.get(field) or {})
        for field in cls.LIST_FIELDS:
            setattr(accumulator, field, data.get(field) or [])
        return accumulator

    def add(self, report: dict[str, Any]) -> None:
        self.record_count += 1

        virus = report.get("virus") or {}
        if not self.virus_name:
            self.virus_name = _clean_string(virus.get("organism_name"))

        _add_counter(self.virus_common_names, virus.get("common_name"))
        _add_counter(self.keywords, virus.get("common_name"))
        _add_counter(self.pangolin, virus.get("pangolin_classification"))
        _add_counter(self.keywords, virus.get("pangolin_classification"))

        virus_names = virus.get("infraspecific_names") or {}
        _add_counter(self.virus_strains, virus_names.get("strain"))

        source_db = _source_database(report)
        _add_counter(self.source_databases, source_db)
        _add_counter(self.completeness, report.get("completeness"))

        accession = _report_accession(report)
        _add_counter(self.accessions, accession, limit=MAX_AGGREGATE_VALUES * 2)

        release = _parse_iso_date(report.get("release_date"))
        if release and (
            self.earliest_release is None
            or _date_key(release) < _date_key(self.earliest_release)
        ):
            self.earliest_release = release

        update = _parse_iso_date(report.get("update_date"))
        if update and (
            self.latest_update is None
            or _date_key(update, is_end=True) > _date_key(self.latest_update, is_end=True)
        ):
            self.latest_update = update

        isolate = report.get("isolate") or {}
        isolate_source = _clean_string(isolate.get("source"))
        _add_counter(self.sample_sources, isolate_source)
        collection_start, collection_end = _collection_date_bounds(isolate.get("collection_date"))
        if collection_start and (
            self.earliest_collection is None
            or _date_key(collection_start) < _date_key(self.earliest_collection)
        ):
            self.earliest_collection = collection_start
        if collection_end and (
            self.latest_collection is None
            or _date_key(collection_end, is_end=True) > _date_key(self.latest_collection, is_end=True)
        ):
            self.latest_collection = collection_end

        _add_counter(self.mol_types, report.get("mol_type"))
        purpose = _clean_enum(report.get("purpose_of_sampling"), "PURPOSE_OF_SAMPLING_")
        if purpose and purpose != "Unknown":
            _add_counter(self.purposes, purpose)
        _add_counter(self.lab_hosts, report.get("lab_host"))

        for place in _sample_location(report):
            key = "|".join(str(place.get(k, "")) for k in ("name", "alternateName"))
            _add_limited_dict(self.spatial_coverage, key, place)
            _add_limited_dict(self.sample_locations, key, place)

        host = report.get("host") or {}
        host_names = host.get("infraspecific_names") or {}
        host_term = _taxonomy_term(
            host.get("tax_id"),
            host.get("organism_name"),
            [host.get("common_name")],
        )
        if host_term:
            _add_limited_dict(self.hosts, host_term.get("identifier") or host_term.get("name"), host_term)
        _add_counter(self.keywords, host.get("common_name"))
        _add_counter(self.sample_sexes, host_names.get("sex"))
        _add_counter(self.sample_alternate_names, host_names.get("isolate"))
        for genotype in _sample_associated_genotype(report):
            _add_counter(self.sample_genotypes, genotype)

        for value in _as_list(report.get("bioprojects")):
            _add_counter(self.bioprojects, value)
        biosample = _biosample_accession(report)
        _add_counter(self.biosamples, biosample)
        for value in _as_list(report.get("sra_accessions")):
            _add_counter(self.sra_accessions, value)

        submitter = report.get("submitter") or {}
        affiliation = _clean_string(submitter.get("affiliation"))
        for name in _as_list(submitter.get("names")):
            clean_name = _clean_string(name)
            if not clean_name or clean_name in self.authors or len(self.authors) >= MAX_AGGREGATE_VALUES:
                continue
            person = {"@type": "Person", "name": clean_name}
            if affiliation:
                person["affiliation"] = affiliation
            self.authors[clean_name] = person

        if len(self.example_works) < MAX_EXAMPLES:
            example = _example_work_from_report(report)
            if example:
                self.example_works.append(example)

        if len(self.sample_examples) < MAX_EXAMPLES:
            sample = _sample_from_report(report)
            if sample:
                self.sample_examples.append(sample)

    def _infectious_agent(self) -> dict[str, Any]:
        alternate_names = []
        alternate_names.extend(_counter_values(self.virus_common_names, 10))
        alternate_names.extend(_counter_values(self.virus_strains, 10))
        return _taxonomy_term(self.tax_id, self.virus_name or f"NCBI Taxonomy {self.tax_id}", alternate_names) or {
            "@type": "DefinedTerm",
            "identifier": f"taxonomy:{self.tax_id}",
            "url": _taxonomy_url(self.tax_id),
        }

    def _sample_collection(self) -> Optional[dict[str, Any]]:
        sample: dict[str, Any] = {
            "@type": "SampleCollection",
            "numberOfItems": {
                "value": self.record_count,
                "unitText": "samples",
            },
        }

        if self.sample_examples:
            sample["itemListElement"] = _dedupe(self.sample_examples)

        aggregate: dict[str, Any] = {}
        sample_types = [
            {"@type": "DefinedTerm", "name": value}
            for value in _counter_values(self.sample_sources, 25)
        ]
        sample_types.extend(
            {"@type": "DefinedTerm", "name": value}
            for value in _counter_values(self.mol_types, 25)
        )
        if sample_types:
            aggregate["sampleType"] = _dedupe(sample_types)

        purposes = [
            {"@type": "DefinedTerm", "name": value}
            for value in _counter_values(self.purposes, 25)
        ]
        if purposes:
            aggregate["experimentalPurpose"] = purposes

        lab_hosts = _counter_values(self.lab_hosts, 25)
        if lab_hosts:
            aggregate["cellType"] = [{"@type": "DefinedTerm", "name": value} for value in lab_hosts]
            aggregate["sampleProcess"] = [f"Passaged in {value}" for value in lab_hosts]

        sexes = _counter_values(self.sample_sexes, 10)
        if sexes:
            aggregate["sex"] = sexes

        genotypes = _counter_values(self.sample_genotypes, 25)
        if genotypes:
            aggregate["associatedGenotype"] = genotypes

        alternate_names = _counter_values(self.sample_alternate_names, 25)
        if alternate_names:
            aggregate["alternateName"] = alternate_names

        if aggregate:
            sample["aggregateElement"] = aggregate

        if self.sample_locations:
            sample["locationOfOrigin"] = list(self.sample_locations.values())

        if self.earliest_collection or self.latest_collection:
            sample["temporalCoverage"] = {
                "@type": "TemporalInterval",
                "temporalType": "sample collection",
            }
            if self.earliest_collection:
                sample["temporalCoverage"]["startDate"] = self.earliest_collection
            if self.latest_collection:
                sample["temporalCoverage"]["endDate"] = self.latest_collection

        return sample if len(sample) > 2 or self.sample_examples else None

    def _is_part_of(self) -> list[dict[str, Any]]:
        items = []
        for accession in _counter_values(self.bioprojects, 25):
            items.append(
                {
                    "@type": "CreativeWork",
                    "identifier": accession,
                    "url": _ncbi_url("bioproject", accession),
                }
            )
        return items

    def _has_part(self) -> list[dict[str, Any]]:
        items = []
        for accession in _counter_values(self.accessions, 25):
            items.append(
                {
                    "@type": "CreativeWork",
                    "identifier": accession,
                    "url": _ncbi_url("nuccore", accession),
                }
            )
        for accession in _counter_values(self.sra_accessions, 25):
            items.append(
                {
                    "@type": "Dataset",
                    "identifier": accession,
                    "url": _ncbi_url("sra", accession),
                }
            )
        return _dedupe(items)

    def _data_collection_additional_properties(self) -> list[dict[str, Any]]:
        properties = []
        completeness_summary = [
            f"{name}: {count}" for name, count in self.completeness.most_common(10)
        ]
        if completeness_summary:
            prop = _make_property("Completeness distribution", completeness_summary)
            if prop:
                properties.append(prop)

        source_db_summary = [
            f"{name}: {count}" for name, count in self.source_databases.most_common(10)
        ]
        if source_db_summary:
            prop = _make_property("Source database distribution", source_db_summary)
            if prop:
                properties.append(prop)
        return properties

    def to_record(self) -> dict[str, Any]:
        name = self.virus_name or f"NCBI Taxonomy {self.tax_id}"
        url = _virus_search_url(self.tax_id, name)
        count_text = f"{self.record_count:,}"
        source_dbs = _counter_values(self.source_databases, 5)
        source_db_text = "/".join(source_dbs) if source_dbs else "NCBI"

        record: dict[str, Any] = {
            "_id": f"ncbi_virus_{self.tax_id}",
            "@type": "DataCollection",
            "identifier": f"taxonomy:{self.tax_id}",
            "name": f"{name} virus sequence records at NCBI Virus",
            "description": (
                f"Viral nucleotide sequence records for {name} available "
                f"through NCBI Virus. This DataCollection aggregates "
                f"{count_text} {source_db_text} sequence records for "
                f"NCBI Taxonomy ID {self.tax_id}. For more details, visit: {url}"
            ),
            "url": url,
            "about": copy.deepcopy(ABOUT),
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": SOURCE_NAME,
                "url": SOURCE_URL,
                "description": SOURCE_DESCRIPTION,
                "archivedAt": url,
                "versionDate": dt.datetime.now(dt.timezone.utc).date().isoformat(),
            },
            "conditionsOfAccess": "Open",
            "isAccessibleForFree": True,
            "license": LICENSE,
            "usageInfo": copy.deepcopy(USAGE_INFO),
            "topicCategory": copy.deepcopy(TOPIC_CATEGORY),
            "measurementTechnique": copy.deepcopy(MEASUREMENT_TECHNIQUE),
            "variableMeasured": copy.deepcopy(VARIABLE_MEASURED),
            "collectionSize": {
                "minValue": self.record_count,
                "unitText": "Viral nucleotide sequences",
            },
            "infectiousAgent": self._infectious_agent(),
            "isBasedOn": self._is_based_on(name),
        }

        if self.latest_update:
            record["date"] = self.latest_update
            record["dateModified"] = self.latest_update
        if self.earliest_release:
            record["dateCreated"] = self.earliest_release
            record["datePublished"] = self.earliest_release
        if self.earliest_release or self.latest_update:
            record["temporalCoverage"] = {"@type": "TemporalInterval"}
            if self.earliest_release:
                record["temporalCoverage"]["startDate"] = self.earliest_release
            if self.latest_update:
                record["temporalCoverage"]["endDate"] = self.latest_update

        if self.hosts:
            record["species"] = list(self.hosts.values())
        if self.spatial_coverage:
            record["spatialCoverage"] = list(self.spatial_coverage.values())
        if self.example_works:
            record["exampleOfWork"] = _dedupe(self.example_works)

        sample = self._sample_collection()
        if sample:
            record["sample"] = sample

        is_part_of = self._is_part_of()
        if is_part_of:
            record["isPartOf"] = is_part_of

        has_part = self._has_part()
        if has_part:
            record["hasPart"] = has_part

        if self.authors:
            record["author"] = list(self.authors.values())

        if self.source_databases:
            record["sdPublisher"] = [
                {"@type": "Organization", "name": name}
                for name in _counter_values(self.source_databases, 10)
            ]

        keywords = []
        keywords.extend(_counter_values(self.keywords, 25))
        keywords.extend(_counter_values(self.virus_common_names, 25))
        keywords.extend(_counter_values(self.virus_strains, 25))
        keywords.extend(_counter_values(self.pangolin, 25))
        if keywords:
            record["keywords"] = _dedupe(keywords)

        additional = self._data_collection_additional_properties()
        if additional:
            record["additionalProperty"] = additional

        return record

    def _is_based_on(self, name: str) -> list[dict[str, Any]]:
        action_obj = {
            "@type": "Action",
            "name": "DataCollection Generation Process in the NIAID Data Ecosystem",
            "description": (
                f"How this NCBI Virus {name} DataCollection Record was "
                "generated for the NIAID Data Ecosystem."
            ),
            "actionProcess": {
                "@type": "HowTo",
                "step": copy.deepcopy(HOWTO_STEPS),
            },
            "target": _report_api_url(self.tax_id),
        }

        source_obj = {
            "@type": "nde:ResourceCatalog",
            "name": SOURCE_NAME,
            "url": SOURCE_CATALOG_RECORD,
        }

        return [action_obj, source_obj]


# ---------------------------------------------------------------------------
# SQLite cache / resume support
# ---------------------------------------------------------------------------


class VirusSQLiteCache:
    def __init__(self, seed_taxa: list[str]) -> None:
        self.seed_taxa = seed_taxa
        self.db_path = os.path.join(CACHE_DIR, CACHE_DB)
        self.signature = self._build_signature()
        self.complete_and_fresh = False

    def _build_signature(self) -> str:
        payload = {
            "schema": CACHE_SCHEMA_VERSION,
            "api_base": API_BASE,
            "seed_taxa": self.seed_taxa,
            "page_size": min(PAGE_SIZE, 1000),
            "max_pages_per_seed": MAX_PAGES_PER_SEED,
            "max_records_per_seed": MAX_RECORDS_PER_SEED,
            "max_examples": MAX_EXAMPLES,
            "max_aggregate_values": MAX_AGGREGATE_VALUES,
            "dedup_accessions": DEDUP_ACCESSIONS,
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()

    def _connect(self) -> sqlite3.Connection:
        con = sqlite3.connect(self.db_path, timeout=60)
        con.execute("PRAGMA journal_mode=WAL")
        con.execute("PRAGMA synchronous=NORMAL")
        return con

    def _init_schema(self) -> None:
        os.makedirs(CACHE_DIR, exist_ok=True)
        with self._connect() as con:
            con.execute(
                """CREATE TABLE IF NOT EXISTS metadata (
                    name text NOT NULL PRIMARY KEY,
                    value text NOT NULL
                )"""
            )
            con.execute(
                """CREATE TABLE IF NOT EXISTS taxon_progress (
                    seed_taxon text NOT NULL PRIMARY KEY,
                    status text NOT NULL,
                    next_page_token text,
                    pages_fetched integer NOT NULL DEFAULT 0,
                    records_seen integer NOT NULL DEFAULT 0,
                    total_count integer,
                    updated_at text NOT NULL
                )"""
            )
            con.execute(
                """CREATE TABLE IF NOT EXISTS taxon_aggregates (
                    tax_id text NOT NULL PRIMARY KEY,
                    data text NOT NULL,
                    updated_at text NOT NULL
                )"""
            )
            con.execute(
                """CREATE TABLE IF NOT EXISTS seen_accessions (
                    accession text NOT NULL PRIMARY KEY
                )"""
            )

    def _metadata_value(self, con: sqlite3.Connection, name: str) -> Optional[str]:
        row = con.execute("SELECT value FROM metadata WHERE name=?", (name,)).fetchone()
        return row[0] if row else None

    def _upsert_metadata(self, con: sqlite3.Connection, name: str, value: str) -> None:
        con.execute(
            """INSERT INTO metadata VALUES(?, ?)
               ON CONFLICT(name) DO UPDATE SET value=excluded.value""",
            (name, value),
        )

    def _today(self) -> str:
        return dt.datetime.now(dt.timezone.utc).date().isoformat()

    def _now(self) -> str:
        return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    def _is_expired(self, con: sqlite3.Connection) -> bool:
        date_created = self._metadata_value(con, "date_created")
        if not date_created:
            return True
        try:
            created = dt.date.fromisoformat(date_created)
        except ValueError:
            return True
        today = dt.datetime.now(dt.timezone.utc).date()
        return today - created >= dt.timedelta(days=CACHE_EXPIRE_DAYS)

    def prepare(self) -> None:
        self._init_schema()
        with self._connect() as con:
            cached_signature = self._metadata_value(con, "cache_signature")
            load_complete = self._metadata_value(con, "load_complete") == "true"
            schema_version = self._metadata_value(con, "schema_version")

            if cached_signature != self.signature or schema_version != CACHE_SCHEMA_VERSION:
                reason = "signature changed" if cached_signature else "cache missing"
                self.new_cache(con, reason)
                return

            if self._is_expired(con):
                self.new_cache(con, "cache expired")
                return

            if load_complete:
                logger.info("Using complete NCBI Virus SQLite cache at %s", self.db_path)
                self.complete_and_fresh = True
                return

            if RESUME_INCOMPLETE_CACHE:
                logger.info("Resuming incomplete NCBI Virus SQLite cache at %s", self.db_path)
                return

            self.new_cache(con, "incomplete cache resume disabled")

    def new_cache(self, con: Optional[sqlite3.Connection] = None, reason: str = "requested") -> None:
        close_con = con is None
        con = con or self._connect()
        try:
            logger.info("Creating new NCBI Virus SQLite cache at %s (%s)", self.db_path, reason)
            con.execute("DELETE FROM taxon_progress")
            con.execute("DELETE FROM taxon_aggregates")
            con.execute("DELETE FROM seen_accessions")
            con.execute("DELETE FROM metadata")
            self._upsert_metadata(con, "date_created", self._today())
            self._upsert_metadata(con, "date_updated", self._today())
            self._upsert_metadata(con, "load_complete", "false")
            self._upsert_metadata(con, "schema_version", CACHE_SCHEMA_VERSION)
            self._upsert_metadata(con, "cache_signature", self.signature)
            self.complete_and_fresh = False
        finally:
            if close_con:
                con.commit()
                con.close()

    def mark_complete(self) -> None:
        with self._connect() as con:
            self._upsert_metadata(con, "load_complete", "true")
            self._upsert_metadata(con, "date_updated", self._today())
        self.complete_and_fresh = True
        logger.info("Marked NCBI Virus SQLite cache complete")

    def seed_progress(self, seed_taxon: str) -> dict[str, Any]:
        with self._connect() as con:
            row = con.execute(
                """SELECT status, next_page_token, pages_fetched, records_seen, total_count
                   FROM taxon_progress WHERE seed_taxon=?""",
                (seed_taxon,),
            ).fetchone()
        if not row:
            return {
                "status": "not_started",
                "next_page_token": None,
                "pages_fetched": 0,
                "records_seen": 0,
                "total_count": None,
            }
        return {
            "status": row[0],
            "next_page_token": row[1],
            "pages_fetched": int(row[2] or 0),
            "records_seen": int(row[3] or 0),
            "total_count": row[4],
        }

    def load_accumulator(self, tax_id: str) -> Optional[TaxonAccumulator]:
        with self._connect() as con:
            row = con.execute("SELECT data FROM taxon_aggregates WHERE tax_id=?", (tax_id,)).fetchone()
        if not row:
            return None
        return TaxonAccumulator.from_cache_dict(json.loads(row[0]))

    def has_seen_accession(self, accession: str) -> bool:
        with self._connect() as con:
            row = con.execute("SELECT accession FROM seen_accessions WHERE accession=?", (accession,)).fetchone()
        return row is not None

    def save_page(
        self,
        seed_taxon: str,
        accumulators: dict[str, TaxonAccumulator],
        touched_taxa: set[str],
        seen_accessions: list[str],
        next_page_token: Optional[str],
        pages_fetched: int,
        records_seen: int,
        total_count: Optional[int],
        status: str,
    ) -> None:
        now = self._now()
        with self._connect() as con:
            for tax_id in sorted(touched_taxa, key=lambda value: int(value) if value.isdigit() else value):
                con.execute(
                    """INSERT INTO taxon_aggregates VALUES(?, ?, ?)
                       ON CONFLICT(tax_id) DO UPDATE SET
                       data=excluded.data,
                       updated_at=excluded.updated_at""",
                    (
                        tax_id,
                        json.dumps(
                            accumulators[tax_id].to_cache_dict(),
                            sort_keys=True,
                            separators=(",", ":"),
                        ),
                        now,
                    ),
                )

            for accession in seen_accessions:
                con.execute(
                    "INSERT OR IGNORE INTO seen_accessions VALUES(?)",
                    (accession,),
                )

            con.execute(
                """INSERT INTO taxon_progress VALUES(?, ?, ?, ?, ?, ?, ?)
                   ON CONFLICT(seed_taxon) DO UPDATE SET
                   status=excluded.status,
                   next_page_token=excluded.next_page_token,
                   pages_fetched=excluded.pages_fetched,
                   records_seen=excluded.records_seen,
                   total_count=excluded.total_count,
                   updated_at=excluded.updated_at""",
                (
                    seed_taxon,
                    status,
                    next_page_token,
                    pages_fetched,
                    records_seen,
                    total_count,
                    now,
                ),
            )

        logger.info(
            "Cached NCBI Virus seed %s page progress: status=%s pages=%s records=%s touched_taxa=%s",
            seed_taxon,
            status,
            pages_fetched,
            records_seen,
            len(touched_taxa),
        )

    def iter_accumulators(self) -> Iterator[TaxonAccumulator]:
        con = self._connect()
        try:
            cursor = con.execute("SELECT data FROM taxon_aggregates ORDER BY CAST(tax_id AS INTEGER)")
            count = 0
            for row in cursor:
                count += 1
                if count % 10000 == 0:
                    logger.info("Loaded %s NCBI Virus cached aggregate records", count)
                yield TaxonAccumulator.from_cache_dict(json.loads(row[0]))
            logger.info("Finished loading %s NCBI Virus cached aggregate records", count)
        finally:
            con.close()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def _tax_id_sort_key(value: str) -> Any:
    return int(value) if value.isdigit() else value


def _add_report_to_accumulators(
    report: dict[str, Any],
    accumulators: dict[str, TaxonAccumulator],
    cache: Optional[VirusSQLiteCache] = None,
) -> Optional[str]:
    accession = _report_accession(report)
    virus = report.get("virus") or {}
    tax_id = _clean_string(virus.get("tax_id"))
    if not tax_id:
        logger.warning("Skipping NCBI Virus report without virus.tax_id: %s", accession)
        return None

    accumulator = accumulators.get(tax_id)
    if accumulator is None and cache is not None:
        accumulator = cache.load_accumulator(tax_id)
    if accumulator is None:
        accumulator = TaxonAccumulator(tax_id)
    accumulators[tax_id] = accumulator
    accumulator.add(report)
    return tax_id


def _parse_without_cache(seed_taxa: list[str]) -> Iterator[dict[str, Any]]:
    accumulators: dict[str, TaxonAccumulator] = {}
    seen_accessions: set[str] = set()
    logger.info("Starting NCBI Virus parse with %d seed taxa", len(seed_taxa))

    for seed_taxon in seed_taxa:
        seed_records = 0
        for report in _fetch_reports_for_seed(seed_taxon):
            accession = _report_accession(report)
            if DEDUP_ACCESSIONS and accession:
                if accession in seen_accessions:
                    continue
                seen_accessions.add(accession)

            if _add_report_to_accumulators(report, accumulators):
                seed_records += 1

        logger.info("Seed taxon %s contributed %d source records", seed_taxon, seed_records)

    logger.info("Built %d NCBI Virus DataCollection groups", len(accumulators))
    for tax_id in sorted(accumulators, key=_tax_id_sort_key):
        yield accumulators[tax_id].to_record()


def _parse_with_cache(seed_taxa: list[str], cache: VirusSQLiteCache) -> Iterator[dict[str, Any]]:
    cache.prepare()
    if cache.complete_and_fresh:
        for accumulator in cache.iter_accumulators():
            yield accumulator.to_record()
        return

    logger.info("Starting cached NCBI Virus parse with %d seed taxa", len(seed_taxa))
    for seed_taxon in seed_taxa:
        progress = cache.seed_progress(seed_taxon)
        if progress["status"] in {"complete", "limited"}:
            logger.info(
                "Skipping cached NCBI Virus seed %s; status=%s pages=%s records=%s",
                seed_taxon,
                progress["status"],
                progress["pages_fetched"],
                progress["records_seen"],
            )
            continue

        page_token = progress["next_page_token"]
        pages_fetched = progress["pages_fetched"]
        records_seen = progress["records_seen"]
        total_count = progress["total_count"]
        accumulators: dict[str, TaxonAccumulator] = {}

        while True:
            logger.info("Fetching NCBI Virus reports for taxon %s page %s", seed_taxon, pages_fetched + 1)
            data = _fetch_report_page(seed_taxon, page_token)
            reports = data.get("reports") or []
            total_count = data.get("total_count", total_count)
            next_page_token = data.get("next_page_token") or data.get("nextPageToken")

            if pages_fetched == 0:
                logger.info("Taxon %s report total_count: %s", seed_taxon, total_count)

            if not reports:
                cache.save_page(
                    seed_taxon=seed_taxon,
                    accumulators=accumulators,
                    touched_taxa=set(),
                    seen_accessions=[],
                    next_page_token=None,
                    pages_fetched=pages_fetched,
                    records_seen=records_seen,
                    total_count=total_count,
                    status="complete",
                )
                break

            touched_taxa: set[str] = set()
            page_seen_accessions: list[str] = []
            page_seen_lookup: set[str] = set()
            accepted_count = 0

            for report in reports:
                if MAX_RECORDS_PER_SEED and records_seen + accepted_count >= MAX_RECORDS_PER_SEED:
                    break

                accession = _report_accession(report)
                if DEDUP_ACCESSIONS and accession:
                    if accession in page_seen_lookup or cache.has_seen_accession(accession):
                        continue
                    page_seen_lookup.add(accession)
                    page_seen_accessions.append(accession)

                tax_id = _add_report_to_accumulators(report, accumulators, cache=cache)
                if tax_id:
                    touched_taxa.add(tax_id)
                    accepted_count += 1

            records_seen += accepted_count
            pages_fetched += 1

            status = "in_progress"
            token_to_store = next_page_token
            if MAX_RECORDS_PER_SEED and records_seen >= MAX_RECORDS_PER_SEED:
                status = "limited"
                token_to_store = None
                logger.info("Reached NCBI_VIRUS_MAX_RECORDS_PER_SEED=%s for %s", MAX_RECORDS_PER_SEED, seed_taxon)
            elif MAX_PAGES_PER_SEED and pages_fetched >= MAX_PAGES_PER_SEED:
                status = "limited"
                token_to_store = None
                logger.info("Reached NCBI_VIRUS_MAX_PAGES_PER_SEED=%s for %s", MAX_PAGES_PER_SEED, seed_taxon)
            elif not next_page_token:
                status = "complete"
                token_to_store = None

            cache.save_page(
                seed_taxon=seed_taxon,
                accumulators=accumulators,
                touched_taxa=touched_taxa,
                seen_accessions=page_seen_accessions,
                next_page_token=token_to_store,
                pages_fetched=pages_fetched,
                records_seen=records_seen,
                total_count=total_count,
                status=status,
            )

            if status in {"complete", "limited"}:
                break

            page_token = next_page_token
            _sleep_between_requests()

    cache.mark_complete()
    for accumulator in cache.iter_accumulators():
        yield accumulator.to_record()


def parse() -> Iterator[dict[str, Any]]:
    """Yield NCBI Virus DataCollection records grouped by virus.tax_id."""
    seed_taxa = _discover_seed_taxa()
    if USE_SQL_CACHE:
        yield from _parse_with_cache(seed_taxa, VirusSQLiteCache(seed_taxa))
    else:
        yield from _parse_without_cache(seed_taxa)


if __name__ == "__main__":
    for record in parse():
        print(json.dumps(record, sort_keys=True))
