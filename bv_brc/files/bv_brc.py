#!/usr/bin/env python3
"""
BV-BRC Genomes DataCollection crawler for the NIAID Data Ecosystem.

Groups Bacterial and Viral Bioinformatics Resource Center genome records
by NCBI Taxonomy ID and generates one DataCollection record per organism.

Data source
-----------
Bulk genome metadata is downloaded from the BV-BRC FTPS server:
    ftps://ftp.bv-brc.org/RELEASE_NOTES/genome_metadata
    ftps://ftp.bv-brc.org/RELEASE_NOTES/genome_summary

Taxon names are batch-resolved from the BV-BRC REST API:
    https://www.bv-brc.org/api/taxonomy/?in(taxon_id,(...))
"""

import copy
import datetime as dt
import json
import logging
import re
import ssl
import time
import urllib.request
from collections import defaultdict
from ftplib import FTP, FTP_TLS
from typing import Any, Optional

from sql_database import NDEDatabase

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")

# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------

API_BASE = "https://www.bv-brc.org/api"

# End-user URL (used as the DataCollection url)
BVBRC_VIEW_URL_TEMPLATE = (
    "https://www.bv-brc.org/view/Taxonomy/{code}#view_tab=genomes"
)

# API version URL
API_VERSION_URL = f"{API_BASE}/"

# ---------------------------------------------------------------------------
# FTP configuration
# ---------------------------------------------------------------------------

FTP_HOST = "ftp.bv-brc.org"
FTP_USER = "anonymous"
FTP_PASS = "guest"
GENOME_METADATA_PATH = "RELEASE_NOTES/genome_metadata"
GENOME_SUMMARY_PATH = "RELEASE_NOTES/genome_summary"
GENOME_LINEAGE_PATH = "RELEASE_NOTES/genome_lineage"

# genome_metadata column indices (0-based)
META_TAXON_ID = 3
META_COMPLETION_DATE = 13
META_COLLECTION_DATE = 37
META_ISOLATION_COUNTRY = 38
META_HOST_NAME = 45
META_DISEASE = 63
META_MIN_COLS = 64

# genome_summary column indices (0-based)
SUMMARY_TAXON_ID = 2
SUMMARY_DATE_MODIFIED = 19
SUMMARY_MIN_COLS = 20

# genome_lineage column indices (0-based)
LINEAGE_TAXON_ID = 2
LINEAGE_SPECIES = 9
LINEAGE_MIN_COLS = 10

# ---------------------------------------------------------------------------
# Curated / static metadata
# ---------------------------------------------------------------------------

SOURCE_METADATA = {
    "name": "BV-BRC (as DataCollection records)",
    "url": "https://www.bv-brc.org/",
    "description": (
        "BV-BRC is a comprehensive resource for bacterial and viral "
        "genomics, bioinformatics, and data analysis. The goal of BV-BRC "
        "is to provide researchers with access to high-quality, curated "
        "data and bioinformatics tools to support their research in the "
        "field of infectious diseases."
    ),
    "license": (
        "https://www.bv-brc.org/docs/system_documentation/"
        "data_management_sharing.html"
    ),
}

_DDE_CURATED_BY = {
    "name": "Data Discovery Engine",
    "url": "https://discovery.biothings.io/",
    "dateModified": "2026-01-03",
}

ABOUT_DEFINED_TERM = {
    "@type": "DefinedTerm",
    "description": "Subclass of BioChemEntity",
    "displayName": "Genome",
    "name": "Genome",
    "url": "http://purl.obolibrary.org/obo/NCIT_C16629",
}

MEASUREMENT_TECHNIQUE = [
    {
        "@type": "DefinedTerm",
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Curation",
        "url": "http://purl.obolibrary.org/obo/NCIT_C48292",
    },
    {
        "@type": "DefinedTerm",
        "inDefinedTermSet": "EDAM",
        "isCurated": True,
        "name": "Annotation",
        "url": "http://edamontology.org/operation_0226",
    },
]

VARIABLE_MEASURED = [
    {
        "@type": "DefinedTerm",
        "inDefinedTermSet": "EFO",
        "isCurated": True,
        "name": "genomic measurement",
        "url": "http://www.ebi.ac.uk/efo/EFO_0004554",
    },
]

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "name": "Genomics",
        "identifier": "topic_0622",
        "url": "http://edamontology.org/topic_0622",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Computational biology",
        "identifier": "topic_2259",
        "url": "http://edamontology.org/topic_2259",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Bioinformatics",
        "identifier": "topic_0091",
        "url": "http://edamontology.org/topic_0091",
        "inDefinedTermSet": "EDAM",
    },
]

ENCODING_FORMAT = [
    {
        "@type": "DefinedTerm",
        "name": "JSON",
        "url": "http://edamontology.org/format_3464",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "TSV",
        "url": "http://edamontology.org/format_3475",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "CSV",
        "url": "http://edamontology.org/format_3752",
        "inDefinedTermSet": "EDAM",
    },
]

HOWTO_STEPS = [
    (
        "Step 1: Download bulk genome data from the BV-BRC FTPS server "
        "(ftps://ftp.bv-brc.org/RELEASE_NOTES/). Stream and aggregate "
        "three tab-separated files: 'genome_metadata' for per-taxon "
        "genome counts, isolation countries, host organisms, diseases, "
        "collection dates, and earliest completion dates; "
        "'genome_summary' for the latest date_modified per taxon; and "
        "'genome_lineage' for species names per taxon_id."
    ),
    (
        "Step 2: Parse the aggregated records to generate an organism-"
        "specific DataCollection per NCBI Taxonomy ID. Use the latest "
        "'date_modified' from genome_summary for 'dateModified' and "
        "the earliest 'completion_date' from genome_metadata for "
        "'dateCreated'. Generate 'temporalCoverage' from the "
        "'collection_date' values across all records for that organism. "
        "Use the count of genome records for 'collectionSize'. Generate "
        "the values for 'species', 'infectiousAgent', 'url', 'name', "
        "and 'description' properties based on the information "
        "retrieved for the NCBI taxonomy code and any templated text."
    ),
    (
        "Step 3. Fill in manually curated fields from the BV-BRC "
        "resource catalog record such as 'measurementTechnique', "
        "'variableMeasured', 'topicCategory', 'conditionsOfAccess', "
        "'usageInfo'."
    ),
]

MAX_RETRIES = 3
RETRY_BACKOFF = 10  # seconds, doubled each retry


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _fetch_json(url: str, timeout: int = 120, retries: int = MAX_RETRIES) -> Any:
    """Fetch a URL and return parsed JSON, with retry on timeout."""
    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "nde-bvbrc-crawler/0.1"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read())
        except Exception as exc:
            last_exc = exc
            if attempt < retries:
                wait = RETRY_BACKOFF * (2 ** (attempt - 1))
                logger.warning(
                    "Attempt %d/%d failed for %s: %s — retrying in %ds",
                    attempt, retries, url[:120], exc, wait,
                )
                time.sleep(wait)
    raise last_exc  # type: ignore[misc]


def _parse_datetime(value: Optional[str]) -> Optional[dt.datetime]:
    if not value:
        return None
    v = value.strip()
    if not v:
        return None
    try:
        if v.endswith("Z"):
            v = v[:-1] + "+00:00"
        parsed = dt.datetime.fromisoformat(v)
        if parsed.tzinfo is None:
            return parsed.replace(tzinfo=dt.timezone.utc)
        return parsed.astimezone(dt.timezone.utc)
    except Exception:
        return None


def _iso_date(value: Optional[dt.datetime]) -> Optional[str]:
    if value is None:
        return None
    return value.date().isoformat()


def _make_record_id(code: str) -> str:
    return f"bvbrc_genome_{code}"


def _fetch_api_version() -> str:
    """Retrieve the BV-BRC API version string from the HTML landing page."""
    try:
        req = urllib.request.Request(
            API_VERSION_URL,
            headers={"User-Agent": "nde-bvbrc-crawler/0.1"},
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            text = resp.read().decode("utf-8", errors="replace")
        # The page contains: <div>Version: 1.9.2</div>
        match = re.search(r"Version:\s*([\d.]+)", text)
        return match.group(1) if match else ""
    except Exception as exc:
        logger.warning("Could not fetch API version: %s", exc)
        return ""


# ---------------------------------------------------------------------------
# FTP bulk download and local aggregation
# ---------------------------------------------------------------------------


def _connect_ftp() -> FTP_TLS:
    """Open an encrypted FTPS connection to BV-BRC.

    Uses a custom FTP_TLS subclass that reuses the control channel's TLS
    session on data connections — required by BV-BRC's FTPS server.
    """

    class _SessionReuseFTP(FTP_TLS):
        """FTP_TLS that reuses the control TLS session on data channels."""

        def ntransfercmd(self, cmd, rest=None):
            conn, size = FTP.ntransfercmd(self, cmd, rest)
            if self._prot_p:
                conn = self.context.wrap_socket(
                    conn, server_hostname=self.host,
                    session=self.sock.session,
                )
            return conn, size

    ftp = _SessionReuseFTP()
    ftp.connect(FTP_HOST, 21, timeout=300)
    ftp.login(FTP_USER, FTP_PASS)
    ftp.prot_p()  # secure data channel
    return ftp


def _stream_ftp_lines(ftp: FTP_TLS, path: str):
    """Yield decoded text lines from an FTPS file (streaming, not buffered)."""
    sock = ftp.transfercmd(f"RETR {path}")
    fp = sock.makefile("rb")
    try:
        for raw_line in fp:
            yield raw_line.decode("utf-8", errors="replace").rstrip("\r\n")
    finally:
        fp.close()
        sock.close()
        try:
            ftp.voidresp()
        except Exception:
            pass


def _aggregate_genome_metadata() -> dict[str, dict[str, Any]]:
    """
    Stream genome_metadata from FTP and aggregate per taxon_id.

    Returns dict of taxon_id -> {
        genome_count, countries, hosts, diseases,
        collection_years, earliest_completion
    }
    """
    taxa: dict[str, dict[str, Any]] = defaultdict(lambda: {
        "genome_count": 0,
        "countries": set(),
        "hosts": set(),
        "diseases": set(),
        "collection_years": set(),
        "earliest_completion": None,
    })

    ftp = _connect_ftp()
    line_count = 0
    start_time = time.time()

    try:
        for line in _stream_ftp_lines(ftp, GENOME_METADATA_PATH):
            if line_count == 0:
                line_count += 1
                continue  # skip header

            fields = line.split("\t")
            if len(fields) <= META_DISEASE:
                line_count += 1
                continue

            taxon_id = fields[META_TAXON_ID].strip()
            if not taxon_id:
                line_count += 1
                continue

            entry = taxa[taxon_id]
            entry["genome_count"] += 1

            country = fields[META_ISOLATION_COUNTRY].strip()
            if country:
                entry["countries"].add(country)

            host = fields[META_HOST_NAME].strip()
            if host:
                entry["hosts"].add(host)

            disease_raw = fields[META_DISEASE].strip()
            if disease_raw:
                for d in disease_raw.split(";"):
                    d = d.strip()
                    if d:
                        entry["diseases"].add(d)

            coll_date = fields[META_COLLECTION_DATE].strip()
            if coll_date:
                try:
                    year = int(coll_date[:4])
                    if 1000 <= year <= 2100:
                        entry["collection_years"].add(year)
                except (ValueError, IndexError):
                    pass

            comp_date = fields[META_COMPLETION_DATE].strip()
            if comp_date and not comp_date.startswith("1900"):
                cur = entry["earliest_completion"]
                if cur is None or comp_date < cur:
                    entry["earliest_completion"] = comp_date

            line_count += 1
            if line_count % 2_000_000 == 0:
                elapsed = time.time() - start_time
                logger.info(
                    "genome_metadata: %dM lines, %d taxa (%.0fs)",
                    line_count // 1_000_000, len(taxa), elapsed,
                )
    finally:
        try:
            ftp.quit()
        except Exception:
            pass

    elapsed = time.time() - start_time
    logger.info(
        "genome_metadata complete: %d lines, %d unique taxa in %.0fs",
        line_count, len(taxa), elapsed,
    )
    return dict(taxa)


def _aggregate_genome_summary() -> dict[str, str]:
    """
    Stream genome_summary from FTP and return the latest date_modified
    per taxon_id.
    """
    latest_modified: dict[str, str] = {}

    ftp = _connect_ftp()
    line_count = 0
    start_time = time.time()

    try:
        for line in _stream_ftp_lines(ftp, GENOME_SUMMARY_PATH):
            if line_count == 0:
                line_count += 1
                continue  # skip header

            fields = line.split("\t")
            if len(fields) <= SUMMARY_DATE_MODIFIED:
                line_count += 1
                continue

            taxon_id = fields[SUMMARY_TAXON_ID].strip()
            date_mod = fields[SUMMARY_DATE_MODIFIED].strip()

            if taxon_id and date_mod:
                cur = latest_modified.get(taxon_id)
                if cur is None or date_mod > cur:
                    latest_modified[taxon_id] = date_mod

            line_count += 1
            if line_count % 2_000_000 == 0:
                elapsed = time.time() - start_time
                logger.info(
                    "genome_summary: %dM lines (%.0fs)",
                    line_count // 1_000_000, elapsed,
                )
    finally:
        try:
            ftp.quit()
        except Exception:
            pass

    elapsed = time.time() - start_time
    logger.info(
        "genome_summary complete: %d lines, %d taxa with dates in %.0fs",
        line_count, len(latest_modified), elapsed,
    )
    return latest_modified


def _extract_taxon_names_ftp() -> dict[str, str]:
    """
    Stream genome_lineage from FTP and extract the 'species' name for each
    taxon_id. Returns {taxon_id: species_name}. Keeps the first (and
    typically only) species name seen per taxon.
    """
    names: dict[str, str] = {}

    ftp = _connect_ftp()
    line_count = 0
    start_time = time.time()

    try:
        for line in _stream_ftp_lines(ftp, GENOME_LINEAGE_PATH):
            if line_count == 0:
                line_count += 1
                continue  # skip header

            fields = line.split("\t")
            if len(fields) <= LINEAGE_SPECIES:
                line_count += 1
                continue

            taxon_id = fields[LINEAGE_TAXON_ID].strip()
            species = fields[LINEAGE_SPECIES].strip()

            if taxon_id and species and taxon_id not in names:
                names[taxon_id] = species

            line_count += 1
            if line_count % 2_000_000 == 0:
                elapsed = time.time() - start_time
                logger.info(
                    "genome_lineage: %dM lines, %d names (%.0fs)",
                    line_count // 1_000_000, len(names), elapsed,
                )
    finally:
        try:
            ftp.quit()
        except Exception:
            pass

    elapsed = time.time() - start_time
    logger.info(
        "genome_lineage complete: %d lines, %d taxon names in %.0fs",
        line_count, len(names), elapsed,
    )
    return names


# ---------------------------------------------------------------------------
# Record building
# ---------------------------------------------------------------------------


def _build_organism_term(code: str, name: str) -> dict[str, Any]:
    return {
        "@type": "DefinedTerm",
        "name": name,
        "identifier": f"taxonomy:{code}",
        "inDefinedTermSet": "NCBI Taxonomy",
        "url": f"https://www.ncbi.nlm.nih.gov/taxonomy/{code}",
    }


def _build_spatial_coverage(countries: set[str]) -> list[dict[str, str]]:
    result = []
    for country in sorted(countries):
        result.append({
            "@type": "AdministrativeArea",
            "administrativeType": "country",
            "locationType": "collection",
            "name": country,
        })
    return result


def _build_species_list(hosts: set[tuple[str, str]]) -> list[dict[str, str]]:
    result = []
    seen = set()
    for host_name, host_sci in sorted(hosts):
        if host_name not in seen:
            result.append({"name": host_name})
            seen.add(host_name)
    return result


def _build_temporal_coverage(
    collection_years: list[int],
) -> Optional[dict[str, Any]]:
    if not collection_years:
        return None
    earliest = min(collection_years)
    latest = max(collection_years)
    return {
        "@type": "TemporalInterval",
        "startDate": str(earliest),
        "endDate": str(latest),
        "temporalType": "collection",
        "name": "Collection Period",
    }


def _build_example_of_work(code: str) -> dict[str, Any]:
    return {
        "@type": "CreativeWork",
        "about": {
            "@type": "DefinedTerm",
            "displayName": "Genome",
            "name": "Genome",
            "url": "http://purl.obolibrary.org/obo/NCIT_C16629",
            "termCode": "NCIT:C16629",
            "inDefinedTermSet": "NCIT",
        },
        "encodingFormat": copy.deepcopy(ENCODING_FORMAT),
        "schemaVersion": f"{API_BASE}/doc/genome",
        "potentialAction": {
            "@type": "Action",
            "name": "Use API call for an example record",
            "target": (
                f"{API_BASE}/genome/"
                f"?eq(taxon_id,{code})"
                "&select(genome_id,genome_name,taxon_id,isolation_country,"
                "host_common_name,host_name,host_scientific_name,"
                "host_taxon_id)"
                "&sort(+genome_id)&http_accept=application/json&limit(1,0)"
            ),
        },
    }


def _build_is_based_on(about_name: str) -> list[dict[str, Any]]:
    action_obj = {
        "@type": "Action",
        "name": (
            "DataCollection Generation Process in the NIAID Data Ecosystem"
        ),
        "description": (
            f"How this BV-BRC {about_name} DataCollection Record was "
            "generated for the NIAID Data Ecosystem."
        ),
        "actionProcess": {
            "@type": "HowTo",
            "step": copy.deepcopy(HOWTO_STEPS),
        },
    }

    source_obj = {
        "@type": "nde:ResourceCatalog",
        "name": "Bacterial and Viral Bioinformatics Resource Center",
        "url": "https://data.niaid.nih.gov/resources?id=dde_42e839db86d4166d",
    }

    return [action_obj, source_obj]


def _build_description(
    taxon_name: str,
    url: str,
    version: str,
    countries: set[str],
    hosts: set[tuple[str, str]],
) -> str:
    parts = []
    parts.append(
        f"{taxon_name} genomes publicly available from BV-BRC"
    )
    if version:
        parts[-1] += f" API version {version}"
    parts[-1] += "."

    if countries:
        country_list = ", ".join(sorted(countries))
        parts.append(
            f"These genomes were collected from {country_list}"
        )

    if hosts:
        host_names = sorted(set(h[0] for h in hosts if h[0]))
        if host_names:
            host_list = ", ".join(host_names)
            parts.append(
                f"from various host organisms including: {host_list}"
            )

    parts.append(f"For more details, visit: {url}")
    return " ".join(parts)


def _build_data_collection(
    code: str,
    taxon_name: str,
    genome_count: int,
    group_data: dict[str, Any],
    version: str,
) -> dict[str, Any]:
    url = BVBRC_VIEW_URL_TEMPLATE.format(code=code)

    earliest = group_data["earliest_inserted"]
    latest = group_data["latest_modified"]
    countries = group_data["countries"]
    hosts = group_data["hosts"]
    collection_years = group_data["collection_years"]

    name = f"BV-BRC {taxon_name} genomes"
    description = _build_description(
        taxon_name, url, version, countries, hosts,
    )

    record: dict[str, Any] = {
        "_id": _make_record_id(code),
        "@type": "DataCollection",
        "about": copy.deepcopy(ABOUT_DEFINED_TERM),
        "includedInDataCatalog": [
            {
                "@type": "DataCatalog",
                "name": "Bacterial and Viral Bioinformatics Resource Center",
                "alternateName": "BV-BRC",
                "url": "https://www.bv-brc.org/",
                "archivedAt": url,
            }
        ],
        "name": name,
        "url": url,
        "description": description,
        "collectionSize": {
            "minValue": genome_count,
            "unitText": "Genomes",
            "@type": "QuantitativeValue",
        },
        "date": _iso_date(latest) or _iso_date(earliest),
        "dateModified": _iso_date(latest),
        "dateCreated": _iso_date(earliest),
        "conditionsOfAccess": "Open",
        "isAccessibleForFree": True,
        "license": (
            "https://www.bv-brc.org/docs/system_documentation/"
            "data_management_sharing.html"
        ),
        "usageInfo": {
            "url": (
                "https://www.bv-brc.org/docs/system_documentation/"
                "data_management_sharing.html"
            ),
        },
        "creditText": (
            "see https://www.bv-brc.org/citation for how to cite "
            "this resource"
        ),
        "measurementTechnique": copy.deepcopy(MEASUREMENT_TECHNIQUE),
        "variableMeasured": copy.deepcopy(VARIABLE_MEASURED),
        "topicCategory": copy.deepcopy(TOPIC_CATEGORY),
        "exampleOfWork": _build_example_of_work(code),
        "isBasedOn": _build_is_based_on(ABOUT_DEFINED_TERM["name"]),
    }

    # Infectious agent (the organism itself)
    record["infectiousAgent"] = _build_organism_term(code, taxon_name)

    # Species (hosts)
    host_list = _build_species_list(hosts)
    if host_list:
        record["species"] = host_list

    # Spatial coverage
    spatial = _build_spatial_coverage(countries)
    if spatial:
        record["spatialCoverage"] = spatial

    # Temporal coverage from collection years
    temporal = _build_temporal_coverage(collection_years)
    if temporal:
        record["temporalCoverage"] = temporal

    # Health conditions (diseases)
    if group_data["diseases"]:
        record["healthCondition"] = [
            {"@type": "DefinedTerm", "name": d}
            for d in sorted(group_data["diseases"])
        ]

    # Version
    # if version:
        # record["version"] = version

    # Funding
    record["funding"] = [
        {
            "@type": "MonetaryGrant",
            "funder": {
                "@type": "Organization",
                "name": (
                    "National Institute of Allergy and Infectious Diseases, "
                    "National Institutes of Health, "
                    "Department of Health and Human Services"
                ),
            },
        }
    ]

    return record


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------


class BvBrc(NDEDatabase):
    SQL_DB = "bv_brc_solar.db"
    EXPIRE = dt.timedelta(days=30)

    def __init__(self):
        super().__init__()
        self._version = ""

    def load_cache(self):
        """Download BV-BRC genome data via FTP bulk files and yield (_id, json_str) tuples."""
        overall_start = time.time()

        # Fetch the API version
        self._version = _fetch_api_version()
        if self._version:
            logger.info("BV-BRC API version: %s", self._version)

        # Step 1 — Download and aggregate genome_metadata from FTP
        logger.info("Step 1/3: Downloading genome_metadata from FTP (%s)...", FTP_HOST)
        meta_agg = _aggregate_genome_metadata()
        if not meta_agg:
            logger.error("No genome records found in FTP genome_metadata file.")
            return

        # Step 2 — Download and aggregate genome_summary for date_modified
        logger.info("Step 2/3: Downloading genome_summary from FTP...")
        date_mods = _aggregate_genome_summary()

        # Step 3 — Download genome_lineage for taxon names (no API needed)
        logger.info("Step 3/3: Downloading genome_lineage from FTP for taxon names...")
        taxon_names = _extract_taxon_names_ftp()

        # Build and yield cache entries
        taxon_ids = sorted(meta_agg.keys())
        total_taxa = len(taxon_ids)
        count = 0
        for taxon_id in taxon_ids:
            agg = meta_agg[taxon_id]
            taxon_name = taxon_names.get(taxon_id, f"Organism {taxon_id}")

            earliest = agg["earliest_completion"]
            if earliest:
                parsed = _parse_datetime(earliest)
                earliest = _iso_date(parsed)

            latest = date_mods.get(taxon_id)
            if latest:
                parsed = _parse_datetime(latest)
                latest = _iso_date(parsed)

            cache_entry = {
                "code": taxon_id,
                "taxon_name": taxon_name,
                "genome_count": agg["genome_count"],
                "version": self._version,
                "earliest_inserted": earliest,
                "latest_modified": latest,
                "countries": sorted(agg["countries"]),
                "hosts": sorted([[h, h] for h in agg["hosts"]]),
                "collection_years": sorted(agg["collection_years"]),
                "diseases": sorted(agg["diseases"]),
                "record_count": agg["genome_count"],
            }

            cache_id = _make_record_id(taxon_id)
            yield (cache_id, json.dumps(cache_entry))
            count += 1

            if count % 5000 == 0:
                logger.info("Built %d / %d cache entries", count, total_taxa)

        elapsed = time.time() - overall_start
        logger.info(
            "Finished loading cache. Total: %d taxa, Time: %.0f min",
            count, elapsed / 60,
        )

    def parse(self, records):
        """Parse cached records into NDE-schema DataCollection documents."""
        count = 0
        errors = 0
        for record in records:
            cache_id = record[0]
            try:
                data = json.loads(record[1])

                # Reconstruct the group_data dict from cached JSON
                group_data = {
                    "earliest_inserted": (
                        dt.datetime.fromisoformat(data["earliest_inserted"])
                        if data.get("earliest_inserted") else None
                    ),
                    "latest_modified": (
                        dt.datetime.fromisoformat(data["latest_modified"])
                        if data.get("latest_modified") else None
                    ),
                    "countries": set(data.get("countries", [])),
                    "hosts": set(tuple(h) for h in data.get("hosts", [])),
                    "collection_years": data.get("collection_years", []),
                    "diseases": set(data.get("diseases", [])),
                    "record_count": data.get("record_count", 0),
                }

                result = _build_data_collection(
                    data["code"],
                    data["taxon_name"],
                    data["genome_count"],
                    group_data,
                    data.get("version", ""),
                )
                if result:
                    yield result
                    count += 1
                    if count % 500 == 0:
                        logger.info("Parsed %s records (errors: %s)", count, errors)
            except Exception as e:
                errors += 1
                logger.error("Error parsing %s: %s", cache_id, e)
                continue

        logger.info("Finished parsing. Total records: %s, Errors: %s", count, errors)
