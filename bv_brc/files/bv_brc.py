#!/usr/bin/env python3
"""
BV-BRC Genomes DataCollection crawler for the NIAID Data Ecosystem.

Groups Bacterial and Viral Bioinformatics Resource Center genome records
by NCBI Taxonomy ID and generates one DataCollection record per organism.

API reference
-------------
Taxonomy table (species with genomes):
    https://www.bv-brc.org/api/taxonomy/
    ?gt(genomes,0)&eq(taxon_rank,species)
    &select(taxon_id,taxon_name,genomes)&limit(1000)

Genomes for a given taxonomy code (paginated):
    https://www.bv-brc.org/api/genome/
    ?in(taxon_lineage_ids,({CODE}))
    &select(genome_id,genome_name,taxon_id,isolation_country,
            host_common_name,host_name,host_scientific_name,
            host_taxon_id,date_inserted,date_modified,collection_year)
    &sort(+genome_id)&http_accept=application/json&limit(1000,0)
"""

import copy
import datetime as dt
import json
import logging
import re
import time
import urllib.parse
import urllib.request
from typing import Any, Iterator, Optional

import requests
from requests.adapters import HTTPAdapter
from sql_database import NDEDatabase
from urllib3.util.retry import Retry

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

TAXONOMY_URL = (
    f"{API_BASE}/taxonomy/"
    "?gt(genomes,0)&eq(taxon_rank,species)"
    "&select(taxon_id,taxon_name,genomes)"
    "&sort(+taxon_id)"
    "&http_accept=application/json"
)

GENOME_URL_TEMPLATE = (
    f"{API_BASE}/genome/"
    "?eq(taxon_id,{code})"
    "&select(genome_id,genome_name,taxon_id,"
    "isolation_country,host_common_name,host_name,"
    "host_scientific_name,host_taxon_id,"
    "date_inserted,date_modified,collection_year,disease)"
    "&sort(+genome_id)"
    "&http_accept=application/json"
    "&limit({limit},{offset})"
)

# Count query (solr format) — returns numFound without fetching all docs
GENOME_COUNT_URL_TEMPLATE = (
    f"{API_BASE}/genome/"
    "?eq(taxon_id,{code})"
    "&select(genome_id)"
    "&limit(1,0)"
    "&http_accept=application/solr+json"
)

# End-user URL (used as the DataCollection url)
BVBRC_VIEW_URL_TEMPLATE = (
    "https://www.bv-brc.org/view/Taxonomy/{code}#view_tab=genomes"
)

# API version URL
API_VERSION_URL = f"{API_BASE}/"

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
        "Step 1: For each NCBI Taxonomy code, pull all the BV-BRC "
        "Genome records associated with that NCBI Taxonomy code "
        "(paginating as needed): "
        "https://www.bv-brc.org/api/genome/"
        "?eq(taxon_id,{NCBI_TAXONOMY_ID})"
        "&select(genome_id,genome_name,taxon_id,isolation_country,"
        "host_common_name,host_name,host_scientific_name,host_taxon_id)"
        "&sort(+genome_id)&http_accept=application/json&limit(1000,0)"
    ),
    (
        "Step 2: Parse the records to generate an organism-specific "
        "DataCollection. Use the latest 'date_modified' value for "
        "'dateModified' and the earliest 'date_inserted' value for "
        "'dateCreated'. Generate the 'temporalCoverage' from the "
        "date_inserted info across all records for that organism. Use "
        "the count of the records for 'collectionSize'. Generate the "
        "values for the 'species', 'infectiousAgent', 'url', 'name', "
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

ROWS_PER_PAGE = 1000
REQUEST_DELAY = 0.5  # seconds between API calls
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
# API fetching
# ---------------------------------------------------------------------------


def _fetch_taxonomy_species() -> list[dict[str, Any]]:
    """
    Page through the BV-BRC taxonomy API to get all species-level taxa
    that have at least one genome.
    Returns list of {taxon_id, taxon_name, genomes, lineage_names}.
    """
    all_taxa: list[dict[str, Any]] = []
    offset = 0
    page_size = 1000

    while True:
        url = f"{TAXONOMY_URL}&limit({page_size},{offset})"
        logger.info(
            "Fetching taxonomy species, offset %d (%d so far)...",
            offset,
            len(all_taxa),
        )

        try:
            data = _fetch_json(url)
        except Exception as exc:
            logger.warning("Error fetching taxonomy page at offset %d: %s (after retries)", offset, exc)
            break

        if not data:
            break

        all_taxa.extend(data)

        if len(data) < page_size:
            break

        offset += len(data)
        time.sleep(REQUEST_DELAY)

    logger.info("Found %d species-level taxa with genomes.", len(all_taxa))
    return all_taxa


def _fetch_genome_count(code: str) -> int:
    """Fetch the total genome count for a taxon using solr response format."""
    url = GENOME_COUNT_URL_TEMPLATE.format(code=code)
    try:
        data = _fetch_json(url, timeout=120)
        return data.get("response", {}).get("numFound", 0)
    except Exception as exc:
        logger.warning("Error fetching genome count for %s: %s", code, exc)
        return 0


def _fetch_genomes_for_taxon(
    code: str,
    expected_count: int,
) -> dict[str, Any]:
    """
    Paginate through all genome records for a given NCBI taxonomy code.

    Returns aggregated data:
      earliest_inserted   – Optional[datetime]
      latest_modified     – Optional[datetime]
      countries           – set of country names
      hosts               – set of (host_name, host_scientific_name) tuples
      collection_years    – list of collection years
      diseases            – set of disease names
      record_count        – int (actual records fetched)
    """
    insertions: list[dt.datetime] = []
    modifications: list[dt.datetime] = []
    countries: set[str] = set()
    hosts: set[tuple[str, str]] = set()
    collection_years: list[int] = []
    diseases: set[str] = set()
    record_count = 0
    offset = 0

    while True:
        url = GENOME_URL_TEMPLATE.format(
            code=code,
            limit=ROWS_PER_PAGE,
            offset=offset,
        )
        logger.info(
            "Fetching genomes for taxon %s, offset %d (%d/%d)...",
            code,
            offset,
            record_count,
            expected_count,
        )

        try:
            records = _fetch_json(url)
        except Exception as exc:
            logger.warning(
                "Error fetching genomes at offset %d for taxon %s: %s (after retries)",
                offset,
                code,
                exc,
            )
            break

        if not records:
            break

        for rec in records:
            record_count += 1

            # Dates
            ins = _parse_datetime(rec.get("date_inserted"))
            mod = _parse_datetime(rec.get("date_modified"))
            if ins:
                insertions.append(ins)
            if mod:
                modifications.append(mod)

            # Country
            country = (rec.get("isolation_country") or "").strip()
            if country:
                countries.add(country)

            # Host
            host_name = (rec.get("host_name") or "").strip()
            host_sci = (rec.get("host_scientific_name") or "").strip()
            if host_name or host_sci:
                hosts.add((host_name or host_sci, host_sci or host_name))

            # Collection year
            year = rec.get("collection_year")
            if year is not None:
                try:
                    collection_years.append(int(year))
                except (ValueError, TypeError):
                    pass

            # Disease
            raw_disease = rec.get("disease")
            if isinstance(raw_disease, list):
                for d in raw_disease:
                    d = (d or "").strip()
                    if d:
                        diseases.add(d)
            elif isinstance(raw_disease, str):
                d = raw_disease.strip()
                if d:
                    diseases.add(d)

        if len(records) < ROWS_PER_PAGE:
            break

        offset += len(records)
        time.sleep(REQUEST_DELAY)

    return {
        "earliest_inserted": min(insertions) if insertions else None,
        "latest_modified": max(modifications) if modifications else None,
        "countries": countries,
        "hosts": hosts,
        "collection_years": collection_years,
        "diseases": diseases,
        "record_count": record_count,
    }


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
    if version:
        record["version"] = version

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
    SQL_DB = "bv_brc.db"
    EXPIRE = dt.timedelta(days=30)

    def __init__(self):
        super().__init__()
        self._version = ""

    def _get_session(self):
        """Create a requests session with retry logic and connection pooling."""
        session = requests.Session()
        retries = Retry(total=5, backoff_factor=2, status_forcelist=[500, 502, 503, 504])
        adapter = HTTPAdapter(max_retries=retries, pool_connections=10, pool_maxsize=10)
        session.mount("https://", adapter)
        return session

    def load_cache(self):
        """Download all BV-BRC taxonomy + genome data and yield (_id, json_str) tuples for the cache."""
        session = self._get_session()

        # Fetch the API version
        self._version = _fetch_api_version()
        if self._version:
            logger.info("BV-BRC API version: %s", self._version)

        # Step 1 — Get all species-level taxa with genomes from taxonomy API
        taxa = _fetch_taxonomy_species()
        if not taxa:
            logger.error("No taxonomy entries returned from BV-BRC taxonomy API.")
            return

        logger.info("Starting processing of %d species-level taxa", len(taxa))

        count = 0
        errors = 0
        for taxon in taxa:
            code = str(taxon.get("taxon_id", "")).strip()
            taxon_name = (taxon.get("taxon_name") or f"Organism {code}").strip()
            genome_count_from_taxonomy = taxon.get("genomes", 0)

            if not code:
                continue

            try:
                logger.info(
                    "Processing taxon %s (%s, ~%d genomes)...",
                    code, taxon_name, genome_count_from_taxonomy,
                )

                # Get the actual count via solr query
                genome_count = _fetch_genome_count(code)
                if genome_count == 0:
                    genome_count = genome_count_from_taxonomy
                time.sleep(REQUEST_DELAY)

                # Fetch detailed genome metadata
                group_data = _fetch_genomes_for_taxon(code, genome_count)

                # Use the higher of the two counts
                final_count = max(genome_count, group_data["record_count"])
                if final_count == 0:
                    logger.warning("No genomes found for taxon %s, skipping.", code)
                    continue

                # Serialize aggregated data for the cache
                # Convert sets and datetimes to JSON-serializable types
                cache_entry = {
                    "code": code,
                    "taxon_name": taxon_name,
                    "genome_count": final_count,
                    "version": self._version,
                    "earliest_inserted": _iso_date(group_data["earliest_inserted"]),
                    "latest_modified": _iso_date(group_data["latest_modified"]),
                    "countries": sorted(group_data["countries"]),
                    "hosts": sorted([list(h) for h in group_data["hosts"]]),
                    "collection_years": sorted(group_data["collection_years"]),
                    "diseases": sorted(group_data["diseases"]),
                    "record_count": group_data["record_count"],
                }

                cache_id = _make_record_id(code)
                yield (cache_id, json.dumps(cache_entry))
                count += 1
                if count % 500 == 0:
                    logger.info("Cached %s / %s taxa (errors: %s)", count, len(taxa), errors)
            except Exception as e:
                errors += 1
                logger.error("Error processing taxon %s: %s", code, e)
                continue

            time.sleep(REQUEST_DELAY)

        logger.info("Finished loading cache. Total cached: %s, Errors: %s", count, errors)

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
