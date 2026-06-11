#!/usr/bin/env python3
"""
UniProt DataCollection crawler for the NIAID Data Ecosystem.

UniProtKB entries are individual protein records, so this crawler does not
emit raw UniProtKB entries. Instead, it creates taxonomy-based
DataCollections for:

* protein sets, represented by UniProt proteomes;
* proteins, represented by UniProtKB records for the organism.
"""

import copy
import datetime as dt
import logging
import os
import re
import time
import urllib.parse
from collections import defaultdict
from typing import Any, Iterator, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")

API_BASE = "https://rest.uniprot.org"
WEB_BASE = "https://www.uniprot.org"
USER_AGENT = "nde-uniprot-datacollection-crawler/0.1"
REQUEST_TIMEOUT = 90
PAGE_SIZE = int(os.environ.get("UNIPROT_PAGE_SIZE", "500"))
TAXONOMY_CHUNK_SIZE = int(os.environ.get("UNIPROT_TAXONOMY_CHUNK_SIZE", "100"))
DEFAULT_PROTEOME_QUERY = os.environ.get(
    "UNIPROT_PROTEOME_QUERY",
    "proteome_type:REFERENCE",
)

SOURCE_ORGANIZATION = {
    "@type": "Organization",
    "name": "UniProt Consortium",
    "url": "https://www.uniprot.org/",
}

CATALOG = {
    "@type": "DataCatalog",
    "name": "UniProt Knowledgebase",
    "alternateName": "UniProtKB",
    "url": "https://www.uniprot.org/",
}

ABOUT_PROTEIN = {
    "@type": "DefinedTerm",
    "description": (
        "A group of complex organic macromolecules composed of one or more "
        "chains of amino acids linked by peptide bonds."
    ),
    "displayName": "Protein",
    "name": "Protein",
    "url": "http://purl.obolibrary.org/obo/NCIT_C17021",
    "inDefinedTermSet": "NCIT",
}

ABOUT_PROTEOME = {
    "@type": "DefinedTerm",
    "description": (
        "The total protein-encoding capability of an organism's genome, or "
        "the actual expressed protein complement of an organism, tissue, "
        "cell type, or compartment."
    ),
    "displayName": "Proteome",
    "name": "Proteome",
    "url": "http://purl.obolibrary.org/obo/NCIT_C18276",
    "inDefinedTermSet": "NCIT",
}

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "name": "Proteomics",
        "identifier": "topic_0121",
        "url": "http://edamontology.org/topic_0121",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Protein properties",
        "identifier": "topic_0123",
        "url": "http://edamontology.org/topic_0123",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Sequence analysis",
        "identifier": "topic_0080",
        "url": "http://edamontology.org/topic_0080",
        "inDefinedTermSet": "EDAM",
    },
]

MEASUREMENT_TECHNIQUE = [
    {
        "@type": "DefinedTerm",
        "name": "Curation",
        "alternateName": ["Curate", "Biocuration"],
        "url": "http://purl.obolibrary.org/obo/NCIT_C48292",
        "inDefinedTermSet": "NCIT",
    },
    {
        "@type": "DefinedTerm",
        "name": "Sequence analysis",
        "url": "http://edamontology.org/operation_2403",
        "inDefinedTermSet": "EDAM",
    },
]

VARIABLE_MEASURED = [
    {
        "@type": "DefinedTerm",
        "name": "Protein Sequence",
        "url": "http://edamontology.org/data_2976",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Annotation",
        "url": "http://purl.obolibrary.org/obo/NCIT_C44272",
        "inDefinedTermSet": "NCIT",
    },
]

ENCODING_FORMATS = [
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
        "name": "FASTA",
        "url": "http://edamontology.org/format_1929",
        "inDefinedTermSet": "EDAM",
    },
]

UNIPROTKB_SCHEMA_VERSION = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/uniprot.xsd"
)

PROTEIN_KEYWORDS = [
    "UniProt",
    "UniProtKB",
    "Proteins",
    "Swiss-Prot",
    "TrEMBL",
]
PROTEIN_SET_KEYWORDS = ["UniProt", "UniProtKB", "Protein Set", "Proteome"]

CITATION = [
    {
        "@type": "ScholarlyArticle",
        "name": "UniProt: the Universal Protein Knowledgebase in 2025",
        "author": [{"@type": "Organization", "name": "UniProt Consortium"}],
        "datePublished": "2025-01-06",
        "journalName": "Nucleic Acids Research",
        "doi": "10.1093/nar/gkae1010",
        "url": "https://doi.org/10.1093/nar/gkae1010",
    }
]

FUNDING = [
    {
        "@type": "MonetaryGrant",
        "identifier": "HG002273",
        "funder": {
            "@type": "Organization",
            "name": "National Human Genome Research Institute",
        },
    },
    {
        "@type": "MonetaryGrant",
        "identifier": "U24HG007822",
        "funder": {
            "@type": "Organization",
            "name": "Office of the Director, National Institutes of Health",
        },
    },
    {
        "@type": "MonetaryGrant",
        "identifier": "R35GM141873",
        "funder": {
            "@type": "Organization",
            "name": "National Institute of General Medical Sciences",
        },
    },
    {
        "@type": "MonetaryGrant",
        "funder": {
            "@type": "Organization",
            "name": "National Institute of Allergy and Infectious Diseases",
        },
    },
    {
        "@type": "MonetaryGrant",
        "identifier": "BB/T015608/1",
        "funder": {
            "@type": "Organization",
            "name": "Biotechnology and Biological Sciences Research Council",
        },
    },
    {
        "@type": "MonetaryGrant",
        "identifier": "825575",
        "funder": {
            "@type": "Organization",
            "name": (
                "Horizon 2020 - Research and Innovation Framework Programme"
            ),
        },
    },
]

USAGE_INFO = {
    "@type": "CreativeWork",
    "name": "License and disclaimer",
    "url": "https://www.uniprot.org/help/license",
}


def _int_env(name: str) -> Optional[int]:
    value = os.environ.get(name)
    if not value:
        return None
    return int(value)


def _taxon_ids_env() -> Optional[set[int]]:
    value = os.environ.get("UNIPROT_TAXON_IDS")
    if not value:
        return None
    ids = set()
    for item in re.split(r"[\s,]+", value):
        if item:
            ids.add(int(item))
    return ids


def _get_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=2,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "HEAD"),
    )
    adapter = HTTPAdapter(max_retries=retries, pool_connections=10, pool_maxsize=10)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": USER_AGENT})
    return session


def _request(
    session: requests.Session,
    url: str,
    params: Optional[dict[str, Any]] = None,
) -> requests.Response:
    response = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    return response


def _next_link(link_header: Optional[str]) -> Optional[str]:
    if not link_header:
        return None
    for part in link_header.split(","):
        if 'rel="next"' in part:
            match = re.search(r"<([^>]+)>", part)
            if match:
                return match.group(1)
    return None


def _iter_search(
    session: requests.Session,
    endpoint: str,
    params: dict[str, Any],
    limit: Optional[int] = None,
) -> Iterator[dict[str, Any]]:
    url = f"{API_BASE}/{endpoint}/search"
    request_params = {"format": "json", "size": PAGE_SIZE, **params}
    seen = 0

    while url:
        response = _request(session, url, request_params)
        request_params = None
        payload = response.json()
        for result in payload.get("results", []):
            yield result
            seen += 1
            if limit is not None and seen >= limit:
                return
        url = _next_link(response.headers.get("Link"))


def _release_info(session: requests.Session) -> dict[str, Optional[str]]:
    response = _request(
        session,
        f"{API_BASE}/uniprotkb/search",
        {"query": "*", "size": 0, "format": "json"},
    )
    release = response.headers.get("X-UniProt-Release")
    release_date = _parse_release_date(response.headers.get("X-UniProt-Release-Date"))
    return {"release": release, "release_date": release_date}


def _parse_release_date(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    for fmt in ("%d-%B-%Y", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(value, fmt).date().isoformat()
        except ValueError:
            pass
    return None


def _today() -> str:
    return dt.datetime.now(dt.timezone.utc).date().isoformat()


def _ui_url(path: str, query: str) -> str:
    return f"{WEB_BASE}/{path}?{urllib.parse.urlencode({'query': query})}"


def _api_search_url(endpoint: str, query: str) -> str:
    return (
        f"{API_BASE}/{endpoint}/search?"
        f"{urllib.parse.urlencode({'query': query, 'format': 'json'})}"
    )


def _proteome_query(selected_taxa: Optional[set[int]]) -> str:
    if selected_taxa:
        return " OR ".join(f"organism_id:{taxon_id}" for taxon_id in sorted(selected_taxa))
    return DEFAULT_PROTEOME_QUERY


def _make_record_id(*parts: Any) -> str:
    raw = "_".join(str(part) for part in parts if part is not None)
    raw = raw.replace(":", "_")
    raw = re.sub(r"[^A-Za-z0-9_]+", "_", raw)
    return re.sub(r"_+", "_", raw).strip("_").lower()


def _is_empty(value: Any) -> bool:
    return value is None or value == "" or value == [] or value == {}


def _clean(value: Any) -> Any:
    if isinstance(value, dict):
        cleaned = {
            key: _clean(item)
            for key, item in value.items()
            if not _is_empty(item)
        }
        return cleaned
    if isinstance(value, list):
        cleaned = [_clean(item) for item in value if not _is_empty(item)]
        return [item for item in cleaned if not _is_empty(item)]
    return value


def _taxonomy_name(taxonomy: dict[str, Any]) -> str:
    return taxonomy.get("scientificName") or taxonomy.get("commonName") or str(taxonomy.get("taxonId"))


def _classify_taxonomy(taxonomy: dict[str, Any]) -> str:
    scientific_names = [
        item.get("scientificName")
        for item in taxonomy.get("lineage", [])
        if item.get("scientificName")
    ]

    if "Deuterostomia" in scientific_names:
        return "host"

    if "Embryophyta" in scientific_names and not any(
        parasite in scientific_names
        for parasite in ["Arceuthobium", "Cuscuta", "Orobanche", "Striga", "Phoradendron"]
    ):
        return "host"

    if "Arthropoda" in scientific_names:
        if "Acari" in scientific_names:
            if "Ixodida" in scientific_names:
                return "host"
            return "infectiousAgent"
        return "host"

    return "infectiousAgent"


def _organism_term(taxonomy: dict[str, Any]) -> dict[str, Any]:
    taxon_id = taxonomy.get("taxonId")
    taxon_id = str(taxon_id) if taxon_id is not None else None
    name = _taxonomy_name(taxonomy)
    term = {
        "@type": "DefinedTerm",
        "identifier": taxon_id,
        "inDefinedTermSet": "UniProt",
        "url": f"https://www.uniprot.org/taxonomy/{taxon_id}" if taxon_id else None,
        "originalName": name,
        "isCurated": False,
        "name": name,
        "displayName": name,
        "classification": _classify_taxonomy(taxonomy),
    }

    alternative_names = []
    if common_name := taxonomy.get("commonName"):
        term["commonName"] = common_name
        term["displayName"] = f"{common_name} | {name}"
        alternative_names.append(common_name)

    if other_names := taxonomy.get("otherNames"):
        alternative_names.extend(other_names)

    if alternative_names:
        term["alternateName"] = list(dict.fromkeys(alternative_names))

    return _clean(term)


def _organism_schema_field(organism: dict[str, Any]) -> str:
    if organism.get("classification") == "infectiousAgent":
        return "infectiousAgent"
    return "species"


def _catalog(url: str, release: dict[str, Optional[str]]) -> dict[str, Any]:
    catalog = copy.deepcopy(CATALOG)
    catalog["archivedAt"] = url
    catalog["version"] = release.get("release")
    catalog["versionDate"] = release.get("release_date") or _today()
    return _clean(catalog)


def _example_of_work(about: dict[str, Any], api_url: str) -> dict[str, Any]:
    return _clean(
        {
            "@type": "CreativeWork",
            "about": copy.deepcopy(about),
            "encodingFormat": copy.deepcopy(ENCODING_FORMATS),
            "schemaVersion": UNIPROTKB_SCHEMA_VERSION,
            "potentialAction": {
                "@type": "SearchAction",
                "target": api_url,
            },
        }
    )


def _is_based_on(
    *,
    collection_kind: str,
    query_url: str,
    release: dict[str, Optional[str]],
) -> list[dict[str, Any]]:
    return [
        {
            "@type": "Action",
            "name": (
                "DataCollection Generation Process in the NIAID Data "
                "Ecosystem"
            ),
            "description": (
                f"How this UniProt {collection_kind} DataCollection record "
                "was generated for the NIAID Data Ecosystem."
            ),
            "actionProcess": {
                "@type": "HowTo",
                "step": [
                    "1. Query the UniProt REST API for UniProt proteomes and taxonomy metadata.",
                    "2. Group records by organism taxonomy identifier.",
                    "3. Create DataCollection records for protein sets and proteins.",
                    "4. Populate counts, organism metadata, release metadata, and source links.",
                ],
            },
        },
        {
            "@type": "ResourceCatalog",
            "name": "UniProt",
            "url": "https://www.uniprot.org/",
            "version": release.get("release"),
            "versionDate": release.get("release_date"),
        },
        {
            "@type": "CreativeWork",
            "name": "UniProt REST API query",
            "url": query_url,
        },
    ]


def _base_record(
    *,
    record_id: str,
    about: dict[str, Any],
    name: str,
    url: str,
    description: str,
    release: dict[str, Optional[str]],
    collection_size: Optional[int],
    unit_text: str,
    taxonomy: Optional[dict[str, Any]],
    query_url: str,
    collection_kind: str,
    keywords: list[str],
    date_modified: Optional[str] = None,
) -> dict[str, Any]:
    organism = _organism_term(taxonomy) if taxonomy else None
    record = {
        "_id": record_id,
        "@type": "DataCollection",
        "about": copy.deepcopy(about),
        "includedInDataCatalog": [_catalog(url, release)],
        "name": name,
        "url": url,
        "description": description,
        "collectionSize": {
            "minValue": collection_size,
            "unitText": unit_text,
        }
        if collection_size is not None
        else None,
        "author": copy.deepcopy(SOURCE_ORGANIZATION),
        "creator": copy.deepcopy(SOURCE_ORGANIZATION),
        "sourceOrganization": copy.deepcopy(SOURCE_ORGANIZATION),
        "citation": copy.deepcopy(CITATION),
        "conditionsOfAccess": "Open",
        "funding": copy.deepcopy(FUNDING),
        "license": "https://creativecommons.org/licenses/by/4.0/",
        "measurementTechnique": copy.deepcopy(MEASUREMENT_TECHNIQUE),
        "topicCategory": copy.deepcopy(TOPIC_CATEGORY),
        "usageInfo": copy.deepcopy(USAGE_INFO),
        "variableMeasured": copy.deepcopy(VARIABLE_MEASURED),
        "date": date_modified or release.get("release_date") or _today(),
        "dateModified": date_modified or release.get("release_date") or _today(),
        # "version": release.get("release"),
        "isAccessibleForFree": True,
        "keywords": copy.deepcopy(keywords),
        "exampleOfWork": _example_of_work(about, query_url),
        "isBasedOn": _is_based_on(
            collection_kind=collection_kind,
            query_url=query_url,
            release=release,
        ),
    }
    if organism:
        record[_organism_schema_field(organism)] = [organism]
    return _clean(record)


def _build_protein_set_record(
    taxonomy: dict[str, Any],
    proteomes: list[dict[str, Any]],
    release: dict[str, Optional[str]],
) -> dict[str, Any]:
    taxon_id = taxonomy.get("taxonId")
    organism = _taxonomy_name(taxonomy)
    stats = taxonomy.get("statistics") or {}
    proteome_count = stats.get("proteomeCount")
    if proteome_count is None:
        proteome_count = len(proteomes)

    query = f"organism_id:{taxon_id}"
    url = _ui_url("proteomes", query)
    api_url = _api_search_url("proteomes", query)

    representative = proteomes[0] if proteomes else {}

    return _base_record(
        record_id=_make_record_id("uniprot", "proteomes", taxon_id),
        about=ABOUT_PROTEOME,
        name=f"{organism} UniProt proteome collection",
        url=url,
        description=(
            f"UniProt proteomes available for {organism} "
            f"(taxonomy:{taxon_id}). These protein set records group "
            "UniProtKB protein records by organism proteome."
        ),
        release=release,
        collection_size=proteome_count,
        unit_text="UniProt proteomes",
        taxonomy=taxonomy,
        query_url=api_url,
        collection_kind="Protein Set",
        keywords=PROTEIN_SET_KEYWORDS,
        date_modified=representative.get("modified") or release.get("release_date"),
    )


def _build_protein_record(
    taxonomy: dict[str, Any],
    release: dict[str, Optional[str]],
) -> dict[str, Any]:
    taxon_id = taxonomy.get("taxonId")
    organism = _taxonomy_name(taxonomy)
    stats = taxonomy.get("statistics") or {}
    reviewed = stats.get("reviewedProteinCount") or 0
    unreviewed = stats.get("unreviewedProteinCount") or 0
    protein_count = reviewed + unreviewed
    query = f"organism_id:{taxon_id}"
    url = _ui_url("uniprotkb", query)
    api_url = _api_search_url("uniprotkb", query)

    return _base_record(
        record_id=_make_record_id("uniprot", "proteins", taxon_id),
        about=ABOUT_PROTEIN,
        name=f"{organism} UniProtKB protein collection",
        url=url,
        description=(
            f"UniProtKB protein records for {organism} "
            f"(taxonomy:{taxon_id}), including reviewed Swiss-Prot records "
            "and unreviewed TrEMBL records."
        ),
        release=release,
        collection_size=protein_count,
        unit_text="UniProtKB protein records",
        taxonomy=taxonomy,
        query_url=api_url,
        collection_kind="Proteins",
        keywords=PROTEIN_KEYWORDS,
    )


def _taxonomy_query_for_ids(ids: list[int]) -> str:
    return " OR ".join(f"tax_id:{taxon_id}" for taxon_id in ids)


def _fetch_taxonomy(
    session: requests.Session,
    taxon_ids: set[int],
) -> dict[int, dict[str, Any]]:
    taxonomy: dict[int, dict[str, Any]] = {}
    sorted_ids = sorted(taxon_ids)

    for start in range(0, len(sorted_ids), TAXONOMY_CHUNK_SIZE):
        chunk = sorted_ids[start : start + TAXONOMY_CHUNK_SIZE]
        query = _taxonomy_query_for_ids(chunk)
        logger.info("Fetching taxonomy metadata for %s taxa", len(chunk))
        for item in _iter_search(
            session,
            "taxonomy",
            {"query": query, "size": len(chunk)},
        ):
            taxon_id = item.get("taxonId")
            if taxon_id is not None:
                taxonomy[taxon_id] = item
        time.sleep(0.1)

    missing = sorted(set(sorted_ids) - set(taxonomy))
    for taxon_id in missing:
        try:
            response = _request(session, f"{API_BASE}/taxonomy/{taxon_id}")
            item = response.json()
            taxonomy[taxon_id] = item
        except Exception as exc:
            logger.warning("Could not fetch taxonomy %s: %s", taxon_id, exc)

    return taxonomy


def _collect_proteomes(
    session: requests.Session,
    selected_taxa: Optional[set[int]],
) -> list[dict[str, Any]]:
    query = _proteome_query(selected_taxa)
    limit = _int_env("UNIPROT_MAX_PROTEOMES")
    proteomes = list(
        _iter_search(
            session,
            "proteomes",
            {"query": query},
            limit=limit,
        )
    )
    logger.info("Collected %s UniProt proteomes", len(proteomes))
    if not proteomes:
        logger.warning("UniProt proteome query returned no records: %s", query)
    return proteomes


def parse() -> Iterator[dict[str, Any]]:
    session = _get_session()
    release = _release_info(session)
    selected_taxa = _taxon_ids_env()
    proteomes = _collect_proteomes(session, selected_taxa)

    taxon_ids = set(selected_taxa or set())
    proteomes_by_taxon: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for proteome in proteomes:
        taxon_id = ((proteome.get("taxonomy") or {}).get("taxonId"))
        if taxon_id is not None:
            taxon_ids.add(taxon_id)
            proteomes_by_taxon[taxon_id].append(proteome)

    taxonomy = _fetch_taxonomy(session, taxon_ids)

    for taxon_id in sorted(proteomes_by_taxon):
        taxon = taxonomy.get(taxon_id)
        if not taxon:
            continue
        yield _build_protein_set_record(
            taxon,
            proteomes_by_taxon[taxon_id],
            release,
        )

    for taxon_id in sorted(taxon_ids):
        taxon = taxonomy.get(taxon_id)
        if not taxon:
            continue
        yield _build_protein_record(taxon, release)
