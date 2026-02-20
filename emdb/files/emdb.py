#!/usr/bin/env python3
"""
EMDB DataCollection crawler for the NIAID Data Ecosystem.

Groups Electron Microscopy Data Bank records by NCBI Taxonomy ID and
generates one DataCollection record per organism.

API reference
-------------
Facet (taxonomy code frequency table):
    https://www.ebi.ac.uk/emdb/api/facet/
    natural_source_ncbi_code%3A%5B*%20TO%20*%5D?field=natural_source_ncbi_code

Records for a given taxonomy code (paginated):
    https://www.ebi.ac.uk/emdb/api/search/
    natural_source_ncbi_code%3A{CODE}?rows=1000&page=1
"""

import copy
import datetime as dt
import json
import logging
import time
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any, Iterator, Optional

# NOTE: species vs infectiousAgent classification is handled downstream
# by the pubtator standardizer in the uploader. All organisms are emitted
# under "species" here.

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------

FACET_URL = (
    "https://www.ebi.ac.uk/emdb/api/facet/"
    "natural_source_ncbi_code%3A%5B*%20TO%20*%5D"
    "?field=natural_source_ncbi_code"
)

SEARCH_URL_TEMPLATE = (
    "https://www.ebi.ac.uk/emdb/api/search/"
    "natural_source_ncbi_code%3A{code}?rows={rows}&page={page}"
)

NCBI_EFETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
)

# End-user search URL (used as the DataCollection url)
EMDB_SEARCH_URL_TEMPLATE = (
    "https://www.ebi.ac.uk/emdb/emsearch/"
    "?q=natural_source_ncbi_code:%22{code}%22"
)

# ---------------------------------------------------------------------------
# Curated / static metadata
# ---------------------------------------------------------------------------

SOURCE_METADATA = {
    "name": "EMDB (as DataCollection records)",
    "url": "https://www.ebi.ac.uk/emdb/",
    "description": (
        "The Electron Microscopy Data Bank is a repository of 3D volume "
        "maps and other processed images from Electron microscopy and "
        "crystallography experiments."
    ),
    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
}

ABOUT_DEFINED_TERM = {
    "@type": "DefinedTerm",
    "description": (
        "For schema, consider ImageObject https://schema.org/ImageObject"
    ),
    "displayName": "Image",
    "name": "Image",
    "url": "http://purl.obolibrary.org/obo/NCIT_C48179",
}

_DDE_CURATED_BY = {
    "name": "Data Discovery Engine",
    "url": "https://discovery.biothings.io/",
    "dateModified": "2026-01-03",
}

MEASUREMENT_TECHNIQUE = [
    {
        "@type": "DefinedTerm",
        "alternateName": ["Biological imaging"],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "EDAM",
        "isCurated": True,
        "name": "Bioimaging",
        "url": "http://edamontology.org/topic_3383",
    },
    {
        "@type": "DefinedTerm",
        "alternateName": ["Electron Microscopy"],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Electron Microscopy",
        "url": "http://purl.obolibrary.org/obo/NCIT_C16854",
    },
    {
        "@type": "DefinedTerm",
        "alternateName": ["Atomic Resolution X-Ray Crystallography"],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Atomic Resolution X-Ray Crystallography",
        "url": "http://purl.obolibrary.org/obo/NCIT_C19588",
    },
    {
        "@type": "DefinedTerm",
        "alternateName": ["Tomography"],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Tomography",
        "url": "http://purl.obolibrary.org/obo/NCIT_C38093",
    },
]

VARIABLE_MEASURED = [
    {
        "@type": "DefinedTerm",
        "alternateName": ["Conformal", "Conformation"],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Conformation",
        "url": "http://purl.obolibrary.org/obo/NCIT_C13738",
    },
    {
        "@type": "DefinedTerm",
        "alternateName": [
            "Mass",
            "Molar Mass",
            "Molecular Mass",
            "Molecular Weight",
            "molecular mass",
            "mw",
        ],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Molecular Mass",
        "url": "http://purl.obolibrary.org/obo/NCIT_C28272",
    },
    {
        "@type": "DefinedTerm",
        "alternateName": ["Resolution", "Resolution Property"],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Resolution",
        "url": "http://purl.obolibrary.org/obo/NCIT_C47876",
    },
]

TOPIC_CATEGORY = [
    "http://edamontology.org/topic_0611",  # Electron microscopy
    "http://edamontology.org/topic_1317",  # Structural biology
    "http://edamontology.org/topic_0781",  # Virology
]

EXAMPLE_OF_WORK_ABOUT = [
    {
        "@type": "DefinedTerm",
        "inDefinedTermSet": "CHMO",
        "name": "electron micrograph",
        "displayName": "Electron Micrograph",
        "url": "http://purl.obolibrary.org/obo/CHMO_0001800",
        "termCode": "CHMO:0001800",
    },
    {
        "@type": "DefinedTerm",
        "inDefinedTermSet": "NCIT",
        "name": "Tomogram",
        "displayName": "Tomogram",
        "url": "http://purl.obolibrary.org/obo/NCIT_C25221",
        "termCode": "NCIT:C25221",
    },
    {
        "@type": "DefinedTerm",
        "inDefinedTermSet": "EDAM",
        "name": "3D EM Map",
        "displayName": "3D EM Map",
        "url": "http://edamontology.org/data_3805",
        "termCode": "EDAM:3805",
    },
]

ENCODING_FORMAT = [
    {
        "@type": "DefinedTerm",
        "name": "PNG",
        "url": "http://edamontology.org/format_3603",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "MAP",
        "url": "http://edamontology.org/format_2060",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "CIF",
        "url": "http://edamontology.org/format_4024",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "XML",
        "url": "http://edamontology.org/format_2332",
        "inDefinedTermSet": "EDAM",
    },
]

HOWTO_STEPS = [
    (
        "Step 1. Use API call to get a frequency table of all records by "
        "NCBI taxonomy ID code: https://www.ebi.ac.uk/emdb/api/facet/"
        "natural_source_ncbi_code%3A%5B*%20TO%20*%5D"
        "?field=natural_source_ncbi_code"
    ),
    (
        "Step 2: For each code, pull all the EMDB records associated with "
        "that NCBI Taxonomy code (paginating as needed): "
        "https://www.ebi.ac.uk/emdb/api/search/"
        "natural_source_ncbi_code%3A{NCBI_CODE}?rows=1000&page=1"
    ),
    (
        "Step 3: Parse the records to generate an organism-specific "
        "DataCollection. Use the latest 'key_dates.update' value for "
        "'dateModified' and the earliest 'key_dates.deposition' value for "
        "'dateCreated'. Generate the 'temporalCoverage' from the "
        "'key_dates' info across all records for that organism. Use the "
        "count of the records for 'collectionSize'. Generate the values "
        "for the 'species', 'infectiousAgent', 'url', 'name', and "
        "'description' properties based on the information retrieved for "
        "the NCBI taxonomy code and any templated text."
    ),
    (
        "Step 4. Fill in manually curated fields from the EMDB resource "
        "catalog record such as 'measurementTechnique', "
        "'variableMeasured', 'topicCategory', 'conditionsOfAccess', "
        "'usageInfo'."
    ),
]

ROWS_PER_PAGE = 1000
REQUEST_DELAY = 0.34  # seconds between API calls (≈3 req/s)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _fetch_json(url: str) -> Any:
    """Fetch a URL and return parsed JSON."""
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "nde-emdb-crawler/0.1"},
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read())


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
    return f"emdb_{code}"


# ---------------------------------------------------------------------------
# API fetching
# ---------------------------------------------------------------------------


def _fetch_facet_codes() -> dict[str, int]:
    """Return {ncbi_taxonomy_code: record_count} from the EMDB facet API."""
    logger.info("Fetching EMDB facet data...")
    data = _fetch_json(FACET_URL)
    return data.get("natural_source_ncbi_code", {})


def _fetch_taxonomy_info(
    codes: list[str],
    batch_size: int = 50,
) -> dict[str, dict[str, str]]:
    """
    Query NCBI taxonomy for a list of codes.
    Returns {code: {"name": ..., "division": ...}}.
    """
    result: dict[str, dict[str, str]] = {}

    for i in range(0, len(codes), batch_size):
        batch = codes[i : i + batch_size]
        ids = ",".join(batch)
        url = f"{NCBI_EFETCH_URL}?db=taxonomy&id={ids}&rettype=xml"
        logger.info(
            "Fetching taxonomy info for %d codes (batch %d)...",
            len(batch),
            i // batch_size + 1,
        )

        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "nde-emdb-crawler/0.1"},
            )
            with urllib.request.urlopen(req, timeout=120) as resp:
                tree = ET.parse(resp)

            for taxon in tree.iter("Taxon"):
                tid = taxon.find("TaxId")
                name = taxon.find("ScientificName")
                div = taxon.find("Division")
                lineage = taxon.find("Lineage")
                # Only process top-level entries (those with a Lineage)
                if tid is not None and name is not None and lineage is not None:
                    result[tid.text] = {
                        "name": name.text or "",
                        "division": div.text if div is not None else "",
                    }
        except Exception as exc:
            logger.warning("Error fetching taxonomy batch %d: %s", i, exc)

        if i + batch_size < len(codes):
            time.sleep(REQUEST_DELAY)

    return result


def _fetch_records_for_code(
    code: str,
    expected_count: int,
) -> dict[str, Any]:
    """
    Paginate through all EMDB records for a given NCBI taxonomy code.

    Returns aggregated data:
      organism_name           – str (from the first EMDB record)
      earliest_deposition     – Optional[datetime]
      latest_update           – Optional[datetime]
      funding                 – list of unique grant dicts
      record_count            – int (actual records fetched)
    """
    organism_name = ""
    depositions: list[dt.datetime] = []
    updates: list[dt.datetime] = []
    funding_set: set[tuple[str, str]] = set()  # (code, funder)
    record_count = 0
    page = 1

    while True:
        url = SEARCH_URL_TEMPLATE.format(
            code=code,
            rows=ROWS_PER_PAGE,
            page=page,
        )
        logger.info(
            "Fetching EMDB records for code %s, page %d (%d/%d)...",
            code,
            page,
            record_count,
            expected_count,
        )

        try:
            records = _fetch_json(url)
        except Exception as exc:
            logger.warning(
                "Error fetching page %d for code %s: %s", page, code, exc
            )
            break

        if not records:
            break

        for rec in records:
            record_count += 1
            admin = rec.get("admin", {})

            # Organism name from the first usable record
            if not organism_name:
                for supra in (
                    rec.get("sample", {})
                    .get("supramolecule_list", {})
                    .get("supramolecule", [])
                ):
                    for ns in supra.get("natural_source", []):
                        org = ns.get("organism", {})
                        if str(org.get("ncbi", "")) == str(code):
                            organism_name = org.get("valueOf_", "")
                            break
                    if organism_name:
                        break

            # Dates
            key_dates = admin.get("key_dates", {})
            dep = _parse_datetime(key_dates.get("deposition"))
            upd = _parse_datetime(key_dates.get("update"))
            if dep:
                depositions.append(dep)
            if upd:
                updates.append(upd)

            # Funding
            for grant in (
                admin.get("grant_support", {}).get("grant_reference", [])
            ):
                grant_code = (grant.get("code") or "").strip()
                funder = (grant.get("funding_body") or "").strip()
                if (
                    grant_code
                    and funder
                    and funder.lower() != "not funded"
                ):
                    funding_set.add((grant_code, funder))

        if len(records) < ROWS_PER_PAGE:
            break

        page += 1
        time.sleep(REQUEST_DELAY)

    # Build de-duplicated funding list
    funding = []
    for g_code, g_funder in sorted(funding_set):
        funding.append(
            {
                "@type": "MonetaryGrant",
                "identifier": g_code,
                "funder": {
                    "@type": "Organization",
                    "name": g_funder,
                },
            }
        )

    return {
        "organism_name": organism_name,
        "earliest_deposition": min(depositions) if depositions else None,
        "latest_update": max(updates) if updates else None,
        "funding": funding,
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


def _build_example_of_work(code: str) -> dict[str, Any]:
    return {
        "@type": "CreativeWork",
        "about": copy.deepcopy(EXAMPLE_OF_WORK_ABOUT),
        "encodingFormat": copy.deepcopy(ENCODING_FORMAT),
        "schemaVersion": "https://www.emdataresource.org/mapformat.html",
        "potentialAction": {
            "@type": "Action",
            "name": "Use API call for an example record",
            "target": (
                "https://www.ebi.ac.uk/emdb/api/search/"
                f"natural_source_ncbi_code%3A{code}?rows=1&page=1"
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
            f"How this EMDB {about_name} DataCollection Record was "
            "generated for the NIAID Data Ecosystem."
        ),
        "actionProcess": {
            "@type": "HowTo",
            "step": copy.deepcopy(HOWTO_STEPS),
        },
    }

    source_obj = {
        "@type": "nde:ResourceCatalog",
        "name": "Electron Microscopy Data Bank",
        "url": "https://data.niaid.nih.gov/resources?id=dde_f6161fcb840b7162",
    }

    return [action_obj, source_obj]


def _build_data_collection(
    code: str,
    facet_count: int,
    taxonomy_info: dict[str, str],
    group_data: dict[str, Any],
) -> dict[str, Any]:
    organism_name = (
        group_data["organism_name"]
        or taxonomy_info.get("name", f"Organism {code}")
    )

    url = EMDB_SEARCH_URL_TEMPLATE.format(code=code)
    identifier = f"taxonomy:{code}"

    earliest = group_data["earliest_deposition"]
    latest = group_data["latest_update"]

    name = f"{organism_name} electron microscopy data records at EMDB"
    description = (
        "3D cryo-EM volume maps and Tomograms available for proteins, "
        "complexes, ligands, nucleic acids and other (macro)molecular "
        f"structures sourced from {organism_name} ({identifier}). "
        f"For more details, visit: {url}"
    )

    record: dict[str, Any] = {
        "_id": _make_record_id(code),
        "@type": "DataCollection",
        "about": copy.deepcopy(ABOUT_DEFINED_TERM),
        "includedInDataCatalog": [
            {
                "@type": "DataCatalog",
                "name": "Electron Microscopy Data Bank",
                "url": "https://www.ebi.ac.uk/emdb/",
                "versionDate": dt.datetime.now(dt.timezone.utc)
                .date()
                .isoformat(),
            }
        ],
        "name": name,
        "url": url,
        "description": description,
        "collectionSize": facet_count,
        "date": _iso_date(latest) or _iso_date(earliest),
        "dateModified": _iso_date(latest),
        "dateCreated": _iso_date(earliest),
        "conditionsOfAccess": "Open",
        "license": "https://creativecommons.org/publicdomain/zero/1.0/",
        "measurementTechnique": copy.deepcopy(MEASUREMENT_TECHNIQUE),
        "variableMeasured": copy.deepcopy(VARIABLE_MEASURED),
        "topicCategory": copy.deepcopy(TOPIC_CATEGORY),
        "usageInfo": (
            "https://www.ebi.ac.uk/emdb/faq#:~:text="
            "Data%20files%20contained%20in%20the,"
            "and%20the%20EMDB%20accession%20id."
        ),
        "creditText": (
            "To cite a record from EMDB, please visit "
            "https://www.ebi.ac.uk/emdb/about."
        ),
        "version": "3.0.9.3",
        "exampleOfWork": _build_example_of_work(code),
        "isBasedOn": _build_is_based_on(ABOUT_DEFINED_TERM["name"]),
    }

    # All organisms stored as species; the pubtator standardizer in the
    # uploader handles reclassification to infectiousAgent when appropriate.
    record["species"] = _build_organism_term(code, organism_name)

    # Temporal coverage across all records in the group
    if earliest or latest:
        record["temporalCoverage"] = {
            "@type": "TemporalInterval",
            "startDate": _iso_date(earliest),
            "endDate": _iso_date(latest),
        }

    # Aggregated / de-duplicated funding from underlying records
    if group_data["funding"]:
        record["funding"] = group_data["funding"]

    return record


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def parse() -> Iterator[dict[str, Any]]:
    """Yield DataCollection records, one per NCBI taxonomy code."""

    # Step 1 – frequency table of taxonomy codes
    facet = _fetch_facet_codes()
    if not facet:
        logger.error("No taxonomy codes returned from EMDB facet API.")
        return

    logger.info("Found %d taxonomy codes.", len(facet))

    # Step 2 – resolve names and divisions from NCBI Taxonomy (batched)
    codes = list(facet.keys())
    taxonomy_info = _fetch_taxonomy_info(codes)

    # Step 3 – for each code, paginate through records and build output
    for code in sorted(facet, key=lambda c: -facet[c]):
        count = facet[code]
        logger.info(
            "Processing taxonomy code %s (~%d records)...",
            code,
            count,
        )

        tax_info = taxonomy_info.get(
            code,
            {"name": f"Organism {code}", "division": ""},
        )

        group_data = _fetch_records_for_code(code, count)

        yield _build_data_collection(code, count, tax_info, group_data)

        time.sleep(REQUEST_DELAY)


if __name__ == "__main__":
    with open("emdb_sample_output.json", "w") as f:
        for record in parse():
            json.dump(record, f)
            f.write("\n")
