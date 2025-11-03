#!/usr/bin/env python3
"""
Fetch ARK Portal dataset metadata using Synapse client and map it to NDE fields.

This cleaner implementation uses the Synapse Python client to query the ARK datasets
table directly, avoiding the need for multiple REST API calls.
"""
from __future__ import annotations

import logging
import os
from typing import Any, Dict, Generator, List, Optional

import pandas as pd
import synapseclient
from synapseclient import Synapse

# Synapse table containing ARK datasets
ARK_DATASETS_TABLE = "syn68554562"
ARK_PORTAL_DATASET_URL = "https://arkportal.synapse.org/Explore/Datasets/DetailsPage?id={syn_id}"
SYNAPSE_OBJECT_URL = "https://www.synapse.org/Synapse:{syn_id}"

# Path to local file containing Synapse auth token
TOKEN_FILE = ".synapse_token"

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

CONDITIONS_OF_ACCESS = {
    "released": "Restricted",
    "under peer review": "Embargoed",
    "unreleased": "Embargoed",
    "test": "Test",
    "deprecated": "Deprecated",
}

DIAGNOSIS_TERMS = {
    "control": {
        "identifier": "http://purl.obolibrary.org/obo/NCIT_C61299",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/ncit.owl",
    },
    "sle": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0007915",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "ra": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0008383",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "vitiligo": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0008661",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "dermatomyositis": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0016367",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "pso": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0005083",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "psa": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0011849",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "scleroderma": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0019340",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "sjd": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0010030",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "ln": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0005556",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "cle": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0005282",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
    "oa": {
        "identifier": "http://purl.obolibrary.org/obo/MONDO_0005178",
        "inDefinedTermSet": "http://purl.obolibrary.org/obo/mondo.owl",
    },
}

IDENTIFIER_BASE_URLS = (
    (
        "dbGapAccession",
        "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={value}",
    ),
    ("ImmPortAccession", "https://www.immport.org/resources/study/{value}"),
)


def clean_value(value: Any) -> Optional[Any]:
    """Clean pandas values, converting NaN to None."""
    if pd.isna(value):
        return None
    return value


def ensure_list(value: Any) -> List[Any]:
    """Convert value to list if it isn't already."""
    if value is None:
        return []
    if isinstance(value, list):
        return [v for v in value if v is not None and not (isinstance(v, float) and pd.isna(v))]
    # For scalar values, check if it's NaN
    try:
        if pd.isna(value):
            return []
    except (TypeError, ValueError):
        # pd.isna might fail on some types, just proceed
        pass
    return [value]


def read_token_from_file(token_file: str = TOKEN_FILE) -> Optional[str]:
    """Read auth token from a local file."""
    # Try to read from the script's directory first
    script_dir = os.path.dirname(os.path.abspath(__file__))
    token_path = os.path.join(script_dir, token_file)

    if os.path.exists(token_path):
        try:
            with open(token_path, 'r') as f:
                token = f.read().strip()
                if token:
                    return token
        except Exception as error:
            logger.warning("Failed to read token from %s: %s", token_path, error)

    # Fallback to current working directory
    if os.path.exists(token_file):
        try:
            with open(token_file, 'r') as f:
                token = f.read().strip()
                if token:
                    return token
        except Exception as error:
            logger.warning("Failed to read token from %s: %s", token_file, error)

    return None


def get_synapse_client(auth_token: Optional[str] = None) -> Synapse:
    """Initialize and login to Synapse client."""
    syn = synapseclient.Synapse()
    token = auth_token or read_token_from_file()
    if token:
        syn.login(authToken=token)
    else:
        syn.login()
    return syn


def fetch_publication_details(syn: Synapse, pub_syn_id: str) -> Optional[Dict]:
    """Fetch publication metadata from Synapse."""
    try:
        annotations = syn.get_annotations(pub_syn_id)
        return {
            "synapseId": pub_syn_id,
            "pmid": clean_value(annotations.get("PMID", [None])[0] if annotations.get("PMID") else None),
            "doi": clean_value(annotations.get("DOI", [None])[0] if annotations.get("DOI") else None),
            "url": clean_value(annotations.get("URL", [None])[0] if annotations.get("URL") else None),
            "title": clean_value(annotations.get("title", [None])[0] if annotations.get("title") else None),
        }
    except Exception as error:
        logger.warning("Failed to load publication %s: %s", pub_syn_id, error)
        return None


def fetch_wiki_content(syn: Synapse, entity_id: str, wiki_id: Optional[str] = None) -> Optional[str]:
    """Fetch wiki markdown content from Synapse."""
    try:
        wiki = syn.getWiki(entity_id, wiki_id) if wiki_id else syn.getWiki(entity_id)
        return wiki.get("markdown")
    except synapseclient.core.exceptions.SynapseHTTPError as error:
        if error.response.status_code == 404:
            return None
        logger.warning("Failed to load wiki %s/%s: %s", entity_id, wiki_id, error)
        return None
    except Exception as error:
        logger.warning("Failed to load wiki %s/%s: %s", entity_id, wiki_id, error)
        return None


def process_dataset_row(syn: Synapse, row: pd.Series) -> Dict:
    """Convert a Synapse table row to NDE format."""
    syn_id = row["id"]
    logger.info("Processing ARK dataset %s", syn_id)

    try:
        # Get entity details for additional metadata
        entity = syn.get(syn_id, downloadFile=False)
    except Exception as error:
        logger.error("Failed to get entity %s: %s", syn_id, error)
        raise

    doc: Dict = {
        "_id": syn_id,
        "url": ARK_PORTAL_DATASET_URL.format(syn_id=syn_id),
        "name": clean_value(row["name"]),
        "description": clean_value(row["description"]),
        "dateCreated": clean_value(row["createdOn"]),
        "dateModified": clean_value(row["modifiedOn"]),
        # Assigned fields from mapping
        "license": "https://arkportal.synapse.org/Data%20Access",
        "usageInfo": "https://help.arkportal.org/help/data-use-certificate#DataUse&Acknowledgement-Acknowledgement",
        "genre": "IID",
    }

    # Add program as author (Organization type)
    program = clean_value(row.get("program"))
    if program is not None:
        doc["author"] = {
            "@type": "Organization",
            "name": program,
        }

    # Conditions of access (mapped from datasetStatus)
    status = clean_value(row.get("datasetStatus"))
    if status is not None:
        mapped = CONDITIONS_OF_ACCESS.get(status.lower())
        if mapped:
            doc["conditionsOfAccess"] = mapped

    # Keywords
    keywords: List[str] = []
    dataset_type = clean_value(row.get("datasetType"))
    if dataset_type is not None:
        keywords.append(dataset_type)
    data_subtype = clean_value(row.get("dataSubtype"))
    if data_subtype is not None:
        keywords.append(data_subtype)

    # Measurement techniques (assay + dataType)
    measurement_values = ensure_list(row.get("assay")) + ensure_list(row.get("dataType"))
    measurement = list(dict.fromkeys(measurement_values))  # Remove duplicates while preserving order
    if len(measurement) > 0:
        doc["measurementTechnique"] = [{"name": value} for value in measurement]

    # Sample/biospecimen information
    biospecimen_type = ensure_list(row.get("biospecimenType"))
    biospecimen_subtype = ensure_list(row.get("biospecimenSubtype"))
    if len(biospecimen_type) > 0 or len(biospecimen_subtype) > 0:
        anatomical_structure: Dict[str, Any] = {}
        if len(biospecimen_type) > 0:
            anatomical_structure["name"] = biospecimen_type
        if len(biospecimen_subtype) > 0:
            anatomical_structure["sampleType"] = biospecimen_subtype
        doc["sample"] = {"anatomicalStructure": anatomical_structure}

    # Health conditions from diagnosis
    diagnosis_values = ensure_list(row.get("diagnosis"))
    health_condition: List[Dict] = []
    for value in diagnosis_values:
        if not value:
            continue
        normalized = str(value).strip()
        # "At-Risk RA" is not an official diagnosis, add to keywords
        if normalized.lower() == "at-risk ra":
            keywords.append(normalized)
            continue
        lookup = DIAGNOSIS_TERMS.get(normalized.lower())
        if lookup:
            entry = {"name": normalized}
            entry.update(lookup)
            health_condition.append(entry)
        else:
            # Unknown diagnosis terms go to keywords
            keywords.append(normalized)
    if len(health_condition) > 0:
        doc["healthCondition"] = health_condition

    # DOI
    doi = clean_value(row.get("doi"))
    if doi is not None:
        doc["doi"] = doi

    # Credit/acknowledgment statement
    acknowledgment_ref = clean_value(row.get("acknowledgmentStatement"))
    if acknowledgment_ref is not None and hasattr(entity, "parentId"):
        credit_text = fetch_wiki_content(syn, entity.parentId, acknowledgment_ref.split("/")[-1])
        if credit_text:
            doc["creditText"] = credit_text
        else:
            doc["creditText"] = SYNAPSE_OBJECT_URL.format(
                syn_id=f"{entity.parentId}/wiki/{acknowledgment_ref.split('/')[-1]}"
            )

    # Publications (citation field)
    pub_syn_ids = ensure_list(row.get("publicationSynID"))
    publications: List[Dict] = []
    for pub_syn_id in pub_syn_ids:
        if not pub_syn_id:
            continue
        pub_details = fetch_publication_details(syn, pub_syn_id)
        if pub_details:
            publications.append(pub_details)
    if len(publications) > 0:
        doc["citation"] = publications

    # Associated code repositories (isRelatedTo with ComputationalTool type)
    code_urls = ensure_list(row.get("associatedCodeURL"))
    if len(code_urls) > 0:
        doc["isRelatedTo"] = [
            {"@type": "ComputationalTool", "codeRepository": url}
            for url in code_urls if url
        ]

    # Identifiers (dbGap, ImmPort)
    identifiers: List[Dict] = []
    for key, template in IDENTIFIER_BASE_URLS:
        value = clean_value(row.get(key))
        if value is not None:
            identifiers.append(
                {
                    "type": key,
                    "value": value,
                    "url": template.format(value=value),
                }
            )
    if len(identifiers) > 0:
        doc["identifier"] = identifiers

    # Add keywords
    if len(keywords) > 0:
        doc["keywords"] = sorted(set(keywords))

    return doc


def parse(
    auth_token: Optional[str] = None,
) -> Generator[Dict, None, None]:
    """
    Fetch and parse ARK datasets from Synapse.

    Args:
        auth_token: Synapse authentication token. If not provided, reads from TOKEN_FILE.

    Yields:
        Parsed dataset documents in NDE format.
    """
    syn = get_synapse_client(auth_token)

    # Query all datasets that are not depricated or test
    query = f"SELECT * FROM {ARK_DATASETS_TABLE} WHERE datasetStatus NOT IN ('deprecated', 'test')"

    logger.info("Querying ARK datasets: %s", query)

    try:
        result = syn.tableQuery(query)
        df = result.asDataFrame()
    except Exception as error:
        logger.error("Failed to query ARK datasets: %s", error)
        return

    logger.info("Found %d datasets", len(df))

    processed = 0
    for _, row in df.iterrows():
        try:
            doc = process_dataset_row(syn, row)
            processed += 1
            yield doc
        except Exception as error:
            logger.error("Failed to process dataset %s: %s", row.get("id", "unknown"), error, exc_info=True)
            continue

    logger.info("Finished parsing %d ARK datasets", processed)


if __name__ == "__main__":
    # Example usage
    for dataset in parse():
        print(f"Processed: {dataset['_id']} - {dataset['name']}")
