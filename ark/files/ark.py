#!/usr/bin/env python3
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
        pub: Dict[str, Any] = {"synapseId": pub_syn_id}

        pmid = annotations.get("PMID", [None])[0] if annotations.get("PMID") else None
        if pmid and not pd.isna(pmid):
            pub["pmid"] = pmid

        doi = annotations.get("DOI", [None])[0] if annotations.get("DOI") else None
        if doi and not pd.isna(doi):
            pub["doi"] = doi

        url = annotations.get("URL", [None])[0] if annotations.get("URL") else None
        if url and not pd.isna(url):
            pub["url"] = url

        title = annotations.get("title", [None])[0] if annotations.get("title") else None
        if title and not pd.isna(title):
            pub["title"] = title

        return pub
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
        # Assigned fields from mapping
        "license": "https://arkportal.synapse.org/Data%20Access",
        "usageInfo": "https://help.arkportal.org/help/data-use-certificate#DataUse&Acknowledgement-Acknowledgement",
        "genre": "IID",
    }

    # Add name if present
    name = row.get("name")
    if name and not pd.isna(name) and str(name).strip():
        doc["name"] = str(name).strip()

    # Add description if present
    description = row.get("description")
    if description and not pd.isna(description) and str(description).strip():
        doc["description"] = str(description).strip()

    # Add dates if present
    date_created = row.get("createdOn")
    if date_created and not pd.isna(date_created):
        doc["dateCreated"] = date_created

    date_modified = row.get("modifiedOn")
    if date_modified and not pd.isna(date_modified):
        doc["dateModified"] = date_modified

    # Add program as author (Organization type)
    program = row.get("program")
    if program and not pd.isna(program) and str(program).strip():
        doc["author"] = {
            "@type": "Organization",
            "name": str(program).strip(),
        }

    # Conditions of access (mapped from datasetStatus)
    status = row.get("datasetStatus")
    if status and not pd.isna(status):
        mapped = CONDITIONS_OF_ACCESS.get(str(status).lower())
        if mapped:
            doc["conditionsOfAccess"] = mapped

    # Keywords
    keywords: List[str] = []
    dataset_type = row.get("datasetType")
    if dataset_type and not pd.isna(dataset_type) and str(dataset_type).strip():
        keywords.append(str(dataset_type).strip())

    data_subtype = row.get("dataSubtype")
    if data_subtype and not pd.isna(data_subtype) and str(data_subtype).strip():
        keywords.append(str(data_subtype).strip())

    # Measurement techniques (assay + dataType)
    measurement_values: List[str] = []

    assay = row.get("assay")
    if assay is not None and not pd.isna(assay):
        if isinstance(assay, (list, tuple)):
            for item in assay:
                if item and not pd.isna(item) and str(item).strip():
                    measurement_values.append(str(item).strip())
        elif str(assay).strip():
            measurement_values.append(str(assay).strip())

    data_type = row.get("dataType")
    if data_type is not None and not pd.isna(data_type):
        if isinstance(data_type, (list, tuple)):
            for item in data_type:
                if item and not pd.isna(item) and str(item).strip():
                    measurement_values.append(str(item).strip())
        elif str(data_type).strip():
            measurement_values.append(str(data_type).strip())

    # Remove duplicates while preserving order
    measurement = list(dict.fromkeys(measurement_values))
    if len(measurement) > 0:
        doc["measurementTechnique"] = [{"name": value} for value in measurement]

    # Sample/biospecimen information
    biospecimen_type_values: List[str] = []
    biospecimen_type = row.get("biospecimenType")
    if biospecimen_type is not None and not pd.isna(biospecimen_type):
        if isinstance(biospecimen_type, (list, tuple)):
            for item in biospecimen_type:
                if item and not pd.isna(item) and str(item).strip():
                    biospecimen_type_values.append(str(item).strip())
        elif str(biospecimen_type).strip():
            biospecimen_type_values.append(str(biospecimen_type).strip())

    biospecimen_subtype_values: List[str] = []
    biospecimen_subtype = row.get("biospecimenSubtype")
    if biospecimen_subtype is not None and not pd.isna(biospecimen_subtype):
        if isinstance(biospecimen_subtype, (list, tuple)):
            for item in biospecimen_subtype:
                if item and not pd.isna(item) and str(item).strip():
                    biospecimen_subtype_values.append(str(item).strip())
        elif str(biospecimen_subtype).strip():
            biospecimen_subtype_values.append(str(biospecimen_subtype).strip())

    if len(biospecimen_type_values) > 0 or len(biospecimen_subtype_values) > 0:
        anatomical_structure: Dict[str, Any] = {}
        if len(biospecimen_type_values) > 0:
            anatomical_structure["name"] = biospecimen_type_values
        if len(biospecimen_subtype_values) > 0:
            anatomical_structure["sampleType"] = biospecimen_subtype_values
        doc["sample"] = {"anatomicalStructure": anatomical_structure}

    # Health conditions from diagnosis
    diagnosis = row.get("diagnosis")
    diagnosis_values: List[str] = []
    if diagnosis is not None and not pd.isna(diagnosis):
        if isinstance(diagnosis, (list, tuple)):
            for item in diagnosis:
                if item and not pd.isna(item) and str(item).strip():
                    diagnosis_values.append(str(item).strip())
        elif str(diagnosis).strip():
            diagnosis_values.append(str(diagnosis).strip())

    health_condition: List[Dict] = []
    for value in diagnosis_values:
        normalized = value.strip()
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
    doi = row.get("doi")
    if doi and not pd.isna(doi) and str(doi).strip():
        doc["doi"] = str(doi).strip()

    # Credit/acknowledgment statement
    acknowledgment_ref = row.get("acknowledgmentStatement")
    if acknowledgment_ref and not pd.isna(acknowledgment_ref) and str(acknowledgment_ref).strip() and hasattr(entity, "parentId"):
        credit_text = fetch_wiki_content(syn, entity.parentId, str(acknowledgment_ref).split("/")[-1])
        if credit_text:
            doc["creditText"] = credit_text
        else:
            doc["creditText"] = SYNAPSE_OBJECT_URL.format(
                syn_id=f"{entity.parentId}/wiki/{str(acknowledgment_ref).split('/')[-1]}"
            )

    # Publications (citation field)
    publication_syn_id = row.get("publicationSynID")
    pub_syn_ids: List[str] = []
    if publication_syn_id is not None and not pd.isna(publication_syn_id):
        if isinstance(publication_syn_id, (list, tuple)):
            for item in publication_syn_id:
                if item and not pd.isna(item) and str(item).strip():
                    pub_syn_ids.append(str(item).strip())
        elif str(publication_syn_id).strip():
            pub_syn_ids.append(str(publication_syn_id).strip())

    publications: List[Dict] = []
    for pub_syn_id in pub_syn_ids:
        pub_details = fetch_publication_details(syn, pub_syn_id)
        if pub_details:
            publications.append(pub_details)
    if len(publications) > 0:
        doc["citation"] = publications

    # Associated code repositories (isRelatedTo with ComputationalTool type)
    code_url = row.get("associatedCodeURL")
    code_urls: List[str] = []
    if code_url is not None and not pd.isna(code_url):
        if isinstance(code_url, (list, tuple)):
            for item in code_url:
                if item and not pd.isna(item) and str(item).strip():
                    code_urls.append(str(item).strip())
        elif str(code_url).strip():
            code_urls.append(str(code_url).strip())

    if len(code_urls) > 0:
        doc["isRelatedTo"] = [
            {"@type": "ComputationalTool", "codeRepository": url}
            for url in code_urls
        ]

    # Identifiers (dbGap, ImmPort)
    identifiers: List[Dict] = []
    for key, template in IDENTIFIER_BASE_URLS:
        value = row.get(key)
        if value and not pd.isna(value) and str(value).strip():
            identifiers.append(
                {
                    "type": key,
                    "value": str(value).strip(),
                    "url": template.format(value=str(value).strip()),
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
