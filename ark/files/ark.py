#!/usr/bin/env python3
from __future__ import annotations

import datetime
import logging
import re
from typing import Any, Dict, Generator, List, Optional

import pandas as pd
import synapseclient
from api_secret import SYNAPSE_TOKEN
from synapseclient import Synapse

# Synapse table containing ARK datasets
ARK_DATASETS_TABLE = "syn68554562"
ARK_PORTAL_DATASET_URL = "https://arkportal.synapse.org/Explore/Datasets/DetailsPage?id={syn_id}"
SYNAPSE_OBJECT_URL = "https://www.synapse.org/Synapse:{syn_id}"

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

CONDITIONS_OF_ACCESS = {
    "released": "Restricted",
    "under peer review": "Embargoed",
    "unreleased": "Embargoed",
    "test": "Test",
    "deprecated": "Deprecated",
}

TOPIC_CATEGORY_MAPPING = {
    "genomics": {
        "name": "Genomics",
        "identifier": "topic_0622",
        "url": "http://edamontology.org/topic_0622",
        "inDefinedTermSet": "EDAM",
    },
    "proteomics": {
        "name": "Proteomics",
        "identifier": "topic_0121",
        "url": "http://edamontology.org/topic_0121",
        "inDefinedTermSet": "EDAM",
    },
    "transcriptomics": {
        "name": "Transcriptomics",
        "identifier": "topic_3308",
        "url": "http://edamontology.org/topic_3308",
        "inDefinedTermSet": "EDAM",
    },
    "epigenomics": {
        "name": "Epigenomics",
        "identifier": "topic_3173",
        "url": "http://edamontology.org/topic_3173",
        "inDefinedTermSet": "EDAM",
    },
}

DIAGNOSIS_TERMS = {
    # "control": {
    #     "identifier": "C61299",
    #     "inDefinedTermSet": "NCIT",
    #     "isCurated": False,
    #     "name": "control",
    #     "url": "http://purl.obolibrary.org/obo/NCIT_C61299",
    #     "alternateName": ["Control Group"],
    # },
    "sle": {
        "identifier": "0007915",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "systemic lupus erythematosus",
        "url": "http://purl.obolibrary.org/obo/MONDO_0007915",
        "alternateName": ["SLE", "lupus erythematosus, systemic"],
    },
    "ra": {
        "identifier": "0008383",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "rheumatoid arthritis",
        "url": "http://purl.obolibrary.org/obo/MONDO_0008383",
        "alternateName": ["RA", "rheumatoid arthritis (disease)"],
    },
    "vitiligo": {
        "identifier": "0008661",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "vitiligo",
        "url": "http://purl.obolibrary.org/obo/MONDO_0008661",
        "alternateName": ["vitiligo (disease)"],
    },
    "dermatomyositis": {
        "identifier": "0016367",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "dermatomyositis",
        "url": "http://purl.obolibrary.org/obo/MONDO_0016367",
        "alternateName": ["dermatomyositis (disease)"],
    },
    "pso": {
        "identifier": "0005083",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "psoriasis",
        "url": "http://purl.obolibrary.org/obo/MONDO_0005083",
        "alternateName": ["PSO", "psoriasis (disease)"],
    },
    "psa": {
        "identifier": "0011849",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "psoriatic arthritis",
        "url": "http://purl.obolibrary.org/obo/MONDO_0011849",
        "alternateName": ["PSA", "arthritis, psoriatic", "psoriatic arthropathy"],
    },
    "scleroderma": {
        "identifier": "0019340",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "scleroderma",
        "url": "http://purl.obolibrary.org/obo/MONDO_0019340",
        "alternateName": ["systemic sclerosis", "systemic scleroderma"],
    },
    "sjd": {
        "identifier": "0010030",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "Sjogren syndrome",
        "url": "http://purl.obolibrary.org/obo/MONDO_0010030",
        "alternateName": ["SJD", "Sjögren's syndrome", "Sjogren's syndrome"],
    },
    "ln": {
        "identifier": "0005556",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "lupus nephritis",
        "url": "http://purl.obolibrary.org/obo/MONDO_0005556",
        "alternateName": ["LN", "nephritis, lupus"],
    },
    "cle": {
        "identifier": "0005282",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "cutaneous lupus erythematosus",
        "url": "http://purl.obolibrary.org/obo/MONDO_0005282",
        "alternateName": ["CLE", "lupus erythematosus, cutaneous"],
    },
    "oa": {
        "identifier": "0005178",
        "inDefinedTermSet": "MONDO",
        "isCurated": False,
        "name": "osteoarthritis",
        "url": "http://purl.obolibrary.org/obo/MONDO_0005178",
        "alternateName": ["OA", "osteoarthritis (disease)", "arthrosis"],
    },
}

IDENTIFIER_BASE_URLS = (
    (
        "dbGapAccession",
        "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={value}",
    ),
    ("ImmPortAccession", "https://www.immport.org/resources/study/{value}"),
)


def extract_grant_ids(credit_text: str) -> List[Dict[str, str]]:
    """
    Extract NIH grant IDs from acknowledgment/credit text.

    Args:
        credit_text: The acknowledgment statement containing grant information

    Returns:
        List of funding objects with identifier field containing grant IDs
    """
    if not credit_text:
        return []

    # Pattern to match NIH grant IDs (e.g., UH2-AR067676, UM2-AR067678)
    # Matches formats like: UH2-AR067676, R01-GM123456, etc.
    grant_pattern = r'\b([A-Z]{1,3}\d{1,2}-[A-Z]{2}\d{6})\b'

    matches = re.findall(grant_pattern, credit_text)

    if not matches:
        return []

    # Remove duplicates while preserving order
    unique_grants = list(dict.fromkeys(matches))

    # Return as list of funding objects
    return [{"identifier": grant_id} for grant_id in unique_grants]


def get_synapse_client(auth_token: Optional[str] = None) -> Synapse:
    """Initialize and login to Synapse client."""
    syn = synapseclient.Synapse()
    token = auth_token or SYNAPSE_TOKEN
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


def build_sample_record(dataset_doc: Dict, biospecimen_type_values: List[str], biospecimen_subtype_values: List[str]) -> Optional[Dict]:
    """Build a separate Sample record from dataset information."""
    if len(biospecimen_type_values) == 0 and len(biospecimen_subtype_values) == 0:
        return None

    dataset_id = dataset_doc.get("_id", "").replace("ark_", "")
    dataset_identifier = dataset_doc.get("identifier") or dataset_id
    dataset_url = dataset_doc.get("url")
    dataset_name = dataset_doc.get("name")

    sample_record = {
        "_id": f"{dataset_doc['_id']}_sample",
        "@type": "Sample",
        "identifier": f"{dataset_identifier}_sample",
        "name": f"{dataset_name} - Sample" if dataset_name else "Sample",
        "url": dataset_url,
        "isPartOf": {
            "@type": "Dataset",
            "identifier": dataset_identifier,
            "url": dataset_url,
        },
    }

    if dataset_name:
        sample_record["isPartOf"]["name"] = dataset_name

    # Copy relevant fields from dataset
    if dataset_doc.get("includedInDataCatalog"):
        sample_record["includedInDataCatalog"] = dataset_doc["includedInDataCatalog"]
    if dataset_doc.get("conditionsOfAccess"):
        sample_record["conditionsOfAccess"] = dataset_doc["conditionsOfAccess"]
    if dataset_doc.get("license"):
        sample_record["license"] = dataset_doc["license"]
    if dataset_doc.get("usageInfo"):
        # usageInfo.url
        sample_record["usageInfo"] = {
            "url": dataset_doc["usageInfo"]
        }
    if dataset_doc.get("measurementTechnique"):
        sample_record["measurementTechnique"] = dataset_doc["measurementTechnique"]
    if dataset_doc.get("healthCondition"):
        sample_record["healthCondition"] = dataset_doc["healthCondition"]

    # Add anatomical structure
    anatomical_structure: Dict[str, Any] = {}
    if len(biospecimen_type_values) > 0:
        anatomical_structure["name"] = biospecimen_type_values
    if len(biospecimen_subtype_values) > 0:
        anatomical_structure["sampleType"] = biospecimen_subtype_values

    sample_record["anatomicalStructure"] = anatomical_structure

    return sample_record


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
        "_id": f'ark_{syn_id}',
        "identifier": syn_id,
        "url": ARK_PORTAL_DATASET_URL.format(syn_id=syn_id),
        # Assigned fields from mapping
        "license": "https://arkportal.synapse.org/Data%20Access",
        "usageInfo": {"url":"https://help.arkportal.org/help/data-use-certificate#DataUse&Acknowledgement-Acknowledgement"},
        "includedInDataCatalog": {
            "@type": "DataCatalog",
            "name": "SAGE ARK Portal",
            "url": "https://arkportal.synapse.org/",
            "versionDate": datetime.date.today().isoformat()
        },
        "@type": "Dataset",
    }

    # Add name if present
    name = row.get("name")
    if name and not pd.isna(name) and str(name).strip():
        doc["name"] = str(name).strip()

    # Add description if present
    description = row.get("description")
    if description and not pd.isna(description) and str(description).strip():
        doc["description"] = str(description).strip()

    # Add dates if present (convert from Unix timestamp in milliseconds to YYYY-MM-DD)
    date_created = row.get("createdOn")
    if date_created and not pd.isna(date_created):
        doc["dateCreated"] = datetime.datetime.fromtimestamp(int(date_created) / 1000).strftime("%Y-%m-%d")

    date_modified = row.get("modifiedOn")
    if date_modified and not pd.isna(date_modified):
        doc["dateModified"] = datetime.datetime.fromtimestamp(int(date_modified) / 1000).strftime("%Y-%m-%d")

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
    if assay is not None:
        if isinstance(assay, (list, tuple, pd.Series)):
            for item in assay:
                if item is not None and not pd.isna(item) and str(item).strip():
                    measurement_values.append(str(item).strip())
        elif not pd.isna(assay) and str(assay).strip():
            measurement_values.append(str(assay).strip())

    data_type = row.get("dataType")
    if data_type is not None:
        if isinstance(data_type, (list, tuple, pd.Series)):
            for item in data_type:
                if item is not None and not pd.isna(item) and str(item).strip():
                    measurement_values.append(str(item).strip())
        elif not pd.isna(data_type) and str(data_type).strip():
            measurement_values.append(str(data_type).strip())

    # Remove duplicates while preserving order
    measurement = list(dict.fromkeys(measurement_values))
    if len(measurement) > 0:
        doc["measurementTechnique"] = [{"name": value} for value in measurement]

    # Topic categories based on measurement techniques
    topic_categories: List[Dict] = []
    for value in measurement_values:
        normalized = value.lower().strip()
        if normalized in TOPIC_CATEGORY_MAPPING:
            topic_term = TOPIC_CATEGORY_MAPPING[normalized].copy()
            topic_term["@type"] = "DefinedTerm"
            topic_categories.append(topic_term)

    # Remove duplicates based on identifier while preserving order
    seen_identifiers = set()
    unique_topics = []
    for topic in topic_categories:
        if topic["identifier"] not in seen_identifiers:
            seen_identifiers.add(topic["identifier"])
            unique_topics.append(topic)

    if len(unique_topics) > 0:
        doc["topicCategory"] = unique_topics

    # Sample/biospecimen information
    biospecimen_type_values: List[str] = []
    biospecimen_type = row.get("biospecimenType")
    if biospecimen_type is not None:
        if isinstance(biospecimen_type, (list, tuple, pd.Series)):
            for item in biospecimen_type:
                if item is not None and not pd.isna(item) and str(item).strip():
                    biospecimen_type_values.append(str(item).strip())
        elif not pd.isna(biospecimen_type) and str(biospecimen_type).strip():
            biospecimen_type_values.append(str(biospecimen_type).strip())

    biospecimen_subtype_values: List[str] = []
    biospecimen_subtype = row.get("biospecimenSubtype")
    if biospecimen_subtype is not None:
        if isinstance(biospecimen_subtype, (list, tuple, pd.Series)):
            for item in biospecimen_subtype:
                if item is not None and not pd.isna(item) and str(item).strip():
                    biospecimen_subtype_values.append(str(item).strip())
        elif not pd.isna(biospecimen_subtype) and str(biospecimen_subtype).strip():
            biospecimen_subtype_values.append(str(biospecimen_subtype).strip())

    # Store sample data for separate Sample record generation
    sample_data = {
        "biospecimen_type_values": biospecimen_type_values,
        "biospecimen_subtype_values": biospecimen_subtype_values,
    }

    # Health conditions from diagnosis
    diagnosis = row.get("diagnosis")
    diagnosis_values: List[str] = []
    if diagnosis is not None:
        if isinstance(diagnosis, (list, tuple, pd.Series)):
            for item in diagnosis:
                if item is not None and not pd.isna(item) and str(item).strip():
                    diagnosis_values.append(str(item).strip())
        elif not pd.isna(diagnosis) and str(diagnosis).strip():
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
            health_condition.append(lookup)
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
            # Extract grant IDs from credit text and add to funding
            funding_list = extract_grant_ids(credit_text)
            if funding_list:
                doc["funding"] = funding_list
        else:
            doc["creditText"] = SYNAPSE_OBJECT_URL.format(
                syn_id=f"{entity.parentId}/wiki/{str(acknowledgment_ref).split('/')[-1]}"
            )

    # Publications (citation field)
    publication_syn_id = row.get("publicationSynID")
    pub_syn_ids: List[str] = []
    if publication_syn_id is not None:
        if isinstance(publication_syn_id, (list, tuple, pd.Series)):
            for item in publication_syn_id:
                if item is not None and not pd.isna(item) and str(item).strip():
                    pub_syn_ids.append(str(item).strip())
        elif not pd.isna(publication_syn_id) and str(publication_syn_id).strip():
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
    if code_url is not None:
        if isinstance(code_url, (list, tuple, pd.Series)):
            for item in code_url:
                if item is not None and not pd.isna(item) and str(item).strip():
                    code_urls.append(str(item).strip())
        elif not pd.isna(code_url) and str(code_url).strip():
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
            identifiers.append(str(value).strip())
    if len(identifiers) > 0:
        doc["identifier"] = identifiers

    # Add keywords
    if len(keywords) > 0:
        doc["keywords"] = sorted(set(keywords))

    return doc, sample_data


def parse(
    auth_token: Optional[str] = None,
) -> Generator[Dict, None, None]:
    """
    Fetch and parse ARK datasets from Synapse.

    Args:
        auth_token: Synapse authentication token. If not provided, uses SYNAPSE_TOKEN from api_secret.

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
            dataset_doc, sample_data = process_dataset_row(syn, row)
            processed += 1

            # Yield the dataset record
            yield dataset_doc

            # Build and yield the sample record if we have sample data
            sample_record = build_sample_record(
                dataset_doc,
                sample_data["biospecimen_type_values"],
                sample_data["biospecimen_subtype_values"]
            )
            if sample_record:
                yield sample_record

        except Exception as error:
            logger.error("Failed to process dataset %s: %s", row.get("id", "unknown"), error, exc_info=True)
            continue

    logger.info("Finished parsing %d ARK datasets", processed)


if __name__ == "__main__":
    # Example usage
    for record in parse():
        record_type = record.get("@type", "Unknown")
        record_id = record.get("_id", "unknown")
        record_name = record.get("name", "")
        print(f"Processed {record_type}: {record_id} - {record_name}")
