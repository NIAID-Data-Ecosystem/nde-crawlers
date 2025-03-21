import logging
from functools import lru_cache

import pandas as pd
import requests

logger = logging.getLogger(__name__)

# CSV file path for mapping
CSV_FILE = "/nvme/nde-hub/standardizers/nctid_lookup/nctid.csv"

@lru_cache(maxsize=128)
def fetch_trial(nctid):
    """
    Query clinicaltrials.gov API for the given NCT identifier.
    Cached to avoid repeated calls for the same trial.
    """
    url = f"https://clinicaltrials.gov/api/int/studies/{nctid}"
    logger.info(f"Fetching trial data for NCT ID: {nctid}")
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Error fetching data for {nctid}: {response.status_code}")

def load_mapping():
    """
    Load the CSV mapping file into a pandas DataFrame.
    """
    try:
        df = pd.read_csv(CSV_FILE)
        logger.info(f"Loaded mapping CSV from {CSV_FILE} with {len(df)} rows")
        return df
    except Exception as e:
        raise Exception(f"Error loading CSV file: {e}")

def extract_trial_info(api_data):
    """
    Extract design-related fields from the API response.
    Returns study_type, intervention_model, allocation, and design_method.
    """
    try:
        design_module = api_data["study"]["protocolSection"]["designModule"]
        study_type = design_module.get("studyType", "").upper()
        design_info = design_module.get("designInfo", {})
        allocation = design_info.get("allocation", "").upper()
        intervention_model = design_info.get("interventionModel", "").upper()
        masking_info = design_info.get("maskingInfo", {})
        design_method = masking_info.get("masking", "NONE").upper()
        if not design_method:
            design_method = "NONE"
        logger.debug(f"Extracted trial info: study_type={study_type}, intervention_model={intervention_model}, allocation={allocation}, design_method={design_method}")
        return study_type, intervention_model, allocation, design_method
    except Exception as e:
        raise Exception("Error extracting trial info: " + str(e))

@lru_cache(maxsize=256)
def get_ncit_name(iri):
    """
    Given an IRI (e.g., http://purl.obolibrary.org/obo/NCIT_C82639),
    extract the NCIT code and query the Ontology Lookup Service (OLS)
    to retrieve the official NCIT term label.
    Cached to avoid repeated lookups.
    """
    try:
        identifier = iri.split("_")[-1]
        obo_id = f"NCIT:{identifier}"
        url = f"https://www.ebi.ac.uk/ols/api/ontologies/ncit/terms?obo_id={obo_id}"
        logger.info(f"Looking up NCIT name for IRI: {iri}")
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            terms = data.get('_embedded', {}).get('terms', [])
            if terms:
                name = terms[0].get('label', f"NCIT Term {identifier}")
                logger.debug(f"NCIT lookup result for {iri}: {name}")
                return name
    except Exception as e:
        logger.error(f"Error retrieving NCIT name for {iri}: {e}")
    return f"NCIT Term {identifier}"

def get_measurement_technique(study_type, intervention_model, allocation, design_method, mapping_df):
    """
    Using the extracted study fields, filter the mapping DataFrame.
    The CSV mapping file is expected to have these columns:
      studytype, studymodel, designmodel, designmethod, Exact match?, IRI, IRI.1, IRI.2
    For the matching row, this function creates a measurementTechnique list
    with one object per non-"NONE" IRI. The object's 'name' is the official NCIT term.
    """
    # Ensure the CSV columns are upper case for matching
    mapping_df['studytype'] = mapping_df['studytype'].str.upper()
    mapping_df['studymodel'] = mapping_df['studymodel'].str.upper()
    mapping_df['designmodel'] = mapping_df['designmodel'].str.upper()
    mapping_df['designmethod'] = mapping_df['designmethod'].str.upper()

    filtered = mapping_df[
        (mapping_df['studytype'] == study_type)
        & (mapping_df['studymodel'] == intervention_model)
        & (mapping_df['designmodel'] == allocation)
        & (mapping_df['designmethod'] == design_method)
    ]
    measurement_techniques = []
    if not filtered.empty:
        row = filtered.iloc[0]
        for col in ['IRI', 'IRI.1', 'IRI.2']:
            iri = row.get(col, "NONE")
            if pd.isna(iri) or str(iri).upper() == "NONE":
                continue
            # Lookup the NCIT name using the cached function
            ncit_name = get_ncit_name(iri)
            identifier = iri.split("_")[-1]
            measurement_obj = {
                "identifier": identifier,
                "inDefinedTermSet": "NCIT",
                "isCurated": False,
                "fromNCT": True,
                "name": ncit_name,
                "url": iri
            }
            measurement_techniques.append(measurement_obj)
        logger.info(f"Found {len(measurement_techniques)} measurementTechnique entries for design: {study_type}, {intervention_model}, {allocation}, {design_method}")
    else:
        logger.info(f"No mapping found for design: {study_type}, {intervention_model}, {allocation}, {design_method}")
    return measurement_techniques

def nctid_helper(docs):
    """
    For each document in the list, if an NCT identifier is present,
    fetch the clinical trial details, extract the design parameters,
    and add a measurementTechnique property to the document.
    Uses caching to avoid repeated API calls.
    """
    mapping_df = load_mapping()
    updated_docs = []
    measurement_added_count = 0
    doc_count = 0
    for doc in docs:
        doc_count += 1
        if doc_count % 100 == 0:
            logger.info(f"Processed {doc_count} documents...")
        nctid = doc.get("nctid")
        if nctid:
            try:
                api_data = fetch_trial(nctid)
                study_type, intervention_model, allocation, design_method = extract_trial_info(api_data)
                measurement_techniques = get_measurement_technique(
                    study_type,
                    intervention_model,
                    allocation,
                    design_method,
                    mapping_df
                )
                if measurement_techniques:
                    doc["measurementTechnique"] = measurement_techniques
                    measurement_added_count += 1
            except Exception as e:
                logger.error(f"Error processing NCT ID {nctid}: {e}")
        updated_docs.append(doc)
    logger.info(f"Processed {doc_count} documents. MeasurementTechnique added to {measurement_added_count} documents.")
    return updated_docs
