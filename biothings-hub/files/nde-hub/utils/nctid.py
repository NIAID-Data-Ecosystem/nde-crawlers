"""measurementTechnique derived from a clinical trial's study design.

Runs for records that carry an `nctid`. The trial's design (study type,
intervention model, allocation, masking) is looked up on clinicaltrials.gov and
matched against the curated NCIT mapping CSV, which yields up to three NCIT
terms per design.
"""

from functools import lru_cache

import pandas as pd
import requests
from config import logger

CSV_FILE = "/nvme/nde-hub/standardizers/nctid_lookup/nctid.csv"

_DESIGN_COLUMNS = ["studytype", "studymodel", "designmodel", "designmethod"]
_IRI_COLUMNS = ["IRI", "IRI.1", "IRI.2"]


def lookup_file():
    return CSV_FILE


@lru_cache(maxsize=None)
def load_mapping():
    """Load the design -> NCIT mapping, upper-cased for matching."""
    df = pd.read_csv(CSV_FILE)
    for column in _DESIGN_COLUMNS:
        df[column] = df[column].str.upper()
    logger.info("Loaded NCT mapping CSV from %s with %s rows", CSV_FILE, len(df))
    return df


@lru_cache(maxsize=4096)
def fetch_trial(nctid):
    """Fetch one trial from clinicaltrials.gov."""
    logger.info("Fetching trial data for NCT ID: %s", nctid)
    response = requests.get(f"https://clinicaltrials.gov/api/int/studies/{nctid}")
    if response.status_code != 200:
        raise Exception(f"Error fetching data for {nctid}: {response.status_code}")
    return response.json()


@lru_cache(maxsize=256)
def get_ncit_name(iri):
    """Resolve an NCIT IRI to its official term label."""
    identifier = iri.split("_")[-1]
    try:
        response = requests.get(f"https://www.ebi.ac.uk/ols/api/ontologies/ncit/terms?obo_id=NCIT:{identifier}")
        if response.status_code == 200:
            terms = response.json().get("_embedded", {}).get("terms", [])
            if terms:
                return terms[0].get("label", f"NCIT Term {identifier}")
    except Exception as e:
        logger.error("Error retrieving NCIT name for %s: %s", iri, e)
    return f"NCIT Term {identifier}"


def extract_trial_info(api_data):
    """Return (study_type, intervention_model, allocation, design_method) for a trial."""
    try:
        design_module = api_data["study"]["protocolSection"]["designModule"]
        design_info = design_module.get("designInfo", {})
        return (
            design_module.get("studyType", "").upper(),
            design_info.get("interventionModel", "").upper(),
            design_info.get("allocation", "").upper(),
            design_info.get("maskingInfo", {}).get("masking", "NONE").upper() or "NONE",
        )
    except Exception as e:
        raise Exception("Error extracting trial info: " + str(e))


def get_measurement_technique(design, mapping_df):
    """Return the NCIT measurementTechnique terms curated for a trial design."""
    study_type, intervention_model, allocation, design_method = design
    filtered = mapping_df[
        (mapping_df["studytype"] == study_type)
        & (mapping_df["studymodel"] == intervention_model)
        & (mapping_df["designmodel"] == allocation)
        & (mapping_df["designmethod"] == design_method)
    ]
    if filtered.empty:
        logger.info("No mapping found for design: %s", design)
        return []

    row = filtered.iloc[0]
    measurement_techniques = []
    for column in _IRI_COLUMNS:
        iri = row.get(column, "NONE")
        if pd.isna(iri) or str(iri).upper() == "NONE":
            continue
        measurement_techniques.append(
            {
                "@type": "DefinedTerm",
                "identifier": iri.split("_")[-1],
                "inDefinedTermSet": "NCIT",
                "isCurated": False,
                "fromNCT": True,
                "name": get_ncit_name(iri),
                "url": iri,
            }
        )
    logger.info("Found %s measurementTechnique entries for design: %s", len(measurement_techniques), design)
    return measurement_techniques


def add_nct_measurement_techniques(docs):
    """Add design-derived measurementTechnique to every record with an nctid."""
    mapping_df = load_mapping()
    added = 0

    for doc in docs:
        if nctid := doc.get("nctid"):
            try:
                design = extract_trial_info(fetch_trial(nctid))
                if measurement_techniques := get_measurement_technique(design, mapping_df):
                    doc["measurementTechnique"] = measurement_techniques
                    added += 1
            except Exception as e:
                logger.error("Error processing NCT ID %s: %s", nctid, e)
        yield doc

    logger.info("MeasurementTechnique added to %s documents from NCT trial designs", added)
