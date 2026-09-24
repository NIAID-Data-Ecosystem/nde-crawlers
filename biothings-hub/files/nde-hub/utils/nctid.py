"""measurementTechnique derived from a clinical trial's study design.

Runs for records that carry an `nctid`. The trial's design (study type,
intervention model, allocation, masking) is looked up on clinicaltrials.gov and
matched against the curated NCIT mapping CSV, which yields up to three NCIT
terms per design.

Designs are fetched a few hundred trials per request and cached in SQLite next
to the mapping CSV. A cached design is refetched once it is `_CACHE_DAYS` old,
since a trial's design can be amended while it is still recruiting; if
clinicaltrials.gov cannot answer, the stale design is used instead.
"""

import os
import re
from datetime import date, timedelta
from functools import cache
from itertools import batched

import pandas as pd
import requests
from config import logger

from .cache import SqliteCache
from .common import as_list, retry

CSV_FILE = "/nvme/nde-hub/standardizers/nctid_lookup/nctid.csv"
DB_PATH = os.path.join(os.path.dirname(CSV_FILE), "nctid_lookup.db")

STUDIES_URL = "https://clinicaltrials.gov/api/v2/studies"
NCIT_TERMS_URL = "https://www.ebi.ac.uk/ols4/api/ontologies/ncit/terms"

# The IDs travel in the query string: 500 is ~7 KB, and 1,000 is rejected as too long.
_REQUEST_SIZE = 400
_CACHE_DAYS = 30
# One malformed ID fails the whole request, so only well-formed ones are sent.
_NCT_ID = re.compile(r"NCT\d{8}")

_DESIGN_COLUMNS = ["studytype", "studymodel", "designmodel", "designmethod"]
_IRI_COLUMNS = ["IRI", "IRI.1", "IRI.2"]

# nctid -> {"design": the trial's designModule, "fetched": ISO date}
DESIGNS = SqliteCache(DB_PATH, "trial_designs", "nctid", "design", normalize=False)
# NCIT IRI -> label, only for labels OLS actually returned
NCIT_NAMES = SqliteCache(DB_PATH, "ncit_names", "iri", "name", preload=True, normalize=False)


def lookup_file():
    return CSV_FILE


def reset_caches():
    """Re-read the mapping CSV and the SQLite caches on the next upload."""
    load_mapping.cache_clear()
    DESIGNS.reset()
    NCIT_NAMES.reset()


@cache
def load_mapping():
    """Load the design -> NCIT mapping, upper-cased for matching."""
    df = pd.read_csv(CSV_FILE)
    for column in _DESIGN_COLUMNS:
        df[column] = df[column].str.upper()
    logger.info("Loaded NCT mapping CSV from %s with %s rows", CSV_FILE, len(df))
    return df


def nct_ids(value):
    """The well-formed NCT IDs in a record's `nctid`, upper-cased."""
    ids = []
    for item in as_list(value):
        nctid = str(item).strip().upper()
        if _NCT_ID.fullmatch(nctid):
            ids.append(nctid)
        else:
            logger.warning("Skipping malformed NCT ID %r", item)
    return ids


@retry(3, 10)
def _fetch_designs(nctids):
    """{nctid: designModule} for up to `_REQUEST_SIZE` trials, in one request.

    Trials clinicaltrials.gov does not know are simply absent from the answer.
    """
    response = requests.get(
        STUDIES_URL,
        params={"filter.ids": ",".join(nctids), "fields": "NCTId,DesignModule", "pageSize": len(nctids)},
        timeout=120,
    )
    response.raise_for_status()
    designs = {}
    for study in response.json().get("studies", []):
        protocol = study.get("protocolSection", {})
        if nctid := protocol.get("identificationModule", {}).get("nctId"):
            designs[nctid.upper()] = protocol.get("designModule", {})
    return designs


def trial_designs(nctids):
    """{nctid: designModule} for `nctids`, from the cache or clinicaltrials.gov."""
    nctids = set(nctids)
    today = date.today()
    fresh_after = (today - timedelta(days=_CACHE_DAYS)).isoformat()
    cached = DESIGNS.get_many(nctids)
    designs = {nctid: entry["design"] for nctid, entry in cached.items() if entry["fetched"] > fresh_after}

    to_fetch = sorted(nctids - designs.keys())
    fetched = {}
    failed = set()
    for chunk in batched(to_fetch, _REQUEST_SIZE):
        try:
            fetched.update(_fetch_designs(chunk))
        except Exception as e:
            failed.update(chunk)
            logger.error("Could not fetch %s trial designs from clinicaltrials.gov: %s", len(chunk), e)
    DESIGNS.put_many({nctid: {"design": design, "fetched": today.isoformat()} for nctid, design in fetched.items()})
    designs.update(fetched)

    unanswered = set(to_fetch) - fetched.keys()
    stale = {nctid: cached[nctid]["design"] for nctid in unanswered if nctid in cached}
    designs.update(stale)

    logger.info(
        "NCT trial designs: trials=%s cache_hits=%s fetched=%s requests=%s not_found=%s failed=%s stale_used=%s",
        len(nctids),
        len(nctids) - len(to_fetch),
        len(fetched),
        -(-len(to_fetch) // _REQUEST_SIZE),
        len(unanswered - failed),
        len(failed),
        len(stale),
    )
    return designs


def design_key(design_module):
    """(study_type, intervention_model, allocation, design_method) of a trial's designModule."""
    design_info = design_module.get("designInfo", {})
    return (
        design_module.get("studyType", "").upper(),
        design_info.get("interventionModel", "").upper(),
        design_info.get("allocation", "").upper(),
        design_info.get("maskingInfo", {}).get("masking", "NONE").upper() or "NONE",
    )


def get_ncit_name(iri):
    """Resolve an NCIT IRI to its official term label."""
    if name := NCIT_NAMES.get(iri):
        return name
    identifier = iri.split("_")[-1]
    try:
        response = requests.get(NCIT_TERMS_URL, params={"obo_id": f"NCIT:{identifier}"}, timeout=60)
        if response.status_code == 200:
            terms = response.json().get("_embedded", {}).get("terms", [])
            if terms and (name := terms[0].get("label")):
                NCIT_NAMES.put(iri, name)
                return name
    except Exception as e:
        logger.error("Error retrieving NCIT name for %s: %s", iri, e)
    return f"NCIT Term {identifier}"


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
        logger.debug("No mapping found for design: %s", design)
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
    logger.debug("Found %s measurementTechnique entries for design: %s", len(measurement_techniques), design)
    return measurement_techniques


def add_nct_measurement_techniques(docs):
    """Add design-derived measurementTechnique to every record with an nctid.

    `docs` is one pipeline batch; every trial in it is looked up at once.
    """
    docs = list(docs)
    mapping_df = load_mapping()
    ids_per_doc = [nct_ids(doc.get("nctid")) for doc in docs]
    designs = trial_designs(nctid for ids in ids_per_doc for nctid in ids)

    techniques_by_design = {}
    added = 0
    for doc, ids in zip(docs, ids_per_doc):
        techniques = []
        for nctid in ids:
            if (design := designs.get(nctid)) is None:
                logger.debug("No trial design for NCT ID %s in %s", nctid, doc.get("_id"))
                continue
            key = design_key(design)
            if key not in techniques_by_design:
                techniques_by_design[key] = get_measurement_technique(key, mapping_df)
            techniques += [t for t in techniques_by_design[key] if t not in techniques]
        if techniques:
            # Each record gets its own copies, so a later stage editing one cannot touch another.
            doc["measurementTechnique"] = [dict(technique) for technique in techniques]
            added += 1

    logger.info("MeasurementTechnique added to %s documents from NCT trial designs", added)
    return docs
