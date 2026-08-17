"""Citations, funding, species and health conditions derived from PubMed IDs.

Runs whenever a batch of records carries `pmids`, `pmcs` or a `citation.doi`:

  * `pmcs` and `citation.doi` are converted to PMIDs,
  * PMIDs are batch-queried against NCBI E-utilities for citation + funding
    metadata (cached in SQLite, so repeat runs skip NCBI entirely),
  * species and diseases PubTator annotated on those PMIDs are added to the
    record when the term also appears in its title / abstract / description.

Helpful links used to write this module:
https://biopython.org/docs/1.76/api/Bio.Entrez.html
https://www.nlm.nih.gov/bsd/mms/medlineelements.html
https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/
"""

import csv
import gzip
import json
import os
import re
import sqlite3
import time
import urllib.error
from datetime import datetime
from email.utils import parsedate_to_datetime
from itertools import batched
from typing import Dict, Iterable, Optional

import orjson
import requests
from Bio import Entrez, Medline
from config import GEO_API_KEY, GEO_EMAIL, logger

from .common import as_list, dict_entries, retry
from .funding import standardize_funder
from .term_matching import (
    is_ambiguous_short_mention as _is_ambiguous_short_mention,
    mentioned_in as _mentioned_in,
    normalize_term_text as _normalize_term_text,
    species_term_matches_mention as _species_term_matches_mention,
    term_matches_mention as _term_matches_mention,
)
from .terms import DB_PATH as PUBTATOR_DB_PATH, get_species_details, query_condition

PMID_DB_PATH = "/data/nde-hub/standardizers/pmid_lookup/pmid_lookup.db"
PUBTATOR_DIR = "/data/nde-hub/standardizers/pmid_lookup/"

# PubTator3 annotation dumps, keyed by the table they populate. Updated monthly.
# Five tab-separated columns: PMID, Type, Concept ID, Mentions, Resource --
# https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/README.txt
PUBTATOR_DUMPS = {
    "species": "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/species2pubtator3.gz",
    "disease": "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/disease2pubtator3.gz",
}

# We keep three of them. They become `{species,disease}_data` rows of
# (pmid, entity_id, names), which `get_data_for_pmids` reads back as
# {concept id: [mentions]}. In the document that finally gets yielded:
# The uniprot lineage of a species is used to classify it as either `species` or `infectiousAgent`.
#   Concept ID -> the term's `identifier` and `url`. A taxon id for
#                 species (resolved against UniProt); a MeSH id for disease, used to
#                 find its MONDO/DOID/HPO/NCIT equivalent and kept only if there
#                 isn't one.
#   Mentions   -> the term's `name`. usually not part of the output. One of them has to appear in the record's
#                 name, abstract or description before we add the term at all, and
#                 the first match becomes the cache key. `name` usually comes from
#                 UniProt or the ontology (MESH ID). mention is only used as `name` when
#                 UniProt has no scientific name.
#
# So a "goats" mention on taxon 9925 yields species: [{name: "Capra hircus",
# identifier: "9925", fromPMID: True, ...}] -- the mention itself is gone.
_PMID_COLUMN, _CONCEPT_ID_COLUMN, _MENTIONS_COLUMN = 0, 2, 3

# Rows per executemany when loading a dump; also how often the load commits.
_DUMP_INSERT_BATCH = 50_000

# NCBI rejects an entire request if any id in it is malformed.
_VALID_PMID_RE = re.compile(r"^[1-9]\d{0,8}$")

# The PMC id converter accepts a couple hundred ids per request.
_PMC_CHUNK_SIZE = 200

_MESH_REQUEST_ATTEMPTS = 7
_MESH_REQUEST_RETRY_SECONDS = 5
_MESH_RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
_MESH_TIMEOUT = (5, 30)

# PubTator sometimes annotates these as diseases for COVID-19; they are not.
_INCORRECT_COVID_TERMS = frozenset(
    {
        "novel tumor",
        "hypervirulent covs",
        "covr-covs",
        "cancer stemness affords novel cancer",
        "2 tumors",
        "novel disease",
        "novel icos deficiency",
        "kucap-2 tumors",
        "hif-1/hif-2",
        "ncp",
    }
)
_COVID_MESH_ID = "MESH:C000657245"
_COVID_MESH_REPLACEMENT = "MESH:D000086382"

# Species names PubTator picks up that are never the study organism.
_SPECIES_BLACKLIST = frozenset({"PERCH", "D-FISH"})

# Real words, but not health conditions. PubTator3 tags them anyway.
_NON_SPECIFIC_DISEASE_MENTIONS = frozenset(
    {
        "acute",
        "antibiotic",
        "blood",
        "cell",
        "child",
        "children",
        "chronic",
        "clinical",
        "dead",
        "death",
        "development",
        "die",
        "died",
        "deaths",
        "dying",
        "faeces",
        "feces",
        "food insecurity",
        "geographic",
        "infant",
        "infected",
        "infection",
        "infections",
        "infectious",
        "infectious disease",
        "inflammatory",
        "low",
        "maternal",
        "mortality",
        "mortality rate",
        "newborn",
        "nutrition",
        "protozoal",
        "systemic",
        "wash",
        "weight",
        "weight gain",
    }
)

# Resolved concepts that are measurements or demographics, never conditions
# Unambiguous cases only: symptoms and broad disease families stay eligible.
_NON_CONDITION_TERM_NAMES = frozenset(
    {
        "dead",
        "death",
        "death domain",
        "feces",
        "food insecurity",
        "infected with sars cov 2",
        "inflammatory",
        "mortality rate",
        "newborn",
        "nutrition",
        "protozoal",
        "weight gain",
    }
)
_NON_CONDITION_TERM_PREFIXES = ("how often experienced ", "obsolete ")

_pmid_conn = None
_pubtator_conn = None
_pubtator_cache = {}
_incompatible_disease_terms = set()
_dumps_checked = False


class PubTatorDumpUnavailable(RuntimeError):
    """A PubTator annotation dump could not be reached at its configured URL."""


# ---------------------------------------------------------------------------
# Connections and one-time setup
# ---------------------------------------------------------------------------
def _get_pmid_conn():
    """Open (once) the PMID lookup DB, creating its tables if needed."""
    global _pmid_conn
    if _pmid_conn is None:
        os.makedirs(os.path.dirname(PMID_DB_PATH), exist_ok=True)
        _pmid_conn = sqlite3.connect(PMID_DB_PATH)
        cur = _pmid_conn.cursor()
        for entity_type in PUBTATOR_DUMPS:
            cur.execute(
                f"""CREATE TABLE IF NOT EXISTS {entity_type}_data
                       (pmid TEXT, entity_id TEXT, names TEXT, PRIMARY KEY (pmid, entity_id))"""
            )
        cur.execute("""CREATE TABLE IF NOT EXISTS eutils_cache (pmid TEXT PRIMARY KEY, data TEXT)""")
        cur.execute("""CREATE TABLE IF NOT EXISTS doi_cache (doi TEXT PRIMARY KEY, pmid TEXT)""")
        _pmid_conn.commit()
    return _pmid_conn


def _get_pubtator_conn():
    global _pubtator_conn
    if _pubtator_conn is None:
        _pubtator_conn = sqlite3.connect(PUBTATOR_DB_PATH)
    return _pubtator_conn


def _refresh_pubtator_dumps():
    """Download and load the PubTator annotation dumps if they changed upstream.

    Runs at most once per successful check per process, and only once a source
    actually has PMIDs to look up -- the dumps are large and the freshness check
    costs a request each. Raises `PubTatorDumpUnavailable` if a dump can't be
    reached, rather than carrying on with whatever was loaded last.
    """
    global _dumps_checked
    if _dumps_checked:
        return
    _get_pmid_conn()
    os.makedirs(PUBTATOR_DIR, exist_ok=True)

    for entity_type, url in PUBTATOR_DUMPS.items():
        filename = os.path.basename(url)
        if not _file_needs_update(url, filename):
            continue
        logger.info("Downloading %s data from %s...", entity_type, url)
        _download_file(url, filename)
        logger.info("Storing %s data...", entity_type)
        _stream_and_store(filename, entity_type)

    # Only after every dump checked out, so a failure is raised again on the
    # next call instead of being skipped as "already checked".
    _dumps_checked = True


def _download_file(url, local_filename):
    full_path = os.path.join(PUBTATOR_DIR, local_filename)
    partial_path = f"{full_path}.part"
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(partial_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        os.replace(partial_path, full_path)
    except Exception:
        try:
            os.remove(partial_path)
        except FileNotFoundError:
            pass
        raise


def _file_needs_update(url, local_filename):
    """True when the remote dump is newer than our copy.

    Raises `PubTatorDumpUnavailable` if the URL doesn't resolve. NCBI has moved
    these files before (PubTatorCentral -> PubTator3) and swallowing that just
    means the annotation tables quietly serve years-old data, so a dead URL
    fails the upload until `PUBTATOR_DUMPS` is corrected.
    """
    full_path = os.path.join(PUBTATOR_DIR, local_filename)
    response = requests.head(url)
    if response.status_code != 200:
        raise PubTatorDumpUnavailable(
            f"PubTator dump unreachable (HTTP {response.status_code}): {url}. "
            f"Check {os.path.dirname(url)}/ and update PUBTATOR_DUMPS in utils/citations.py."
        )

    # No local copy: take whatever is there, headers or not.
    if not os.path.exists(full_path):
        return True

    last_modified = response.headers.get("Last-Modified")
    if not last_modified:
        logger.warning("No Last-Modified for %s; keeping the copy already on disk", url)
        return False

    # parsedate_to_datetime handles every date form HTTP allows and returns an
    # aware datetime, so compare POSIX timestamps and stay clear of naive/aware.
    try:
        remote_timestamp = parsedate_to_datetime(last_modified).timestamp()
    except (TypeError, ValueError) as e:
        logger.warning("Unparseable Last-Modified %r for %s: %s", last_modified, url, e)
        return False
    return remote_timestamp > os.path.getmtime(full_path)


def _dump_rows(filename):
    """Yield (pmid, concept id, mentions) from a PubTator dump, skipping short rows."""
    with gzip.open(os.path.join(PUBTATOR_DIR, filename), "rt") as file:
        for row in csv.reader(file, delimiter="\t"):
            if len(row) <= _MENTIONS_COLUMN:
                continue
            yield row[_PMID_COLUMN], row[_CONCEPT_ID_COLUMN], row[_MENTIONS_COLUMN]


def _stream_and_store(filename, entity_type):
    """Replace one PubTator annotation table from a complete dump.

    The staging table is committed in batches because these dumps contain
    hundreds of millions of rows. The active table remains untouched until the
    load completes, then the two are swapped in one transaction. This both
    preserves the old cache after an interrupted load and removes annotations
    that disappeared upstream instead of retaining them through upserts.
    """
    if entity_type not in PUBTATOR_DUMPS:
        raise ValueError(f"Unknown PubTator entity type: {entity_type}")

    conn = _get_pmid_conn()
    active_table = f"{entity_type}_data"
    staging_table = f"{active_table}_staging"
    conn.execute(f"DROP TABLE IF EXISTS {staging_table}")
    conn.execute(
        f"""CREATE TABLE {staging_table}
               (pmid TEXT, entity_id TEXT, names TEXT, PRIMARY KEY (pmid, entity_id))"""
    )
    conn.commit()

    total = 0
    try:
        insert = f"""INSERT INTO {staging_table} (pmid, entity_id, names) VALUES (?, ?, ?)
                     ON CONFLICT(pmid, entity_id) DO UPDATE SET names=excluded.names"""
        for rows in batched(_dump_rows(filename), _DUMP_INSERT_BATCH):
            conn.executemany(insert, rows)
            conn.commit()
            total += len(rows)
            if total % (_DUMP_INSERT_BATCH * 20) == 0:
                logger.info("Loaded %s %s annotations", f"{total:,}", entity_type)

        conn.execute("BEGIN IMMEDIATE")
        conn.execute(f"DROP TABLE {active_table}")
        conn.execute(f"ALTER TABLE {staging_table} RENAME TO {active_table}")
        conn.commit()
    except Exception:
        conn.rollback()
        conn.execute(f"DROP TABLE IF EXISTS {staging_table}")
        conn.commit()
        raise

    logger.info("Loaded %s %s annotations from %s", f"{total:,}", entity_type, filename)


def get_data_for_pmids(pmids, entity_type):
    """Return {entity_id: [names]} PubTator annotated on any of `pmids`."""
    cur = _get_pmid_conn().cursor()
    placeholders = ", ".join("?" for _ in pmids)
    cur.execute(f"SELECT entity_id, names FROM {entity_type}_data WHERE pmid IN ({placeholders})", pmids)

    data = {}
    for entity_id, names in cur.fetchall():
        data.setdefault(entity_id, []).extend(names.split("|"))
    return data


# ---------------------------------------------------------------------------
# Standardized term lookups
# ---------------------------------------------------------------------------
def pubtator_lookup(name, table):
    """Return a fresh copy of the cached standardized term for `name`, if any."""
    cache_key = (name.lower().strip(), table)
    if cache_key in _pubtator_cache:
        return json.loads(_pubtator_cache[cache_key])
    c = _get_pubtator_conn().cursor()
    c.execute(f"SELECT * FROM {table} WHERE original_name=?", (cache_key[0],))
    result = c.fetchone()
    if result:
        _pubtator_cache[cache_key] = result[1]
        return json.loads(result[1])
    return None


def pubtator_add(name, table, standard_dict):
    conn = _get_pubtator_conn()
    c = conn.cursor()
    normalized_name = name.lower().strip()
    # Keep one authoritative mapping for each mention. This also lets a newly
    # validated result replace a stale incompatible entry already on disk.
    c.execute(f"DELETE FROM {table} WHERE original_name=?", (normalized_name,))
    c.execute(f"INSERT INTO {table} VALUES (?, ?)", (normalized_name, standard_dict))
    conn.commit()
    _pubtator_cache[(normalized_name, table)] = standard_dict


def _as_augmented(term):
    """Mark a standardized term as PMID-derived rather than curated."""
    term["fromPMID"] = True
    term["isCurated"] = False
    term.pop("curatedBy", None)
    term.pop("originalName", None)
    return term


def _is_non_condition_term(term):
    """True for resolved concepts that cannot represent a health condition."""
    normalized_name = _normalize_term_text(term.get("name"))
    return normalized_name in _NON_CONDITION_TERM_NAMES or normalized_name.startswith(_NON_CONDITION_TERM_PREFIXES)


def _mesh_literal(value, preferred_language="en"):
    """Return one string from a JSON-LD literal or list of literals."""
    if isinstance(value, str):
        return value

    literals = []
    for item in as_list(value):
        if not isinstance(item, dict) or not item.get("@value"):
            continue
        literals.append(item)

    if not literals:
        return None
    for item in literals:
        if item.get("@language") == preferred_language:
            return item["@value"]
    for item in literals:
        if not item.get("@language"):
            return item["@value"]
    return literals[0]["@value"]


def _fetch_mesh_record(identifier):
    """Fetch MeSH JSON, retrying only transient request failures."""
    url = f"https://id.nlm.nih.gov/mesh/{identifier}.json"
    for attempt in range(_MESH_REQUEST_ATTEMPTS):
        try:
            response = requests.get(url, timeout=_MESH_TIMEOUT)
            response.raise_for_status()
        except requests.HTTPError as error:
            status_code = error.response.status_code if error.response is not None else None
            if status_code not in _MESH_RETRY_STATUSES or attempt + 1 == _MESH_REQUEST_ATTEMPTS:
                raise
            logger.warning(
                "MeSH returned HTTP %s for %s; retrying in %ss",
                status_code,
                identifier,
                _MESH_REQUEST_RETRY_SECONDS,
            )
        except (requests.ConnectionError, requests.Timeout) as error:
            if attempt + 1 == _MESH_REQUEST_ATTEMPTS:
                raise
            logger.warning(
                "MeSH request failed for %s (%s); retrying in %ss",
                identifier,
                type(error).__name__,
                _MESH_REQUEST_RETRY_SECONDS,
            )
        else:
            # JSON decoding and schema errors are deterministic for a response;
            # keep them outside the retry handlers so they fail immediately.
            return response.json()

        time.sleep(_MESH_REQUEST_RETRY_SECONDS)

    raise RuntimeError("MeSH request retry loop ended unexpectedly")


def get_disease_details(identifier, original_name):
    """Standardize a disease from its MeSH id, preferring an existing ontology term."""
    identifier = identifier.split(":")[-1]
    cache_key = (identifier, original_name.casefold().strip())

    if cache_key in _incompatible_disease_terms:
        return None

    if lookup_result := pubtator_lookup(original_name, "health_conditions"):
        if not _is_non_condition_term(lookup_result) and _term_matches_mention(lookup_result, original_name):
            return _as_augmented(lookup_result)
        logger.debug(
            "Ignoring incompatible cached disease mapping for %s: %s",
            original_name,
            lookup_result.get("name"),
        )

    logger.debug("Converting %s from MeSH %s to standard format", original_name, identifier)
    if non_mesh_result := query_condition(original_name, identifier):
        if _is_non_condition_term(non_mesh_result) or not _term_matches_mention(non_mesh_result, original_name):
            logger.debug(
                "Ignoring incompatible ontology mapping for %s: %s",
                original_name,
                non_mesh_result.get("name"),
            )
        else:
            pubtator_add(original_name, "health_conditions", json.dumps(non_mesh_result))
            return _as_augmented(non_mesh_result)

    logger.debug("Fetching details for %s with ID %s", original_name, identifier)
    disease_info = _fetch_mesh_record(identifier)

    standard_dict = {
        "@type": "DefinedTerm",
        "identifier": identifier,
        "inDefinedTermSet": "MeSH",
        "url": f"https://id.nlm.nih.gov/mesh/{identifier}.html",
        "isCurated": False,
    }
    if terms := disease_info.get("terms"):
        alternative_names = []
        for term in terms:
            if not isinstance(term, dict) or not (term_label := _mesh_literal(term.get("label"))):
                continue
            if term.get("preferred"):
                standard_dict["name"] = term_label
            else:
                alternative_names.append(term_label)
        if alternative_names:
            standard_dict["alternateName"] = alternative_names
    if label := _mesh_literal(disease_info.get("label")):
        standard_dict["name"] = label
    if "name" not in standard_dict:
        raise Exception(f"No name found for {identifier}")

    if _is_non_condition_term(standard_dict) or not _term_matches_mention(standard_dict, original_name):
        logger.debug("Ignoring incompatible MeSH mapping for %s: %s", original_name, standard_dict.get("name"))
        _incompatible_disease_terms.add(cache_key)
        return None

    pubtator_add(original_name, "health_conditions", json.dumps(standard_dict))
    standard_dict["fromPMID"] = True
    return standard_dict


# ---------------------------------------------------------------------------
# Adding PubTator species / diseases to a record
# ---------------------------------------------------------------------------
def remove_first_by_name(lst, target):
    target_lower = target.lower()
    for i, item in enumerate(lst):
        if item.get("name", "").lower() == target_lower:
            lst.pop(i)
            break


def is_bacdive_record(rec):
    for catalog in as_list(rec.get("includedInDataCatalog")):
        if isinstance(catalog, dict) and str(catalog.get("name", "")).lower() == "bacdive":
            return True
        if isinstance(catalog, str) and catalog.lower() == "bacdive":
            return True
    return False


def _record_text(rec):
    """The record text a PubTator term must appear in before we trust it."""
    return (
        rec.get("abstract", "").lower(),
        rec.get("description", "").lower(),
        rec.get("name", "").lower(),
    )


def update_record_disease(rec, disease_data):
    """Add PubTator diseases mentioned in the record's own text."""
    if isinstance(rec.get("healthCondition"), dict):
        rec["healthCondition"] = [rec["healthCondition"]]

    haystacks = _record_text(rec)

    for mesh_id, diseases in disease_data.items():
        if "MESH" not in mesh_id:
            logger.debug("Invalid MeSH ID %s", mesh_id)
            continue
        for disease in diseases:
            name = disease.strip()
            if not name or not mesh_id.strip():
                logger.debug("Empty disease name or MeSH ID: %r / %r", name, mesh_id)
                continue
            if name.lower() in _INCORRECT_COVID_TERMS and mesh_id == _COVID_MESH_ID:
                logger.debug("Incorrect Covid-19 mapping found for %s", name)
                continue
            if mesh_id == _COVID_MESH_ID:
                mesh_id = _COVID_MESH_REPLACEMENT
            if not _mentioned_in(name, haystacks):
                continue
            if _is_ambiguous_short_mention(name):
                logger.debug("Skipping ambiguous short disease mention %r in %s", name, rec.get("_id"))
                continue
            if _normalize_term_text(name) in _NON_SPECIFIC_DISEASE_MENTIONS:
                logger.debug("Skipping non-specific disease mention %r in %s", name, rec.get("_id"))
                continue

            logger.debug("Adding %s to record %s", name, rec["_id"])
            try:
                standardized_dict = get_disease_details(mesh_id, name)
            except Exception as e:
                logger.warning("Could not get details for %s with ID %s: %s", name, mesh_id, e)
                continue
            if not standardized_dict:
                continue

            if any(d.get("name", "").lower() == name.lower() for d in rec.get("healthCondition", [])):
                remove_first_by_name(rec["healthCondition"], name)
            rec["healthCondition"] = rec.get("healthCondition", []) + [standardized_dict]
            break


def update_record_species(rec, species_data):
    """Add PubTator species mentioned in the record's own text."""
    if is_bacdive_record(rec):
        logger.debug("Skipping PMID species augmentation for BacDive record %s", rec.get("_id"))
        return

    if isinstance(rec.get("species"), dict):
        rec["species"] = [rec["species"]]
    if isinstance(rec.get("infectiousAgent"), dict):
        rec["infectiousAgent"] = [rec["infectiousAgent"]]

    haystacks = _record_text(rec)

    for taxonomy_id, species_names in species_data.items():
        for name in species_names:
            name = name.strip()
            if not name or not taxonomy_id.strip():
                logger.debug("Empty species name or taxonomy ID: %r / %r", name, taxonomy_id)
                continue
            if not _mentioned_in(name, haystacks):
                continue
            if _is_ambiguous_short_mention(name):
                logger.debug("Skipping ambiguous short species mention %r in %s", name, rec.get("_id"))
                continue
            if name.upper() in _SPECIES_BLACKLIST:
                logger.debug("Blacklisted: %s in record: %s, skipping", name, rec["_id"])
                continue

            if lookup_result := pubtator_lookup(name, "species"):
                if not _species_term_matches_mention(lookup_result, name):
                    logger.debug(
                        "Ignoring incompatible cached species mapping for %s: %s",
                        name,
                        lookup_result.get("name"),
                    )
                    continue
                standardized_dict = _as_augmented(lookup_result)
            else:
                try:
                    standardized_dict = get_species_details(name, taxonomy_id)
                except Exception as e:
                    logger.warning("Could not get details for %s with ID %s: %s", name, taxonomy_id, e)
                    continue
                if standardized_dict is None:
                    logger.debug("Skipping %s with ID %s: filtered by drop list", name, taxonomy_id)
                    continue
                if not _species_term_matches_mention(standardized_dict, name):
                    logger.debug(
                        "Ignoring incompatible species mapping for %s: %s",
                        name,
                        standardized_dict.get("name"),
                    )
                    continue
                pubtator_add(name, "species", json.dumps(standardized_dict))
                standardized_dict = _as_augmented(standardized_dict)

            if any(spec.get("name", "").lower() == name.lower() for spec in rec.get("species", [])):
                remove_first_by_name(rec["species"], name)
            elif any(spec.get("name", "").lower() == name.lower() for spec in rec.get("infectiousAgent", [])):
                remove_first_by_name(rec["infectiousAgent"], name)

            if standardized_dict.get("classification") == "infectiousAgent":
                rec["infectiousAgent"] = rec.get("infectiousAgent", []) + [standardized_dict]
            else:
                if "classification" not in standardized_dict:
                    logger.debug("Could not classify %s with ID %s", name, taxonomy_id)
                rec["species"] = rec.get("species", []) + [standardized_dict]
            break


# ---------------------------------------------------------------------------
# Identifier conversion
# ---------------------------------------------------------------------------
@retry(7, 5)
def _convert_pmc_chunk(pmc_ids, pmc_pmid):
    base_url = (
        "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=my_email@example.com&format=json&"
    )
    request = requests.get(base_url + "ids=" + ",".join(pmc_ids)).json()
    for record in request.get("records"):
        if pmid := record.get("pmid"):
            pmc_pmid[record.get("pmcid")] = str(pmid)
    time.sleep(0.5)


def _convert_pmcs(pmc_ids):
    """Convert PMC ids to PMIDs, a couple hundred at a time."""
    pmc_pmid = {}
    for chunk in batched(sorted(set(pmc_ids)), _PMC_CHUNK_SIZE):
        _convert_pmc_chunk(chunk, pmc_pmid)
    return pmc_pmid


@retry(7, 5)
def _convert_doi(doi, doi_pmid):
    """Resolve a DOI to a PMID via ESearch, caching the answer (including misses)."""
    if doi in doi_pmid:
        return

    api_key = GEO_API_KEY
    Entrez.email = GEO_EMAIL
    if api_key:
        Entrez.api_key = api_key

    conn = _get_pmid_conn()
    cur = conn.cursor()
    cur.execute("SELECT pmid FROM doi_cache WHERE doi = ?", (doi,))
    row = cur.fetchone()
    if row is not None:
        doi_pmid[doi] = row[0]  # may be None if previously resolved to nothing
        return

    handle = Entrez.esearch(db="pubmed", term=f"{doi}[DOI]", retmode="json")
    data = json.loads(handle.read())
    handle.close()
    pmids = data.get("esearchresult", {}).get("idlist", [])
    result_pmid = pmids[0] if pmids else None
    doi_pmid[doi] = result_pmid

    cur.execute("INSERT OR REPLACE INTO doi_cache (doi, pmid) VALUES (?, ?)", (doi, result_pmid))
    conn.commit()
    time.sleep(0.05 if api_key else 0.35)


# ---------------------------------------------------------------------------
# E-utilities
# ---------------------------------------------------------------------------
def _get_pub_date(date: str):
    """Turn a MedLine publication date such as "2000 Spring" into an ISO date.

    https://www.nlm.nih.gov/bsd/mms/medlineelements.html#dp

    Seasons use the metrological start (Winter: Dec 1, Spring: Mar 1,
    Summer: Jun 1, Fall: Sep 1). Y/M/D-D takes the first day, Y/M or Y/M-M the
    first day of that month, Y or Y-Y the first day of that year.
    TODO: Not important since only one instance so far but fix edge case "2016 11-12"
    """
    months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    seasons = {"spring": " mar 1", "summer": " jun 1", "fall": " sep 1", "winter": " dec 1"}

    s_date = date.lower().split()
    date_len = len(s_date)
    # if length is 1 can either be year or year-year
    if date_len == 1:
        return datetime.strptime(s_date[0].split("-")[0], "%Y").date().isoformat()
    # if length is 2 can either be year season or year month or year month-month
    if date_len == 2:
        if s_date[1][:3] in months:
            return datetime.strptime(s_date[0] + " " + s_date[1][:3], "%Y %b").date().isoformat()
        if season := seasons.get(s_date[1]):
            return datetime.strptime(s_date[0] + season, "%Y %b %d").date().isoformat()
        logger.warning("Need to update isoformat transformation: %s", date)
        return None
    # if length is 3 should be year month day or year month day-day or year month-month day
    if date_len == 3:
        year = s_date[0]
        og_month = s_date[1].split("-")[0]
        day = s_date[2].split("-")[0]
        month = next((month for month in months if month in og_month), None)
        if month and day.isdigit():
            return datetime.strptime(year + " " + month + " " + day, "%Y %b %d").date().isoformat()
        if month and day[:3] in months:
            # malformed month-range like "2005 Feb Nov" -- use the starting month
            return datetime.strptime(year + " " + month, "%Y %b").date().isoformat()
        logger.warning("Need to update isoformat transformation: %s", date)
        return None
    # exception case there are quite a few entries with this case "2020 Jan - Feb"
    if date_len == 4 and s_date[1] in months and s_date[3] in months and s_date[2] == "-":
        return datetime.strptime(s_date[0] + " " + s_date[1], "%Y %b").date().isoformat()
    logger.warning("Need to update isoformat transformation: %s", date)
    return None


def _parse_citation(record):
    citation = {}
    if name := record.get("TI"):
        citation["name"] = name
    if pmid := record.get("PMID"):
        citation["pmid"] = pmid
        citation["identifier"] = "PMID:" + pmid
        citation["url"] = "https://pubmed.ncbi.nlm.nih.gov/" + pmid + "/"
    for aid in record.get("AID") or []:
        if aid.endswith(" [doi]"):
            citation["doi"] = aid[: -len(" [doi]")].strip()
            break
    if journal_name := record.get("JT"):
        citation["journalName"] = journal_name
    if date_published := record.get("DP"):
        if date := _get_pub_date(date_published):
            citation["datePublished"] = date

    # make an empty list if there is some kind of author
    if record.get("AU") or record.get("CN"):
        citation["author"] = []
    for author in record.get("AU") or []:
        citation["author"].append({"@type": "Person", "name": author})
    for corp_author in record.get("CN") or []:
        citation["author"].append({"@type": "Organization", "name": corp_author})

    if citation:
        citation["@type"] = "ScholarlyArticle"
        citation["fromPMID"] = True
    return citation


@retry(7, 30)
def batch_get_pmid_eutils(pmids: Iterable[str], email: str, api_key: Optional[str] = None) -> Dict:
    """Fetch citation and funding metadata for `pmids` in one pair of requests."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    ct_fd = {}
    try:
        handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
    except urllib.error.HTTPError as err:
        logger.error("This is the length of the pmids %s", len(pmids))
        logger.error("The list of pmids %s", pmids)
        logger.error("HTTP url: %s", err.url)
        raise err

    for record in Medline.parse(handle):
        ct_fd[record.get("PMID")] = {"citation": _parse_citation(record)}

    # throttle request rates, NCBI says up to 10 requests per second with API Key, 3/s without.
    time.sleep(0.1 if api_key else 0.35)

    # get the funding using the xml file because of problems parsing the medline file
    # https://www.nlm.nih.gov/bsd/mms/medlineelements.html#gr
    handle = Entrez.efetch(db="pubmed", id=pmids, retmode="xml")
    # Have to use Entrez.read() instead of Entrez.parse(): https://github.com/biopython/biopython/issues/1027
    # This can get an incompleteread error, which the retry above covers.
    records = Entrez.read(handle)["PubmedArticle"]

    for record in records:
        funding = []
        for grant in record["MedlineCitation"]["Article"].get("GrantList") or []:
            fund = {"@type": "MonetaryGrant"}
            if grant_id := grant.get("GrantID"):
                fund["identifier"] = str(grant_id)
            if agency := grant.get("Agency"):
                agency = str(agency)
                fund["funder"] = standardize_funder(agency) or {"@type": "Organization", "name": agency}
            if grant.get("Agency") or grant.get("GrantID"):
                fund["fromPMID"] = True
            funding.append(fund)

        pmid = record["MedlineCitation"].get("PMID")
        if pmid and funding:
            ct_fd.setdefault(pmid, {})["funding"] = funding

    return ct_fd


def _filter_valid_pmids(pmid_list):
    """Drop malformed PMIDs, which would otherwise fail the whole NCBI request."""
    valid = []
    invalid = []
    for raw in pmid_list:
        pmid = str(raw).strip().lstrip("0")
        if _VALID_PMID_RE.match(pmid):
            valid.append(pmid)
        else:
            invalid.append(raw)
    if invalid:
        logger.warning("Dropping %d malformed PMID(s) before NCBI request: %s", len(invalid), invalid[:20])
    return valid


def cached_batch_get_pmid_eutils(pmid_list, email, api_key):
    """`batch_get_pmid_eutils` backed by a SQLite cache of previous results."""
    conn = _get_pmid_conn()
    cur = conn.cursor()

    pmid_list = _filter_valid_pmids(pmid_list)
    if not pmid_list:
        return {}

    placeholders = ", ".join("?" for _ in pmid_list)
    cur.execute(f"SELECT pmid, data FROM eutils_cache WHERE pmid IN ({placeholders})", pmid_list)
    cached_results = {row[0]: orjson.loads(row[1]) for row in cur.fetchall()}
    uncached_pmids = [p for p in pmid_list if p not in cached_results]

    if cached_results:
        logger.debug("PMID eutils cache: %d cached, %d to fetch", len(cached_results), len(uncached_pmids))

    fresh_results = {}
    if uncached_pmids:
        fresh_results = batch_get_pmid_eutils(uncached_pmids, email, api_key)
        for pmid, data in fresh_results.items():
            cur.execute(
                "INSERT OR REPLACE INTO eutils_cache (pmid, data) VALUES (?, ?)",
                (str(pmid), orjson.dumps(data).decode("utf-8")),
            )
        conn.commit()

    results = {**cached_results, **fresh_results}
    for info in results.values():
        citation = info.get("citation")
        if isinstance(citation, dict) and citation:
            citation.setdefault("@type", "ScholarlyArticle")
        for grant in info.get("funding") or []:
            if not isinstance(grant, dict):
                continue
            grant.setdefault("@type", "MonetaryGrant")
            for funder in as_list(grant.get("funder")):
                if isinstance(funder, dict):
                    funder.setdefault("@type", "Organization")
    return results


# ---------------------------------------------------------------------------
# The pipeline stage
# ---------------------------------------------------------------------------
def _is_doi_stub(entry, enriched_doi):
    """True when `entry` is the bare {"doi": ...} the parser left for this paper.

    DOIs are case-insensitive, but registrars mint them in mixed case, so the
    stub and the DOI E-utilities returns can differ only in casing.
    """
    return bool(enriched_doi) and isinstance(entry.get("doi"), str) and entry["doi"].lower() == enriched_doi


def _attach_citation(rec, citation):
    # citation entries have to be objects, and the record's own DOI stub for this
    # paper is now redundant -- the enriched citation carries that DOI.
    enriched_doi = (citation.get("doi") or "").lower()
    kept = []
    for entry in as_list(rec.get("citation")):
        if not isinstance(entry, dict):
            logger.warning("Dropping non-object citation %r from %s", entry, rec.get("_id"))
            continue
        if _is_doi_stub(entry, enriched_doi):
            continue
        kept.append(entry)

    kept.append(citation)
    rec["citation"] = kept


def _attach_funding(rec, funding):
    rec["funding"] = as_list(rec.get("funding")) + list(funding)


def _collect_pmids(docs):
    """Resolve every doc's pmcs / citation DOIs into its `pmids` field."""
    pmc_ids = []
    doi_pmid = {}
    for doc in docs:
        if pmcs := doc.get("pmcs"):
            pmc_ids += [pmc.strip() for pmc in pmcs.split(",")]
        for citation in dict_entries(doc, "citation"):
            if doi := citation.get("doi"):
                _convert_doi(doi, doi_pmid)

    pmc_pmid = _convert_pmcs(pmc_ids) if pmc_ids else {}

    pmid_list = set()
    for doc in docs:
        for pmc in [pmc.strip() for pmc in doc.pop("pmcs", "").split(",") if pmc.strip()]:
            if pmid := pmc_pmid.get(pmc):
                doc["pmids"] = doc.get("pmids") + "," + pmid if doc.get("pmids") else pmid
            else:
                logger.debug("There is an issue with this PMCID. PMCID: %s, rec_id: %s", pmc, doc["_id"])
        for citation in dict_entries(doc, "citation"):
            if pmid := doi_pmid.get(citation.get("doi")):
                doc["pmids"] = doc.get("pmids") + "," + pmid if doc.get("pmids") else pmid
        if pmids := doc.get("pmids"):
            pmid_list.update(pmid.strip() for pmid in pmids.split(","))

    return sorted(pmid_list)


def add_citations(docs):
    """Add citations, funding, species and diseases from each record's PMIDs.

    `docs` is one batch of records; PMIDs across the whole batch are queried in
    a single request.
    """
    docs = list(docs)
    _refresh_pubtator_dumps()

    pmid_list = _collect_pmids(docs)
    eutils_info = cached_batch_get_pmid_eutils(pmid_list, GEO_EMAIL, GEO_API_KEY) if pmid_list else {}

    for rec in docs:
        if pmids := rec.pop("pmids", None):
            # fixes issue where pmid numbers under 10 are read as 04 instead of 4
            pmids = sorted({pmid.strip().lstrip("0") for pmid in pmids.split(",")})
            if species := get_data_for_pmids(pmids, "species"):
                update_record_species(rec, species)
            if diseases := get_data_for_pmids(pmids, "disease"):
                update_record_disease(rec, diseases)

            for pmid in pmids:
                info = eutils_info.get(pmid)
                if not info:
                    # this covers records whose pmid does not resolve, e.g.
                    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41964
                    logger.debug("There is an issue with this pmid. PMID: %s, rec_id: %s", pmid, rec["_id"])
                    continue
                if citation := info.get("citation"):
                    _attach_citation(rec, citation)
                if funding := info.get("funding"):
                    _attach_funding(rec, funding)
        yield rec


def standardize_fields(docs, batch_size=1000):
    """Replace ScholarlyArticle PMID stubs with full citations (used by DDE).

    Applies to `isBasedOn`, `isBasisFor`, `citedBy`, `isPartOf` and `hasPart`.
    """
    for batch in batched(docs, batch_size):
        yield from _standardize_fields_batch(batch)


def _standardize_fields_batch(docs):
    fields = ["isBasedOn", "isBasisFor", "citedBy", "isPartOf", "hasPart"]

    pmid_list = set()
    for doc in docs:
        for field in fields:
            for entry in dict_entries(doc, field):
                if entry.get("@type") == "ScholarlyArticle" and entry.get("pmid"):
                    pmid_list.add(str(entry["pmid"]).strip())

    eutils_info = cached_batch_get_pmid_eutils(sorted(pmid_list), GEO_EMAIL, GEO_API_KEY) if pmid_list else {}

    for doc in docs:
        for field in fields:
            entries = as_list(doc.get(field))
            if not entries:
                continue
            field_pmids = [
                str(entry.get("pmid"))
                for entry in entries
                if isinstance(entry, dict) and entry.get("pmid") and entry.get("@type") == "ScholarlyArticle"
            ]
            doc[field] = [
                entry
                for entry in entries
                if not (isinstance(entry, dict) and entry.get("pmid") and entry.get("@type") == "ScholarlyArticle")
            ]
            for pmid in field_pmids:
                if citation := (eutils_info.get(pmid) or {}).get("citation"):
                    citation["@type"] = "ScholarlyArticle"
                    doc[field].append(citation)
        yield doc
