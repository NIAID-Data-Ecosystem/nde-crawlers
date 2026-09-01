"""Species and health conditions mined from a record's description.

Runs for resource records and Samples whose `additionalType` is `BioSample`
when they have a description but are missing taxonomy or health conditions.
Descriptions go to the EXTRACT tagger (tagger.jensenlab.org), whose response is
cached in SQLite per record id, so only records seen for the first time cost a
request. Request failures are cached briefly and repeated service failures open
a worker-local circuit breaker. The tagged names are then validated against
EXTRACT's candidate NCBI taxon ids and standardized with UniProt taxonomy.
Successful taxa and deterministic rejections are cached, and accepted terms
are marked `fromEXTRACT`.
"""

import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from email.utils import parsedate_to_datetime
from itertools import batched

import requests
from config import logger

from .cache import SqliteCache, SqliteKeySet
from .common import as_list, sqlite, supports_description_enrichment
from .taxonomy import classify_from_lineage
from .term_matching import (
    compact_term_text,
    is_ambiguous_short_mention,
    mentioned_in,
    species_term_matches_mention,
    term_expansion_mentioned,
    term_labels,
    term_matches_mention,
)
from .terms import DB_PATH as PUBTATOR_DB_PATH, SPECIES_CACHE_DB_PATH as DB_PATH, fetch_taxon, normalize_taxon_id, query_condition

_NEGATIVE_DISEASE_TABLE = "health_conditions_negative"
_NEGATIVE_SPECIES_TABLE = "description_species_negative"
_TAXON_DETAILS_TABLE = "extract_taxon_details"
_TAXON_FAILURE_TABLE = "extract_taxon_failures"

# Include a version in rejection keys so a future matching-rule change can
# invalidate old rejections without a data migration.
_SPECIES_REJECTION_CACHE_VERSION = "v1"
_DISEASE_IDENTIFIER_CACHE_VERSION = "v1"
_UNIPROT_WORKERS = 10

# These acronym/candidate pairs were confirmed to be meaning collisions in
# Figshare. Keep a match only when the resolved label is also written out in
# the title or description. Pairing the acronym with its identifier avoids
# suppressing the same acronym if EXTRACT selects a different, correct term.
_CONTEXT_REQUIRED_TAXA = frozenset(
    {
        ("cb1", "45142"),  # cannabinoid receptor -> Cucurbitaria berberidis
        ("cmv", "12305"),  # cytomegalovirus -> Cucumber mosaic virus
        ("sag", "59303"),  # sagopilone -> Sagiyama virus
        ("scv", "135656"),  # small colony variant -> Strawberry crinkle virus
        ("tsv", "12317"),  # trichodysplasia polyomavirus -> Tobacco streak virus
    }
)
_CONTEXT_REQUIRED_DISEASES = frozenset(
    {
        ("fdh", "DOID:2120"),  # FIDDLEHEAD/formate dehydrogenase -> focal dermal hypoplasia
        ("hnpp", "DOID:0060843"),  # phosphoester notation -> hereditary neuropathy
    }
)

# EXTRACT entity type codes.
_SPECIES_TYPE = "-2"
_DISEASE_TYPE = "-26"
_CACHE_TABLES = {_SPECIES_TYPE: "species", _DISEASE_TYPE: "disease"}

# SQLite allows 999 bound variables by default; keep headroom.
_SQL_CHUNK_SIZE = 900
_EXTRACT_CHUNK_SIZE = 5000
_EXTRACT_URL = "https://tagger.jensenlab.org/GetEntities"
_EXTRACT_TIMEOUT = (5, 60)
_EXTRACT_MAX_ATTEMPTS = 4
_EXTRACT_RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
_EXTRACT_FAILURE_TTL_SECONDS = 60 * 60
_EXTRACT_CIRCUIT_FAILURE_THRESHOLD = 3
_EXTRACT_CIRCUIT_COOLDOWN_SECONDS = 60 * 60


class ExtractCircuitOpen(RuntimeError):
    """Raised when EXTRACT requests are paused after repeated service failures."""


_extract_circuit_failures = 0
_extract_circuit_open_until = 0.0

# Reuse HTTP connections for EXTRACT.
_EXTRACT_SESSION = requests.Session()

# Place names and terms EXTRACT reliably gets wrong.
BASIC_DROP_LIST = frozenset(
    {
        "tonga",
        "alabama",
        "argentina",
        "namibia",
        "panama",
        "virginia",
        "bulgaria",
        "arge",
        "togo",
        "serendip",
        "vector",
        "metagenome",
        "gut microbiome",
        "arizona",
        "california",
        "omicron",
        "sonoma",
    }
)

# Drop rules with NCBI Taxon IDs. `ignore_children` also drops anything whose
# lineage contains the id.
ADVANCED_DROP_RULES = {
    "other sequences": {
        "id": "28384",
        "ignore_children": True,
        "rationale": "EXTRACT does not do this well. PubTator seems to do better",
    },
    "collection": {
        "id": "1768868",
        "ignore_children": False,
        "rationale": "too easy for EXTRACT to get this wrong",
    },
    "omicron": {
        "id": "2613138",
        "ignore_children": True,
        "rationale": "COVID-19 confusion, EXTRACT will get it wrong",
    },
    "sonoma": {
        "id": "1535511",
        "ignore_children": False,
        "rationale": "Likelihood it's a place rather than organism is very high",
    },
    "china": {
        "id": "3034371",
        "ignore_children": False,
        "rationale": "Likelihood it's a place rather than organism is very high",
    },
    "nevada": {
        "id": "359889",
        "ignore_children": False,
        "rationale": "Likelihood it's a place rather than organism is very high",
    },
    "montana": {
        "id": "441235",
        "ignore_children": False,
        "rationale": "Likelihood it's a place rather than organism is very high",
    },
    "vector": {
        "id": "2971083",
        "ignore_children": False,
        "rationale": "generic study language rather than an organism",
    },
    "arge": {
        "id": "95269",
        "ignore_children": False,
        "rationale": "HTML-split substring of 'Large-scale' rather than an organism",
    },
    "metagenome": {
        "id": "256318",
        "ignore_children": False,
        "rationale": "environmental sequence collection rather than an organism",
    },
    "gut microbiome": {
        "id": "749906",
        "ignore_children": False,
        "rationale": "environmental sample rather than an infectious agent",
    },
    "age strata": {
        "id": "1208515",
        "ignore_children": False,
        "rationale": "demographic phrase matched to the Agestrata taxon",
    },
    "matara": {
        "id": "2612617",
        "ignore_children": False,
        "rationale": "place name matched to an unrelated taxon",
    },
    "kerala": {
        "id": "2249684",
        "ignore_children": False,
        "rationale": "place name matched to an unrelated taxon",
    },
    "transformation": {
        "id": "2839062",
        "ignore_children": False,
        "rationale": "generic biomedical process matched to an unrelated taxon",
    },
    "scleroderma": {
        "id": "68787",
        "ignore_children": False,
        "rationale": "disease mention matched to the Scleroderma fungus taxon",
    },
    "syncope": {
        "id": "1271638",
        "ignore_children": False,
        "rationale": "clinical condition matched to an unrelated taxon",
    },
    "latina": {
        "id": "1325907",
        "ignore_children": False,
        "rationale": "demographic term matched to an unrelated taxon",
    },
    "human microbiome": {
        "id": "646099",
        "ignore_children": False,
        "rationale": "microbiome study phrase rather than an infectious agent",
    },
    "venus": {
        "id": "55714",
        "ignore_children": False,
        "rationale": "fluorescent reporter name matched to the Venus taxon",
    },
    "napo": {
        "id": "706958",
        "ignore_children": False,
        "rationale": "river/place name matched to an unrelated taxon",
    },
}

ONTOLOGY_PRIORITY = {"MONDO": 0, "HPO": 1, "DOID": 2, "NCIT": 3}

_RESPONSE_CACHE_DDL = (
    "CREATE TABLE IF NOT EXISTS species (ndeid TEXT PRIMARY KEY, text_response TEXT)",
    "CREATE TABLE IF NOT EXISTS disease (ndeid TEXT PRIMARY KEY, text_response TEXT)",
)
_FAILURE_CACHE_DDL = """
CREATE TABLE IF NOT EXISTS extract_failures (
    ndeid TEXT NOT NULL,
    entity_type TEXT NOT NULL,
    retry_after REAL NOT NULL,
    error TEXT,
    PRIMARY KEY (ndeid, entity_type)
)
"""
_EXTRACT_CACHE_DDL = (*_RESPONSE_CACHE_DDL, _FAILURE_CACHE_DDL)

# All but TAXON_DETAILS are preloaded whole; they're read for nearly every term.
SPECIES_DETAILS = SqliteCache(DB_PATH, "species_details", preload=True)
TAXON_DETAILS = SqliteCache(
    DB_PATH,
    _TAXON_DETAILS_TABLE,
    key_column="identifier",
    value_column="standard_dict",
    normalize=False,
)
FAILED_TAXA = SqliteKeySet(DB_PATH, _TAXON_FAILURE_TABLE, key_column="identifier", preload=True, normalize=False)
REJECTED_SPECIES = SqliteKeySet(
    DB_PATH,
    _NEGATIVE_SPECIES_TABLE,
    key_column="lookup_key",
    preload=True,
    normalize=False,
)
HEALTH_CONDITIONS = SqliteCache(PUBTATOR_DB_PATH, "health_conditions", preload=True)
NEGATIVE_DISEASES = SqliteKeySet(PUBTATOR_DB_PATH, _NEGATIVE_DISEASE_TABLE, preload=True)


def reset_caches():
    """Drop the process-wide caches so a new upload sees fresh lookup data."""
    SPECIES_DETAILS.reset()
    TAXON_DETAILS.reset()
    FAILED_TAXA.reset()
    REJECTED_SPECIES.reset()
    HEALTH_CONDITIONS.reset()
    NEGATIVE_DISEASES.reset()


# ---------------------------------------------------------------------------
# Drop rules
# ---------------------------------------------------------------------------
def _matches_drop_rule(name, identifier="", lineage_ids=()):
    """Return the drop rule this term violates, if any."""
    name = (name or "").lower()
    identifier = str(identifier or "")
    for term_name, rule in ADVANCED_DROP_RULES.items():
        if name == term_name.lower() or identifier == rule["id"] or identifier.endswith(rule["id"]):
            return term_name, rule
        if rule["ignore_children"] and rule["id"] in lineage_ids:
            return term_name, rule
    return None


def _lineage_ids(details):
    return {str(item.get("taxId", "")) for item in (details or {}).get("lineage", []) if item.get("taxId")}


def should_filter_species_entry(entry, species_mapping=None):
    """True when a species entry matches a drop rule (checking its lineage too)."""
    name = entry.get("name", "").lower()
    lineage_ids = _lineage_ids((species_mapping or {}).get(name))
    return bool(_matches_drop_rule(name, entry.get("identifier", ""), lineage_ids))


def filter_species_by_advanced_rules(species_mapping, species_list):
    """Drop names from `species_list` that match a drop rule."""
    to_remove = set()
    for species in species_list:
        details = species_mapping.get(species.lower())
        if not details:
            continue
        match = _matches_drop_rule(species, details.get("identifier", ""), _lineage_ids(details))
        if match:
            term_name, rule = match
            logger.debug("Filtering species '%s' - drop rule for '%s': %s", species, term_name, rule["rationale"])
            to_remove.add(species)

    if to_remove:
        logger.debug("Filtered out %s species based on advanced drop rules: %s", len(to_remove), to_remove)
    return [species for species in species_list if species not in to_remove]


# ---------------------------------------------------------------------------
# The EXTRACT tagger
# ---------------------------------------------------------------------------
def _entity_types_param(entity_types):
    """Format one or more EXTRACT entity type codes for the API."""
    if isinstance(entity_types, str):
        return entity_types
    return " ".join(entity_types)


def _retry_delay(response, attempt):
    """Honor Retry-After when possible, otherwise use exponential backoff."""
    retry_after = response.headers.get("Retry-After")
    if retry_after:
        try:
            return max(0.0, float(retry_after))
        except ValueError:
            try:
                retry_at = parsedate_to_datetime(retry_after)
                return max(0.0, retry_at.timestamp() - time.time())
            except (TypeError, ValueError, OverflowError):
                pass
    return float(2**attempt)


def _reset_extract_circuit():
    """Reset process-local EXTRACT service failure state."""
    global _extract_circuit_failures, _extract_circuit_open_until
    _extract_circuit_failures = 0
    _extract_circuit_open_until = 0.0


def _extract_circuit_is_open():
    """True while this worker is cooling down after repeated EXTRACT failures."""
    global _extract_circuit_failures, _extract_circuit_open_until
    if not _extract_circuit_open_until:
        return False
    if time.monotonic() < _extract_circuit_open_until:
        return True

    logger.info("EXTRACT circuit-breaker cooldown ended; requests will resume")
    _extract_circuit_failures = 0
    _extract_circuit_open_until = 0.0
    return False


def _record_extract_service_failure(reason, open_immediately=False):
    """Open the process-local circuit after repeated service-level failures."""
    global _extract_circuit_failures, _extract_circuit_open_until
    _extract_circuit_failures += 1
    if not open_immediately and _extract_circuit_failures < _EXTRACT_CIRCUIT_FAILURE_THRESHOLD:
        return

    _extract_circuit_open_until = time.monotonic() + _EXTRACT_CIRCUIT_COOLDOWN_SECONDS
    logger.warning(
        "EXTRACT circuit breaker opened for %.0fs after %s consecutive service failures; last failure: %s",
        _EXTRACT_CIRCUIT_COOLDOWN_SECONDS,
        _extract_circuit_failures,
        reason,
    )


def _extract_request(session, params):
    """Keep descriptions in the request body and avoid URL-length limits."""
    return session.post(_EXTRACT_URL, data=params, timeout=_EXTRACT_TIMEOUT)


def _extract_error_summary(error):
    """Describe a request failure without copying the record text into logs or SQLite."""
    if isinstance(error, requests.HTTPError) and error.response is not None:
        return f"HTTP {error.response.status_code}"
    if isinstance(error, requests.RequestException):
        return type(error).__name__
    return f"{type(error).__name__}: {error}"[:500]


def query_extract_api(description, entity_types, session=None):
    """Query EXTRACT for one or more entity types, retrying transient failures."""
    if _extract_circuit_is_open():
        raise ExtractCircuitOpen("EXTRACT circuit breaker is open")

    entity_types = _entity_types_param(entity_types)
    session = session or _EXTRACT_SESSION
    params = {"document": description, "entity_types": entity_types, "format": "tsv"}

    for attempt in range(_EXTRACT_MAX_ATTEMPTS):
        try:
            response = _extract_request(session, params)
        except (requests.ConnectionError, requests.Timeout) as e:
            if attempt + 1 == _EXTRACT_MAX_ATTEMPTS:
                _record_extract_service_failure(e)
                raise
            delay = float(2**attempt)
            logger.warning(
                "EXTRACT request failed for entity types %s (%s); retrying in %.1fs",
                entity_types,
                e,
                delay,
            )
            time.sleep(delay)
            continue

        if response.status_code in _EXTRACT_RETRY_STATUSES:
            if attempt + 1 == _EXTRACT_MAX_ATTEMPTS:
                _record_extract_service_failure(f"HTTP {response.status_code}")
                response.raise_for_status()
            delay = _retry_delay(response, attempt)
            logger.warning(
                "EXTRACT returned HTTP %s for entity types %s; retrying in %.1fs",
                response.status_code,
                entity_types,
                delay,
            )
            time.sleep(delay)
            continue

        try:
            response.raise_for_status()
        except requests.HTTPError:
            if response.status_code == 403:
                _record_extract_service_failure("HTTP 403", open_immediately=True)
            else:
                _reset_extract_circuit()
            raise

        _reset_extract_circuit()
        return response.text

    raise RuntimeError("EXTRACT retry loop ended unexpectedly")


def _iter_extract_tsv_lines(text_response):
    """Yield (extracted_text, entity_type, onto_id) from an EXTRACT TSV response."""
    for line in (text_response or "").splitlines():
        if not line:
            continue
        # We only care about the first 3 columns.
        parts = line.split("\t", 2)
        if len(parts) < 3:
            continue
        yield parts[0], parts[1], parts[2]


def _response_for_entity_type(text_response, entity_type):
    """Return only the TSV rows for one type from a combined EXTRACT response."""
    lines = []
    for line in (text_response or "").splitlines():
        parts = line.split("\t", 2)
        if len(parts) >= 3 and parts[1] == entity_type:
            lines.append(line)
    return "\n".join(lines)


def _fetch_cached_responses(cursor, entity_type, ndeids):
    """Fetch cached EXTRACT responses for many record ids: {ndeid: text_response}."""
    table = _CACHE_TABLES.get(entity_type)
    if not table or not ndeids:
        return {}

    out = {}
    for chunk in batched(ndeids, _SQL_CHUNK_SIZE):
        placeholders = ",".join("?" for _ in chunk)
        cursor.execute(f"SELECT ndeid, text_response FROM {table} WHERE ndeid IN ({placeholders})", chunk)
        out.update(cursor.fetchall())
    return out


def _fetch_cached_failures(cursor, entity_type, ndeids, now):
    """Fetch record ids whose latest EXTRACT failure is still cooling down."""
    if not ndeids:
        return set()

    out = set()
    for chunk in batched(ndeids, _SQL_CHUNK_SIZE):
        placeholders = ",".join("?" for _ in chunk)
        cursor.execute(
            f"SELECT ndeid FROM extract_failures "
            f"WHERE entity_type = ? AND retry_after > ? AND ndeid IN ({placeholders})",
            (entity_type, now, *chunk),
        )
        out.update(row[0] for row in cursor.fetchall())
    return out


def _cache_extract_results(rows_by_type, failed_rows):
    """Write fetched responses and temporary failures after network I/O."""
    if not any(rows_by_type.values()) and not failed_rows:
        return

    with sqlite(DB_PATH, *_EXTRACT_CACHE_DDL) as conn:
        if failed_rows:
            conn.executemany(
                "INSERT OR REPLACE INTO extract_failures "
                "(ndeid, entity_type, retry_after, error) VALUES (?, ?, ?, ?)",
                failed_rows,
            )
        for entity_type, rows in rows_by_type.items():
            table = _CACHE_TABLES.get(entity_type)
            if table and rows:
                conn.executemany(f"INSERT OR REPLACE INTO {table} VALUES (?, ?)", rows)
                conn.executemany(
                    "DELETE FROM extract_failures WHERE ndeid = ? AND entity_type = ?",
                    ((ndeid, entity_type) for ndeid, _ in rows),
                )


def _tagged_entities(response_text, entity_type):
    """Return the entity names of `entity_type` from an EXTRACT response."""
    if not response_text:
        return []

    return [
        (extracted_text, onto_id)
        for extracted_text, tagged_type, onto_id in _iter_extract_tsv_lines(response_text)
        if tagged_type == entity_type
        and extracted_text.lower() not in BASIC_DROP_LIST
        and not _matches_drop_rule(extracted_text, onto_id)
    ]


def _extract_entities(doc_list):
    """Add `fromEXTRACT` species / healthCondition stubs from each description."""
    started = time.monotonic()
    count = 0
    requested = {_SPECIES_TYPE: 0, _DISEASE_TYPE: 0}
    cache_hits = 0
    cache_misses = 0
    failure_cache_hits = 0
    api_calls = 0
    combined_calls = 0
    circuit_skips = 0
    api_seconds = 0.0

    for chunk_docs in batched(doc_list, _EXTRACT_CHUNK_SIZE):
        species_ids = []
        disease_ids = []
        for doc in chunk_docs:
            if not doc.get("description"):
                continue
            if "species" not in doc and "infectiousAgent" not in doc:
                species_ids.append(doc["_id"].lower())
            if "healthCondition" not in doc:
                disease_ids.append(doc["_id"].lower())

        with sqlite(DB_PATH, *_EXTRACT_CACHE_DDL) as conn:
            c = conn.cursor()
            cached_species = _fetch_cached_responses(c, _SPECIES_TYPE, species_ids)
            cached_disease = _fetch_cached_responses(c, _DISEASE_TYPE, disease_ids)
            now = time.time()
            failed_species = _fetch_cached_failures(c, _SPECIES_TYPE, species_ids, now)
            failed_disease = _fetch_cached_failures(c, _DISEASE_TYPE, disease_ids, now)

        cached_by_type = {_SPECIES_TYPE: cached_species, _DISEASE_TYPE: cached_disease}
        failed_by_type = {_SPECIES_TYPE: failed_species, _DISEASE_TYPE: failed_disease}
        pending_writes = {_SPECIES_TYPE: [], _DISEASE_TYPE: []}
        pending_failures = []

        for doc in chunk_docs:
            count += 1
            if count % 1000 == 0:
                logger.info("EXTRACT: processed %s documents", count)

            description = doc.get("description")
            if not description:
                continue
            ndeid = doc["_id"].lower()
            needed_types = []
            if "species" not in doc and "infectiousAgent" not in doc:
                needed_types.append(_SPECIES_TYPE)
            if "healthCondition" not in doc:
                needed_types.append(_DISEASE_TYPE)

            responses = {}
            missing_types = []
            for entity_type in needed_types:
                requested[entity_type] += 1
                response_text = cached_by_type[entity_type].get(ndeid)
                if response_text is None:
                    cache_misses += 1
                    if ndeid in failed_by_type[entity_type]:
                        failure_cache_hits += 1
                    else:
                        missing_types.append(entity_type)
                else:
                    cache_hits += 1
                    responses[entity_type] = response_text

            if missing_types:
                if _extract_circuit_is_open():
                    circuit_skips += 1
                else:
                    api_calls += 1
                    if len(missing_types) > 1:
                        combined_calls += 1
                    api_started = time.monotonic()
                    try:
                        fetched_response = query_extract_api(description, missing_types)
                    except ExtractCircuitOpen:
                        circuit_skips += 1
                    except Exception as e:
                        error = _extract_error_summary(e)
                        logger.error("Error querying EXTRACT for document %s: %s", doc.get("_id"), error)
                        retry_after = time.time() + _EXTRACT_FAILURE_TTL_SECONDS
                        pending_failures.extend(
                            (ndeid, entity_type, retry_after, error) for entity_type in missing_types
                        )
                    else:
                        for entity_type in missing_types:
                            # A combined request is split before caching, preserving
                            # the existing one-entity-type-per-table cache contents.
                            response_text = (
                                _response_for_entity_type(fetched_response, entity_type)
                                if len(missing_types) > 1
                                else fetched_response
                            )
                            responses[entity_type] = response_text
                            cached_by_type[entity_type][ndeid] = response_text
                            pending_writes[entity_type].append((ndeid, response_text))
                    finally:
                        api_seconds += time.monotonic() - api_started

            try:
                if _SPECIES_TYPE in responses:
                    for name, onto_id in _tagged_entities(responses[_SPECIES_TYPE], _SPECIES_TYPE):
                        if not mentioned_in(name, (description,)):
                            logger.debug("Skipping EXTRACT species substring %r in %s", name, doc.get("_id"))
                            continue
                        species = doc.setdefault("species", [])
                        # EXTRACT can return several taxonomy candidates for one
                        # mention. Keep each candidate until UniProt labels let
                        # us select the one that actually agrees with the text.
                        if not any(
                            entry.get("name") == name and str(entry.get("identifier")) == str(onto_id)
                            for entry in species
                        ):
                            species.append(
                                {
                                    "@type": "DefinedTerm",
                                    "name": name,
                                    "identifier": onto_id,
                                    "fromEXTRACT": True,
                                }
                            )

                if _DISEASE_TYPE in responses:
                    for name, onto_id in _tagged_entities(responses[_DISEASE_TYPE], _DISEASE_TYPE):
                        if not mentioned_in(name, (description,)):
                            logger.debug("Skipping EXTRACT disease substring %r in %s", name, doc.get("_id"))
                            continue
                        conditions = doc.setdefault("healthCondition", [])
                        if not any(
                            entry.get("name") == name and str(entry.get("identifier")) == str(onto_id)
                            for entry in conditions
                        ):
                            conditions.append(
                                {
                                    "@type": "DefinedTerm",
                                    "name": name,
                                    "identifier": onto_id,
                                    "fromEXTRACT": True,
                                }
                            )
            except Exception as e:
                logger.error("Error processing EXTRACT response for document %s: %s", doc.get("_id"), e)

        _cache_extract_results(pending_writes, pending_failures)

    logger.info(
        "EXTRACT: docs=%s species_requested=%s disease_requested=%s cache_hits=%s "
        "cache_misses=%s failure_cache_hits=%s api_calls=%s combined_calls=%s "
        "circuit_skips=%s api_seconds=%.1fs total_seconds=%.1fs",
        count,
        requested[_SPECIES_TYPE],
        requested[_DISEASE_TYPE],
        cache_hits,
        cache_misses,
        failure_cache_hits,
        api_calls,
        combined_calls,
        circuit_skips,
        api_seconds,
        time.monotonic() - started,
    )
    return doc_list


def _already_named(entries, name):
    return any(entry.get("name") == name or name in entry.get("alternateName", []) for entry in entries)


# ---------------------------------------------------------------------------
# UniProt taxonomy for extracted species
# ---------------------------------------------------------------------------
def get_species_details(original_name, identifier):
    """Standardize an extracted species name from its UniProt taxonomy entry."""
    term = fetch_taxon(original_name, identifier, classify=classify_from_lineage)
    # Nothing to classify from means we cannot call it a host.
    term.setdefault("classification", "infectiousAgent")
    term["isCurated"] = False
    term["fromEXTRACT"] = True
    return term


def _species_candidate_matches_mention(term, mention):
    """True when a taxon label safely supports an EXTRACT mention.

    Short acronyms are especially collision-prone. Require their original
    casing to occur in an authoritative UniProt label, while retaining the
    shared allowlist for well-known biomedical acronyms such as HCV and H1N1.
    """
    if not species_term_matches_mention(term, mention):
        return False
    if not is_ambiguous_short_mention(mention):
        return True

    normalized_mention = " ".join(str(mention).split())
    return any(" ".join(str(label).split()) == normalized_mention for label in term_labels(term))


# ---------------------------------------------------------------------------
# Standardizing extracted species
# ---------------------------------------------------------------------------
def _normalize_term_entries(doc, field):
    """Coerce a term field into a list of dicts."""
    if field not in doc:
        return
    entries = as_list(doc[field])
    doc[field] = [{"@type": "DefinedTerm", "name": entry} if isinstance(entry, str) else entry for entry in entries]


def _build_species_lineage_info(species_list, species_mapping):
    """Map each species' scientific name to the set of its ancestors."""
    lineage_info = {}
    for species in species_list:
        details = species_mapping.get(species.lower())
        if details and details.get("lineage"):
            lineage_info[details["name"].lower()] = {item["scientificName"].lower() for item in details["lineage"]}
    return lineage_info


def _filter_species_terms_for_ancestors(species_mapping, species_list, lineage_info):
    """Drop species that are an ancestor of another species in the same record."""
    to_remove = set()
    if "mus" in species_list and "mus sp." in species_list:
        to_remove.add("mus sp.")

    for species in species_list:
        details = species_mapping.get(species.lower())
        if not details:
            continue
        scientific_name = details["name"].lower()
        for other_species in species_list:
            other_details = species_mapping.get(other_species.lower())
            if not other_details:
                continue
            if other_species != species and scientific_name in lineage_info.get(other_details["name"].lower(), set()):
                to_remove.add(scientific_name)

    return [
        species
        for species in species_list
        if species_mapping.get(species.lower()) and species_mapping[species.lower()]["name"].lower() not in to_remove
    ]


def _insert_species(doc_list, species_mapping):
    """Replace extracted species stubs with their standardized terms."""
    for doc in doc_list:
        combined_entries = as_list(doc.get("species")) + as_list(doc.get("infectiousAgent"))
        if not combined_entries:
            continue

        species_names = [entry.get("name") for entry in combined_entries if "name" in entry]
        lineage_info = _build_species_lineage_info(species_names, species_mapping)
        keep_names = _filter_species_terms_for_ancestors(species_mapping, species_names, lineage_info)
        keep_names = filter_species_by_advanced_rules(species_mapping, keep_names)

        hosts = []
        infectious_agents = []
        seen_host_names = set()
        seen_agent_names = set()
        seen_identifiers = set()

        for entry in combined_entries:
            name = entry.get("name")
            if not name:
                continue
            # Uncurated names the filters rejected are dropped entirely.
            if not entry.get("curatedBy") and name not in keep_names:
                continue

            lower_name = name.lower()
            if not entry.get("fromEXTRACT", False):
                # Curated entries keep whatever classification they came with.
                classification = entry.get("classification")
                if classification == "host":
                    hosts.append(entry)
                    seen_host_names.add(lower_name)
                elif classification == "infectiousAgent":
                    infectious_agents.append(entry)
                    seen_agent_names.add(lower_name)
                seen_identifiers.add(entry.get("identifier"))
                continue

            new_obj = species_mapping.get(lower_name)
            if not new_obj:
                continue
            needs_context = (
                compact_term_text(name),
                str(new_obj.get("identifier") or ""),
            ) in _CONTEXT_REQUIRED_TAXA
            if needs_context and not term_expansion_mentioned(
                new_obj,
                name,
                (str(doc.get("name") or ""), str(doc.get("description") or "")),
            ):
                logger.debug(
                    "Ignoring unexpanded EXTRACT species acronym %r in %s",
                    name,
                    doc.get("_id"),
                )
                continue
            # `lineage` is only used to filter; it never belongs in the record.
            new_obj = {key: value for key, value in new_obj.items() if key != "lineage"}
            identifier = new_obj.get("identifier")
            classification = new_obj.get("classification")
            if classification == "infectiousAgent":
                if lower_name not in seen_agent_names and identifier not in seen_identifiers:
                    infectious_agents.append(new_obj)
                    seen_agent_names.add(lower_name)
                    seen_identifiers.add(identifier)
            elif classification == "host":
                if lower_name not in seen_host_names and identifier not in seen_identifiers:
                    hosts.append(new_obj)
                    seen_host_names.add(lower_name)
                    seen_identifiers.add(identifier)

        if hosts:
            doc["species"] = hosts
        else:
            doc.pop("species", None)
        if infectious_agents:
            doc["infectiousAgent"] = infectious_agents
        else:
            doc.pop("infectiousAgent", None)
    return doc_list


def _standardize_extracted_species(doc_list):
    """Standardize species/infectiousAgent for records whose terms all came from EXTRACT."""
    for doc in doc_list:
        _normalize_term_entries(doc, "species")
        _normalize_term_entries(doc, "infectiousAgent")

    # Only standardize a record if every one of its terms came from EXTRACT;
    # a record with curated terms is left alone.
    docs_to_standardize = []
    for doc in doc_list:
        combined = as_list(doc.get("species")) + as_list(doc.get("infectiousAgent"))
        if combined and all(term.get("fromEXTRACT", False) for term in combined):
            docs_to_standardize.append(doc)

    term_names = {
        term["name"]
        for doc in docs_to_standardize
        for field in ("species", "infectiousAgent")
        for term in as_list(doc.get(field))
        if term.get("fromEXTRACT", False) and term.get("name")
    }
    candidate_identifiers = {}
    for doc in docs_to_standardize:
        for field in ("species", "infectiousAgent"):
            for term in as_list(doc.get(field)):
                if not term.get("fromEXTRACT", False) or not term.get("name") or not term.get("identifier"):
                    continue
                candidates = candidate_identifiers.setdefault(term["name"], [])
                identifier = str(term["identifier"])
                if identifier not in candidates:
                    candidates.append(identifier)
    logger.info(
        "Species standardization: total_docs=%s docs_to_standardize=%s unique_terms=%s",
        len(doc_list),
        len(docs_to_standardize),
        len(term_names),
    )
    if not term_names:
        return doc_list

    formatted_species = []
    try:
        missing_terms = []
        species_cache_hits = 0
        incompatible_cache_hits = 0
        for original_name in term_names:
            standardized = SPECIES_DETAILS.get(original_name)
            if standardized and _species_candidate_matches_mention(standardized, original_name):
                formatted_species.append(dict(standardized, fromEXTRACT=True, originalName=original_name))
                species_cache_hits += 1
            else:
                if standardized:
                    incompatible_cache_hits += 1
                    logger.debug(
                        "Ignoring incompatible cached EXTRACT species mapping for %s: %s",
                        original_name,
                        standardized.get("name"),
                    )
                missing_terms.append(original_name)

        logger.info(
            "Species standardization cache: hits=%s incompatible_hits=%s misses=%s",
            species_cache_hits,
            incompatible_cache_hits,
            len(missing_terms) - incompatible_cache_hits,
        )
        if missing_terms:
            formatted_species.extend(_resolve_missing_species(missing_terms, candidate_identifiers))
    except Exception as e:
        logger.error("Error during species standardization: %s", e)
        return doc_list

    species_mapping = {(sp.get("originalName") or sp["name"]).lower(): sp for sp in formatted_species}
    # Runs even with no candidates: this is also what removes unstandardizable EXTRACT stubs
    _insert_species(docs_to_standardize, species_mapping)
    return doc_list


def _candidate_taxon_ids(identifiers):
    """Return unique, normalized NCBI taxon ids in their original order."""
    normalized = []
    seen = set()
    for identifier in identifiers:
        taxon_id = normalize_taxon_id(identifier)
        if taxon_id and taxon_id not in seen:
            normalized.append(taxon_id)
            seen.add(taxon_id)
    return normalized


def _species_rejection_key(original_name, identifiers):
    """Key a rejection by mention and candidate set, preserving meaningful case."""
    mention = " ".join(str(original_name or "").split())
    candidates = sorted(
        {
            normalize_taxon_id(identifier) or f"raw:{str(identifier).strip()}"
            for identifier in identifiers
            if str(identifier).strip()
        }
    )
    return json.dumps(
        [_SPECIES_REJECTION_CACHE_VERSION, mention, candidates],
        ensure_ascii=False,
        separators=(",", ":"),
    )


def _cacheable_taxon_details(details):
    """Strip mention-specific provenance before caching a taxon by identifier."""
    return {key: value for key, value in details.items() if key not in {"originalName", "fromEXTRACT"}}


def _taxon_details_for_mention(details, original_name):
    """Stamp cached identifier-level details for one EXTRACT mention."""
    candidate = dict(details)
    candidate["originalName"] = original_name
    candidate["isCurated"] = False
    candidate["fromEXTRACT"] = True
    return candidate


def _permanent_taxon_failure(error):
    """True when retrying the same identifier cannot reasonably help."""
    if isinstance(error, ValueError):
        return True
    if isinstance(error, requests.HTTPError) and error.response is not None:
        return error.response.status_code in {400, 404, 410}
    return False


def _fetch_candidate_taxon(identifier):
    """Fetch identifier-level taxonomy details for concurrent resolution."""
    try:
        details = get_species_details(identifier, identifier)
    except Exception as error:
        return identifier, None, _permanent_taxon_failure(error), error
    if not details:
        return identifier, None, True, ValueError("empty UniProt taxonomy response")
    return identifier, _cacheable_taxon_details(details), False, None


def _load_candidate_taxa(identifiers):
    """Load candidate taxa by id, fetching uncached ids concurrently.

    Returns ``(details_by_id, transient_failure_ids)``. Permanent failures are
    cached separately, while transient failures deliberately remain retryable.
    """
    identifiers = list(dict.fromkeys(identifiers))
    if not identifiers:
        return {}, set()

    started = time.monotonic()
    failed_ids = FAILED_TAXA.known(identifiers)
    details_by_id = TAXON_DETAILS.get_many(set(identifiers) - failed_ids)
    missing = [
        identifier for identifier in identifiers if identifier not in failed_ids and identifier not in details_by_id
    ]
    cache_hits = len(details_by_id)
    negative_hits = len(failed_ids)

    permanent_failures = []
    transient_failures = set()
    fetched_details = {}
    if missing:
        with ThreadPoolExecutor(max_workers=min(_UNIPROT_WORKERS, len(missing))) as executor:
            futures = {executor.submit(_fetch_candidate_taxon, identifier): identifier for identifier in missing}
            for future in as_completed(futures):
                identifier, details, permanent, error = future.result()
                if details:
                    details_by_id[identifier] = details
                    fetched_details[identifier] = details
                elif permanent:
                    permanent_failures.append(identifier)
                    logger.debug("Permanent UniProt taxonomy failure for %s: %s", identifier, error)
                else:
                    transient_failures.add(identifier)
                    logger.warning("Transient UniProt taxonomy failure for %s: %s", identifier, error)

    TAXON_DETAILS.put_many(fetched_details)
    FAILED_TAXA.add_many(permanent_failures)

    logger.info(
        "Species candidate taxonomy: ids=%s cache_hits=%s negative_hits=%s api_calls=%s "
        "permanent_failures=%s transient_failures=%s total_seconds=%.1fs",
        len(identifiers),
        cache_hits,
        negative_hits,
        len(missing),
        len(permanent_failures),
        len(transient_failures),
        time.monotonic() - started,
    )
    return details_by_id, transient_failures


def _resolve_missing_species(missing_terms, candidate_identifiers=None):
    """Resolve EXTRACT mentions from its candidate taxon ids.

    Description enrichment intentionally does not fall back to text2term. In a
    six-hour production sample, rebuilding text2term's NCBITaxon TF-IDF index
    Deliberately no text2term fallback: in a production sample it cost ~2.9 hours
    of index rebuilds for 3 safe mappings out of 6,686 terms. EXTRACT already
    supplies candidate taxon ids.
    """
    candidate_identifiers = candidate_identifiers or {}
    resolved = []
    pending = {}
    rejection_cache_hits = 0

    for original_name in missing_terms:
        raw_identifiers = candidate_identifiers.get(original_name, ())
        rejection_key = _species_rejection_key(original_name, raw_identifiers)
        if rejection_key in REJECTED_SPECIES:
            rejection_cache_hits += 1
            continue
        pending[original_name] = (_candidate_taxon_ids(raw_identifiers), rejection_key)

    candidate_ids = [identifier for identifiers, _ in pending.values() for identifier in identifiers]
    details_by_id, transient_failures = _load_candidate_taxa(candidate_ids)

    rejected_keys = []
    resolved_cache = {}
    for original_name, (identifiers, rejection_key) in pending.items():
        species_details = None
        for identifier in identifiers:
            details = details_by_id.get(identifier)
            if not details:
                continue
            candidate = _taxon_details_for_mention(details, original_name)
            if _species_candidate_matches_mention(candidate, original_name):
                species_details = candidate
                break

        if species_details:
            resolved.append(species_details)
            resolved_cache[original_name] = {
                key: value for key, value in species_details.items() if key != "fromEXTRACT"
            }
            continue

        # Do not turn a service outage into a permanent semantic rejection.
        if not any(identifier in transient_failures for identifier in identifiers):
            rejected_keys.append(rejection_key)

    SPECIES_DETAILS.put_many(resolved_cache)
    REJECTED_SPECIES.add_many(rejected_keys)

    failed = len(missing_terms) - len(resolved)
    logger.info(
        "Species standardization: resolved=%s failed=%s rejection_cache_hits=%s new_rejections=%s",
        len(resolved),
        failed,
        rejection_cache_hits,
        len(rejected_keys),
    )
    return resolved


def _dedupe_species(doc_list):
    """Deduplicate species by identifier and apply the drop rules once more."""
    for doc in doc_list:
        for field in ("species", "infectiousAgent"):
            if field not in doc:
                continue
            seen = set()
            kept = []
            for entry in doc[field]:
                # Description-specific deny rules must never remove a curated
                # or source-provided taxon that merely shares the same label.
                # Filter before recording the identifier so a rejected EXTRACT
                # entry cannot hide a later source-provided duplicate.
                if entry.get("fromEXTRACT", False) and should_filter_species_entry(entry):
                    continue
                identifier = entry.get("identifier")
                if identifier is not None:
                    if identifier in seen:
                        continue
                    seen.add(identifier)
                kept.append(entry)
            if kept:
                doc[field] = kept
            else:
                doc.pop(field, None)
    return doc_list


def _remove_redundant_species(doc_list):
    """Drop species that duplicate a curated infectiousAgent in the same record."""
    for doc in doc_list:
        if "infectiousAgent" not in doc or "species" not in doc:
            continue
        curated_names = {
            agent["name"].strip().lower()
            for agent in doc["infectiousAgent"]
            if agent.get("isCurated", False) and agent.get("name")
        }
        if not curated_names:
            continue

        kept = [sp for sp in doc["species"] if sp.get("name", "").strip().lower() not in curated_names]
        if kept:
            doc["species"] = kept
        else:
            # Every species duplicated a curated agent. Drop the field rather than
            # leave an empty list for later stages to trip over.
            doc.pop("species", None)
    return doc_list


# ---------------------------------------------------------------------------
# Standardizing extracted health conditions
# ---------------------------------------------------------------------------
def _normalize_disease_ontology_id(identifier):
    """Normalize ontology ids supplied by EXTRACT, ignoring other identifiers."""
    identifier = str(identifier or "").strip()
    if ":" not in identifier:
        return None
    prefix, value = identifier.split(":", 1)
    prefix = prefix.upper()
    if prefix == "HPO":
        prefix = "HP"
    if prefix not in {"MONDO", "HP", "DOID", "NCIT"} or not value:
        return None
    return f"{prefix}:{value}"


def _disease_term_ontology_id(term):
    """Return a standardized term's prefixed ontology identifier."""
    identifier = str(term.get("identifier") or "").strip()
    if ":" in identifier:
        return _normalize_disease_ontology_id(identifier)
    term_set = str(term.get("inDefinedTermSet") or "").strip()
    if not term_set or not identifier:
        return None
    return _normalize_disease_ontology_id(f"{term_set}:{identifier}")


def _disease_lookup_key(name, identifier=None):
    """Key text searches by mention and direct lookups by ontology identifier."""
    name = str(name or "").lower().strip()
    identifier = _normalize_disease_ontology_id(identifier)
    if not identifier:
        return name
    return json.dumps(
        [_DISEASE_IDENTIFIER_CACHE_VERSION, identifier],
        ensure_ascii=False,
        separators=(",", ":"),
    )


def _disease_mapping_key(disease):
    """Return the mention/candidate pair used to map one extracted stub."""
    name = str(disease.get("name") or "").lower().strip()
    identifier = None
    if disease.get("fromEXTRACT"):
        identifier = _normalize_disease_ontology_id(disease.get("identifier"))
    return name, identifier


def _insert_disease(doc_list, disease_mapping):
    """Replace extracted health condition stubs with their standardized terms."""
    for doc in doc_list:
        if "healthCondition" not in doc:
            continue

        # PMID-derived conditions are already standardized; keep them as they are.
        preserved = [disease for disease in doc["healthCondition"] if disease.get("fromPMID", False)]
        preserved_names = {disease["name"].lower() for disease in preserved if disease.get("name")}
        preserved_identifiers = {disease.get("identifier") for disease in preserved if disease.get("identifier")}

        updated = []
        for disease in doc["healthCondition"]:
            original_name = (disease.get("name") or "").lower()
            identifier = disease.get("identifier")
            if original_name in preserved_names or (identifier and identifier in preserved_identifiers):
                continue
            standardized = disease_mapping.get(_disease_mapping_key(disease))
            if not standardized:
                updated.append(disease)
                continue
            needs_context = (
                compact_term_text(disease.get("name")),
                _disease_term_ontology_id(standardized),
            ) in _CONTEXT_REQUIRED_DISEASES
            if needs_context and not term_expansion_mentioned(
                standardized,
                disease.get("name"),
                (str(doc.get("name") or ""), str(doc.get("description") or "")),
            ):
                logger.debug(
                    "Ignoring unexpanded EXTRACT disease acronym %r in %s",
                    disease.get("name"),
                    doc.get("_id"),
                )
                continue
            updated.append(standardized)

        if conditions := preserved + updated:
            doc["healthCondition"] = conditions
        else:
            doc.pop("healthCondition", None)
    return doc_list


def _standardize_extracted_diseases(doc_list):
    """Standardize the health conditions of every record that has uncurated ones."""

    disease_lookups = set()
    for doc in doc_list:
        if "healthCondition" not in doc:
            continue
        _normalize_term_entries(doc, "healthCondition")
        for term in doc["healthCondition"]:
            if "isCurated" in term or not term.get("name"):
                continue
            identifier = None
            if term.get("fromEXTRACT"):
                identifier = _normalize_disease_ontology_id(term.get("identifier"))
            disease_lookups.add((term["name"], identifier))

    disease_mapping = {}
    incompatible_lookups = set()
    for disease_name, candidate_id in disease_lookups:
        disease_key = disease_name.lower().strip()
        mapping_key = (disease_key, candidate_id)
        cache_key = _disease_lookup_key(disease_name, candidate_id)
        if cache_key in NEGATIVE_DISEASES:
            incompatible_lookups.add(mapping_key)
            continue

        if standardized := HEALTH_CONDITIONS.get(cache_key):
            identifier_matches = not candidate_id or _disease_term_ontology_id(standardized) == candidate_id
            mention_matches = term_matches_mention(standardized, disease_name) or bool(
                candidate_id and is_ambiguous_short_mention(disease_name)
            )
            if identifier_matches and mention_matches:
                disease_mapping[mapping_key] = dict(
                    standardized,
                    fromEXTRACT=True,
                    originalName=disease_name,
                )
                continue
            logger.debug(
                "Ignoring incompatible cached EXTRACT disease mapping for %s: %s",
                disease_name,
                standardized.get("name"),
            )

        try:
            if candidate_id:
                disease_details = query_condition(disease_name, ontology_id=candidate_id)
            else:
                disease_details = query_condition(disease_name)
        except Exception as e:
            logger.debug("An error occurred while processing %s: %s", disease_name, e)
            continue

        if not disease_details:
            NEGATIVE_DISEASES.add(cache_key)
            incompatible_lookups.add(mapping_key)
            continue
        identifier_matches = not candidate_id or _disease_term_ontology_id(disease_details) == candidate_id
        mention_matches = term_matches_mention(disease_details, disease_name) or bool(
            candidate_id and is_ambiguous_short_mention(disease_name)
        )
        if not identifier_matches or not mention_matches:
            logger.debug(
                "Ignoring incompatible EXTRACT disease mapping for %s: %s",
                disease_name,
                disease_details.get("name"),
            )
            incompatible_lookups.add(mapping_key)
            continue

        disease_details.setdefault("originalName", disease_name)
        disease_details.pop("curatedBy", None)
        disease_mapping[mapping_key] = dict(disease_details, fromEXTRACT=True)
        HEALTH_CONDITIONS.put(cache_key, disease_details)

    if incompatible_lookups:
        for doc in doc_list:
            if "healthCondition" not in doc:
                continue
            kept = [
                disease
                for disease in doc["healthCondition"]
                if disease.get("isCurated") is not None or _disease_mapping_key(disease) not in incompatible_lookups
            ]
            if kept:
                doc["healthCondition"] = kept
            else:
                doc.pop("healthCondition", None)

    if not disease_mapping:
        return doc_list
    return _insert_disease(doc_list, disease_mapping)


def _dedupe_diseases(doc_list):
    """Unify duplicate health conditions, preferring PMID-derived terms then ontology rank."""
    for doc in doc_list:
        if "healthCondition" not in doc:
            continue

        unique_map = {}
        for disease in doc["healthCondition"]:
            key = disease.get("identifier") or (disease.get("name") or "").lower()
            existing = unique_map.get(key)
            if existing is None:
                unique_map[key] = disease
                continue

            if existing.get("fromPMID") is True:
                continue
            if disease.get("fromPMID") is True:
                unique_map[key] = disease
                continue

            new_priority = ONTOLOGY_PRIORITY.get(disease.get("inDefinedTermSet"), 999)
            existing_priority = ONTOLOGY_PRIORITY.get(existing.get("inDefinedTermSet"), 999)
            if new_priority < existing_priority:
                unique_map[key] = disease

        doc["healthCondition"] = list(unique_map.values())
    return doc_list


# ---------------------------------------------------------------------------
# The pipeline stage
# ---------------------------------------------------------------------------
def augment_from_descriptions(docs, *, filter_supported_types=True):
    """Mine species and health conditions out of each record's description.

    Normal pipeline callers retain the record-type guard. Cache prewarming can
    disable it so sources that are temporarily skipping this stage can still
    populate their caches offline.
    """
    # Keep this guard at the augmentation boundary as well as in the pipeline.
    doc_list = list(docs)
    if filter_supported_types:
        doc_list = [doc for doc in doc_list if supports_description_enrichment(doc)]
    started = time.monotonic()
    timings = {}
    steps = (
        ("extract", _extract_entities),
        ("standardize_species", _standardize_extracted_species),
        ("dedupe_species", _dedupe_species),
        ("remove_redundant_species", _remove_redundant_species),
        ("standardize_diseases", _standardize_extracted_diseases),
        ("dedupe_diseases", _dedupe_diseases),
    )
    for name, step in steps:
        step_started = time.monotonic()
        step(doc_list)
        timings[name] = time.monotonic() - step_started

    logger.info(
        "Descriptions: docs=%s extract=%.1fs standardize_species=%.1fs dedupe_species=%.1fs "
        "remove_redundant_species=%.1fs standardize_diseases=%.1fs dedupe_diseases=%.1fs total=%.1fs",
        len(doc_list),
        timings["extract"],
        timings["standardize_species"],
        timings["dedupe_species"],
        timings["remove_redundant_species"],
        timings["standardize_diseases"],
        timings["dedupe_diseases"],
        time.monotonic() - started,
    )
    return doc_list
