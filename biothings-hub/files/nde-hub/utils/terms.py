"""Species, infectiousAgent and healthCondition standardization.

Runs whenever a batch of records carries any of those three fields. Each term
is resolved in this order:

  1. the PubTator lookup DB (`pubtator_lookup.db`), matched on name, original
     name or any alternate name,
  2. the resolved-species cache (`extract_lookup.db`),
  3. UniProt taxonomy, when the record supplies a numeric NCBI taxon id,
  4. text2term against ncbitaxon, for names with no usable identifier.

Names that resolve to nothing are remembered in a negative cache so later runs
don't retry them. Species are split into `species` (hosts) and
`infectiousAgent` from their UniProt lineage.
"""

import datetime
import json
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from config import logger

from .cache import SqliteCache, SqliteKeySet
from .common import as_list, sqlite
from .taxonomy import classify_from_lineage

DB_PATH = "/data/nde-hub/standardizers/pubtator_lookup/pubtator_lookup.db"
SPECIES_CACHE_DB_PATH = "/data/nde-hub/standardizers/extract_lookup/extract_lookup.db"

_NEGATIVE_SPECIES_TABLE = "species_negative"
_TERM_FIELDS = ("species", "infectiousAgent", "healthCondition")
_SPECIES_FIELDS = ("species", "infectiousAgent")

# Flags that mean a term has already been standardized by someone else.
_STANDARDIZED_FLAGS = ("curatedBy", "fromPMID", "fromEXTRACT")

_TAXID_RE = re.compile(r"^\d+$")

# Place names commonly mistaken for organisms.
DROP_LIST_TERMS = {
    "sonoma": {"id": "1535511"},
    "china": {"id": "3034371"},
    "nevada": {"id": "359889"},
    "montana": {"id": "441235"},
}

# Reuse HTTP connections for UniProt lookups.
_UNIPROT_SESSION = requests.Session()

# Species we resolved ourselves, and names that resolved to nothing. memoize=False
# because a big source can meet hundreds of thousands of distinct species and each
# batch only needs its own; writes still reach later batches through the table.
SPECIES_DETAILS = SqliteCache(SPECIES_CACHE_DB_PATH, "species_details", memoize=False)
NEGATIVE_SPECIES = SqliteKeySet(SPECIES_CACHE_DB_PATH, _NEGATIVE_SPECIES_TABLE)

_lookup_cache = None
_identifier_resolutions = {}


def reset_caches():
    """Drop the process-wide caches so a new upload sees fresh lookup data."""
    global _lookup_cache
    _lookup_cache = None
    _identifier_resolutions.clear()
    SPECIES_DETAILS.reset()
    NEGATIVE_SPECIES.reset()


# ---------------------------------------------------------------------------
# Drop list
# ---------------------------------------------------------------------------
def should_filter_term(term_name, identifier=None):
    """True when a term is a known place name rather than an organism."""
    if not term_name:
        return False

    if term_name.lower().strip() in DROP_LIST_TERMS:
        logger.info("Filtering term '%s': place name", term_name)
        return True

    if identifier:
        clean_id = str(identifier).split("*")[-1].strip()
        for term_config in DROP_LIST_TERMS.values():
            if term_config["id"] == clean_id:
                logger.info("Filtering term '%s' by ID %s: place name", term_name, clean_id)
                return True
    return False


# ---------------------------------------------------------------------------
# The PubTator lookup dictionaries
# ---------------------------------------------------------------------------
class AliasLookupDict(dict):
    """A lookup dict that already indexes every alias, so no scan is needed."""

    aliases_indexed = True


_LOOKUP_DDL = (
    "CREATE TABLE IF NOT EXISTS health_conditions (original_name text, standard_dict text)",
    "CREATE TABLE IF NOT EXISTS species (original_name text, standard_dict text)",
)


def lookup_db():
    """Open the PubTator lookup DB, creating it and its tables if needed."""
    db_dir = os.path.dirname(DB_PATH)
    if db_dir:
        os.makedirs(db_dir, exist_ok=True)
    return sqlite(DB_PATH, *_LOOKUP_DDL)


def _add_lookup_alias(lookup_dict, alias, item_data):
    if not alias:
        return
    normalized_alias = str(alias).lower().strip()
    if normalized_alias:
        lookup_dict.setdefault(normalized_alias, item_data)


def _build_lookup_dict(rows):
    lookup_dict = AliasLookupDict()
    loaded_items = []

    for original_name, standard_dict in rows:
        if not standard_dict:
            continue
        item_data = json.loads(standard_dict)
        _add_lookup_alias(lookup_dict, original_name, item_data)
        loaded_items.append(item_data)

    for item_data in loaded_items:
        _add_lookup_alias(lookup_dict, item_data.get("name"), item_data)
        _add_lookup_alias(lookup_dict, item_data.get("originalName"), item_data)
        for alternate_name in as_list(item_data.get("alternateName")):
            _add_lookup_alias(lookup_dict, alternate_name, item_data)

    return lookup_dict


def _lookup_dicts():
    """Load (once per upload) the health condition and species lookup dictionaries."""
    global _lookup_cache
    if _lookup_cache is None:
        with lookup_db() as conn:
            hc_rows = conn.execute("SELECT original_name, standard_dict FROM health_conditions").fetchall()
            species_rows = conn.execute("SELECT original_name, standard_dict FROM species").fetchall()
        _lookup_cache = (_build_lookup_dict(hc_rows), _build_lookup_dict(species_rows))
        logger.info("Term lookup loaded: %s health conditions, %s species", len(_lookup_cache[0]), len(_lookup_cache[1]))
    return _lookup_cache


def lookup_item(original_name, data_dict):
    """Find a standardized term for `original_name` in a lookup dictionary."""
    if not original_name or not data_dict:
        return None

    original_name_lower = original_name.lower().strip()
    if original_name_lower in data_dict:
        return data_dict[original_name_lower]

    if getattr(data_dict, "aliases_indexed", False):
        return None

    for item_data in data_dict.values():
        if "name" in item_data and original_name_lower == item_data["name"].lower().strip():
            return item_data
        for alternate_name in as_list(item_data.get("alternateName")):
            if original_name_lower == alternate_name.lower().strip():
                return item_data
    return None


# ---------------------------------------------------------------------------
# Ontology lookups for health conditions
# ---------------------------------------------------------------------------
ONTOLOGY_URLS = [
    "https://biothings.transltr.io/mondo",
    "https://biothings.transltr.io/hpo",
    "https://biothings.transltr.io/doid",
    "https://biothings.transltr.io/ncit",
]


def _retry_request(url, retries=7):
    for _ in range(retries):
        response = requests.get(url)
        try:
            return response.json()
        except json.decoder.JSONDecodeError:
            logger.info("Retrying...")
    logger.info("Failed to decode JSON")
    return None


def process_synonyms(synonym_field):
    if isinstance(synonym_field, dict):
        return synonym_field.get("exact", [])
    return [syn.split('"')[1] for syn in synonym_field if "EXACT" in syn]


def create_return_object(hit, alternate_names, original_name):
    ontology, identifier = hit["_id"].split(":")[:2]
    standard_dict = {
        "@type": "DefinedTerm",
        "identifier": identifier,
        "inDefinedTermSet": ontology,
        "isCurated": True,
        "name": hit.get("label") or hit.get("name"),
        "originalName": original_name,
        "url": f"http://purl.obolibrary.org/obo/{ontology}_{identifier}",
        "curatedBy": {
            "name": "Biothings API",
            "url": "https://biothings.io/",
            "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
        },
    }
    if alternate_names:
        standard_dict["alternateName"] = list(set(alternate_names))
    return standard_dict


def _handle_response(data, condition, base_url, match_condition=True):
    for hit in (data or {}).get("hits", []):
        alternate_names = process_synonyms(hit.get("synonym", {}))
        term_name = hit.get("label") or hit.get("name")
        if not term_name:
            continue
        if not match_condition:
            logger.info("Found %s via xrefs.mesh in ontology: %s", condition, base_url.split("/")[-1])
            return create_return_object(hit, alternate_names, condition)
        condition_lower = condition.lower().strip()
        if term_name.lower().strip() == condition_lower or any(name.lower().strip() == condition_lower for name in alternate_names):
            logger.info("Found %s in ontology: %s", condition, base_url.split("/")[-1])
            return create_return_object(hit, alternate_names, condition)
    return None


def query_condition(health_condition, mesh_id=None):
    """Look a health condition up in MONDO, HPO, DOID then NCIT."""
    logger.info('Querying for "%s"...', health_condition)
    for base_url in ONTOLOGY_URLS:
        try:
            if mesh_id:
                data = _retry_request(f'{base_url}/query?q=xrefs.mesh:"{mesh_id}"&limit=1000')
                if result := _handle_response(data, health_condition, base_url, match_condition=False):
                    return result

            for query in (
                f'label:("{health_condition}")',
                f'name:("{health_condition}")',
                f'synonym.exact:"{health_condition}"',
            ):
                data = _retry_request(f"{base_url}/query?q={query}&limit=1000")
                if result := _handle_response(data, health_condition, base_url):
                    return result
        except Exception as e:
            logger.info("An error occurred while querying %s: %s", base_url, e)
    logger.info("Unable to find %s", health_condition)
    return None


# ---------------------------------------------------------------------------
# UniProt taxonomy
# ---------------------------------------------------------------------------
def normalize_taxon_id(identifier):
    """Return a numeric NCBI taxonomy id, or None for non-taxon placeholders."""
    if not identifier:
        return None
    candidate = str(identifier).split("*")[-1].strip()
    if candidate.isdigit():
        return candidate
    last_part = candidate.split(":")[-1].strip()
    return last_part if last_part.isdigit() else None


def fetch_taxon(original_name, identifier, classify=classify_from_lineage, max_retries=3):
    """Build a DefinedTerm for a taxon from its UniProt taxonomy entry.

    Provenance is the caller's to stamp: this sets no `isCurated`, `curatedBy` or
    `fromEXTRACT`, and leaves `classification` unset when UniProt reports no
    lineage. `lineage` is included so callers can filter on it; drop it before
    the term reaches a record.

    Raises ValueError for an identifier that isn't a numeric NCBI taxon id, since
    UniProt taxonomy would reject it anyway.
    """
    taxon_id = normalize_taxon_id(identifier)
    if not taxon_id:
        raise ValueError(f"Invalid NCBI Taxonomy ID: {identifier}")

    for attempt in range(max_retries):
        response = _UNIPROT_SESSION.get(f"https://rest.uniprot.org/taxonomy/{taxon_id}", timeout=30)
        if response.status_code == 429:
            retry_after = int(response.headers.get("Retry-After", 2**attempt))
            logger.warning("UniProt 429 for %s, retrying in %ss", taxon_id, retry_after)
            time.sleep(retry_after)
            continue
        response.raise_for_status()
        break
    else:
        raise requests.exceptions.HTTPError(f"UniProt rate limit exceeded after {max_retries} retries for {taxon_id}")

    species_info = response.json()
    standard_dict = {
        "@type": "DefinedTerm",
        "identifier": taxon_id,
        "inDefinedTermSet": "UniProt",
        "url": f"https://www.uniprot.org/taxonomy/{taxon_id}",
        "originalName": original_name,
    }
    _add_uniprot_names(standard_dict, species_info, original_name)

    if lineage := species_info.get("lineage"):
        standard_dict["classification"] = classify(standard_dict["name"], lineage)
        standard_dict["lineage"] = lineage
    else:
        logger.warning("No lineage found for %s", taxon_id)
    return standard_dict


def get_species_details(original_name, identifier):
    """Standardize a species from UniProt, curated by PubTator.

    Used when a PubTator annotation supplies the taxon id.
    """
    logger.info("Getting details for %s", original_name)
    if should_filter_term(original_name, identifier):
        logger.info("Skipping %s: filtered by drop list", original_name)
        return None

    term = fetch_taxon(original_name, identifier)
    term.pop("lineage", None)
    term["isCurated"] = True
    term["curatedBy"] = {
        "name": "PubTator",
        "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
        "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
    }
    return term


def _get_uniprot_details(original_name, identifier):
    """Fetch species details from UniProt for our own resolution (not curated)."""
    term = fetch_taxon(original_name, identifier)
    term["isCurated"] = False
    # Nothing to classify from means we cannot call it a host.
    term.setdefault("classification", "infectiousAgent")
    return term


def _add_uniprot_names(standard_dict, species_info, original_name):
    """Fill name, commonName, displayName and alternateName from a UniProt taxon."""
    standard_dict["name"] = species_info.get("scientificName") or original_name

    alternative_names = []
    if common_name := species_info.get("commonName"):
        standard_dict["commonName"] = common_name
        alternative_names.append(common_name)
        standard_dict["displayName"] = f"{common_name} | {standard_dict['name']}"
    else:
        standard_dict["displayName"] = standard_dict["name"]

    alternative_names.extend(species_info.get("otherNames") or [])
    if alternative_names:
        standard_dict["alternateName"] = list(set(alternative_names))


# ---------------------------------------------------------------------------
# Scanning a batch for terms that still need resolving
# ---------------------------------------------------------------------------
def _resolve_name_via_identifier(item, field):
    """Name a species entry that has only an identifier, via UniProt.

    Mutates `item` in place on success and returns the resolved name, or None
    when the entry should be skipped. Resolutions are cached per upload.
    """
    if field not in _SPECIES_FIELDS:
        return None
    identifier = item.get("identifier")
    if not identifier:
        return None

    cache_key = str(identifier).split("*")[-1].strip()
    # UniProt taxonomy only accepts numeric NCBI tax IDs; skip anything else
    # (e.g. IPC###, internal codes) to avoid guaranteed 400s.
    if not _TAXID_RE.match(cache_key):
        logger.debug("Skipping non-numeric %s identifier: %s", field, identifier)
        _identifier_resolutions[cache_key] = None
        return None

    if cache_key in _identifier_resolutions:
        details = _identifier_resolutions[cache_key]
    else:
        try:
            details = _get_uniprot_details(item.get("originalName") or str(identifier), identifier)
        except Exception as e:
            logger.warning("Could not resolve %s entry via identifier %s: %s", field, identifier, e)
            _identifier_resolutions[cache_key] = None
            return None
        _identifier_resolutions[cache_key] = details

    name = details.get("name") if details else None
    if not name:
        return None
    for field_name, value in details.items():
        if field_name not in item or not item[field_name]:
            item[field_name] = value
    return name


def _scan_doc(doc, species_dict, unstandardized):
    """Name any identifier-only entries, and collect species we still have to resolve."""
    for field in _TERM_FIELDS:
        for entry in as_list(doc.get(field)):
            if not isinstance(entry, dict):
                continue
            if "curatedBy" in entry:
                continue

            name = entry.get("name")
            if not name:
                name = _resolve_name_via_identifier(entry, field)
                if not name:
                    logger.info("Skipping %s entry without resolvable name: %s", field, entry)
                    continue

            if field not in _SPECIES_FIELDS or "inDefinedTermSet" in entry:
                continue
            if should_filter_term(name, entry.get("identifier")):
                continue
            if lookup_item(name, species_dict):
                continue

            key = name.lower().strip()
            if key:
                unstandardized.setdefault(key, (name, entry.get("identifier")))


# ---------------------------------------------------------------------------
# Resolving species the lookup DB doesn't know
# ---------------------------------------------------------------------------
def _resolve_one(original_name, taxon_id):
    """Resolve one species via UniProt. Returns (key, details) with details None on failure."""
    key = original_name.lower().strip()
    try:
        return key, _get_uniprot_details(original_name, taxon_id)
    except ValueError as e:
        logger.info("Skipping UniProt lookup for %s (ID %s): %s", original_name, taxon_id, e)
        return key, None
    except Exception as e:
        logger.warning("UniProt lookup failed for %s (ID %s): %s", original_name, taxon_id, e)
        return key, None


def _resolve_via_uniprot(lookups, resolved):
    """Resolve (name, taxon_id) pairs concurrently, caching results. Returns the failures."""
    failed = []
    if not lookups:
        return failed
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(_resolve_one, name, taxon_id): name for name, taxon_id in lookups}
        for future in as_completed(futures):
            key, details = future.result()
            if details:
                resolved[key] = details
                SPECIES_DETAILS.put(details["originalName"], details)
            else:
                NEGATIVE_SPECIES.add(key)
                failed.append(futures[future])
    return failed


def _resolve_via_text2term(names, resolved):
    """Map names to ncbitaxon with text2term, then resolve the hits via UniProt."""
    import text2term

    if not os.path.exists("cache/ncbitaxon"):
        logger.info("Building text2term ncbitaxon ontology cache...")
        text2term.cache_ontology("https://purl.obolibrary.org/obo/ncbitaxon.owl", "ncbitaxon")

    results = text2term.map_terms(names, "ncbitaxon", use_cache=True)
    results.sort_values(["Source Term", "Mapping Score"], ascending=[True, False], inplace=True)
    results = results.drop_duplicates(subset=["Source Term"], keep="first")
    logger.info("text2term mapped %s of %s terms", len(results), len(names))

    lookups = []
    mapped = set()
    for _, row in results.iterrows():
        original_name = row["Source Term"]
        key = original_name.lower().strip()
        mapped.add(key)
        identifier = normalize_taxon_id(row["Mapped Term CURIE"])
        if identifier:
            lookups.append((original_name, identifier))
        else:
            NEGATIVE_SPECIES.add(key)

    # Remember the names text2term couldn't map at all.
    for name in names:
        if name.lower().strip() not in mapped:
            NEGATIVE_SPECIES.add(name)

    _resolve_via_uniprot(lookups, resolved)


def _resolve_species(unstandardized):
    """Resolve `{key: (name, identifier)}` into `{key: standardized term}`."""
    if not unstandardized:
        return {}

    negative_cache = NEGATIVE_SPECIES.known(unstandardized)
    resolved = SPECIES_DETAILS.get_many(set(unstandardized) - negative_cache)

    need_uniprot = []
    need_text2term = []
    for key, (name, identifier) in unstandardized.items():
        if key in negative_cache or key in resolved:
            continue
        if taxon_id := normalize_taxon_id(identifier):
            need_uniprot.append((name, taxon_id))
        else:
            need_text2term.append(name)

    logger.info(
        "Species resolution: %s unstandardized, %s cached, %s via UniProt, %s via text2term",
        len(unstandardized),
        len(resolved),
        len(need_uniprot),
        len(need_text2term),
    )

    need_text2term.extend(_resolve_via_uniprot(need_uniprot, resolved))
    if need_text2term:
        try:
            _resolve_via_text2term(need_text2term, resolved)
        except Exception as e:
            logger.error("Error during text2term resolution: %s", e)

    logger.info("Species resolution: %s resolved", len(resolved))
    return resolved


def _apply_resolved_species(doc, resolved):
    """Replace a record's unstandardized species/infectiousAgent with resolved terms."""
    entries = [(field, entry) for field in _SPECIES_FIELDS for entry in as_list(doc.get(field)) if isinstance(entry, dict)]
    if not any("inDefinedTermSet" not in entry and "curatedBy" not in entry and "name" in entry for _, entry in entries):
        return

    new_species = []
    new_infectious_agents = []
    seen_identifiers = set()

    for original_field, entry in entries:
        name = entry.get("name")
        if not name:
            continue

        term = entry
        if "inDefinedTermSet" not in entry and "curatedBy" not in entry:
            standardized = resolved.get(name.lower().strip())
            if standardized is None:
                # Nothing resolved this name; leave it where the source put it.
                (new_infectious_agents if original_field == "infectiousAgent" else new_species).append(entry)
                continue
            term = standardized.copy()
            term.pop("lineage", None)
            term["originalName"] = name

        identifier = term.get("identifier")
        if identifier and identifier in seen_identifiers:
            continue
        if identifier:
            seen_identifiers.add(identifier)

        if term is entry:
            is_agent = original_field == "infectiousAgent" or term.get("classification") == "infectiousAgent"
        else:
            is_agent = term.get("classification") == "infectiousAgent"
        (new_infectious_agents if is_agent else new_species).append(term)

    if new_species:
        doc["species"] = new_species
    else:
        doc.pop("species", None)
    if new_infectious_agents:
        doc["infectiousAgent"] = new_infectious_agents


# ---------------------------------------------------------------------------
# Applying the lookup dictionaries to a record
# ---------------------------------------------------------------------------
def _standardize_section(section, lookup_dict, is_species_section=False):
    new_section = []
    for original_obj in section:
        if not isinstance(original_obj, dict):
            logger.error("Invalid object: %s", original_obj)
            continue
        if any(key in original_obj for key in _STANDARDIZED_FLAGS):
            new_section.append(original_obj)
            continue

        original_name = original_obj.get("name")
        if not original_name:
            if is_species_section:
                original_name = _resolve_name_via_identifier(original_obj, "species")
            if not original_name:
                new_section.append(original_obj)
                continue

        if is_species_section and should_filter_term(original_name, original_obj.get("identifier")):
            logger.info("Filtering out '%s' from species section", original_name)
            continue

        new_obj = lookup_item(original_name, lookup_dict)
        if not new_obj:
            new_section.append(original_obj)
            continue

        if is_species_section:
            if original_obj.get("classification") and not new_obj.get("classification"):
                new_obj = new_obj.copy()
                new_obj["classification"] = original_obj["classification"]
            if should_filter_term(new_obj.get("name"), new_obj.get("identifier")):
                logger.info("Filtering out retrieved '%s'", new_obj.get("name"))
                continue
        new_section.append(new_obj)
    return new_section


def _dedupe_by_identifier(entries):
    seen = set()
    unique = []
    for entry in entries:
        if not entry:
            continue
        identifier = entry.get("identifier")
        if identifier not in seen:
            seen.add(identifier)
            unique.append(entry)
    return unique


def standardize_doc_terms(doc, hc_dict, species_dict):
    """Standardize one record's species, infectiousAgent and healthCondition."""
    health_conditions = _standardize_section(as_list(doc.get("healthCondition")), hc_dict)
    organisms = _standardize_section(
        as_list(doc.get("species")) + as_list(doc.get("infectiousAgent")),
        species_dict,
        is_species_section=True,
    )

    infectious_agents = [item for item in organisms if item.get("classification") == "infectiousAgent"]
    species = [item for item in organisms if item.get("classification") != "infectiousAgent"]

    # A species reclassified as an infectiousAgent must not stay in both fields.
    converted_names = {
        name.lower().strip()
        for agent in infectious_agents
        for name in (agent.get("originalName"), agent.get("name"))
        if name
    }
    species = [item for item in species if (item.get("name") or "").lower().strip() not in converted_names]

    if health_conditions := _dedupe_by_identifier(health_conditions):
        doc["healthCondition"] = health_conditions
    if species := _dedupe_by_identifier(species):
        doc["species"] = species
    else:
        doc.pop("species", None)
    if infectious_agents := _dedupe_by_identifier(infectious_agents):
        doc["infectiousAgent"] = infectious_agents

    return doc


# ---------------------------------------------------------------------------
# The pipeline stage
# ---------------------------------------------------------------------------
def standardize_terms(docs):
    """Standardize the species, infectiousAgent and healthCondition of one batch.

    Applying the lookup dictionaries costs a few microseconds per record, so
    this runs in-process: a worker pool spent more time pickling records than
    the work itself, and each fork held its own copy of the lookup tables.
    """
    docs = list(docs)
    hc_dict, species_dict = _lookup_dicts()

    unstandardized = {}
    for doc in docs:
        _scan_doc(doc, species_dict, unstandardized)

    resolved = _resolve_species(unstandardized)

    for doc in docs:
        standardize_doc_terms(doc, hc_dict, species_dict)
        if resolved:
            _apply_resolved_species(doc, resolved)
        yield doc
