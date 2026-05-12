import datetime
import json
import logging
import os
import re
import sqlite3
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool

import orjson
import requests
from biothings.utils.dataload import tab2dict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

DB_PATH = "/data/nde-hub/standardizers/pubtator_lookup/pubtator_lookup.db"
EXTRACT_DB_PATH = "/data/nde-hub/standardizers/extract_lookup/extract_lookup.db"
_NEGATIVE_SPECIES_TABLE = "species_negative"

# Reuse HTTP connections for UniProt lookups.
_UNIPROT_SESSION = requests.Session()

# Regular expression for matching taxon IDs
_TAXID_RE = re.compile(r"^\d+$")

# Drop list for filtering out place names commonly confused with organisms
DROP_LIST_TERMS = {
    # Term names to filter (case-insensitive)
    "sonoma": {
        "id": "1535511"
    },
    "china": {
        "id": "3034371"
    },
    "nevada": {
        "id": "359889"
    },
    "montana": {
        "id": "441235"
    }
}


def should_filter_term(term_name, identifier=None):
    """
    Check if a term should be filtered based on the drop list.

    Args:
        term_name (str): The name of the term to check
        identifier (str): The identifier of the term (optional)

    Returns:
        bool: True if the term should be filtered, False otherwise
    """
    if not term_name:
        return False

    term_lower = term_name.lower().strip()

    # Check by term name
    if term_lower in DROP_LIST_TERMS:
        logger.info(f"Filtering term '{term_name}': place name")
        return True

    # Check by identifier if provided
    if identifier:
        # Clean identifier (remove any prefixes/suffixes)
        clean_id = identifier.split("*")[-1].strip()
        for term_config in DROP_LIST_TERMS.values():
            if term_config["id"] == clean_id:
                logger.info(f"Filtering term '{term_name}' by ID {clean_id}: place name")
                return True

    return False


_IDENTIFIER_RESOLUTION_CACHE = {}
_TRANSFORM_HC_DICT = None
_TRANSFORM_SPECIES_DICT = None


class AliasLookupDict(dict):
    aliases_indexed = True


def _resolve_name_via_identifier(item, key):
    """For species/infectiousAgent items lacking a name, look up the taxon via UniProt.

    Mutates the item in place with the resolved fields when successful.
    Returns the resolved name, or None when the item should be skipped.
    Resolutions are cached per-run by normalized identifier.
    """
    if key not in ("species", "infectiousAgent"):
        return None
    identifier = item.get("identifier")
    if not identifier:
        return None

    cache_key = str(identifier).split("*")[-1].strip()

    # UniProt taxonomy only accepts numeric NCBI tax IDs; skip anything else
    # (e.g. IPC###, internal codes) to avoid guaranteed 400s.
    if not _TAXID_RE.match(cache_key):
        logger.debug(f"Skipping non-numeric {key} identifier: {identifier}")
        _IDENTIFIER_RESOLUTION_CACHE[cache_key] = None
        return None

    if cache_key in _IDENTIFIER_RESOLUTION_CACHE:
        details = _IDENTIFIER_RESOLUTION_CACHE[cache_key]
    else:
        try:
            details = _get_uniprot_details(item.get("originalName") or str(identifier), identifier)
        except Exception as e:
            logger.warning(f"Could not resolve {key} entry via identifier {identifier}: {e}")
            _IDENTIFIER_RESOLUTION_CACHE[cache_key] = None
            return None
        _IDENTIFIER_RESOLUTION_CACHE[cache_key] = details

    name = details.get("name") if details else None
    if not name:
        return None
    for field, value in details.items():
        if field not in item or not item[field]:
            item[field] = value
    return name


def _iter_data_file(data_folder):
    with open(os.path.join(os.fspath(data_folder), "data.ndjson"), "rb") as f:
        for line in f:
            yield orjson.loads(line)


def _iter_ndjson_file(file_path):
    with open(file_path, "rb") as f:
        for line in f:
            yield orjson.loads(line)


def _is_reiterable(data):
    try:
        return iter(data) is not data
    except TypeError:
        return False


def _add_unique_value(values, name):
    if not name:
        return
    normalized_name = name.lower().strip()
    if normalized_name:
        values.setdefault(normalized_name, None)


def _iter_section_items(doc, key):
    value = doc.get(key)
    if isinstance(value, list):
        return value
    if isinstance(value, dict) and value:
        return [value]
    return []


def _collect_values_from_doc(doc, key, values):
    for item in _iter_section_items(doc, key):
        if not isinstance(item, dict):
            continue
        if "curatedBy" in item:
            logger.info(f"{item.get('name', '(unnamed)')} has already been curated, skipping...")
            continue
        name = item.get("name")
        if not name:
            name = _resolve_name_via_identifier(item, key)
            if not name:
                logger.info(f"Skipping {key} entry without resolvable name: {item}")
                continue
        _add_unique_value(values, name)


def extract_values(doc_list, key):
    values = {}
    for doc in doc_list:
        _collect_values_from_doc(doc, key, values)
    return list(values.keys())


def _collect_unstandardized_species_from_doc(doc, species_dict, unstandardized_species):
    for field in ("species", "infectiousAgent"):
        for entry in _iter_section_items(doc, field):
            if not isinstance(entry, dict):
                continue
            if "inDefinedTermSet" in entry or "curatedBy" in entry:
                continue

            name = entry.get("name")
            if not name:
                name = _resolve_name_via_identifier(entry, field)
                if not name:
                    continue
                if "inDefinedTermSet" in entry or "curatedBy" in entry:
                    continue

            if should_filter_term(name, entry.get("identifier")):
                continue
            if species_dict and lookup_item(name, species_dict):
                continue

            key = name.lower().strip()
            if key:
                unstandardized_species.setdefault(key, (name, entry.get("identifier")))


def _scan_documents(docs, species_dict=None, spool_file=None):
    health_conditions = {}
    species = {}
    infectious_agents = {}
    unstandardized_species = {}
    count = 0

    for doc in docs:
        count += 1
        if count % 1000 == 0:
            logger.info(f"Processed {count} lines")

        _collect_values_from_doc(doc, "healthCondition", health_conditions)
        _collect_values_from_doc(doc, "species", species)
        _collect_values_from_doc(doc, "infectiousAgent", infectious_agents)

        if species_dict is not None:
            _collect_unstandardized_species_from_doc(doc, species_dict, unstandardized_species)

        if spool_file is not None:
            spool_file.write(orjson.dumps(doc))
            spool_file.write(b"\n")

    return (
        list(health_conditions.keys()),
        list(species.keys()),
        list(infectious_agents.keys()),
        unstandardized_species,
        count,
    )


def standardize_data(data):
    return _standardize_data(data)


def _standardize_data(data):
    logger.info("Standardizing data...")
    temp_path = None

    try:
        _ensure_lookup_db()
        hc_dict, species_dict = fetch_data_from_db()

        if isinstance(data, (str, os.PathLike)):
            logger.info("Reading data from file...")
            source_factory = lambda: _iter_data_file(data)
            scan_result = _scan_documents(source_factory(), species_dict=species_dict)
        elif _is_reiterable(data):
            logger.info("Reading data from reusable iterable...")
            source_factory = lambda: iter(data)
            scan_result = _scan_documents(source_factory(), species_dict=species_dict)
        else:
            logger.info("Reading data from one-pass iterable; spooling documents to disk for a streaming second pass...")
            with tempfile.NamedTemporaryFile(prefix="pubtator-standardize-", suffix=".ndjson", delete=False) as f:
                temp_path = f.name
                scan_result = _scan_documents(data, species_dict=species_dict, spool_file=f)
            source_factory = lambda: _iter_ndjson_file(temp_path)

        health_conditions_list, species_list, infectious_agents_list, unstandardized_species, count = scan_result
        logger.info(f"Scanned {count} documents")
        logger.info(
            f"Found {len(health_conditions_list)} health conditions, {len(species_list)} species, {len(infectious_agents_list)} infectious agents"
        )

        docs = update_lookup_dict(
            health_conditions_list,
            species_list,
            infectious_agents_list,
            source_factory,
            unstandardized_species=unstandardized_species,
            hc_dict=hc_dict,
            species_dict=species_dict,
        )
        for doc in docs:
            yield doc
    finally:
        if temp_path:
            try:
                os.remove(temp_path)
            except FileNotFoundError:
                pass


def classify_as_host_or_agent(lineage):
    """
    Classifies a species as either host or infectious agent based on its lineage.

    Parameters:
    - lineage (list): The lineage of the species.

    Returns:
    - new_classification (str): The classification of the species.
    """
    # Extracting scientific names for easy processing
    scientific_names = [item["scientificName"] for item in lineage]

    # Check for host species conditions
    if "Deuterostomia" in scientific_names:
        logger.info(f"Found Deuterostomia in {scientific_names}, classifying as host")
        new_classification = "host"
    elif "Embryophyta" in scientific_names and not any(
        parasite in scientific_names for parasite in ["Arceuthobium", "Cuscuta", "Orobanche", "Striga", "Phoradendron"]
    ):
        logger.info(f"Found Embryophyta in {scientific_names}, classifying as host")
        new_classification = "host"
    elif "Arthropoda" in scientific_names:
        if "Acari" in scientific_names:
            if "Ixodida" in scientific_names:
                logger.info(f"Found Ixodida in {scientific_names}, classifying as host")
                new_classification = "host"
            else:
                logger.info(f"Found Acari in {scientific_names}, classifying as infectiousAgent")
                new_classification = "infectiousAgent"
        else:
            logger.info(f"Found Arthropoda in {scientific_names}, classifying as host")
            new_classification = "host"
    else:
        # If not falling under the above host conditions, classify as infectiousAgent
        logger.info(f"Found {scientific_names}, classifying as infectiousAgent")
        new_classification = "infectiousAgent"

    return new_classification


def get_species_details(original_name, identifier):
    logger.info(f"Getting details for {original_name}")

    # Check if this term should be filtered out
    should_filter = should_filter_term(original_name, identifier)
    if should_filter:
        logger.info(f"Skipping {original_name}: filtered by drop list")
        return None

    # Fetch details from the UniProt API
    identifier = identifier.split("*")[-1]
    # try:
    species_info = requests.get(f"https://rest.uniprot.org/taxonomy/{identifier}")
    species_info.raise_for_status()
    species_info = species_info.json()
    standard_dict = {
        "@type": "DefinedTerm",
        "identifier": identifier,
        "inDefinedTermSet": "UniProt",
        "url": f"https://www.uniprot.org/taxonomy/{identifier}",
        "originalName": original_name,
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
        },
    }
    if scientific_name := species_info.get("scientificName"):
        standard_dict["name"] = scientific_name
    else:
        standard_dict["name"] = original_name

    alternative_names = []
    if common_name := species_info.get("commonName"):
        standard_dict["commonName"] = common_name
        alternative_names.append(common_name)

        standard_dict["displayName"] = f"{common_name} | {scientific_name}"
    else:
        standard_dict["displayName"] = scientific_name

    if other_names := species_info.get("otherNames"):
        alternative_names.extend(other_names)

    if alternative_names:
        standard_dict["alternateName"] = list(set(alternative_names))

    if lineage := species_info.get("lineage"):
        standard_dict["classification"] = classify_as_host_or_agent(lineage)
    else:
        logger.warning(f"No lineage found for {identifier}")

    return standard_dict


def load_json_response(response):
    try:
        return response.json()
    except json.decoder.JSONDecodeError:
        return None


def retry_request(url, retries=7):
    for _ in range(retries):
        response = requests.get(url)
        data = load_json_response(response)
        if data is not None:
            return data
        logger.info("Retrying...")
    logger.info("Failed to decode JSON")
    return None


def handle_response(data, condition, base_url, match_condition=True):
    if "hits" in data:
        for hit in data["hits"]:
            alternate_names = process_synonyms(hit.get("synonym", {}))
            term_name = hit.get("label") or hit.get("name")
            if term_name:
                if match_condition:
                    if term_name.lower().strip() == condition.lower().strip() or any(
                        name.lower().strip() == condition.lower().strip() for name in alternate_names
                    ):
                        logger.info(f"Found {condition} in ontology: {base_url.split('/')[-1]}")
                        return create_return_object(hit, alternate_names, condition)
                else:
                    logger.info(f"Found {condition} via xrefs.mesh in ontology: {base_url.split('/')[-1]}")
                    return create_return_object(hit, alternate_names, condition)
    return None


def query_condition(health_condition, mesh_id=None):
    BASE_URLS = [
        "https://biothings.transltr.io/mondo",
        "https://biothings.transltr.io/hpo",
        "https://biothings.transltr.io/doid",
        "https://biothings.transltr.io/ncit",
    ]
    logger.info(f'Querying for "{health_condition}"...')
    for base_url in BASE_URLS:
        try:
            if mesh_id:
                # Query using xrefs.mesh
                url = f'{base_url}/query?q=xrefs.mesh:"{mesh_id}"&limit=1000'
                data = retry_request(url)
                result = handle_response(data, health_condition, base_url, match_condition=False)
                if result is not None:
                    return result

            url = f'{base_url}/query?q=label:("{health_condition}")&limit=1000'
            data = retry_request(url)
            result = handle_response(data, health_condition, base_url)
            if result is not None:
                return result

            url = f'{base_url}/query?q=name:("{health_condition}")&limit=1000'
            data = retry_request(url)
            result = handle_response(data, health_condition, base_url)
            if result is not None:
                return result

            url = f'{base_url}/query?q=synonym.exact:"{health_condition}"&limit=1000'
            data = retry_request(url)
            result = handle_response(data, health_condition, base_url)
            if result is not None:
                return result
        except Exception as e:
            logger.info(f"An error occurred while querying {base_url}: {e}")
    logger.info(f"Unable to find {health_condition}")
    return None


def process_synonyms(synonym_field):
    if isinstance(synonym_field, dict):
        # If the synonym field is a dictionary, return the exact synonyms.
        return synonym_field.get("exact", [])
    else:
        # If the synonym field is a list, filter it by 'EXACT'.
        return [syn.split('"')[1] for syn in synonym_field if "EXACT" in syn]


def create_return_object(hit, alternate_names, original_name):
    ontology = hit["_id"].split(":")[0]
    identifier = hit["_id"].split(":")[1]
    term_name = hit.get("label") or hit.get("name")

    standard_dict = {
        "@type": "DefinedTerm",
        "identifier": identifier,
        "inDefinedTermSet": ontology,
        "isCurated": True,
        "name": term_name,
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


def get_new_health_conditions(health_conditions):
    logger.info("Getting new health conditions...")
    found_health_conditions = []
    not_found_health_conditions = []
    for health_condition in health_conditions:
        logger.info(
            f"Processed {health_conditions.index(health_condition) + 1} / {len(health_conditions)} health conditions"
        )
        result = query_condition(health_condition)
        if result is not None:
            logger.info(f"Found {health_condition} details")
            found_health_conditions.append(health_condition)
            logger.info(f"Adding {health_condition} to database...")
            # result_dict = get_health_condition_details(result, health_condition)
            conn = sqlite3.connect(DB_PATH)
            c = conn.cursor()
            c.execute(
                "INSERT INTO health_conditions VALUES (?, ?)",
                (health_condition.lower().strip(), json.dumps(result)),
            )
            logger.info(f"Added {health_condition}")
            conn.commit()
            conn.close()
        else:
            logger.info(f"{health_condition} not found")
            not_found_health_conditions.append(health_condition)
    return (found_health_conditions, not_found_health_conditions)


def get_new_species(species):
    return None
    # Split the species list into chunks of 1000
    chunks = [species[x : x + 1000] for x in range(0, len(species), 1000)]
    logger.info(f"Splitting into {len(chunks)} chunks")

    chunk_count = 0
    for chunk in chunks:
        while True:
            chunk_count += 1
            logger.info(f"Processing chunk {chunk_count} ")

            # Convert the chunk to a JSON string
            data = json.dumps(".    ".join(chunk))

            # Submit the data to the PubTator API for annotation
            try:
                submit_response = requests.post(
                    "https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/submit/species",
                    data=data,
                )
            except requests.exceptions.ConnectionError as e:
                logger.info(f"Connection error: {e}")
                logger.info("Retrying...")
                continue

            logger.info(f"Waiting for response, {submit_response.text}")
            timeout = 0
            retries = 0
            response_recieved = False

            while True:
                # Wait for 10 seconds before checking the response
                time.sleep(10)
                timeout += 10

                # Retry the request if the timeout exceeds 100 seconds
                if timeout > 120:
                    retries += 1
                    # if retries > 5:
                    #     raise Exception("Attempted 5 times, giving up")
                    logger.info("Timeout, retrying...")
                    try:
                        os.remove(f"{submit_response.text}_response.csv")
                    except FileNotFoundError:
                        logger.info("Issue removing file: %s_response.csv", submit_response.text)
                    # get_new_species(species)
                    break

                # Check the response from the PubTator API
                try:
                    retrieve_response = requests.get(
                        f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/retrieve/{submit_response.text}"
                    )
                except requests.exceptions.ConnectionError as e:
                    logger.info(f"Connection error: {e}")
                    logger.info("Retrying...")
                    continue

                # If the response is successful, break the loop
                if retrieve_response.status_code == 200:
                    logger.info("Got response")
                    response_recieved = True
                    break
            if response_recieved:
                break

        # Process the retrieved result
        result = retrieve_response.text.split("00000|a|-NoAbstract-")[1].strip()

        if len(result.split("\n")) == 1:
            result = result + "\n" + result

        with open(f"{submit_response.text}_response.csv", "a") as f:
            for line in result.split("\n"):
                if "Species" in line:
                    f.write(line + "\n")
                else:
                    logger.info(f"Removed {line}")

        with open(f"{submit_response.text}_response.csv", "r") as f:
            lines = f.readlines()

        with open(f"{submit_response.text}_response.csv", "w") as f:
            for line in lines:
                if line.strip():
                    f.write(line)

            f.seek(0, 0)
            f.write("\n")

    if os.stat(f"{submit_response.text}_response.csv").st_size != 0:
        result = tab2dict(
            f"{submit_response.text}_response.csv",
            cols=[3, 4, 5],
            key=0,
            sep="\t",
            alwayslist=True,
        )

        remove_dupes = {}

        for key, value in result.items():
            if key.strip().lower() not in remove_dupes:
                remove_dupes[key.strip().lower()] = value
        for pubtator_species_name in remove_dupes:
            # check if the pubtator species name is already in the database
            conn = sqlite3.connect(DB_PATH)
            c = conn.cursor()
            c.execute("SELECT * FROM species WHERE original_name=?", (pubtator_species_name.lower().strip(),))
            result = c.fetchone()
            conn.close()
            if result:
                logger.info(f"Skipping {pubtator_species_name}, already in database")
                continue
            result_dict = get_species_details(
                # name of the species
                pubtator_species_name,
                # identifier
                remove_dupes[pubtator_species_name][0][1],
            )
            # Skip if the species was filtered out
            if result_dict is None:
                logger.info(f"Skipping {pubtator_species_name}, filtered by drop list")
                continue
            conn = sqlite3.connect(DB_PATH)
            c = conn.cursor()
            c.execute(
                "INSERT INTO species VALUES (?, ?)",
                (pubtator_species_name.lower(), json.dumps(result_dict)),
            )
            logger.info(f"Added {pubtator_species_name}")
            conn.commit()
            conn.close()

        difference = list(set([x.lower().strip() for x in species]) - set(remove_dupes.keys()))
        try:
            os.remove(f"{submit_response.text}_response.csv")
        except FileNotFoundError:
            logger.info("Issue removing file: %s_response.csv", submit_response.text)
        logger.info(f"New species: {difference}")

        not_found = []
        for submitted_species in species:
            if submitted_species.lower() not in remove_dupes.keys():
                not_found.append(submitted_species)

        return (difference, not_found)
    else:
        try:
            os.remove(f"{submit_response.text}_response.csv")
        except FileNotFoundError:
            logger.info("Issue removing file: %s_response.csv", submit_response.text)
        logger.info("No new species found")
        return None


def lookup_item(original_name, data_dict):
    if not original_name or not data_dict:
        return None

    original_name_lower = original_name.lower().strip()
    if original_name_lower in data_dict:
        return data_dict[original_name_lower]

    if getattr(data_dict, "aliases_indexed", False) or not hasattr(data_dict, "values"):
        return None

    for item_data in data_dict.values():
        if "name" in item_data and original_name_lower == item_data["name"].lower().strip():
            return item_data
        if "alternateName" in item_data:
            for alternate_name in item_data["alternateName"]:
                if original_name_lower == alternate_name.lower().strip():
                    return item_data
    return None


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
        alternate_names = item_data.get("alternateName", [])
        if isinstance(alternate_names, str):
            alternate_names = [alternate_names]
        for alternate_name in alternate_names:
            _add_lookup_alias(lookup_dict, alternate_name, item_data)

    return lookup_dict


def fetch_data_from_db():
    with sqlite3.connect(DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT original_name, standard_dict FROM health_conditions")
        hc_cursor = c.fetchall()
        c.execute("SELECT original_name, standard_dict FROM species")
        species_cursor = c.fetchall()

    hc_dict = _build_lookup_dict(hc_cursor)
    species_dict = _build_lookup_dict(species_cursor)

    return hc_dict, species_dict


def process_section(section, cursor_dict, is_species_section=False):
    new_section_list = []
    for original_obj in section:
        if isinstance(original_obj, str):
            logger.error(f"Invalid object: {original_obj}")
            continue
        # Check if the item has already been curated or standardized
        if any(key in original_obj for key in ["curatedBy", "fromPMID", "fromEXTRACT"]):
            logger.info(f"{original_obj.get('name')} has already been processed, skipping...")
            new_section_list.append(original_obj)
            continue
        original_name = original_obj.get("name")
        if not original_name:
            if is_species_section:
                original_name = _resolve_name_via_identifier(original_obj, "species")
            if not original_name:
                new_section_list.append(original_obj)
                continue

        # Apply drop list filtering for species-related sections
        if is_species_section:
            should_filter = should_filter_term(original_name,
                                               original_obj.get("identifier"))
            if should_filter:
                logger.info(f"Filtering out '{original_name}' from species section")
                continue

        new_obj = lookup_item(original_name, cursor_dict)
        if new_obj:
            # Check if the retrieved object should be filtered as well
            if is_species_section:
                should_filter = should_filter_term(new_obj.get("name"),
                                                   new_obj.get("identifier"))
                if should_filter:
                    logger.info(f"Filtering out retrieved '{new_obj.get('name')}'")
                    continue
            new_section_list.append(new_obj)
        else:
            new_section_list.append(original_obj)
    return new_section_list


def _init_transform_worker(hc_dict, species_dict, identifier_resolution_cache):
    global _TRANSFORM_HC_DICT, _TRANSFORM_SPECIES_DICT
    _TRANSFORM_HC_DICT = hc_dict
    _TRANSFORM_SPECIES_DICT = species_dict
    _IDENTIFIER_RESOLUTION_CACHE.update(identifier_resolution_cache or {})


def process_document(args):
    if len(args) == 4:
        doc, hc_dict, species_dict, doc_index = args
    else:
        doc, doc_index = args
        hc_dict = _TRANSFORM_HC_DICT or {}
        species_dict = _TRANSFORM_SPECIES_DICT or {}

    # Get sections from the document
    health_conditions_list = doc.get("healthCondition", {})
    species_list = doc.get("species", {})
    infectious_agent_list = doc.get("infectiousAgent", {})

    # Ensure each section is a list
    if isinstance(species_list, dict):
        species_list = [species_list]
    if isinstance(infectious_agent_list, dict):
        infectious_agent_list = [infectious_agent_list]
    if isinstance(health_conditions_list, dict):
        health_conditions_list = [health_conditions_list]

    # Process sections using lookup dictionaries.
    new_health_conditions_list = process_section(health_conditions_list, hc_dict)
    new_species_and_infectious_agents_list = process_section(species_list + infectious_agent_list, species_dict, is_species_section=True)

    def remove_duplicates_from_list(data_list):
        seen_identifiers = set()
        unique_list = []
        for entry in data_list:
            if entry:
                identifier = entry.get("identifier")
                if identifier not in seen_identifiers:
                    seen_identifiers.add(identifier)
                    unique_list.append(entry)
        return unique_list

    # Split processed items into species and infectious agents.
    new_species_list = [
        item for item in new_species_and_infectious_agents_list if item.get("classification") != "infectiousAgent"
    ]
    new_infectious_agent_list = [
        item for item in new_species_and_infectious_agents_list if item.get("classification") == "infectiousAgent"
    ]

    # Build a set of normalized names for converted species (using both originalName and name).
    converted_names = set()
    for agent in new_infectious_agent_list:
        if agent.get("originalName"):
            converted_names.add(agent.get("originalName").lower().strip())
        if agent.get("name"):
            converted_names.add(agent.get("name").lower().strip())

    # Filter out any species entry that matches any converted name.
    filtered_species_list = []
    for item in new_species_list:
        sp_name = None
        if isinstance(item, dict):
            sp_name = item.get("name")
        else:
            sp_name = item

        if sp_name and sp_name.lower().strip() in converted_names:
            # Species was converted to an infectiousAgent; skip it.
            continue
        else:
            filtered_species_list.append(item)

    # Remove duplicates from each section.
    no_hc_dupes = remove_duplicates_from_list(new_health_conditions_list) if new_health_conditions_list else []
    no_ia_dupes = remove_duplicates_from_list(new_infectious_agent_list) if new_infectious_agent_list else []
    no_sp_dupes = remove_duplicates_from_list(filtered_species_list) if filtered_species_list else []

    # Update document sections.
    if no_hc_dupes:
        doc["healthCondition"] = no_hc_dupes
    if no_sp_dupes:
        doc["species"] = no_sp_dupes
    else:
        # Remove species property if list is empty.
        if "species" in doc:
            del doc["species"]
    if no_ia_dupes:
        doc["infectiousAgent"] = no_ia_dupes

    # Log progress every 1000 documents.
    if doc_index % 1000 == 0:
        logger.info(f"Processed {doc_index} documents")

    return doc


def transform(doc_list):
    hc_dict, species_dict = fetch_data_from_db()
    return transform_with_lookup(doc_list, hc_dict, species_dict)


def transform_with_lookup(doc_list, hc_dict, species_dict):
    processes = int(os.environ.get("PUBTATOR_TRANSFORM_PROCESSES", "15"))
    chunksize = int(os.environ.get("PUBTATOR_TRANSFORM_CHUNKSIZE", "500"))
    identifier_resolution_cache = dict(_IDENTIFIER_RESOLUTION_CACHE)

    with Pool(
        processes,
        initializer=_init_transform_worker,
        initargs=(hc_dict, species_dict, identifier_resolution_cache),
    ) as pool:
        args = ((doc, index) for index, doc in enumerate(doc_list, start=1))
        for doc in pool.imap(process_document, args, chunksize=chunksize):
            yield doc


def update_table(conn, table_name, item_list, get_new_items_func):
    # Get the cursor for the SQLite connection and fetch all the data from the given table
    c = conn.cursor()
    c.execute(f"SELECT * FROM {table_name}")
    data = c.fetchall()

    # Initialize a list to store items not found in the data
    no_matches = []
    for original_name in item_list:
        # Try to find a match in the data using the lookup function
        found_match = lookup_item(original_name, data)
        # If no match is found and the item is not already in the no_matches list, add it
        if not found_match and original_name not in no_matches:
            continue
            logger.info(f"No {original_name} in lookup dictionary, saving for batch...")
            no_matches.append(original_name)

    # If there are any items in the no_matches list
    if no_matches:
        logger.info(f"Getting new {table_name}")
        # Call the get_new_items_func to fetch new data for the missing items
        no_results = get_new_items_func(no_matches)
        # If no_results is None, skip adding items to the table
        if no_results is None:
            logger.info(f"No results for {no_matches}, skipping adding to lookup dictionary")


def _extract_taxon_id(identifier):
    """Extract a numeric taxon ID from identifier string (e.g., 'taxonomy:NEWT:9606')."""
    if not identifier:
        return None
    parts = identifier.split(":")
    last_part = parts[-1].strip()
    if last_part.isdigit():
        return last_part
    return None


def _load_species_cache():
    """Load the species_details cache from the extract lookup DB."""
    try:
        with sqlite3.connect(EXTRACT_DB_PATH) as conn:
            c = conn.cursor()
            c.execute("CREATE TABLE IF NOT EXISTS species_details (original_name TEXT PRIMARY KEY, standard_dict TEXT)")
            c.execute("SELECT original_name, standard_dict FROM species_details")
            return {row[0].lower().strip(): json.loads(row[1]) for row in c.fetchall() if row[1]}
    except Exception as e:
        logger.error(f"text2term: error loading species cache: {e}")
        return {}


def _cache_species_details(species_details):
    """Cache a species_details entry in the extract lookup DB."""
    try:
        with sqlite3.connect(EXTRACT_DB_PATH) as conn:
            c = conn.cursor()
            c.execute("CREATE TABLE IF NOT EXISTS species_details (original_name TEXT PRIMARY KEY, standard_dict TEXT)")
            c.execute(
                "INSERT OR REPLACE INTO species_details VALUES (?, ?)",
                (species_details["originalName"].lower().strip(), json.dumps(species_details)),
            )
            conn.commit()
    except Exception as e:
        logger.error(f"text2term: error caching species details: {e}")


def _fetch_negative_species_cache():
    """Load the set of species names known to be unresolvable."""
    try:
        with sqlite3.connect(EXTRACT_DB_PATH) as conn:
            c = conn.cursor()
            c.execute(f"CREATE TABLE IF NOT EXISTS {_NEGATIVE_SPECIES_TABLE} (original_name TEXT PRIMARY KEY)")
            c.execute(f"SELECT original_name FROM {_NEGATIVE_SPECIES_TABLE}")
            return {row[0] for row in c.fetchall()}
    except Exception as e:
        logger.error(f"text2term: error loading negative species cache: {e}")
        return set()


def _cache_negative_species(original_name):
    """Record a species name that could not be resolved, to avoid retrying."""
    try:
        with sqlite3.connect(EXTRACT_DB_PATH) as conn:
            c = conn.cursor()
            c.execute(f"CREATE TABLE IF NOT EXISTS {_NEGATIVE_SPECIES_TABLE} (original_name TEXT PRIMARY KEY)")
            c.execute(
                f"INSERT OR IGNORE INTO {_NEGATIVE_SPECIES_TABLE} (original_name) VALUES (?)",
                (original_name.lower().strip(),),
            )
            conn.commit()
    except Exception as e:
        logger.error(f"text2term: error caching negative species: {e}")


def _get_uniprot_details(original_name, identifier, max_retries=3):
    """Fetch species details from UniProt and classify as host or infectiousAgent."""
    identifier = str(identifier).split("*")[-1]
    for attempt in range(max_retries):
        response = _UNIPROT_SESSION.get(f"https://rest.uniprot.org/taxonomy/{identifier}", timeout=30)
        if response.status_code == 429:
            retry_after = int(response.headers.get("Retry-After", 2 ** attempt))
            logger.warning(f"text2term: UniProt 429 for {identifier}, retrying in {retry_after}s")
            time.sleep(retry_after)
            continue
        response.raise_for_status()
        break
    else:
        raise requests.exceptions.HTTPError(f"UniProt rate limit exceeded after {max_retries} retries for {identifier}")
    species_info = response.json()

    standard_dict = {
        "@type": "DefinedTerm",
        "identifier": identifier,
        "inDefinedTermSet": "UniProt",
        "url": f"https://www.uniprot.org/taxonomy/{identifier}",
        "originalName": original_name,
        "isCurated": False,
    }

    if scientific_name := species_info.get("scientificName"):
        standard_dict["name"] = scientific_name
    else:
        standard_dict["name"] = original_name

    alternative_names = []
    if common_name := species_info.get("commonName"):
        standard_dict["commonName"] = common_name
        alternative_names.append(common_name)
        standard_dict["displayName"] = f"{common_name} | {standard_dict['name']}"
    else:
        standard_dict["displayName"] = standard_dict["name"]

    if other_names := species_info.get("otherNames"):
        alternative_names.extend(other_names)
    if alternative_names:
        standard_dict["alternateName"] = list(set(alternative_names))

    if lineage := species_info.get("lineage"):
        standard_dict["classification"] = classify_as_host_or_agent(lineage)
        standard_dict["lineage"] = lineage
    else:
        logger.warning(f"No lineage found for {identifier}")
        standard_dict["classification"] = "infectiousAgent"

    return standard_dict


def _apply_resolved_species(doc, resolved):
    """Replace unstandardized species/infectiousAgent entries with resolved ones."""
    has_unstandardized = False
    for field in ("species", "infectiousAgent"):
        entries = doc.get(field, [])
        if isinstance(entries, dict):
            entries = [entries]
        if isinstance(entries, list):
            for entry in entries:
                if isinstance(entry, dict) and "name" in entry and "inDefinedTermSet" not in entry and "curatedBy" not in entry:
                    has_unstandardized = True
                    break
        if has_unstandardized:
            break

    if not has_unstandardized:
        return

    new_species = []
    new_infectious_agents = []
    seen_identifiers = set()

    all_entries = []
    for field in ("species", "infectiousAgent"):
        entries = doc.get(field, [])
        if isinstance(entries, dict):
            entries = [entries]
        if isinstance(entries, list):
            for entry in entries:
                if isinstance(entry, dict):
                    all_entries.append((field, entry))

    for original_field, entry in all_entries:
        name = entry.get("name")
        if not name:
            continue

        identifier = entry.get("identifier")

        if "inDefinedTermSet" in entry or "curatedBy" in entry:
            if identifier and identifier in seen_identifiers:
                continue
            if identifier:
                seen_identifiers.add(identifier)
            if original_field == "infectiousAgent" or entry.get("classification") == "infectiousAgent":
                new_infectious_agents.append(entry)
            else:
                new_species.append(entry)
            continue

        key = name.lower().strip()
        if key in resolved:
            std = resolved[key].copy()
            std.pop("lineage", None)
            std["originalName"] = name
            new_id = std.get("identifier")
            if new_id and new_id in seen_identifiers:
                continue
            if new_id:
                seen_identifiers.add(new_id)
            if std.get("classification") == "infectiousAgent":
                new_infectious_agents.append(std)
            else:
                new_species.append(std)
        else:
            if original_field == "infectiousAgent":
                new_infectious_agents.append(entry)
            else:
                new_species.append(entry)

    if new_species:
        doc["species"] = new_species
    else:
        doc.pop("species", None)

    if new_infectious_agents:
        doc["infectiousAgent"] = new_infectious_agents


def _normalize_unstandardized_species(unstandardized):
    normalized = {}
    for key, value in unstandardized.items():
        if isinstance(value, tuple):
            name, identifier = value
        else:
            name, identifier = key, value
        normalized_key = str(name).lower().strip()
        if normalized_key:
            normalized.setdefault(normalized_key, (name, identifier))
    return normalized


def resolve_species_terms(unstandardized):
    """Resolve unique species names using the local cache, UniProt IDs, then text2term."""
    unstandardized = _normalize_unstandardized_species(unstandardized)

    if not unstandardized:
        logger.info("text2term: all species already standardized")
        return {}

    logger.info(f"text2term: {len(unstandardized)} unstandardized species to resolve")

    species_cache = _load_species_cache()
    negative_cache = _fetch_negative_species_cache()

    resolved = {}
    need_uniprot = []  # (name, taxon_id) tuples for direct UniProt lookup
    need_text2term = []
    negative_count = 0

    for key, (name, identifier) in unstandardized.items():
        # Skip if previously confirmed unresolvable
        if key in negative_cache:
            negative_count += 1
            continue

        # Check positive cache first
        if key in species_cache:
            resolved[key] = species_cache[key]
            continue

        # Try to extract taxon ID from identifier for direct UniProt lookup
        taxon_id = _extract_taxon_id(identifier)
        if taxon_id:
            need_uniprot.append((name, taxon_id))
        else:
            need_text2term.append(name)

    logger.info(
        "text2term: %s cached, %s negative-cached, %s need UniProt, %s need text2term",
        len(resolved), negative_count,
        len(need_uniprot), len(need_text2term),
    )

    def _resolve_one(original_name, taxon_id):
        """Resolve a single species via UniProt. Returns (key, details) or (key, None)."""
        key = original_name.lower().strip()
        try:
            details = _get_uniprot_details(original_name, taxon_id)
            return key, details
        except Exception as e:
            logger.warning(f"text2term: UniProt lookup failed for {original_name} (ID {taxon_id}): {e}")
            return key, None

    # Concurrent UniProt lookups for entries with known taxon IDs
    if need_uniprot:
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = {
                executor.submit(_resolve_one, name, tid): name
                for name, tid in need_uniprot
            }
            for future in as_completed(futures):
                key, details = future.result()
                if details:
                    resolved[key] = details
                    _cache_species_details(details)
                    species_cache[key] = details
                else:
                    _cache_negative_species(key)
                    negative_cache.add(key)
                    # Fall back to text2term for this term
                    need_text2term.append(futures[future])

    logger.info(f"text2term: {len(resolved)} resolved from cache/ID, {len(need_text2term)} need text2term")

    if need_text2term:
        try:
            import text2term

            if not os.path.exists("cache/ncbitaxon"):
                logger.info("text2term: building ncbitaxon ontology cache...")
                text2term.cache_ontology("https://purl.obolibrary.org/obo/ncbitaxon.owl", "ncbitaxon")

            t2t_results = text2term.map_terms(need_text2term, "ncbitaxon", use_cache=True)
            t2t_results.sort_values(["Source Term", "Mapping Score"], ascending=[True, False], inplace=True)
            t2t_results = t2t_results.drop_duplicates(subset=["Source Term"], keep="first")

            logger.info(f"text2term: mapped {len(t2t_results)} terms")

            # Collect text2term results for concurrent UniProt lookups
            t2t_lookups = []
            mapped_source_terms = set()
            for _, row in t2t_results.iterrows():
                identifier = row["Mapped Term CURIE"].split(":")[1]
                original_name = row["Source Term"]
                mapped_source_terms.add(original_name.lower().strip())
                t2t_lookups.append((original_name, identifier))

            # Cache negatives for terms text2term couldn't map at all
            for name in need_text2term:
                key = name.lower().strip()
                if key not in mapped_source_terms and key not in negative_cache:
                    _cache_negative_species(key)
                    negative_cache.add(key)

            with ThreadPoolExecutor(max_workers=10) as executor:
                futures = {
                    executor.submit(_resolve_one, name, tid): name
                    for name, tid in t2t_lookups
                }
                for future in as_completed(futures):
                    key, details = future.result()
                    if details:
                        resolved[key] = details
                        _cache_species_details(details)
                    else:
                        _cache_negative_species(key)
                        negative_cache.add(key)
        except Exception as e:
            logger.error(f"text2term: error during text2term resolution: {e}")

    logger.info(f"text2term: total resolved={len(resolved)}")
    return resolved


def _collect_unstandardized_species(docs, spool_file=None):
    unstandardized_species = {}
    for doc in docs:
        _collect_unstandardized_species_from_doc(doc, None, unstandardized_species)
        if spool_file is not None:
            spool_file.write(orjson.dumps(doc))
            spool_file.write(b"\n")
    return unstandardized_species


def _apply_resolved_species_to_stream(docs, resolved):
    for doc in docs:
        if resolved:
            _apply_resolved_species(doc, resolved)
        yield doc


def resolve_unstandardized_species(doc_list):
    """Resolve species not standardized by PubTator DB lookup using text2term + UniProt."""
    temp_path = None
    try:
        if _is_reiterable(doc_list):
            source_factory = lambda: iter(doc_list)
            unstandardized = _collect_unstandardized_species(source_factory())
        else:
            with tempfile.NamedTemporaryFile(prefix="pubtator-resolve-species-", suffix=".ndjson", delete=False) as f:
                temp_path = f.name
                unstandardized = _collect_unstandardized_species(doc_list, spool_file=f)
            source_factory = lambda: _iter_ndjson_file(temp_path)

        resolved = resolve_species_terms(unstandardized)
        yield from _apply_resolved_species_to_stream(source_factory(), resolved)
    finally:
        if temp_path:
            try:
                os.remove(temp_path)
            except FileNotFoundError:
                pass


def _ensure_lookup_db():
    if not os.path.exists(DB_PATH) or os.stat(DB_PATH).st_size == 0:
        logger.info("No lookup dictionary found, creating new one")
    db_dir = os.path.dirname(DB_PATH)
    if db_dir:
        os.makedirs(db_dir, exist_ok=True)
    with sqlite3.connect(DB_PATH) as conn:
        c = conn.cursor()
        c.execute(
            """CREATE TABLE IF NOT EXISTS health_conditions
                        (original_name text, standard_dict text)"""
        )
        c.execute(
            """CREATE TABLE IF NOT EXISTS species
                        (original_name text, standard_dict text)"""
        )


def update_lookup_dict(
    health_conditions_list,
    species_list,
    infectious_agents_list,
    doc_list,
    unstandardized_species=None,
    hc_dict=None,
    species_dict=None,
):
    _ensure_lookup_db()
    if hc_dict is None or species_dict is None:
        hc_dict, species_dict = fetch_data_from_db()

    source_factory = doc_list if callable(doc_list) else lambda: iter(doc_list)

    try:
        # Update the health_conditions and species tables with new data
        # update_table(conn, "health_conditions", health_conditions_list, get_new_health_conditions)
        # update_table(conn, "species", species_list + infectious_agents_list, get_new_species)
        if unstandardized_species is None:
            unstandardized_species = _collect_unstandardized_species(source_factory())

        # Log the completion of updating the lookup dictionary and standardize the documents
        logger.info("Finished updating lookup dictionary")
        resolved_species = resolve_species_terms(unstandardized_species)
        logger.info("Standardizing document")
        transformed_docs = transform_with_lookup(source_factory(), hc_dict, species_dict)
        return _apply_resolved_species_to_stream(transformed_docs, resolved_species)
    except Exception as e:
        logger.error(f"An error occurred while updating lookup dictionary: {e}")
        raise


def add_to_drop_list(term_name, ncbi_tax_id):
    """
    Add a new term to the drop list for filtering.

    Args:
        term_name (str): The name of the term to filter
        ncbi_tax_id (str): The NCBI Taxonomy ID
    """
    DROP_LIST_TERMS[term_name.lower().strip()] = {
        "id": str(ncbi_tax_id)
    }
    logger.info(f"Added '{term_name}' to drop list")


def get_drop_list_summary():
    """
    Get a summary of all terms in the drop list.

    Returns:
        dict: Summary of drop list terms and their configurations
    """
    summary = {
        "total_terms": len(DROP_LIST_TERMS),
        "terms": DROP_LIST_TERMS
    }
    return summary
