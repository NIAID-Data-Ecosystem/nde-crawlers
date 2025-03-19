import json
import os
import re
import sqlite3

import orjson
import pandas as pd
import requests
import text2term
from config import logger

from .pubtator import DB_PATH as PUBTATOR_DB_PATH, query_condition

DB_PATH = "/data/nde-hub/standardizers/extract_lookup/extract_lookup.db"


# Function to query the EXTRACT API
def query_extract_api(description, entity_type):
    base_url = "http://tagger.jensenlab.org/GetEntities"
    params = {"document": description,
              "entity_types": entity_type, "format": "tsv"}
    response = requests.get(base_url, params=params)
    return response.text

# Function to parse the TSV response from the EXTRACT API


def parse_tsv(ndeid, text_response):
    dictlist = []
    records = text_response.strip().split("\n")
    for record in records:
        if record:
            parts = record.split("\t")
            dictlist.append({
                "_id": ndeid,
                "extracted_text": parts[0],
                "entity_type": parts[1],
                "onto_id": parts[2]
            })
    return dictlist


def extract(doc_list):
    logger.info("Processing descriptions...")
    count = 0
    drop_list = {
        "tonga", "alabama", "nevada", "argentina", "namibia", "panama",
        "virginia", "bulgaria", "togo", "china", "serendip", "arizona",
        "california", "montana",
    }

    with sqlite3.connect(DB_PATH) as conn:
        c = conn.cursor()
        for doc in doc_list:
            count += 1
            if count % 100 == 0:
                logger.info(f"Extracted {count} documents")
            logger.info(f"Processing document {doc['_id']}...")

            for entity_type in ["-2", "-26"]:
                # For species (entity type "-2") skip if the record already contains species or infectiousAgent.
                if entity_type == "-2" and ("species" in doc or "infectiousAgent" in doc):
                    logger.info(
                        f"Species or infectiousAgent already present in document {doc['_id']}. Skipping species extraction..."
                    )
                    continue
                # For diseases (entity type "-26") skip if healthCondition is present.
                if entity_type == "-26" and "healthCondition" in doc:
                    logger.info(
                        f"HealthCondition already present in document {doc['_id']}. Skipping disease extraction..."
                    )
                    continue

                try:
                    response_text = get_cached_description(
                        doc["_id"].lower(), entity_type, c)
                    if response_text is None:
                        logger.info(f"No cached entry for {doc['_id']}. Pinging API...")
                        response_text = query_extract_api(doc["description"], entity_type)
                        cache_description(doc["_id"].lower(), response_text, entity_type, c)
                    elif response_text == "":
                        logger.info(f"Cached entry for {doc['_id']} is empty. Skipping...")
                        continue
                    extracted_entities = parse_tsv(
                        doc["_id"].lower(), response_text)
                    if extracted_entities:
                        logger.info(
                            f"Found {len(extracted_entities)} entities of type {entity_type}")
                    for entity in extracted_entities:
                        extracted_text = entity["extracted_text"].lower()
                        if extracted_text in drop_list:
                            logger.info(
                                f"Found {extracted_text} in drop list. Skipping...")
                            continue
                        if entity["entity_type"] == "-2":
                            # Add species from extraction and tag with fromEXTRACT=True
                            doc.setdefault("species", [])
                            if not any(
                                s.get("name") == entity["extracted_text"]
                                or entity["extracted_text"] in s.get("alternateName", [])
                                for s in doc["species"]
                            ):
                                logger.info(f"Adding species {entity['extracted_text']} to document {doc['_id']}")
                                doc["species"].append(
                                    {
                                        "name": entity["extracted_text"],
                                        "identifier": entity["onto_id"],
                                        "fromEXTRACT": True,
                                    }
                                )
                        elif entity["entity_type"] == "-26":
                            doc.setdefault("healthCondition", [])
                            if not any(
                                d.get("name") == entity["extracted_text"]
                                or entity["extracted_text"] in d.get("alternateName", [])
                                for d in doc["healthCondition"]
                            ):
                                logger.info(
                                    f"Adding disease {entity['extracted_text']} to document {doc['_id']}")
                                doc["healthCondition"].append(
                                    {"name": entity["extracted_text"]})
                except Exception as e:
                    logger.error(
                        f"Error processing document {doc['_id']}: {e}")
    return doc_list


def get_cached_description(ndeid, entity_type, cursor):
    logger.info(f"Checking cache for {ndeid}...")
    table = "species" if entity_type == "-2" else "disease" if entity_type == "-26" else None
    if table:
        cursor.execute(
            f"SELECT text_response FROM {table} WHERE ndeid = ?", (ndeid,))
        result = cursor.fetchone()
        return result[0] if result else None


def cache_description(ndeid, text_response, entity_type, cursor):
    logger.info(f"Caching description for {ndeid}...")
    table = "species" if entity_type == "-2" else "disease" if entity_type == "-26" else None
    if table:
        cursor.execute(
            f"CREATE TABLE IF NOT EXISTS {table} (ndeid TEXT PRIMARY KEY, text_response TEXT)")
        cursor.execute(
            f"INSERT OR REPLACE INTO {table} VALUES (?, ?)", (ndeid, text_response))


def apply_heuristics(df):
    # Exclude terms with 3 or less characters
    df = df[df["Source Term"].str.len() > 3]

    filtered_df = pd.DataFrame()
    grouped = df.groupby("Source Term")

    for name, group in grouped:
        if len(name) == 4:
            group = group[group["Mapping Score"] > 0.91]
        elif len(name) == 5:
            group = group[group["Mapping Score"] > 0.95]
        else:
            group = group.sort_values(
                by="Mapping Score", ascending=False).head(1)

        group = group[
            group.apply(
                lambda row: not re.match(
                    r"\b[A-Z].\s[^\s]+\b", row["Source Term"])
                or (row["Source Term"].split(".")[0] == row["Mapped Term Label"].split(" ")[0][0]),
                axis=1,
            )
        ]
        group = group.sort_values(by="Mapping Score", ascending=False).head(1)
        filtered_df = pd.concat([filtered_df, group])

    return filtered_df


def build_species_lineage_info(species_list, species_mapping):
    lineage_info = {}
    for species in species_list:
        logger.info(f"Building lineage info for {species}")
        details = species_mapping.get(species.lower())
        if details and "lineage" in details:
            scientific_name = details["name"].lower()
            ancestors = set(item["scientificName"].lower()
                            for item in details.pop("lineage"))
            lineage_info[scientific_name] = ancestors
    return lineage_info


def filter_species_terms_for_ancestors(species_mapping, species_list, lineage_info):
    to_remove = set()

    if "mus" in species_list and "mus sp." in species_list:
        logger.info(
            "Found both 'mus' and 'mus sp.' in species list. Removing 'mus sp.'...")
        to_remove.add("mus sp.")

    for species in species_list:
        details = species_mapping.get(species.lower())
        if not details:
            logger.info(f"Unable to find species details for {species}")
            continue
        logger.info(f"Building lineage info for {species}")
        scientific_name = details["name"].lower()

        for other_species in species_list:
            other_details = species_mapping.get(other_species.lower())
            if not other_details:
                logger.info(
                    f"Unable to find other species details for {other_species}")
                continue
            other_scientific_name = other_details["name"].lower()

            if other_species != species and scientific_name in lineage_info.get(other_scientific_name, set()):
                logger.info(
                    f"Found {scientific_name} as ancestor of {other_scientific_name}. Removing...")
                to_remove.add(scientific_name)

    filtered_species = [
        species
        for species in species_list
        if species_mapping.get(species.lower())
        and species_mapping.get(species.lower())["name"].lower() not in to_remove
    ]
    return filtered_species


def insert_species(doc_list, species_mapping):
    for doc in doc_list:
        # Lists for preserving curated entries and for adding updated (extracted) entries.
        preserved_hosts = []  # for curated host entries
        preserved_infectious_agents = []  # for curated infectiousAgent entries
        updated_hosts = []  # for extracted host entries (after mapping)
        updated_infectious_agents = []  # for extracted infectiousAgent entries (after mapping)

        # Sets for tracking duplicates.
        existing_host_names = set()
        existing_infectious_agent_names = set()
        existing_identifiers = set()

        # Combine entries from both "species" and "infectiousAgent" fields.
        combined_entries = []
        if "species" in doc and doc["species"]:
            if isinstance(doc["species"], list):
                combined_entries.extend(doc["species"])
            else:
                combined_entries.append(doc["species"])
        if "infectiousAgent" in doc and doc["infectiousAgent"]:
            if isinstance(doc["infectiousAgent"], list):
                combined_entries.extend(doc["infectiousAgent"])
            else:
                combined_entries.append(doc["infectiousAgent"])

        if combined_entries:
            # Build a list of names for lineage filtering.
            species_names = [entry.get("name") for entry in combined_entries if "name" in entry]
            lineage_info = build_species_lineage_info(species_names, species_mapping)
            filtered_species_names = filter_species_terms_for_ancestors(
                species_mapping, species_names, lineage_info)

            for entry in combined_entries:
                name = entry.get("name")
                from_extract = entry.get("fromEXTRACT", False)
                logger.info(f"Processing entry: name='{name}', fromEXTRACT={from_extract}")
                if not name:
                    logger.info("Skipping entry without a name.")
                    continue

                # If the entry is not curated and its name is not in the filtered list, skip it.
                if not entry.get("curatedBy") and name not in filtered_species_names:
                    logger.info(f"Skipping '{name}' because it is not in the filtered species list: {filtered_species_names}")
                    continue

                # Process curated entries (those that did NOT come from extraction)
                if not entry.get("fromEXTRACT", False):
                    classification = entry.get("classification")
                    lower_name = name.lower()
                    identifier = entry.get("identifier")
                    if classification == "host":
                        logger.info(f"Preserving host: {name}")
                        preserved_hosts.append(entry)
                        existing_host_names.add(lower_name)
                    elif classification == "infectiousAgent":
                        logger.info(f"Preserving infectiousAgent: {name}")
                        preserved_infectious_agents.append(entry)
                        existing_infectious_agent_names.add(lower_name)
                    existing_identifiers.add(identifier)
                else:
                    # Process extracted entries by looking up details in species_mapping.
                    original_name = name
                    lower_name = original_name.lower()
                    new_obj = species_mapping.get(lower_name)
                    if not new_obj:
                        logger.info(f"Unable to find species details for {original_name}")
                        continue
                    logger.info(f"Adding extracted entry for '{original_name}' with classification '{new_obj.get('classification')}'")
                    classification = new_obj.get("classification")
                    identifier = new_obj.get("identifier")
                    if classification == "infectiousAgent":
                        if lower_name not in existing_infectious_agent_names and identifier not in existing_identifiers:
                            updated_infectious_agents.append(new_obj)
                            existing_infectious_agent_names.add(lower_name)
                            existing_identifiers.add(identifier)
                        else:
                            logger.info(f"Duplicate or conflict for infectiousAgent {original_name} in document {doc['_id']}")
                    elif classification == "host":
                        if lower_name not in existing_host_names and identifier not in existing_identifiers:
                            updated_hosts.append(new_obj)
                            existing_host_names.add(lower_name)
                            existing_identifiers.add(identifier)
                        else:
                            logger.info(f"Duplicate or conflict for host {original_name} in document {doc['_id']}")
                    else:
                        logger.info(f"Unknown classification for {original_name} in document {doc['_id']}")

            # Reassign fields based on classification.
            all_hosts = preserved_hosts + updated_hosts
            all_infectious_agents = preserved_infectious_agents + updated_infectious_agents

            if all_hosts:
                doc["species"] = all_hosts
            else:
                doc.pop("species", None)

            if all_infectious_agents:
                doc["infectiousAgent"] = all_infectious_agents
            else:
                doc.pop("infectiousAgent", None)

            logger.info(f"Updated species in document {doc['_id']}")
    logger.info(f"After reassigning, doc[{doc['_id']}]: species={doc.get('species')}, infectiousAgent={doc.get('infectiousAgent')}")
    return doc_list

def insert_disease(doc_list, disease_mapping):
    for doc in doc_list:
        if "healthCondition" in doc:
            updated_disease = []

            # Preserve curated diseases.
            preserved_diseases = [disease for disease in doc["healthCondition"] if disease.get("fromPMID", False)]
            preserved_disease_names = {disease["name"].lower() for disease in preserved_diseases}
            preserved_disease_identifiers = {
                disease.get("identifier") for disease in preserved_diseases if disease.get("identifier")
            }

            # Process non-curated diseases.
            for disease_obj in doc["healthCondition"]:
                original_name = disease_obj["name"].lower()
                identifier = disease_obj.get("identifier")
                if original_name in preserved_disease_names or (
                    identifier and identifier in preserved_disease_identifiers
                ):
                    continue
                new_obj = disease_mapping.get(original_name, None)
                if new_obj:
                    updated_disease.append(new_obj)
                else:
                    updated_disease.append(disease_obj)

            doc["healthCondition"] = preserved_diseases + updated_disease
            logger.info(f"Updated diseases in document {doc['_id']}")

    return doc_list


def classify_as_host_or_agent(scientific_name, lineage):
    hosts = ["Deuterostomia", "Embryophyta", "Arthropoda", "Archaea", "Mollusca"]
    if scientific_name in hosts:
        logger.info(f"Found {scientific_name}, classifying as host")
        return "host"
    scientific_names = [item["scientificName"] for item in lineage]
    if "Viruses" in scientific_names:
        logger.info(
            f"Found Viruses in {scientific_names}, classifying as infectiousAgent")
        return "infectiousAgent"
    elif "Archaea" in scientific_names:
        logger.info(
            f"Found Archaea in {scientific_names}, classifying as host")
        return "host"
    elif "Mollusca" in scientific_names:
        logger.info(
            f"Found Mollusca in {scientific_names}, classifying as host")
        return "host"
    if "Deuterostomia" in scientific_names:
        logger.info(
            f"Found Deuterostomia in {scientific_names}, classifying as host")
        return "host"
    elif "Embryophyta" in scientific_names and not any(
        parasite in scientific_names for parasite in ["Arceuthobium", "Cuscuta", "Orobanche", "Striga", "Phoradendron"]
    ):
        logger.info(
            f"Found Embryophyta in {scientific_names}, classifying as host")
        return "host"
    elif "Arthropoda" in scientific_names:
        if "Acari" in scientific_names:
            if "Ixodida" in scientific_names or "Ixodes" in scientific_names:
                logger.info(
                    f"Found Ixodida or Ixodes in {scientific_names}, classifying as host")
                return "host"
            else:
                logger.info(
                    f"Found Acari in {scientific_names}, classifying as infectiousAgent")
                return "infectiousAgent"
        else:
            logger.info(
                f"Found Arthropoda in {scientific_names}, classifying as host")
            return "host"
    else:
        logger.info(f"Found {scientific_names}, classifying as infectiousAgent")
        return "infectiousAgent"


def get_species_details(original_name, identifier):
    logger.info(f"Getting details for {original_name}")
    identifier = str(identifier)
    identifier = identifier.split("*")[-1]
    species_info = requests.get(f"https://rest.uniprot.org/taxonomy/{identifier}")
    species_info.raise_for_status()
    species_info = species_info.json()
    standard_dict = {
        "identifier": identifier,
        "inDefinedTermSet": "UniProt",
        "url": f"https://www.uniprot.org/taxonomy/{identifier}",
        "originalName": original_name,
        "isCurated": False,
        "fromEXTRACT": True,
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
        standard_dict["displayName"] = scientific_name if scientific_name else original_name

    if other_names := species_info.get("otherNames"):
        alternative_names.extend(other_names)
    if alternative_names:
        standard_dict["alternateName"] = list(set(alternative_names))
    if lineage := species_info.get("lineage"):
        standard_dict["classification"] = classify_as_host_or_agent(
            standard_dict["name"], lineage)
        standard_dict["lineage"] = lineage
    else:
        logger.warning(f"No lineage found for {identifier}")
        standard_dict["classification"] = "infectiousAgent"
    return standard_dict


def process_species(doc_list):
    logger.info("Processing species and infectiousAgent...")
    # Fetch species from the SQLite database.
    species_dict = fetch_species_from_db()


    for doc in doc_list:
        # Normalize "species"
        if "species" in doc:
            if isinstance(doc["species"], dict):
                doc["species"] = [doc["species"]]
            elif isinstance(doc["species"], list):
                normalized = []
                for s in doc["species"]:
                    if isinstance(s, str):
                        normalized.append({"name": s})
                    else:
                        normalized.append(s)
                doc["species"] = normalized

        # Normalize "infectiousAgent"
        if "infectiousAgent" in doc:
            if isinstance(doc["infectiousAgent"], dict):
                doc["infectiousAgent"] = [doc["infectiousAgent"]]
            elif isinstance(doc["infectiousAgent"], list):
                normalized = []
                for agent in doc["infectiousAgent"]:
                    if isinstance(agent, str):
                        normalized.append({"name": agent})
                    else:
                        normalized.append(agent)
                doc["infectiousAgent"] = normalized

    # We want to process a document only if all entries in its species/infectiousAgent fields
    # came solely from EXTRACT (i.e. they all have "fromEXTRACT": True). If any entry is curated,
    # we skip standardization for that document.
    docs_to_standardize = []
    for doc in doc_list:
        combined = []
        if "species" in doc:
            combined.extend(doc["species"])
        if "infectiousAgent" in doc:
            combined.extend(doc["infectiousAgent"])
        if not combined:
            logger.info(f"Document {doc['_id']} has no species or infectiousAgent. Skipping standardization.")
            continue

        if any(not term.get("fromEXTRACT", False) for term in combined):
            logger.info(
                f"Document {doc['_id']} already contains curated species/infectiousAgent. Skipping standardization."
            )
            continue
        else:
            docs_to_standardize.append(doc)

    term_names = set()
    for doc in docs_to_standardize:
        if "species" in doc:
            names = [s["name"] for s in doc["species"] if s.get("fromEXTRACT", False)]
            term_names.update(names)
            logger.info(f"Document {doc['_id']} has {len(names)} species from EXTRACT.")
        if "infectiousAgent" in doc:
            names = [a["name"] for a in doc["infectiousAgent"] if a.get("fromEXTRACT", False)]
            term_names.update(names)
            logger.info(f"Document {doc['_id']} has {len(names)} infectiousAgent from EXTRACT.")

    unique_term_names = list(term_names)
    logger.info(f"Total unique species/infectiousAgent to process: {len(unique_term_names)}")

    formatted_species = []
    if unique_term_names:
        try:
            logger.info("Processing species/infectiousAgent with Text2Term...")
            if not os.path.exists("cache/ncbitaxon"):
                text2term.cache_ontology("https://purl.obolibrary.org/obo/ncbitaxon.owl", "ncbitaxon")
            t2t_results = text2term.map_terms(unique_term_names, "ncbitaxon", use_cache=True)
            logger.info(f"Found {len(t2t_results)} results from Text2Term")
            t2t_results.sort_values(["Source Term", "Mapping Score"], ascending=[True, False], inplace=True)

            for index, row in t2t_results.iterrows():
                identifier = row["Mapped Term CURIE"].split(":")[1]  # e.g., 'NCBITAXON:ID'
                original_name = row["Source Term"]
                logger.info(f"Processing species/infectiousAgent: {original_name}")

                # Look up in the in-memory species cache first.
                standardized_species = lookup_species_in_db(original_name, species_dict)
                if standardized_species:
                    logger.info(f"Found {original_name} in the database.")
                    standardized_species["fromEXTRACT"] = True
                    standardized_species["originalName"] = original_name
                    formatted_species.append(standardized_species)
                else:
                    logger.info(f"Unable to find {original_name} in the database, retrieving from external API...")
                    try:
                        species_details = get_species_details(original_name, identifier)
                        if species_details:
                            species_details["fromEXTRACT"] = True
                            formatted_species.append(species_details)
                            # Remove the flag before caching
                            species_details.pop("fromEXTRACT")
                            cache_species_in_db(species_details)
                            species_dict[original_name.lower().strip()] = species_details
                    except Exception as e:
                        logger.info(f"An error occurred while processing {original_name}: {e}")
                        continue
        except Exception as e:
            logger.error(f"Error during Text2Term processing: {e}")
            return doc_list


    species_mapping = {
        (sp["originalName"].lower() if "originalName" in sp else sp["name"].lower()): sp for sp in formatted_species
    }

    return insert_species(doc_list, species_mapping) if formatted_species else doc_list


def fetch_species_from_db():
    # conn = sqlite3.connect(DB_PATH)
    # c = conn.cursor()
    with sqlite3.connect(DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT original_name, standard_dict FROM species_details")
        species_cursor = c.fetchall()
        species_dict = {item[0].lower().strip(): json.loads(item[1])
                        for item in species_cursor if item[1]}
        return species_dict
    # conn.close()


def lookup_species_in_db(original_name, species_dict):
    original_name_lower = original_name.lower().strip()
    return species_dict.get(original_name_lower)


def cache_species_in_db(species_details):
    # conn = sqlite3.connect(DB_PATH)
    with sqlite3.connect(DB_PATH) as conn:
        c = conn.cursor()
        c.execute(
            "INSERT OR REPLACE INTO species_details VALUES (?, ?)",
            (species_details["originalName"].lower().strip(),
             json.dumps(species_details)),
        )
        conn.commit()
    # conn.close()


def deduplicate_species(species_list):
    """
    Deduplicate a list of species dictionaries based on the 'identifier' key.
    Keeps the first occurrence of each identifier.
    """
    seen = set()
    deduped = []
    for sp in species_list:
        identifier = sp.get('identifier')
        # If identifier exists and we haven't seen it, add it.
        if identifier is not None:
            if identifier not in seen:
                seen.add(identifier)
                deduped.append(sp)
        else:
            # If no identifier is found, simply include the record.
            deduped.append(sp)
    return deduped

def deduplicate_species_in_docs(doc_list):
    """
    Iterate over documents and deduplicate the 'species' field.
    """
    for doc in doc_list:
        if "species" in doc:
            doc["species"] = deduplicate_species(doc["species"])
    return doc_list


def remove_redundant_species(doc_list):
    """
    Remove species entries that are redundant when a curated infectiousAgent exists.

    For each document, if there is at least one curated infectiousAgent (i.e. where "isCurated" is True),
    then remove any species entry whose "name" matches the infectiousAgent name.
    """
    logger.info("Removing redundant species entries based on curated infectiousAgent...")
    for doc in doc_list:
        if "infectiousAgent" in doc and "species" in doc:
            curated_names = {
                ia["name"].strip().lower()
                for ia in doc["infectiousAgent"]
                if ia.get("isCurated", False)
            }
            if curated_names:
                original_species = doc["species"]
                # Filter out species that match any curated infectiousAgent name.
                doc["species"] = [
                    sp for sp in original_species
                    if sp.get("name", "").strip().lower() not in curated_names
                ]
                if len(original_species) != len(doc["species"]):
                    logger.info(f"Removed redundant species from document {doc['_id']}")
    return doc_list


ontology_priority = {"MONDO": 0, "HPO": 1, "DOID": 2, "NCIT": 3}


def deduplicate_diseases(doc_list):
    """
    Final pass to unify duplicate diseases. If multiple diseases share
    the same identifier (preferred) or the same lowercased name, only keep one.

    Priority rules:
      1. If one has fromPMID=True, keep that one.
      2. Otherwise compare ontology priorities (inDefinedTermSet).
      3. If still tied, keep the first one encountered.
    """
    logger.info("Deduplicating diseases...")
    for doc in doc_list:
        if "healthCondition" not in doc:
            continue

        unique_map = {}

        for disease in doc["healthCondition"]:
            key = disease.get("identifier") or disease["name"].lower()

            if key not in unique_map:
                unique_map[key] = disease
            else:
                existing = unique_map[key]

                existing_pmid = existing.get("fromPMID") == True
                new_pmid = disease.get("fromPMID") == True

                if existing_pmid and not new_pmid:
                    continue
                elif new_pmid and not existing_pmid:
                    unique_map[key] = disease
                else:
                    existing_ont = existing.get("inDefinedTermSet")
                    new_ont = disease.get("inDefinedTermSet")

                    existing_priority = ontology_priority.get(
                        existing_ont, 999)
                    new_priority = ontology_priority.get(new_ont, 999)

                    if new_priority < existing_priority:
                        unique_map[key] = disease
                    else:
                        continue

        # Rebuild healthCondition from unique_map values
        doc["healthCondition"] = list(unique_map.values())

    return doc_list


def process_diseases(doc_list):
    logger.info("Processing diseases...")
    disease_names = set()  # Using a set to avoid duplicates

    # Fetch diseases from the SQLite database
    disease_dict = fetch_diseases_from_db()

    for doc in doc_list:
        logger.info(f"Processing document {doc['_id']}...")
        if "healthCondition" in doc:
            if isinstance(doc["healthCondition"], dict):
                doc["healthCondition"] = [doc["healthCondition"]]
            # Filter out curated diseases and add to the set
            disease_in_doc = [s["name"]
                              for s in doc["healthCondition"] if "isCurated" not in s]
            disease_names.update(disease_in_doc)
            logger.info(f"Found {len(disease_in_doc)} new diseases to process")

    unique_disease_names = list(disease_names)
    logger.info(
        f"Total unique diseases to process: {len(unique_disease_names)}")

    formatted_diseases = []
    for disease_name in unique_disease_names:
        # Check if the disease is in the SQLite database
        standardized_disease = lookup_disease_in_db(disease_name, disease_dict)
        if standardized_disease:
            standardized_disease["fromEXTRACT"] = True
            standardized_disease["originalName"] = disease_name
            formatted_diseases.append(standardized_disease)
        else:
            try:
                # If not in the database, query the external API
                disease_details = query_condition(disease_name)
                if disease_details:
                    disease_details["fromEXTRACT"] = True
                    disease_details.pop("curatedBy")
                    formatted_diseases.append(disease_details)
                    disease_details.pop("fromEXTRACT")
                    # Cache the result in the SQLite database for future use
                    cache_disease_in_db(disease_details)
                    # Update the in-memory disease dictionary
                    disease_dict[disease_name.lower().strip()
                                 ] = disease_details
            except Exception as e:
                logger.info(
                    f"An error occurred while processing {disease_name}: {e}")
                continue

    disease_mapping = {
        sp["originalName"].lower() if "originalName" in sp else sp["name"].lower(): sp for sp in formatted_diseases
    }
    return insert_disease(doc_list, disease_mapping) if formatted_diseases else doc_list


def fetch_diseases_from_db():
    conn = sqlite3.connect(PUBTATOR_DB_PATH)
    c = conn.cursor()
    c.execute("SELECT original_name, standard_dict FROM health_conditions")
    disease_cursor = c.fetchall()
    conn.close()

    disease_dict = {item[0].lower().strip(): json.loads(item[1])
                    for item in disease_cursor if item[1]}
    return disease_dict


def lookup_disease_in_db(original_name, disease_dict):
    original_name_lower = original_name.lower().strip()
    return disease_dict.get(original_name_lower)


def cache_disease_in_db(disease_details):
    conn = sqlite3.connect(PUBTATOR_DB_PATH)
    c = conn.cursor()
    c.execute(
        "INSERT OR REPLACE INTO health_conditions VALUES (?, ?)",
        (disease_details["originalName"].lower().strip(),
         json.dumps(disease_details)),
    )
    conn.commit()
    conn.close()


def process_descriptions(data):
    logger.info("Standardizing data...")
    # Check if data is a file path (str)
    if isinstance(data, str):
        # Read data from the file and process it
        logger.info("Reading data from file...")
        with open(os.path.join(data, "data.ndjson"), "rb") as f:
            doc_list = []
            count = 0
            for line in f:
                count += 1
                if count % 1000 == 0:
                    logger.info(f"Processed {count} lines")
                doc = orjson.loads(line)
                doc_list.append(doc)
    else:
        # If data is a list, process it
        logger.info("Reading data from list...")
        doc_list = list(data)
        logger.info(f"Total documents to process: {len(doc_list)}")

    updated_docs = extract(doc_list)
    updated_docs = process_species(updated_docs)
    updated_docs = deduplicate_species_in_docs(updated_docs)
    updated_docs = remove_redundant_species(updated_docs)
    updated_docs = process_diseases(updated_docs)
    updated_docs = deduplicate_diseases(updated_docs)
    for doc in updated_docs:
        yield doc
