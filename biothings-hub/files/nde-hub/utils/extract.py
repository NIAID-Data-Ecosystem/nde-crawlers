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
    params = {"document": description, "entity_types": entity_type, "format": "tsv"}
    response = requests.get(base_url, params=params)
    return response.text


# Function to parse the TSV response from the EXTRACT API


def parse_tsv(ndeid, text_response):
    dictlist = []
    records = text_response.strip().split("\n")
    for record in records:
        if record:
            parts = record.split("\t")
            dictlist.append({"_id": ndeid, "extracted_text": parts[0], "entity_type": parts[1], "onto_id": parts[2]})
    return dictlist


def extract(doc_list):
    logger.info("Processing descriptions...")
    count = 0
    drop_list = {
        "tonga",
        "alabama",
        "nevada",
        "argentina",
        "namibia",
        "panama",
        "virginia",
        "bulgaria",
        "togo",
        "china",
        "serendip",
        "arizona",
        "california",
        "montana",
    }

    # conn = sqlite3.connect(DB_PATH)
    # c = conn.cursor()
    with sqlite3.connect(DB_PATH) as conn:
        c = conn.cursor()
        for doc in doc_list:
            count += 1
            if count % 100 == 0:
                logger.info(f"Extracted {count} documents")
            logger.info(f"Processing document {doc['_id']}...")

            for entity_type in ["-2", "-26"]:
                if entity_type == "-2" and ("species" in doc or "infectiousAgent" in doc):
                    logger.info(f"Species or infectiousAgent already present in document {doc['_id']}. Skipping...")
                    continue
                if entity_type == "-26" and "healthCondition" in doc:
                    logger.info(f"HealthCondition already present in document {doc['_id']}. Skipping...")
                    continue

                try:
                    response_text = get_cached_description(doc["_id"].lower(), entity_type, c)
                    if response_text is None:
                        # No row in the DB
                        logger.info(f"No cached entry for {doc['_id']}. Pinging API...")
                        response_text = query_extract_api(doc["description"], entity_type)
                        cache_description(doc["_id"].lower(), response_text, entity_type, c)
                    elif response_text == "":
                        # There is a row, but it's empty
                        logger.info(f"Cached entry for {doc['_id']} is empty. Skipping...")
                        continue
                    extracted_entities = parse_tsv(doc["_id"].lower(), response_text)
                    if extracted_entities:
                        logger.info(f"Found {len(extracted_entities)} entities of type {entity_type}")
                    for entity in extracted_entities:
                        extracted_text = entity["extracted_text"].lower()
                        if extracted_text in drop_list:
                            logger.info(f"Found {extracted_text} in drop list. Skipping...")
                            continue
                        if entity["entity_type"] == "-2":
                            doc.setdefault("species", [])
                            if not any(
                                s.get("name") == entity["extracted_text"]
                                or entity["extracted_text"] in s.get("alternateName", [])
                                for s in doc["species"]
                            ):
                                logger.info(f"Adding species {entity['extracted_text']} to document {doc['_id']}")
                                doc["species"].append(
                                    {"name": entity["extracted_text"], "identifier": entity["onto_id"]}
                                )
                        elif entity["entity_type"] == "-26":
                            doc.setdefault("healthCondition", [])
                            if not any(
                                d.get("name") == entity["extracted_text"]
                                or entity["extracted_text"] in d.get("alternateName", [])
                                for d in doc["healthCondition"]
                            ):
                                logger.info(f"Adding disease {entity['extracted_text']} to document {doc['_id']}")
                                doc["healthCondition"].append({"name": entity["extracted_text"]})
                except Exception as e:
                    logger.error(f"Error processing document {doc['_id']}: {e}")
    return doc_list


def get_cached_description(ndeid, entity_type, cursor):
    logger.info(f"Checking cache for {ndeid}...")
    table = "species" if entity_type == "-2" else "disease" if entity_type == "-26" else None
    if table:
        cursor.execute(f"SELECT text_response FROM {table} WHERE ndeid = ?", (ndeid,))
        result = cursor.fetchone()
        return result[0] if result else None


def cache_description(ndeid, text_response, entity_type, cursor):
    logger.info(f"Caching description for {ndeid}...")
    table = "species" if entity_type == "-2" else "disease" if entity_type == "-26" else None
    if table:
        cursor.execute(f"CREATE TABLE IF NOT EXISTS {table} (ndeid TEXT PRIMARY KEY, text_response TEXT)")
        cursor.execute(f"INSERT OR REPLACE INTO {table} VALUES (?, ?)", (ndeid, text_response))


def apply_heuristics(df):
    # Exclude terms with 3 or less characters
    # save any terms with 3 or less characters to csv for debugging
    # debug_df_3 = df[df["Source Term"].str.len() == 3]
    # debug_df_3.to_csv("3_char_terms.csv", mode="a", columns=["Source Term", "Mapping Score"])
    # debug_df_4 = df[df["Source Term"].str.len() == 4]
    # debug_df_4.to_csv("4_char_terms.csv", mode="a", columns=["Source Term", "Mapping Score"])

    df = df[df["Source Term"].str.len() > 3]

    filtered_df = pd.DataFrame()

    grouped = df.groupby("Source Term")

    # Rules for 4 and 5+ letter terms
    for name, group in grouped:
        if len(name) == 4:
            group = group[group["Mapping Score"] > 0.91]
        elif len(name) == 5:
            group = group[group["Mapping Score"] > 0.95]
        else:
            group = group.sort_values(by="Mapping Score", ascending=False).head(1)

        # Handling terms formatted as g. species
        group = group[
            group.apply(
                lambda row: not re.match(r"\b[A-Z].\s[^\s]+\b", row["Source Term"])
                or (row["Source Term"].split(".")[0] == row["Mapped Term Label"].split(" ")[0][0]),
                axis=1,
            )
        ]

        # Sort by 'Mapping Score' in descending order and keep the top entry
        group = group.sort_values(by="Mapping Score", ascending=False).head(1)

        filtered_df = pd.concat([filtered_df, group])

    return filtered_df


def build_species_lineage_info(species_list, species_mapping):
    lineage_info = {}
    for species in species_list:
        logger.info(f"Building lineage info for {species}")
        details = species_mapping.get(species.lower())
        if details and "lineage" in details:
            # Create a set of all ancestors for this species
            scientific_name = details["name"].lower()
            ancestors = set(item["scientificName"].lower() for item in details.pop("lineage"))
            lineage_info[scientific_name] = ancestors
    return lineage_info


def filter_species_terms_for_ancestors(species_mapping, species_list, lineage_info):
    to_remove = set()

    if "mus" in species_list and "mus sp." in species_list:
        logger.info("Found both 'mus' and 'mus sp.' in species list. Removing 'mus sp.'...")
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
                logger.info(f"Unable to find other species details for {other_species}")
                continue
            other_scientific_name = other_details["name"].lower()

            if other_species != species and scientific_name in lineage_info.get(other_scientific_name, set()):
                logger.info(f"Found {scientific_name} as ancestor of {other_scientific_name}. Removing...")
                to_remove.add(scientific_name)

    # Filter out the ancestor terms
    filtered_species = [
        species
        for species in species_list
        if species_mapping.get(species.lower())
        and species_mapping.get(species.lower())["name"].lower() not in to_remove
    ]
    return filtered_species


def insert_species(doc_list, species_mapping):
    for doc in doc_list:
        updated_species = []
        updated_infectious_agents = []

        # Collect names and identifiers from preserved_species to avoid duplicates
        existing_species_names = set()
        existing_infectious_agent_names = set()
        existing_identifiers = set()

        if "species" in doc:
            species_names = [species_obj["name"] for species_obj in doc["species"]]
            lineage_info = build_species_lineage_info(species_names, species_mapping)
            filtered_species_names = filter_species_terms_for_ancestors(species_mapping, species_names, lineage_info)

            # Preserve curated species
            preserved_species = [species for species in doc["species"] if species.get("fromPMID", False)]

            # Initialize existing_* sets with preserved_species data
            for species in preserved_species:
                classification = species.get("classification")
                name_lower = species["name"].lower()
                identifier = species.get("identifier")
                if classification == "host":
                    existing_species_names.add(name_lower)
                    existing_identifiers.add(identifier)
                elif classification == "infectiousAgent":
                    existing_infectious_agent_names.add(name_lower)
                    existing_identifiers.add(identifier)

            # Process species objects that are not curated
            for species_obj in doc["species"]:
                if species_obj["name"] not in filtered_species_names:
                    continue

                if species_obj in preserved_species:
                    updated_species.append(species_obj)  # Keep curated species unchanged
                    continue

                original_name = species_obj["name"]
                lower_name = original_name.lower()  # Used for case-insensitive comparison
                new_obj = species_mapping.get(lower_name, None)  # Get new data or keep existing

                if not new_obj:
                    logger.info(f"Unable to find species details for {original_name}")
                    continue

                # Check classification and categorize accordingly
                if (
                    new_obj.get("classification") == "infectiousAgent"
                    and lower_name not in existing_infectious_agent_names
                    and new_obj.get("identifier") not in existing_identifiers
                ):
                    updated_infectious_agents.append(new_obj)
                    existing_infectious_agent_names.add(lower_name)
                    existing_identifiers.add(new_obj.get("identifier"))
                elif (
                    new_obj.get("classification") == "host"
                    and lower_name not in existing_species_names
                    and new_obj.get("identifier") not in existing_identifiers
                ):
                    updated_species.append(new_obj)
                    existing_species_names.add(lower_name)
                    existing_identifiers.add(new_obj.get("identifier"))
                else:
                    logger.info(f"Unable to add {original_name} to document {doc['_id']}")
                    logger.info(f"Classification: {new_obj.get('classification')}")
                    logger.info(f"Lower name: {lower_name}")
                    logger.info(f"Identifiers: {existing_identifiers}")
                    logger.info(f"identifier: {new_obj.get('identifier')}")

            # Combine curated species with updated species
            doc["species"] = preserved_species + updated_species

            if updated_infectious_agents:
                doc["infectiousAgent"] = updated_infectious_agents
            else:
                doc.pop("infectiousAgent", None)

            if not doc["species"]:
                doc.pop("species", None)

            logger.info(f"Updated species in document {doc['_id']}")

    return doc_list


def insert_disease(doc_list, disease_mapping):
    for doc in doc_list:
        if "healthCondition" in doc:
            updated_disease = []

            # Preserve curated diseases
            preserved_diseases = [disease for disease in doc["healthCondition"] if disease.get("fromPMID", False)]

            # Track names and identifiers of preserved diseases to avoid duplication
            preserved_disease_names = {disease["name"].lower() for disease in preserved_diseases}
            preserved_disease_identifiers = {
                disease.get("identifier") for disease in preserved_diseases if disease.get("identifier")
            }

            # Process diseases that are not curated
            for disease_obj in doc["healthCondition"]:
                original_name = disease_obj["name"].lower()
                identifier = disease_obj.get("identifier")

                # Skip if disease is already preserved by name or identifier
                if original_name in preserved_disease_names or (
                    identifier and identifier in preserved_disease_identifiers
                ):
                    continue  # Skip adding this disease since it's already preserved

                new_obj = disease_mapping.get(original_name, None)
                if new_obj:
                    updated_disease.append(new_obj)
                else:
                    updated_disease.append(disease_obj)  # Keep existing disease if no new data

            # Combine curated diseases with updated diseases
            doc["healthCondition"] = preserved_diseases + updated_disease
            logger.info(f"Updated diseases in document {doc['_id']}")

    return doc_list


def classify_as_host_or_agent(scientific_name, lineage):
    """
    Classifies a species as either host or infectious agent based on its lineage.

    Parameters:
    - lineage (list): The lineage of the species.

    Returns:
    - new_classification (str): The classification of the species.
    """
    hosts = ["Deuterostomia", "Embryophyta", "Arthropoda", "Archaea", "Mollusca"]
    if scientific_name in hosts:
        logger.info(f"Found {scientific_name}, classifying as host")
        return "host"
    # Extracting scientific names for easy processing
    scientific_names = [item["scientificName"] for item in lineage]

    if "Viruses" in scientific_names:
        logger.info(f"Found Viruses in {scientific_names}, classifying as infectiousAgent")
        return "infectiousAgent"
    elif "Archaea" in scientific_names:
        logger.info(f"Found Archaea in {scientific_names}, classifying as host")
        return "host"
    elif "Mollusca" in scientific_names:
        logger.info(f"Found Mollusca in {scientific_names}, classifying as host")
        return "host"

    # Check for host species conditions
    if "Deuterostomia" in scientific_names:
        logger.info(f"Found Deuterostomia in {scientific_names}, classifying as host")
        return "host"
    elif "Embryophyta" in scientific_names and not any(
        parasite in scientific_names for parasite in ["Arceuthobium", "Cuscuta", "Orobanche", "Striga", "Phoradendron"]
    ):
        logger.info(f"Found Embryophyta in {scientific_names}, classifying as host")
        return "host"
    elif "Arthropoda" in scientific_names:
        if "Acari" in scientific_names:
            if "Ixodida" in scientific_names or "Ixodes" in scientific_names:
                logger.info(f"Found Ixodida or Ixodes in {scientific_names}, classifying as host")
                return "host"
            else:
                logger.info(f"Found Acari in {scientific_names}, classifying as infectiousAgent")
                return "infectiousAgent"
        else:
            logger.info(f"Found Arthropoda in {scientific_names}, classifying as host")
            return "host"
    else:
        # If not falling under the above host conditions, classify as infectiousAgent
        logger.info(f"Found {scientific_names}, classifying as infectiousAgent")
        return "infectiousAgent"


def get_species_details(original_name, identifier):
    logger.info(f"Getting details for {original_name}")

    # Fetch details from the UniProt API
    # make it a string to avoid errors
    identifier = str(identifier)
    identifier = identifier.split("*")[-1]

    # try:
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
        standard_dict["classification"] = classify_as_host_or_agent(standard_dict["name"], lineage)
        standard_dict["lineage"] = lineage
    else:
        logger.warning(f"No lineage found for {identifier}")
        standard_dict["classification"] = "infectiousAgent"

    return standard_dict


def process_species(doc_list):
    logger.info("Processing species...")
    species_names = set()  # Using a set to avoid duplicates

    # Fetch species from the SQLite database
    species_dict = fetch_species_from_db()

    for doc in doc_list:
        logger.info(f"Processing document {doc['_id']}...")
        if "species" in doc:
            if isinstance(doc["species"], dict):
                doc["species"] = [doc["species"]]
            species_in_doc = [s["name"] for s in doc["species"] if "isCurated" not in s]
            species_names.update(species_in_doc)
            logger.info(f"Found {len(species_in_doc)} new species to process")

    unique_species_names = list(species_names)
    logger.info(f"Total unique species to process: {len(unique_species_names)}")

    formatted_species = []

    # Process all species names with Text2Term at once
    if unique_species_names:
        try:
            # get current working directory
            logger.info("Processing species with Text2Term...")
            if not os.path.exists("cache/ncbitaxon"):
                text2term.cache_ontology("https://purl.obolibrary.org/obo/ncbitaxon.owl", "ncbitaxon")
            t2t_results = text2term.map_terms(unique_species_names, "ncbitaxon", use_cache=True)
            logger.info(f"Found {len(t2t_results)} results from Text2Term")
            t2t_results.sort_values(["Source Term", "Mapping Score"], ascending=[True, False], inplace=True)

            for index, row in t2t_results.iterrows():
                identifier = row["Mapped Term CURIE"].split(":")[1]  # CURIE format is 'NCBITAXON:ID'
                original_name = row["Source Term"]
                logger.info(f"Processing {original_name}")

                # Check if the species is in the SQLite database
                standardized_species = lookup_species_in_db(original_name, species_dict)
                if standardized_species:
                    logger.info(f"Found {original_name} in the database")
                    standardized_species["fromEXTRACT"] = True
                    standardized_species["originalName"] = original_name
                    formatted_species.append(standardized_species)
                else:
                    logger.info(f"Unable to find {original_name} in the database")
                    try:
                        # Get species details
                        species_details = get_species_details(original_name, identifier)
                        if species_details:
                            species_details["fromEXTRACT"] = True
                            formatted_species.append(species_details)
                            # Remove 'fromEXTRACT' before caching
                            species_details.pop("fromEXTRACT")
                            # Cache the result in the SQLite database for future use
                            cache_species_in_db(species_details)
                            # Update the in-memory species dictionary
                            species_dict[original_name.lower().strip()] = species_details
                    except Exception as e:
                        logger.info(f"An error occurred while processing {original_name}: {e}")
                        continue

        except Exception as e:
            logger.error(f"Error during Text2Term processing: {e}")
            return doc_list

    species_mapping = {
        sp["originalName"].lower() if "originalName" in sp else sp["name"].lower(): sp for sp in formatted_species
    }
    return insert_species(doc_list, species_mapping) if formatted_species else doc_list


def fetch_species_from_db():
    # conn = sqlite3.connect(DB_PATH)
    # c = conn.cursor()
    with sqlite3.connect(DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT original_name, standard_dict FROM species_details")
        species_cursor = c.fetchall()
        species_dict = {item[0].lower().strip(): json.loads(item[1]) for item in species_cursor if item[1]}
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
            (species_details["originalName"].lower().strip(), json.dumps(species_details)),
        )
        conn.commit()
    # conn.close()


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

                    existing_priority = ontology_priority.get(existing_ont, 999)
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
            disease_in_doc = [s["name"] for s in doc["healthCondition"] if "isCurated" not in s]
            disease_names.update(disease_in_doc)
            logger.info(f"Found {len(disease_in_doc)} new diseases to process")

    unique_disease_names = list(disease_names)
    logger.info(f"Total unique diseases to process: {len(unique_disease_names)}")

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
                    disease_dict[disease_name.lower().strip()] = disease_details
            except Exception as e:
                logger.info(f"An error occurred while processing {disease_name}: {e}")
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

    disease_dict = {item[0].lower().strip(): json.loads(item[1]) for item in disease_cursor if item[1]}
    return disease_dict


def lookup_disease_in_db(original_name, disease_dict):
    original_name_lower = original_name.lower().strip()
    return disease_dict.get(original_name_lower)


def cache_disease_in_db(disease_details):
    conn = sqlite3.connect(PUBTATOR_DB_PATH)
    c = conn.cursor()
    c.execute(
        "INSERT OR REPLACE INTO health_conditions VALUES (?, ?)",
        (disease_details["originalName"].lower().strip(), json.dumps(disease_details)),
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
    updated_docs = process_diseases(updated_docs)
    updated_docs = deduplicate_diseases(updated_docs)
    for doc in updated_docs:
        yield doc
