import csv
import datetime
import json
import logging
import os
import re
import sqlite3

import orjson
import pandas as pd
import requests
import text2term

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

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


def cache_description(ndeid, text_response, entity_type):
    logger.info(f"Caching description for {ndeid}...")
    if entity_type == "-2":
        table = "species"
    elif entity_type == "-26":
        table = "disease"
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute(f"CREATE TABLE IF NOT EXISTS {table} (ndeid text, text_response text)")
    c.execute(f"INSERT INTO {table} VALUES (?, ?)", (ndeid, text_response))
    conn.commit()
    conn.close()


def get_cached_description(ndeid, entity_type):
    logger.info(f"Checking cache for {ndeid}...")
    # if db path doesn't exist, create it
    if not os.path.exists(DB_PATH):
        logger.info(f"Creating cache database at {DB_PATH}")
        os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute(f"CREATE TABLE IF NOT EXISTS species (ndeid text, text_response text)")
        c.execute(f"CREATE TABLE IF NOT EXISTS disease (ndeid text, text_response text)")
        conn.commit()
        conn.close()
    if entity_type == "-2":
        table = "species"
    elif entity_type == "-26":
        table = "disease"
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute(f"SELECT text_response FROM {table} WHERE ndeid = ?", (ndeid,))
    result = c.fetchone()
    conn.close()
    if result:
        logger.info(f"Found cached description for {ndeid}")
    else:
        logger.info(f"No cached description found for {ndeid}")
    return result[0] if result else None


def extract(doc_list):
    logger.info("Processing descriptions...")
    count = 0
    drop_list = [
        "Tonga",
        "Alabama",
        "Nevada",
        "Argentina",
        "Namibia",
        "Panama",
        "Virginia",
        "Bulgaria",
        "Togo",
        "China",
        "Serendip",
        "Arizona",
        "California",
        "Montana",
    ]
    # Normalize the drop list for case-insensitive comparison
    drop_list = set(map(str.lower, drop_list))
    for doc in doc_list:
        count += 1
        if count % 100 == 0:
            logger.info(f"Extracted {count} documents")
        logger.info(f"Processing document {doc['_id']}...")
        for entity_type in ["-2", "-26"]:  # -2 for species, -26 for diseases
            # Check conditions based on entity_type
            if entity_type == "-2" and ("species" in doc or "infectiousAgent" in doc):
                logger.info(f"Species or infectiousAgent already present in document {doc['_id']}. Skipping...")
                continue
            if entity_type == "-26" and "healthCondition" in doc:
                logger.info(f"HealthCondition already present in document {doc['_id']}. Skipping...")
                continue
            try:
                response_text = get_cached_description(doc["_id"], entity_type)
                if not response_text and response_text != "":
                    response_text = query_extract_api(doc["description"], entity_type)
                    # save response_text to sqlite3 cache
                    cache_description(doc["_id"], response_text, entity_type)
                extracted_entities = parse_tsv(doc["_id"], response_text)
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
                            doc["species"].append({"name": entity["extracted_text"], "identifier": entity["onto_id"]})
                    # For diseases
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


def apply_heuristics(df):
    # Exclude terms with 3 or less characters
    # save any terms with 3 or less characters to csv for debugging
    debug_df_3 = df[df["Source Term"].str.len() == 3]
    debug_df_3.to_csv("3_char_terms.csv", mode="a", columns=["Source Term", "Mapping Score"])
    debug_df_4 = df[df["Source Term"].str.len() == 4]
    debug_df_4.to_csv("4_char_terms.csv", mode="a", columns=["Source Term", "Mapping Score"])

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

        # Collect names for comparison to avoid duplicates
        existing_species_names = set()
        existing_infectious_agent_names = set()
        existing_identifiers = set()

        if "species" in doc:
            species_names = [species_obj["name"] for species_obj in doc["species"]]
            lineage_info = build_species_lineage_info(species_names, species_mapping)
            filtered_species_names = filter_species_terms_for_ancestors(species_mapping, species_names, lineage_info)
            doc["species"] = [species for species in doc["species"] if species["name"] in filtered_species_names]
            for species_obj in doc["species"]:
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
                    logger.info(f"Unable add {original_name} to document {doc['_id']}")
                    logger.info(f"Classification: {new_obj.get('classification')}")
                    logger.info(f"Lower name: {lower_name}")
                    logger.info(f"Identifiers: {existing_identifiers}")
                    logger.info(f"identifier: {new_obj.get('identifier')}")

            # Update the document with new species and infectious agents
            if updated_species:
                doc["species"] = updated_species
            else:
                doc.pop("species", None)

            if updated_infectious_agents:
                doc["infectiousAgent"] = updated_infectious_agents
            else:
                doc.pop("infectiousAgent", None)

            logger.info(f"Updated species in document {doc['_id']}")

    return doc_list


def insert_disease(doc_list, disease_mapping):
    for doc in doc_list:
        if "healthCondition" in doc:
            updated_disease = []
            for disease_obj in doc["healthCondition"]:
                original_name = disease_obj["name"]
                if original_name in disease_mapping:
                    updated_disease.append(disease_mapping[original_name])
                else:
                    updated_disease.append(disease_obj)  # keep existing disease if no new data
            doc["healthCondition"] = updated_disease
            logger.info(f"Updated disease in document {doc['_id']}")

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
        "isCurated": True,
        "fromEXTRACT": True,
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

    # if not os.path.exists("cache/ncbitaxon"):
    #     text2term.cache_ontology("https://purl.obolibrary.org/obo/ncbitaxon.owl", "ncbitaxon")

    for doc in doc_list:
        logger.info(f"Processing document {doc['_id']}...")
        if "species" in doc:
            # Filter out curated species and add to the set
            if isinstance(doc["species"], dict):
                doc["species"] = [doc["species"]]
            species_in_doc = [s["name"] for s in doc["species"] if not s.get("isCurated", False)]
            species_names.update(species_in_doc)
            logger.info(f"Found {len(species_in_doc)} new species to process")

    # Convert set to list for Text2Term processing
    unique_species_names = list(species_names)

    logger.info(f"Total unique species to process: {len(unique_species_names)}")

    # Process species names with Text2Term
    if unique_species_names:
        t2t_result = text2term.map_terms(unique_species_names, "ncbitaxon", use_cache=False)
        t2t_result.sort_values(["Source Term", "Mapping Score"], ascending=[True, False], inplace=True)

        # Apply heuristics
        t2t_result = apply_heuristics(t2t_result)
        logger.info(f"Total species after applying heuristics: {len(t2t_result)}")

        formatted_species = []
        for index, row in t2t_result.iterrows():
            identifier = row["Mapped Term CURIE"].split(":")[1]  # CURIE format is 'NCBITAXON:ID'
            original_name = row["Source Term"]
            try:
                species_details = get_species_details(original_name, identifier)
            except Exception as e:
                logger.info(f"An error occurred while processing {original_name}: {e}")
                continue
            if species_details:
                formatted_species.append(species_details)
            else:
                logger.info(f"Unable to find species details for {original_name} ({identifier})")

        species_mapping = {}
        for sp in formatted_species:
            lower_name = sp["originalName"].lower()
            species_mapping[lower_name] = sp

        return insert_species(doc_list, species_mapping)
    else:
        logger.info("No new species to process")
        return doc_list


def create_return_object(hit, alternate_names, original_name):
    ontology = hit["_id"].split(":")[0]
    identifier = hit["_id"].split(":")[1]

    standard_dict = {
        "identifier": hit["_id"].split(":")[1],
        "inDefinedTermSet": ontology,
        "isCurated": True,
        "fromEXTRACT": True,
        "name": hit["name"],
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


def process_synonyms(synonym_field):
    if isinstance(synonym_field, dict):
        # If the synonym field is a dictionary, return the exact synonyms.
        return synonym_field.get("exact", [])
    else:
        # If the synonym field is a list, filter it by 'EXACT'.
        return [syn.split('"')[1] for syn in synonym_field if "EXACT" in syn]


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


def handle_response(data, condition, base_url):
    if "hits" in data:
        for hit in data["hits"]:
            alternate_names = process_synonyms(hit.get("synonym", {}))
            if hit["name"].lower().strip() == condition.lower().strip() or any(
                name.lower().strip() == condition.lower().strip() for name in alternate_names
            ):
                logger.info(f"Found {condition} in ontology: {base_url.split('/')[-1]}")
                return create_return_object(hit, alternate_names, condition)
    return None


def query_condition(health_condition):
    BASE_URLS = [
        "https://biothings.ncats.io/mondo",
        "https://biothings.ncats.io/hpo",
        "https://biothings.ncats.io/doid",
        "https://biothings.ncats.io/ncit",
    ]
    logger.info(f'Querying for "{health_condition}"...')
    for base_url in BASE_URLS:
        try:
            url = f'{base_url}/query?q=name:("{health_condition}")&limit=1000'
            data = retry_request(url)
            result = handle_response(data, health_condition, base_url)
            if result is not None:
                return result

            # if not found in name, search in synonyms
            url = (
                f'{base_url}/query?q=synonym.exact:"{health_condition}"&limit=1000'
                if "hpo" in base_url
                else f'{base_url}/query?q=synonym:"{health_condition}"&limit=1000'
            )
            data = retry_request(url)
            result = handle_response(data, health_condition, base_url)
            if result is not None:
                return result
        except Exception as e:
            logger.info(f"An error occurred while querying {base_url}: {e}")
    logger.info(f"Unable to find {health_condition}")
    return None


def process_diseases(doc_list):
    logger.info("Processing diseases...")
    disease_names = set()  # Using a set to avoid duplicates

    for doc in doc_list:
        logger.info(f"Processing document {doc['_id']}...")
        if "healthCondition" in doc:
            if isinstance(doc["healthCondition"], dict):
                doc["healthCondition"] = [doc["healthCondition"]]
            # Filter out curated disease and add to the set
            disease_in_doc = [s["name"] for s in doc["healthCondition"] if not s.get("isCurated", False)]
            disease_names.update(disease_in_doc)
            logger.info(f"Found {len(disease_in_doc)} new disease to process")

    # Convert set to list for Text2Term processing
    unique_disease_names = list(disease_names)

    logger.info(f"Total unique disease to process: {len(unique_disease_names)}")

    # Process disease names with Translator KPIs
    if unique_disease_names:
        formatted_disease = []
        for disease_name in unique_disease_names:
            disease_details = query_condition(disease_name)
            if disease_details:
                formatted_disease.append(disease_details)
            else:
                logger.info(f"Unable to find disease details for {disease_name}")

        disease_mapping = {sp["originalName"]: sp for sp in formatted_disease}
        return insert_disease(doc_list, disease_mapping)
    else:
        logger.info("No new disease to process")
        return doc_list


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
    for doc in updated_docs:
        yield doc
