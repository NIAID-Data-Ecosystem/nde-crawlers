import datetime
import json
import logging
import os
import sqlite3
import time

import orjson
import requests
from Bio import Entrez
from biothings.utils.dataload import tab2dict
from config import GEO_EMAIL, logger

from .date import add_date

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

DB_PATH = "/data/nde-hub/standardizers/pubtator_lookup/pubtator_lookup.db"


def extract_values(doc_list, key):
    # Initialize an empty list to store the extracted values
    values_list = []

    # Iterate through each document in the doc_list
    for doc in doc_list:
        # If the value for the given key is a list, iterate through the items in the list
        if isinstance(doc.get(key, {}), list):
            for item in doc[key]:
                # Check if the item has already been curated
                if "curatedBy" in item:
                    logger.info(f'{item["name"]} has already been curated, skipping...')
                    continue
                # Append the "name" field of the item to the values_list
                values_list.append(item["name"])

        # If the value for the given key is a dictionary, extract the "name" field and append it to the values_list
        elif item_name := doc.get(key, {}).get("name"):
            # Check if the item has already been curated
            if "curatedBy" in doc[key]:
                logger.info(f"{item_name} has already been curated, skipping...")
                continue
            values_list.append(item_name)

    # Remove duplicates and return a list of unique values with whitespace stripped
    return list(dict.fromkeys([x.lower().strip() for x in values_list]))


@add_date
def standardize_data(data):
    logger.info("Standardizing data...")
    # Check if data is a file path (str)
    if isinstance(data, str):
        # Read data from the file and process it
        logger.info("Reading data from file...")
        with open(os.path.join(data, "data.ndjson"), "rb") as f:
            health_conditions_list = []
            species_list = []
            infectious_agents_list = []
            doc_list = []
            count = 0
            for line in f:
                count += 1
                if count % 1000 == 0:
                    logger.info(f"Processed {count} lines")
                doc = orjson.loads(line)
                doc_list.append(doc)

            health_conditions_list = extract_values(doc_list, "healthCondition")
            species_list = extract_values(doc_list, "species")
            infectious_agents_list = extract_values(doc_list, "infectiousAgent")

            logger.info(
                f"Found {len(health_conditions_list)} health conditions, {len(species_list)} species, {len(infectious_agents_list)} infectious agents"
            )
            return update_lookup_dict(health_conditions_list, species_list, infectious_agents_list, doc_list)

    else:
        # If data is a list, process it
        logger.info("Reading data from list...")
        health_conditions_list = []
        species_list = []
        infectious_agents_list = []
        doc_list = list(data)

        health_conditions_list = extract_values(doc_list, "healthCondition")
        species_list = extract_values(doc_list, "species")
        infectious_agents_list = extract_values(doc_list, "infectiousAgent")

        logger.info(
            f"Found {len(health_conditions_list)} health conditions, {len(species_list)} species, {len(infectious_agents_list)} infectious agents"
        )
        # Call the update_lookup_dict function with the extracted data
        return update_lookup_dict(health_conditions_list, species_list, infectious_agents_list, doc_list)


def get_species_details(original_name, identifier):
    logger.info(f"Getting details for {original_name}")

    # Fetch details from the UniProt API or the NCBI Taxonomy API
    identifier = identifier.split("*")[-1]

    try:
        species_info = requests.get(f"https://rest.uniprot.org/taxonomy/{identifier}")
        species_info.raise_for_status()
        species_info = species_info.json()
        standard_dict = {
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

        if other_names := species_info.get("otherNames"):
            alternative_names.extend(other_names)

        if alternative_names:
            standard_dict["alternateName"] = alternative_names

        return standard_dict

    except requests.exceptions.HTTPError as e:
        logger.info(f"No Uniprot information found for {original_name}, {identifier}. Trying NCBI...")
    except requests.exceptions.SSLError as e:
        logger.info(f"SSL Error: {e}")
        logger.info(f"Retrying in 5 seconds...")
        time.sleep(5)
        get_species_details(original_name, identifier)

    Entrez.email = GEO_EMAIL

    while True:
        try:
            handle = Entrez.efetch(db="taxonomy", id=identifier, retmode="xml", max_tries=10)
            record = Entrez.read(handle, validate=False)
            handle.close()
            break
        except Exception as e:
            logger.info(f"Error: {e}")
            logger.info(f"Retrying in 5 seconds...")
            time.sleep(5)
            get_species_details(original_name, identifier)

    standard_dict = {
        "identifier": identifier,
        "inDefinedTermSet": "NCBI Taxonomy",
        "url": f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={identifier}",
        "originalName": original_name,
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
        },
    }
    if scientific_name := record[0].get("ScientificName"):
        standard_dict["name"] = scientific_name
    else:
        standard_dict["name"] = original_name

    alternative_names = []
    if common_name := record[0].get("OtherNames", {}).get("GenbankCommonName"):
        standard_dict["commonName"] = common_name
        alternative_names.append(common_name)

        standard_dict["displayName"] = f"{common_name} | {scientific_name}"

    if other_names := record[0].get("OtherNames", {}).get("Name"):
        for name_obj in other_names:
            if name_obj["ClassCDE"] == "authority":
                if name_obj["DispName"] not in alternative_names:
                    alternative_names.append(name_obj["DispName"])
    if alternative_names:
        standard_dict["alternateName"] = alternative_names

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


def process_synonyms(synonym_field):
    if isinstance(synonym_field, dict):
        # If the synonym field is a dictionary, return the exact synonyms.
        return synonym_field.get("exact", [])
    else:
        # If the synonym field is a list, filter it by 'EXACT'.
        return [syn.split('"')[1] for syn in synonym_field if "EXACT" in syn]


def get_xref_name(xref_ontology, xref_identifier):
    # Dictionary mapping each ontology to its API endpoint
    api_endpoints = {
        "ICD9": "https://clinicaltables.nlm.nih.gov/api/icd9cm_dx/v3/search?terms=",
        "MESH": "https://id.nlm.nih.gov/mesh/lookup/label?resource=",
        "SCTID": "https://www.ebi.ac.uk/ols4/api/v2/ontologies/snomed/classes/http%253A%252F%252Fsnomed.info%252Fid%252F",
        "EFO": "https://www.ebi.ac.uk/ols4/api/v2/ontologies/efo/classes/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_",
        "ICD10CM": "https://clinicaltables.nlm.nih.gov/api/icd10cm/v3/search?terms=",
    }

    # Check if the provided ontology is in our list of endpoints
    if xref_ontology not in api_endpoints:
        return None

    # Form the URL for the request
    url = api_endpoints[xref_ontology] + xref_identifier

    response = requests.get(url)

    response.raise_for_status()

    # Parse the JSON response and return the name
    try:
        if xref_ontology == "ICD9":
            return response.json()[4][1]
        elif xref_ontology == "MESH":
            return response.json()[0]
        elif xref_ontology == "SCTID" or xref_ontology == "EFO" or xref_ontology == "ICD10CM":
            return response.json()["label"]
    except KeyError:
        logger.info(f"Check the response from {url}")
        return None

    return response.json()["name"]  # This assumes that the response is a JSON object with a 'name' field


def create_return_object(hit, alternate_names, original_name):
    ontology = hit["_id"].split(":")[0]
    identifier = hit["_id"].split(":")[1]
    sameas_list = []
    if hit.get("xrefs"):
        for xref_ontology, xref_identifier in hit["xrefs"].items():
            sameas = {}
            sameas["identifier"] = f"{xref_ontology}:{xref_identifier}"
            sameas["url"] = f"http://purl.obolibrary.org/obo/{xref_ontology}_{xref_identifier}"
            # TODO - add name through api call
            sameas_list.append(sameas)
    elif hit.get("xref"):
        for xref in hit["xref"]:
            xref_ontology = xref.split(":")[0]
            xref_identifier = xref.split(":")[1]
            sameas = {}
            sameas["identifier"] = xref
            sameas["url"] = f"http://purl.obolibrary.org/obo/{xref_ontology}_{xref_identifier}"
            sameas_name = get_xref_name(xref_ontology, xref_identifier)
            if sameas_name:
                sameas["name"] = sameas_name
            sameas_list.append(sameas)

    standard_dict = {
        "identifier": hit["_id"],
        "inDefinedTermSet": ontology,
        "isCurated": True,
        "name": hit["name"],
        "originalName": original_name,
        "url": f"http://purl.obolibrary.org/obo/{ontology}_{identifier}",
        "curatedBy": {
            "name": "NCATS Biomedical Data Translator Program",
            "url": "https://biothings.ncats.io/",
            "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
        },
    }
    if sameas_list:
        standard_dict["sameas"] = sameas_list
    if alternate_names:
        standard_dict["alternateName"] = alternate_names
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
    # Split the species list into chunks of 1000
    chunks = [species[x : x + 1000] for x in range(0, len(species), 1000)]
    logger.info(f"Splitting into {len(chunks)} chunks")

    chunk_count = 0
    for chunk in chunks:
        while True:
            chunk_count += 1
            logger.info(f"Processing chunk {chunk_count}")

            # Convert the chunk to a JSON string
            data = json.dumps(".    ".join(chunk))

            # Submit the data to the PubTator API for annotation
            submit_response = requests.post(
                "https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/submit/species",
                data=data,
            )

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
                    if retries > 5:
                        raise Exception("Attempted 5 times, giving up")
                    logger.info("Timeout, retrying...")
                    try:
                        os.remove(f"{submit_response.text}_response.csv")
                    except FileNotFoundError:
                        logger.info("Issue removing file: %s_response.csv", submit_response.text)
                    # get_new_species(species)
                    break

                # Check the response from the PubTator API
                retrieve_response = requests.get(
                    f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/retrieve/{submit_response.text}"
                )

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
            conn = sqlite3.connect(DB_PATH)
            c = conn.cursor()
            # Update the database with new species
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


def lookup_item(original_name, data):
    # Iterate through the data
    for item in data:
        # Compare the original_name with the item name, ignoring case
        if original_name.lower().strip() == item[0].lower().strip():
            # If the item data is not None, return the JSON parsed data
            if item[1] is not None:
                return json.loads(item[1])

        # If the item data is not None, proceed with the comparison
        if item[1] is not None:
            # Load the JSON data for the item
            item_data = json.loads(item[1])

            # Compare the original_name with the item name in the JSON data, ignoring case
            if "name" in item_data and original_name.lower().strip() == item_data["name"].lower().strip():
                return item_data

            # If the item has alternate names in the JSON data
            if "alternateName" in item_data:
                # Iterate through the alternate names
                for alternate_name in item_data["alternateName"]:
                    # Compare the original_name with the alternate name, ignoring case
                    if original_name.lower().strip() == alternate_name.lower().strip():
                        return item_data

    # If no matching item is found, return None
    return None


def transform(doc_list):
    # Connect to the SQLite database and fetch health conditions and species data
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT * FROM health_conditions")
    hc_cursor = c.fetchall()
    c.execute("SELECT * FROM species")
    species_cursor = c.fetchall()
    conn.close()

    def process_section(section, cursor, lookup_fn):
        new_section_list = []

        # If the section is a list, iterate through its items
        if isinstance(section, list):
            count = 0
            for original_obj in section:
                count += 1
                if original_name := original_obj.get("name"):
                    # Look up the object in the database using the provided lookup function
                    new_obj = lookup_fn(original_name, cursor)
                    # If the object is found, append it to the new_section_list
                    if new_obj:
                        logger.info(f"Found {original_name} in database")
                        new_section_list.append(new_obj)
                    # Otherwise, append the original object to the new_section_list
                    else:
                        logger.info(f"Could not find {original_name} in database")
                        new_section_list.append(original_obj)
        # If the section is not a list, process it as a single object
        elif original_name := section.get("name"):
            new_obj = lookup_fn(original_name, cursor)
            if new_obj:
                new_section_list.append(new_obj)
            else:
                new_section_list.append(section)

        return new_section_list

    # Iterate through the documents in the doc_list
    for doc in doc_list:
        health_conditions_list = doc.get("healthCondition", {})
        species_list = doc.get("species", {})
        infectious_agent_list = doc.get("infectiousAgent", {})

        # Process each section in the document using the appropriate lookup function
        new_health_conditions_list = process_section(health_conditions_list, hc_cursor, lookup_item)
        new_species_list = process_section(species_list, species_cursor, lookup_item)
        new_infectious_agent_list = process_section(infectious_agent_list, species_cursor, lookup_item)

        # Update the document with the new lists
        if new_health_conditions_list:
            doc["healthCondition"] = new_health_conditions_list
        if new_species_list:
            doc["species"] = new_species_list
        if new_infectious_agent_list:
            doc["infectiousAgent"] = new_infectious_agent_list

        # Yield the transformed document
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


def update_lookup_dict(health_conditions_list, species_list, infectious_agents_list, doc_list):
    # If the database file does not exist or is empty, create a new one
    if not os.path.exists(DB_PATH) or os.stat(DB_PATH).st_size == 0:
        logger.info("No lookup dictionary found, creating new one")
        os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        # Create the health_conditions table and insert manually provided health conditions
        c.execute(
            """CREATE TABLE health_conditions
                        (original_name text, standard_dict text)"""
        )
        # Create the species table
        c.execute(
            """CREATE TABLE species
                        (original_name text, standard_dict text)"""
        )
    else:
        # If the database file exists, connect to it
        conn = sqlite3.connect(DB_PATH)

    try:
        # Update the health_conditions and species tables with new data
        update_table(conn, "health_conditions", health_conditions_list, get_new_health_conditions)
        update_table(conn, "species", species_list + infectious_agents_list, get_new_species)
        conn.close()

        # Log the completion of updating the lookup dictionary and standardize the documents
        logger.info("Finished updating lookup dictionary")
        logger.info("Standardizing document")
        transformed_docs = transform(doc_list)
        return transformed_docs
    except Exception as e:
        # Log the error, close the connection, and raise the exception
        logger.error(f"An error occurred while updating lookup dictionary: {e}")
        conn.close()
        raise
