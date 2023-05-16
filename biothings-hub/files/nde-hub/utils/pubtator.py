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


MANUAL_HEALTH_CONDITIONS = [
    {
        "identifier": "D000086382",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D000086382",
        "originalName": "covid-19",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02",
        },
        "name": "COVID-19",
        "alternateName": [
            "2019 Novel Coronavirus Disease",
            "2019 Novel Coronavirus Infection",
            "2019-nCoV Disease",
            "2019-nCoV Infection",
            "COVID-19 Pandemic",
            "COVID-19 Pandemics",
            "COVID-19 Virus Disease",
            "COVID-19 Virus Infection",
            "COVID19",
            "Coronavirus Disease 2019",
            "Coronavirus Disease-19",
            "SARS Coronavirus 2 Infection",
            "SARS-CoV-2 Infection",
            "Severe Acute Respiratory Syndrome Coronavirus 2 Infection",
        ],
    },
    {
        "identifier": "D000086382",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D000086382",
        "originalName": "sars-cov-2",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02",
        },
        "name": "COVID-19",
        "alternateName": [
            "2019 Novel Coronavirus Disease",
            "2019 Novel Coronavirus Infection",
            "2019-nCoV Disease",
            "2019-nCoV Infection",
            "COVID-19 Pandemic",
            "COVID-19 Pandemics",
            "COVID-19 Virus Disease",
            "COVID-19 Virus Infection",
            "COVID19",
            "Coronavirus Disease 2019",
            "Coronavirus Disease-19",
            "SARS Coronavirus 2 Infection",
            "SARS-CoV-2 Infection",
            "Severe Acute Respiratory Syndrome Coronavirus 2 Infection",
        ],
    },
    {
        "identifier": "D006262",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D006262",
        "originalName": "healthy",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02",
        },
        "name": "Health",
        "alternateName": ["Individual Health", "Normalcy", "Normality"],
    },
    {
        "identifier": "D000086002",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D000086002",
        "originalName": "malignant mesothelioma",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02",
        },
        "name": "Mesothelioma, Malignant",
        "alternateName": [
            "Malignant Mesothelioma",
            "Malignant Pleural Mesothelioma",
            "Mesothelioma, Malignant Pleural",
        ],
    },
    {
        "identifier": "D003863",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D003863",
        "originalName": "depression",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02",
        },
        "name": "Depression",
        "alternateName": ["Depressive Symptoms", "Emotional Depression"],
    },
    {
        "identifier": "D000074285",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D000074285",
        "originalName": "Smoker",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02",
        },
        "name": "Smokers",
        "alternateName": [
            "Smokers, Non-Tobacco Products",
            "Smokers, Tobacco",
            "Vapers",
        ],
    },
    {
        "identifier": "D003924",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D003924",
        "originalName": "Type 2 Diabetes Mellitus",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-12",
        },
        "name": "Diabetes Mellitus, Type 2",
        "alternateName": [
            "Diabetes Mellitus, Adult-Onset",
            "Diabetes Mellitus, Ketosis-Resistant",
            "Diabetes Mellitus, Maturity-Onset",
            "Diabetes Mellitus, Non Insulin Dependent",
            "Diabetes Mellitus, Non-Insulin-Dependent",
            "Diabetes Mellitus, Noninsulin Dependent",
            "Diabetes Mellitus, Noninsulin-Dependent",
            "Diabetes Mellitus, Slow-Onset",
            "Diabetes Mellitus, Stable",
            "Diabetes Mellitus, Type II",
            "MODY",
            "Maturity-Onset Diabetes Mellitus",
            "Maturity-Onset Diabetes",
            "NIDDM",
            "Noninsulin-Dependent Diabetes Mellitus",
            "Type 2 Diabetes",
            "Type 2 Diabetes Mellitus",
        ],
    },
    {
        "identifier": "D000078822",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D000078822",
        "originalName": "Ex-smoker",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-12",
        },
        "name": "Ex-Smokers",
        "alternateName": ["Exsmokers"],
    },
    {
        "identifier": "D006967",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D006967",
        "originalName": "allergy",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-12",
        },
        "name": "Hypersensitivity",
        "alternateName": ["Allergic Reaction", "Allergy"],
    },
]


def extract_values(doc_list, key):
    # Initialize an empty list to store the extracted values
    values_list = []

    # Iterate through each document in the doc_list
    for doc in doc_list:
        # If the value for the given key is a list, iterate through the items in the list
        if isinstance(doc.get(key, {}), list):
            for item in doc[key]:
                # Append the "name" field of the item to the values_list
                values_list.append(item["name"])

        # If the value for the given key is a dictionary, extract the "name" field and append it to the values_list
        elif item_name := doc.get(key, {}).get("name"):
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


def get_details(original_name, type, identifier):
    logger.info(f"Getting details for {original_name}")

    # Check if the entity is a disease
    if type == "Disease":
        # Fetch details from the MeSH API and create a standardized dictionary
        info = requests.get(f'https://id.nlm.nih.gov/mesh/lookup/details?descriptor={identifier.split(":")[1]}')
        try:
            info.raise_for_status()
            info = info.json()

            standard_dict = {
                "identifier": identifier.split(":")[1],
                "inDefinedTermSet": "MESH",
                "url": info["descriptor"],
                "originalName": original_name,
                "isCurated": True,
                "curatedBy": {
                    "name": "PubTator",
                    "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
                    "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
                },
            }
            if len(info["terms"]) > 0:
                standard_dict["name"] = info["terms"][0]["label"]
            else:
                logger.info(f"No offical name for {original_name}")
                standard_dict["name"] = original_name

            if len(info["terms"]) > 1:
                standard_dict["alternateName"] = [x["label"] for x in info["terms"][1:]]

            return standard_dict
        except requests.exceptions.HTTPError as e:
            logger.info(f"No information found for {original_name}, {identifier}")
            return None

    # Check if the entity is a species
    elif type == "Species":
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
                get_details(original_name, type, identifier)

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


def get_new_health_conditions(health_conditions, sql_data):
    # Split the health_conditions list into chunks of 1000
    chunks = [health_conditions[x : x + 1000] for x in range(0, len(health_conditions), 1000)]
    logger.info(f"Splitting into {len(chunks)} chunks")

    chunk_count = 0
    for chunk in chunks:
        chunk_count += 1
        logger.info(f"Processing chunk {chunk_count}")

        # Convert the chunk to a JSON string
        data = json.dumps(".    ".join(chunk))

        # Submit the data to the PubTator API for annotation
        submit_response = requests.post(
            "https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/submit/disease",
            data=data,
        )

        logger.info(f"Waiting for response, {submit_response.text}")
        timeout = 0
        retries = 0

        while True:
            # Wait for 10 seconds before checking the response
            time.sleep(10)
            timeout += 10

            # Retry the request if the timeout exceeds 100 seconds
            if timeout > 100:
                retries += 1
                # if retries > 3:
                #     raise Exception("Attempted 3 times, giving up")
                logger.info("Timeout, retrying...")
                try:
                    os.remove(f"{submit_response.text}_response.csv")
                except FileNotFoundError:
                    logger.info("Issue removing file: %s_response.csv", submit_response.text)
                get_new_health_conditions(health_conditions, sql_data)

            # Check the response from the PubTator API
            retrieve_response = requests.get(
                f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/retrieve/{submit_response.text}"
            )

            # If the response is successful, break the loop
            if retrieve_response.status_code == 200:
                logger.info("Got response")
                break

        # Process the retrieved result
        result = retrieve_response.text.split("00000|a|-NoAbstract-")[1].strip()

        if len(result.split("\n")) == 1:
            result = result + "\n" + result

        with open(f"{submit_response.text}_response.csv", "a") as f:
            for line in result.split("\n"):
                if "MESH" in line:
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

    # Update the database with new health conditions
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

        for original_name in remove_dupes:
            if original_name not in [x[0].lower().strip() for x in sql_data]:
                result_dict = get_details(
                    # name of the health condition / species
                    original_name,
                    # type (disease, species)
                    remove_dupes[original_name][0][0],
                    # identifier
                    remove_dupes[original_name][0][1],
                )
                conn = sqlite3.connect(DB_PATH)
                c = conn.cursor()
                c.execute(
                    "INSERT INTO health_conditions VALUES (?, ?)",
                    (original_name.lower(), json.dumps(result_dict)),
                )
                logger.info(f"Added {original_name}")
                conn.commit()
                conn.close()

        # compare what was submitted to what was added
        difference = list(set([x.lower().strip() for x in health_conditions]) - set(remove_dupes.keys()))
        logger.info(f"New health conditions: {difference}")
        try:
            os.remove(f"{submit_response.text}_response.csv")
        except FileNotFoundError:
            logger.info("Issue removing file: %s_response.csv", submit_response.text)
        not_found = []
        for submitted_health_condition in health_conditions:
            if submitted_health_condition.lower() not in remove_dupes.keys():
                not_found.append(submitted_health_condition)

        return (difference, not_found)
    else:
        try:
            os.remove(f"{submit_response.text}_response.csv")
        except FileNotFoundError:
            logger.info("Issue removing file: %s_response.csv", submit_response.text)

        logger.info("No new health conditions found")
        return None


def get_new_species(species, sql_data):
    # Split the species list into chunks of 1000
    chunks = [species[x : x + 1000] for x in range(0, len(species), 1000)]
    logger.info(f"Splitting into {len(chunks)} chunks")

    chunk_count = 0
    for chunk in chunks:
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

        while True:
            # Wait for 10 seconds before checking the response
            time.sleep(10)
            timeout += 10

            # Retry the request if the timeout exceeds 100 seconds
            if timeout > 100:
                retries += 1
                # if retries > 3:
                #     raise Exception("Attempted 3 times, giving up")
                logger.info("Timeout, retrying...")
                try:
                    os.remove(f"{submit_response.text}_response.csv")
                except FileNotFoundError:
                    logger.info("Issue removing file: %s_response.csv", submit_response.text)
                get_new_species(species, sql_data)

            # Check the response from the PubTator API
            retrieve_response = requests.get(
                f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/retrieve/{submit_response.text}"
            )

            # If the response is successful, break the loop
            if retrieve_response.status_code == 200:
                logger.info("Got response")
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
        for original_name in remove_dupes:
            if original_name not in [x[0].lower().strip() for x in sql_data]:
                result_dict = get_details(
                    # name of the health condition / species
                    original_name,
                    # type (disease, species)
                    remove_dupes[original_name][0][0],
                    # identifier
                    remove_dupes[original_name][0][1],
                )
                conn = sqlite3.connect(DB_PATH)
                c = conn.cursor()
                # Update the database with new species
                c.execute(
                    "INSERT INTO species VALUES (?, ?)",
                    (original_name.lower(), json.dumps(result_dict)),
                )
                logger.info(f"Added {original_name}")
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


def lookup_health_condition(original_name, health_conditions_data):
    # Iterate through the health conditions data
    for health_condition in health_conditions_data:
        # Compare the original_name with the health_condition name, ignoring case
        if original_name.lower().strip() == health_condition[0].lower().strip():
            # If the health_condition data is not None, return the JSON parsed data
            if health_condition[1] is not None:
                return json.loads(health_condition[1])

        # If the health_condition data is not None, proceed with the comparison
        if health_condition[1] is not None:
            # Load the JSON data for the health_condition
            hc_data = json.loads(health_condition[1])

            # Compare the original_name with the health_condition name in the JSON data, ignoring case
            if "name" in hc_data and original_name.lower().strip() == hc_data["name"].lower().strip():
                return hc_data

            # If the health_condition has alternate names in the JSON data
            if "alternateName" in hc_data:
                # Iterate through the alternate names
                for alternate_name in hc_data["alternateName"]:
                    # Compare the original_name with the alternate name, ignoring case
                    if original_name.lower().strip() == alternate_name.lower().strip():
                        return hc_data

    # If no matching health condition is found, return None
    return None


def lookup_species(original_name, species_data):
    # Iterate through the species data
    for species in species_data:
        # Compare the original_name with the species name, ignoring case
        if original_name.lower().strip() == species[0].lower().strip():
            # If the species data is not None, return the JSON parsed data
            if species[1] is not None:
                return json.loads(species[1])

        # If the species data is not None, proceed with the comparison
        if species[1] is not None:
            # Load the JSON data for the species
            sp_data = json.loads(species[1])

            # Compare the original_name with the species name in the JSON data, ignoring case
            if "name" in sp_data and original_name.lower().strip() == sp_data["name"].lower().strip():
                return sp_data

            # If the species has alternate names in the JSON data
            if "alternateName" in sp_data:
                # Iterate through the alternate names
                for alternate_name in sp_data["alternateName"]:
                    # Compare the original_name with the alternate name, ignoring case
                    if original_name.lower().strip() == alternate_name.lower().strip():
                        return sp_data

    # If no matching species is found, return None
    return None


def transform(doc_list):
    # Connect to the SQLite database and fetch health conditions and species data
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT * FROM health_conditions")
    health_conditions_data = c.fetchall()
    c.execute("SELECT * FROM species")
    species_data = c.fetchall()
    conn.close()

    def process_section(section, lookup_data, lookup_fn):
        new_section_list = []

        # If the section is a list, iterate through its items
        if isinstance(section, list):
            count = 0
            for original_obj in section:
                count += 1
                if original_name := original_obj.get("name"):
                    # Look up the object in the database using the provided lookup function
                    new_obj = lookup_fn(original_name, lookup_data)
                    # If the object is found, append it to the new_section_list
                    if new_obj:
                        new_section_list.append(new_obj)
                    # Otherwise, append the original object to the new_section_list
                    else:
                        new_section_list.append(original_obj)
        # If the section is not a list, process it as a single object
        elif original_name := section.get("name"):
            new_obj = lookup_fn(original_name, lookup_data)
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
        new_health_conditions_list = process_section(
            health_conditions_list, health_conditions_data, lookup_health_condition
        )
        new_species_list = process_section(species_list, species_data, lookup_species)
        new_infectious_agent_list = process_section(infectious_agent_list, species_data, lookup_species)

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
    # TODO iterate over cursor:
    data = c.fetchall()
    # [orignal_name, json_data]

    # Initialize a list to store items not found in the data
    no_matches = []
    for original_name in item_list:
        # Try to find a match in the data using the lookup function
        # TODO replace data with c
        found_match = lookup(original_name, data)
        # If no match is found and the item is not already in the no_matches list, add it
        if not found_match and original_name not in no_matches:
            logger.info(f"No {original_name} in lookup dictionary, saving for batch...")
            no_matches.append(original_name)

    # If there are any items in the no_matches list
    if no_matches:
        logger.info(f"Getting new {table_name}")
        # Call the get_new_items_func to fetch new data for the missing items
        no_results = get_new_items_func(no_matches, data)
        # If no_results is None, add the items to the table with no detailed information
        if no_results is None:
            logger.info(f"No results for {no_matches}, adding to lookup dictionary with no detailed information")
            for item in no_matches:
                c.execute(f"INSERT INTO {table_name} VALUES (?, ?)", (item.lower(), None))
            conn.commit()
        # If there are items in no_results[1], add them to the table with no detailed information
        elif len(no_results[1]):
            for item in no_results[1]:
                logger.info(f"No results for {item}, adding to lookup dictionary with no detailed information")
                c.execute(f"INSERT INTO {table_name} VALUES (?, ?)", (item.lower(), None))
                conn.commit()


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
        for health_condition in MANUAL_HEALTH_CONDITIONS:
            c.execute(
                f"INSERT INTO health_conditions VALUES ('{health_condition['originalName']}', '{json.dumps(health_condition)}')"
            )
            conn.commit()

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


def lookup(original_name, data):
    # Iterate through the provided data
    for item in data:
        # If the original_name matches the item's original_name, check if the item has details
        if original_name.lower().strip() == item[0].lower().strip():
            if item[1] is not None:
                logger.info(f"Found {original_name} in lookup dictionary")
                return True
            else:
                logger.info(f"{original_name} exists in lookup dictionary, but no detailed information")
                return True
        # If the item has details and alternate names, check if original_name matches any alternate names
        elif item[1] is not None and "alternateName" in item[1]:
            for alternate_name in json.loads(item[1])["alternateName"]:
                if original_name.lower().strip() == alternate_name.lower().strip():
                    logger.info(f"Found {original_name} in alternate names of {item[0]}")
                    return True
    # If no match is found, return False
    return False
