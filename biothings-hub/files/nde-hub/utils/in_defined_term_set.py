import csv
import datetime
import logging
import os
from io import StringIO

import orjson
import requests
from utils.pubtator import get_species_details, query_condition

curated_by = {
    "curatedBy": {
        "name": "Data Discovery Engine",
        "url": "https://discovery.biothings.io/",
        "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
    }
}


def process_csv_data(csv_content):
    # Prepare the CSV data for reading
    csv_file = StringIO(csv_content)
    csv_reader = csv.DictReader(csv_file)

    # Initialize a dictionary to hold the processed data
    properties = {}

    # Variable to keep track of the current property being processed
    current_property = None

    for row in csv_reader:
        if ".inDefinedTermSet" in row["Old button text"]:
            # Start of a new property section
            current_property = row["Old button text"].split(".")[0]
            properties[current_property] = {}
        elif row["base url (for mapping purposes)"] and current_property:
            if "http" not in row["base url (for mapping purposes)"]:
                continue
            # Row belongs to the current property
            base_url = row["base url (for mapping purposes)"]
            # The script assumes the 'InDefinedTermSet' field to hold the intended value
            value = row["InDefinedTermSet"]
            properties[current_property][base_url] = value

    return properties


def species_func(species, term_set):
    if term_set == "UniProt":
        identifer = species.get("url").split("_")[-1]
        species.update(get_species_details(species.get("name"), identifer))
        species.update(curated_by)
    return species


def health_condition_func(health_condition, term_set):
    term_sets = ["MeSH", "DOID", "NCIT", "MONDO"]
    if term_set in term_sets:
        result = query_condition(health_condition.get("name"))
        if result:
            health_condition.update(result)
            health_condition.update(curated_by)
    else:
        health_condition["inDefinedTermSet"] = "Other"
    return health_condition


def get_term_set_helper(url, property_dict):
    for prop, term_set in property_dict.items():
        if prop in url:
            return term_set
    else:
        return None


def de_duplicate_dicts(dict_list):
    unique_dicts = []
    seen = set()
    for d in dict_list:
        if isinstance(d, dict):  # Check if the item is a dictionary
            identifier = d.get("identifier")
            if identifier:
                if identifier not in seen:
                    seen.add(identifier)
                    unique_dicts.append(d)
        else:
            for obj in d:
                identifier = obj.get("identifier")
                if identifier:
                    if identifier not in seen:
                        seen.add(identifier)
                        unique_dicts.append(obj)    
    return unique_dicts


def get_in_defined_term_set(doc, properties_dict):
    nde_properties = {
        "species": species_func,
        "infectiousAgent": species_func,
        "healthCondition": health_condition_func,
        # TODO add other properties
    }

    # To handle species and infectiousAgent deliniation
    temp_results = {"species": [], "infectiousAgent": []}
    for nde_property, func in nde_properties.items():
        if doc_property := doc.get(nde_property):
            if isinstance(doc_property, list):
                for item in doc_property:
                    # if we find a url DDE gets curation
                    if url := item.get("url"):
                        if term_set := get_term_set_helper(url, properties_dict[nde_property]):
                            result = func(item, term_set)

                            # Store only species and infectiousAgent results temporarily
                            if nde_property in ["species", "infectiousAgent"]:
                                classification = result.get("classification")
                                if classification == "infectiousAgent":
                                    temp_results["infectiousAgent"].append(result)
                                else:
                                    temp_results["species"].append(result)

                            else:
                                # Update doc directly for other properties
                                index = doc_property.index(item)
                                doc_property[index] = result
                    else:
                        # if there is no url then keep as is for pubtator to curate
                        if nde_property in ["species", "infectiousAgent"]:
                            temp_results[nde_property].append(doc_property)

            else:
                # if we find a url DDE gets curation
                if url := doc_property.get("url"):
                    if term_set := get_term_set_helper(url, properties_dict[nde_property]):
                        result = func(doc_property, term_set)

                        # Store only species and infectiousAgent results temporarily
                        if nde_property in ["species", "infectiousAgent"]:
                            classification = result.get("classification")
                            if classification == "infectiousAgent":
                                temp_results["infectiousAgent"].append(result)
                            else:
                                temp_results["species"].append(result)

                        else:
                            # Update doc directly for other properties
                            doc[nde_property] = result
                else:
                    # if there is no url then keep as is for pubtator to curate
                    if nde_property in ["species", "infectiousAgent"]:
                        temp_results[nde_property].append(doc_property)

    # Update doc with species and infectiousAgent results
    for key, value in temp_results.items():
        if value:
            de_duped = de_duplicate_dicts(value)
            doc[key] = de_duped

    return doc


# def handle_dde_docs(data_folder):
#     csv_url = "https://docs.google.com/spreadsheets/d/107WVX39r_a6xBGZ_gCku0LBNRWWBmk9x7Dg53wj1SiI/export?format=csv"
#     response = requests.get(csv_url).text

#     properties = process_csv_data(response)

#     docs = []
#     if isinstance(data_folder, str):
#         # Read data from the file and process it
#         logging.info("Reading data from file...")
#         with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
#             count = 0
#             for line in f:
#                 count += 1
#                 if count % 1000 == 0:
#                     logging.info(f"Processed {count} lines")
#                 doc = orjson.loads(line)
#                 docs.append(doc)

#     for doc in docs:
#         get_in_defined_term_set(doc, properties)

#     return docs


def handle_dde_docs(data_folder):
    csv_url = "https://docs.google.com/spreadsheets/d/107WVX39r_a6xBGZ_gCku0LBNRWWBmk9x7Dg53wj1SiI/export?format=csv"
    response = requests.get(csv_url).text
    properties = process_csv_data(response)

    if isinstance(data_folder, str):
        # Read data from the file and process it
        logging.info("Reading data from file...")
        with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
            count = 0
            for line in f:
                count += 1
                if count % 1000 == 0:
                    logging.info(f"Processed {count} lines")
                doc = orjson.loads(line)
                get_in_defined_term_set(doc, properties)
                yield doc


# def get_in_defined_term_set_wrapper(func):
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         csv_url = "https://docs.google.com/spreadsheets/d/107WVX39r_a6xBGZ_gCku0LBNRWWBmk9x7Dg53wj1SiI/export?format=csv"
#         response = requests.get(csv_url).text

#         properties = process_csv_data(response)

#         gen = func(*args, **kwargs)
#         for doc in gen:
#             get_in_defined_term_set(doc, properties)
#             yield doc

# return wrapper
