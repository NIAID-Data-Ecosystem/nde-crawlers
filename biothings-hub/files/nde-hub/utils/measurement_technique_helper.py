import csv
import os

import orjson


def get_identifier(url):
    """
    Extracts the identifier from the given URL.
    """
    # Split the URL by '_'
    parts = url.split("_")
    # The identifier is the last part of the URL
    if parts[-1]:
        return parts[-1]
    else:
        return "0000"


# Read the mapping from the CSV file
def load_mapping(name):
    # csv_file = f'mt mappings - {name}.csv'
    csv_file = f"/data/nde-hub/standardizers/measurement_technique_lookup/{name}.csv"
    mapping = {}
    with open(csv_file, "r", newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            repo_technique = row["Repository Technique"]
            mapping[repo_technique] = {
                "name": row["Manually Mapped Term"],
                "inDefinedTermSet": row["Ontology"],
                "url": row["URL"],
                "identifier": get_identifier(row["URL"]),
                "isCurated": True,
            }
    return mapping


def process_measurement_technique(data, name):
    """
    Updates the measurementTechnique field in each document based on the provided mapping.
    """
    mapping = load_mapping(name)

    if isinstance(data, str):
        with open(os.path.join(data, "data.ndjson"), "rb") as f:
            doc_list = [orjson.loads(line) for line in f]
    else:
        doc_list = list(data)

    for doc in doc_list:
        if "measurementTechnique" in doc:
            mt = doc["measurementTechnique"]
            # Handle the case where measurementTechnique is a list
            if isinstance(mt, list):
                for idx, item in enumerate(mt):
                    original_name = item.get("name")
                    if original_name in mapping:
                        mt[idx] = mapping[original_name]
                        mt[idx]["originalName"] = original_name
            # Handle the case where measurementTechnique is an object
            elif isinstance(mt, dict):
                original_name = mt.get("name")
                if original_name in mapping:
                    mt = mapping[original_name]
                    mt["originalName"] = original_name

            doc["measurementTechnique"] = mt
    return doc_list
