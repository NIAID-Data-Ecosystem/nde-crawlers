import csv
import json
import os


def read_ndjson(file_path):
    docs = []
    with open(file_path, "r") as file:
        for line in file:
            docs.append(json.loads(line.strip()))
    return docs


def add_disambiguating_description(docs, source_name):
    """
    Adds 'disambiguatingDescription' to documents from a CSV file located in the /data/nde-hub/disambiguating_descriptions/ directory.

    :param docs: List of dictionaries or string specifying the path to a ndjson file containing documents.
    :param source_name: String specifying the name of the CSV file containing disambiguating descriptions.
    :return: List of documents, potentially updated with separate 'disambiguatingDescription' objects.
    """

    if isinstance(docs, str):
        file_path = os.path.join(docs, "data.ndjson")
        docs = read_ndjson(file_path)

    file_path = f"/data/nde-hub/disambiguating_descriptions/{source_name}.csv"
    with open(file_path, "r") as file:
        reader = csv.DictReader(file)
        disambiguating_descriptions = {row["_id"].lower(): row["Processed Summary"] for row in reader}

    updated_docs = []
    for doc in docs:
        doc_id = doc["_id"].lower()
        if doc_id in disambiguating_descriptions:
            doc["disambiguatingDescription"] = disambiguating_descriptions[doc_id]
        updated_docs.append(doc)

    return updated_docs
