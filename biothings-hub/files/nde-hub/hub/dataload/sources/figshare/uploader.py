import csv
import logging
import os

import orjson
from hub.dataload.nde import NDESourceUploader
from utils.utils import nde_upload_wrapper

logging.basicConfig(level=logging.INFO)


def load_mapping_sheet_from_csv(csv_file):
    """Load mapping sheet from a CSV file and return as a list of dictionaries."""
    mappings = []
    with open(csv_file, "r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            mappings.append(row)
    return mappings


def process_documents(documents, mappings):
    count = 0
    """Process each document and update fields based on mappings."""
    processed_docs = []

    for doc in documents:
        count = 1
        if count % 1000 == 0:
            logging.info(f"Processed {count} documents")
        keywords = doc.get("keywords", [])
        topic_categories = []  # List of DefinedTerm objects
        remaining_keywords = []  # Keywords not mapped

        for keyword in keywords:
            mapping = next((m for m in mappings if m["Source Term"] == keyword), None)

            if mapping:
                # Check if the term is deprecated and GOOD
                if mapping.get("Decision") == "Good":
                    if "obsolete" in mapping.get("Tags", "").lower():
                        replacement = mapping.get("Consider")  # Replacement term
                        if replacement:
                            topic_categories.append(create_defined_term(replacement, "Plant biology"))
                    else:
                        topic_categories.append(
                            create_defined_term(mapping["Mapped Term Label"], mapping["Mapped Term CURIE"])
                        )
                elif mapping.get("Decision") != "Good":
                    better_mapping = mapping.get("Better mapping")
                    if better_mapping == "Ignore":
                        remaining_keywords.append(keyword)
                    else:
                        topic_categories.append(create_defined_term(better_mapping, mapping["Mapped Term CURIE"]))
            else:
                # Handle unmapped terms
                if not contains_special_characters(keyword) and "years" not in keyword.lower():
                    remaining_keywords.append(keyword)

        # Update document fields
        if topic_categories:
            doc["topicCategory"] = topic_categories
        if remaining_keywords:
            doc["keywords"] = remaining_keywords
        processed_docs.append(doc)

    return processed_docs


def create_defined_term(label, iri):
    """Create a DefinedTerm object."""
    return {
        "@type": "DefinedTerm",
        "name": label,
        "identifier": iri,
        "curatedBy": {"name": "Text2Term-assisted manual mapping"},
        "isCurated": True,
    }


def contains_special_characters(term):
    """Check if the term contains special characters."""
    return any(char in term for char in "!@#$%^&*()[]{};:,<>?/|\\~`")


class FigshareUploader(NDESourceUploader):
    name = "figshare"


@nde_upload_wrapper
def load_data(self, data_folder):
    mapping_file = "mappings.csv"
    doc_list = []
    count = 0
    if isinstance(data_folder, str):
        ndjson_file = os.path.join(data_folder, "data.ndjson")
        with open(ndjson_file, "rb") as f:
            for line in f:
                doc = orjson.loads(line)
                doc_list.append(doc)
                count = 1
                if count % 1000 == 0:
                    logging.info(f"Processed {count} documents")
    else:
        doc_list = list(data_folder)

    mappings = load_mapping_sheet_from_csv(mapping_file)
    processed_documents = process_documents(doc_list, mappings)
    for doc in processed_documents:
        yield doc
