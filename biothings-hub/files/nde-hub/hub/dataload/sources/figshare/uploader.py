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
    processed_docs = []

    for doc in documents:
        keywords = doc.get("keywords", [])
        unique_terms = set()  # To store tuples of (label, curie)
        remaining_keywords = []

        for keyword in keywords:
            mapping = next((m for m in mappings if m["Source Term"] == keyword), None)

            if mapping:
                # Check if the term is deprecated and GOOD
                if mapping.get("Decision") == "Good":
                    if "obsolete" in mapping.get("Tags", "").lower():
                        replacement = mapping.get("Consider")  # Replacement term
                        if replacement:
                            unique_terms.add((replacement, "Plant biology"))
                    else:
                        unique_terms.add((mapping["Mapped Term Label"], mapping["Mapped Term CURIE"]))
                else:
                    better_mapping = mapping.get("Better mapping")
                    if better_mapping == "Ignore":
                        remaining_keywords.append(keyword)
                    else:
                        unique_terms.add((better_mapping, mapping["Mapped Term CURIE"]))
            else:
                # Handle unmapped terms
                if not contains_special_characters(keyword) and "years" not in keyword.lower():
                    remaining_keywords.append(keyword)

        # Convert the unique terms into DefinedTerm objects
        topic_categories = [
            create_defined_term(label, curie)
            for (label, curie) in unique_terms
        ]

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
                count += 1
                if count % 1000 == 0:
                    logging.info(f"Processed {count} documents")
    else:
        doc_list = list(data_folder)

    mappings = load_mapping_sheet_from_csv(mapping_file)
    processed_documents = process_documents(doc_list, mappings)
    for doc in processed_documents:
        # Yield only if `topicCategory` exists and has a `name` field
        if "topicCategory" in doc and any("name" in category for category in doc["topicCategory"]):
            yield doc
