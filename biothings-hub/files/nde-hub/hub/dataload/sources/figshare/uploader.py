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
        unique_names = set()  # Handle duplicates
        topic_terms = []
        remaining_keywords = []

        for keyword in keywords:
            # Normalize case for matching
            mapping = next((m for m in mappings if m["Source Term"].lower() == keyword.lower()), None)

            if mapping:
                # Check if the term is deprecated and GOOD
                if mapping.get("Decision") == "Good":
                    if "obsolete" in (mapping.get("Tags") or "").lower():
                        replacement = mapping.get("Consider")  # Replacement term
                        if replacement and replacement.lower() != "ignored":
                            if replacement not in unique_names:
                                unique_names.add(replacement)
                                topic_terms.append((replacement, "Plant biology"))
                            logging.debug(f"Keyword '{keyword}' is obsolete. Using replacement '{replacement}'.")
                    else:
                        label = mapping["Mapped Term Label"]
                        if label.lower() != "ignored" and label not in unique_names:
                            unique_names.add(label)
                            topic_terms.append((label, mapping["Mapped Term CURIE"]))
                else:
                    better_mapping = mapping.get("Better mapping")
                    if better_mapping == "Ignore":
                        remaining_keywords.append(keyword)
                        logging.debug(f"Keyword '{keyword}' ignored due to mapping decision.")
                    else:
                        # Only add if not ignored and not duplicated
                        if better_mapping.lower() != "ignored" and better_mapping not in unique_names:
                            unique_names.add(better_mapping)
                            topic_terms.append((better_mapping, mapping["Mapped Term CURIE"]))
                            logging.debug(f"Keyword '{keyword}' mapped to better mapping: {better_mapping}.")
            else:
                # Handle unmapped terms
                if not contains_special_characters(keyword) and "years" not in keyword.lower():
                    remaining_keywords.append(keyword)

        # Convert the collected terms into DefinedTerm objects
        topic_categories = [create_defined_term(label, curie) for (label, curie) in topic_terms]

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
        logging.debug(f"Loaded mappings: {mappings[:5]}")  # Print first 5 mappings for debugging

        processed_documents = process_documents(doc_list, mappings)

        for i, doc in enumerate(processed_documents, start=1):
            # Yield only if `topicCategory` exists and at least one category has a `name` field
            if "topicCategory" in doc and any("name" in category for category in doc["topicCategory"]):
                yield doc
            else:
                logging.debug(f"Document {i} has no valid topicCategory with 'name'. Not yielding.")
