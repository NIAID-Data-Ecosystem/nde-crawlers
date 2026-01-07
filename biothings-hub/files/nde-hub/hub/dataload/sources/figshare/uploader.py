import csv
import logging
import os
import re

import orjson
from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.extract import process_descriptions
from utils.utils import nde_upload_wrapper

logging.basicConfig(level=logging.INFO)

_SPECIAL_CHARS_RE = re.compile(r"[!@#$%^&*()\[\]{};:,<>?/|\\~`]" )


def load_mapping_sheet_from_csv(csv_file):
    """Load mapping sheet from a CSV file and return as a list of dictionaries."""
    mappings = []
    with open(csv_file, "r", newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        for row in reader:
            mappings.append(row)
    return mappings


def process_documents(documents, mapping_index):
    for doc in documents:
        keywords = doc.get("keywords", [])
        unique_names = set()  # Handle duplicates
        topic_terms = []
        remaining_keywords = []

        for keyword in keywords:
            keyword_lc = keyword.lower()
            mapping = mapping_index.get(keyword_lc)

            if mapping:
                # Check if the term is deprecated and GOOD
                decision = (mapping.get("Decision") or "").strip().lower()
                if decision == "good":
                    if "obsolete" in (mapping.get("Tags") or "").lower():
                        replacement = mapping.get("Consider")  # Replacement term
                        if replacement and replacement.lower() != "ignored":
                            if replacement not in unique_names:
                                unique_names.add(replacement)
                                topic_terms.append((replacement, "Plant biology"))
                            # logging.debug(f"Keyword '{keyword}' is obsolete. Using replacement '{replacement}'.")
                    else:
                        label = mapping["Mapped Term Label"]
                        if label.lower() != "ignored" and label not in unique_names:
                            unique_names.add(label)
                            topic_terms.append((label, mapping["Mapped Term CURIE"]))
                else:
                    better_mapping = mapping.get("Better mapping")
                    if better_mapping == "Ignore":
                        remaining_keywords.append(keyword)
                        # logging.debug(f"Keyword '{keyword}' ignored due to mapping decision.")
                    else:
                        # Only add if not ignored and not duplicated
                        if better_mapping and better_mapping.lower() != "ignored" and better_mapping not in unique_names:
                            unique_names.add(better_mapping)
                            topic_terms.append((better_mapping, mapping["Mapped Term CURIE"]))
                            # logging.debug(f"Keyword '{keyword}' mapped to better mapping: {better_mapping}.")
            else:
                # Handle unmapped terms
                if not contains_special_characters(keyword) and "years" not in keyword_lc:
                    remaining_keywords.append(keyword)

        # Convert the collected terms into DefinedTerm objects
        topic_categories = [create_defined_term(label, curie) for (label, curie) in topic_terms]

        # Update document fields
        if topic_categories:
            doc["topicCategory"] = topic_categories
        if remaining_keywords:
            doc["keywords"] = remaining_keywords

        yield doc


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
    return _SPECIAL_CHARS_RE.search(term) is not None


class FigshareUploader(NDESourceUploader):
    name = "figshare"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        mapping_file = "mappings.csv"

        def _iter_docs(data):
            if isinstance(data, str):
                ndjson_file = os.path.join(data, "data.ndjson")
                count = 0
                with open(ndjson_file, "rb") as f:
                    for line in f:
                        doc = orjson.loads(line)
                        count += 1
                        if count % 1000 == 0:
                            logging.info(f"Read {count} documents")
                        yield doc
            else:
                yield from data

        mappings = load_mapping_sheet_from_csv(mapping_file)
        mapping_index = {}
        for m in mappings:
            source_term = (m.get("Source Term") or "").strip().lower()
            if source_term and source_term not in mapping_index:
                mapping_index[source_term] = m

        # logging.debug(f"Loaded mappings: {mappings[:5]}")

        processed_documents = process_documents(_iter_docs(data_folder), mapping_index)

        def _has_valid_topic_category(doc):
            return "topicCategory" in doc and any("name" in category for category in doc["topicCategory"])

        filtered_documents = (doc for doc in processed_documents if _has_valid_topic_category(doc))
        enriched_documents = process_descriptions(filtered_documents)
        corrected_documents = corrections(enriched_documents)

        for i, doc in enumerate(corrected_documents, start=1):
            yield doc
