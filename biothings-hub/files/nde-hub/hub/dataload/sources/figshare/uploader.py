import csv
import re
from pathlib import Path

from config import logger
from hub.dataload.nde import NDESourceUploader
from utils import iter_ndjson, nde_upload_wrapper

# Figshare keyword -> EDAM topic sheet (Text2term suggestions, manually reviewed).
MAPPING_FILE = Path(__file__).resolve().parent / "topic_mappings.tsv"
ACCEPTED_DECISIONS = {"good", "ok"}

_SPECIAL_CHARS_RE = re.compile(r"[!@#$%^&*()\[\]{};:,<>?/|\\~`]")


def row_topics(row):
    """Return the (name, identifier) topics a sheet row maps its keyword to.

    A Better mapping CURIE overrides the Text2term match; multiple terms are
    '|'-separated in both Better mapping columns. An empty result means the
    keyword isn't a topic (ignored, or rejected without a replacement) and stays a keyword.
    """
    better = (row.get("Better mapping") or "").strip()
    curies = (row.get("Better mapping CURIE") or "").strip()
    if better.lower() == "ignore":
        return []
    if curies:
        names = [name.strip() for name in better.split("|")]
        identifiers = [curie.strip() for curie in curies.split("|")]
        if len(names) != len(identifiers):
            logger.warning(
                "Figshare topic mapping for %r: Better mapping %r doesn't pair up with CURIEs %r",
                row.get("Source Term"),
                better,
                curies,
            )
            return []
        return list(zip(names, identifiers))
    if (row.get("Decision") or "").strip().lower() in ACCEPTED_DECISIONS:
        return [(row["Mapped Term Label"].strip(), row["Mapped Term CURIE"].strip())]
    return []


def load_mapping_index(mapping_file=MAPPING_FILE):
    """Map each lowercased source term to its topics; the first row for a term wins."""
    mapping_index = {}
    with open(mapping_file, "r", newline="", encoding="utf-8") as file:
        for row in csv.DictReader(file, delimiter="\t"):
            source_term = (row.get("Source Term") or "").strip().lower()
            if source_term and source_term not in mapping_index:
                mapping_index[source_term] = row_topics(row)
    return mapping_index


def process_documents(documents, mapping_index):
    for doc in documents:
        topic_categories = []
        seen_identifiers = set()  # Handle duplicates
        remaining_keywords = []

        for keyword in doc.get("keywords", []):
            keyword_lc = keyword.strip().lower()
            topics = mapping_index.get(keyword_lc)

            if topics is None:
                # Handle unmapped terms
                if not contains_special_characters(keyword) and "years" not in keyword_lc:
                    remaining_keywords.append(keyword)
            elif not topics:
                remaining_keywords.append(keyword)
            else:
                for name, identifier in topics:
                    if identifier not in seen_identifiers:
                        seen_identifiers.add(identifier)
                        topic_categories.append(create_defined_term(name, identifier))

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
        "curatedBy": {"@type": "SoftwareApplication", "name": "Text2Term-assisted manual mapping"},
        "isCurated": True,
    }


def contains_special_characters(term):
    """Check if the term contains special characters."""
    return _SPECIAL_CHARS_RE.search(term) is not None


class FigshareUploader(NDESourceUploader):
    name = "figshare"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        mapping_index = load_mapping_index()
        processed_documents = process_documents(iter_ndjson(data_folder), mapping_index)

        def _has_valid_topic_category(doc):
            return "topicCategory" in doc and any("name" in category for category in doc["topicCategory"])

        yield from (doc for doc in processed_documents if _has_valid_topic_category(doc))
