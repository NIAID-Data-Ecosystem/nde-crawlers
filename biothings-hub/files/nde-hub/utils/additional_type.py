"""Curated additionalType terms for related records.

Runs for records whose `isBasedOn` or `hasPart` entries carry an
`additionalType`, when the source has a lookup TSV in `additional_type_lookup/`.
Each row maps a label the repository gives to an ontology term, which replaces
the label as a DefinedTerm. A row whose term name is `DELETE` drops a
placeholder label, and a row without a URL leaves the label as the repository
gave it, so only ontology-linked terms become DefinedTerms.
"""

import csv
from collections import Counter
from functools import cache
from pathlib import Path

from config import logger

from .common import dict_entries

LOOKUP_DIR = Path(__file__).resolve().parent / "additional_type_lookup"
FIELDS = ("isBasedOn", "hasPart")

# The same marker as the measurementTechnique lookups.
DELETE = "DELETE"


def lookup_file(source):
    return str(LOOKUP_DIR / f"{source}.tsv")


@cache
def load_mapping(source):
    """{repository label: DefinedTerm, or None to drop the label}. Labels without a term are left out."""
    mapping = {}
    with open(lookup_file(source), newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            label = row["Repository Label"].strip()
            name = (row.get("Term Name") or "").strip()
            url = (row.get("URL") or "").strip()
            if name.upper() == DELETE:
                mapping[label] = None
            elif url:
                mapping[label] = {
                    "@type": "DefinedTerm",
                    "name": name,
                    "url": url,
                    "identifier": url.rsplit("_", 1)[-1],
                    "inDefinedTermSet": row["Ontology"].strip(),
                    "isCurated": True,
                }
    logger.info("Loaded %s additionalType mappings for %s", len(mapping), source)
    return mapping


def reset_caches():
    load_mapping.cache_clear()


def curate_additional_types(docs, source):
    """Replace the additionalType labels of every record's isBasedOn and hasPart entries in one batch."""
    mapping = load_mapping(source)
    unmapped = Counter()
    for doc in docs:
        for field in FIELDS:
            for entry in dict_entries(doc, field):
                value = entry.get("additionalType")
                label = value.get("name") if isinstance(value, dict) else value
                if not isinstance(label, str) or (isinstance(value, dict) and value.get("@type") == "DefinedTerm"):
                    continue
                if label not in mapping:
                    unmapped[label] += 1
                elif mapping[label] is None:
                    del entry["additionalType"]
                else:
                    entry["additionalType"] = {**mapping[label], "originalName": label}
        yield doc
    if unmapped:
        logger.info("additionalType labels with no term for %s: %s", source, dict(unmapped.most_common(10)))
