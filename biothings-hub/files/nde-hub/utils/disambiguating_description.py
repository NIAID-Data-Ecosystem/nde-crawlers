"""Curated disambiguatingDescription summaries.

Runs when the source has a CSV in `/data/nde-hub/disambiguating_descriptions/`
mapping record ids to a processed summary.
"""

import csv
from functools import cache

from config import logger

LOOKUP_DIR = "/data/nde-hub/disambiguating_descriptions"


def lookup_file(source):
    return f"{LOOKUP_DIR}/{source}.csv"


@cache
def load_descriptions(source):
    """Load a source's summaries: {record id: disambiguating description}."""
    with open(lookup_file(source), "r") as file:
        descriptions = {row["_id"].lower(): row["Processed Summary"] for row in csv.DictReader(file)}
    logger.info("Loaded %s disambiguating descriptions for %s", len(descriptions), source)
    return descriptions


def add_disambiguating_description(docs, source):
    """Add the curated disambiguatingDescription to every record in one batch."""
    descriptions = load_descriptions(source)
    for doc in docs:
        if description := descriptions.get(doc["_id"].lower()):
            doc["disambiguatingDescription"] = description
        yield doc
