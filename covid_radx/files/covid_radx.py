import logging
import os

import orjson

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

# The RADx Data Hub site is offline; records come from a 2025-03-26 crawl.
NDJSON_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "covid_radx.ndjson")
CATALOG_NAME = "COVID Rapid Acceleration of Diagnostics (RADx) Data Hub"

# @type the current schema requires on nested objects the archived crawl predates.
NESTED_TYPES = {
    "sdPublisher": "DataCatalog",
    "healthCondition": "DefinedTerm",
    "measurementTechnique": "DefinedTerm",
    "variableMeasured": "DefinedTerm",
    "funding": "MonetaryGrant",
    "isBasedOn": "CreativeWork",
    "usageInfo": "CreativeWork",
}


def _set_type(value, type_name):
    for item in value if isinstance(value, list) else [value]:
        if isinstance(item, dict):
            item.setdefault("@type", type_name)


def parse(doc):
    """Bring an archived record up to the current parser's output."""
    doc["_id"] = doc["identifier"]

    catalog = doc["includedInDataCatalog"]
    catalog["name"] = CATALOG_NAME
    if archived_at := catalog.pop("dataset", None):
        catalog["archivedAt"] = archived_at

    for field, type_name in NESTED_TYPES.items():
        if field in doc:
            _set_type(doc[field], type_name)

    if funding := doc.get("funding"):
        for grant in funding if isinstance(funding, list) else [funding]:
            _set_type(grant.get("funder"), "Organization")

    if url := doc.get("isBasedOn", {}).get("url"):
        doc["isBasedOn"]["url"] = url.strip("{}")

    if interval := doc.get("temporalCoverage", {}).get("temporalInterval"):
        doc["temporalCoverage"] = {**interval, "@type": "TemporalInterval"}

    return doc


def make_requests():
    count = 0
    with open(NDJSON_PATH, "rb") as f:
        for line in f:
            if line.strip():
                count += 1
                yield parse(orjson.loads(line))
    logger.info("Total number of documents parsed: %s", count)
