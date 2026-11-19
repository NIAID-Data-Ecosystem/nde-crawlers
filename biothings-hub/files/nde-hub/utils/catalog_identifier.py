"""Short catalog identifiers on every record's `includedInDataCatalog`.

Each catalog entry gets the short, lower-cased identifier of its NDE
ResourceCatalog (the RepoMetaCuration sheet's identifier, e.g. "immport"), so a
source's records can be queried as `includedInDataCatalog.identifier:immport`
instead of by the catalog's full name.

`catalog_identifier_lookup.tsv` maps the catalog names the crawlers emit to those
identifiers. The lookup wins over an identifier a crawler set itself. A catalog
the lookup doesn't know is left as it is, with a warning, so a new source needs
a row there.
"""

import csv
from functools import cache
from pathlib import Path

from config import logger

from .common import dict_entries

LOOKUP_FILE = Path(__file__).resolve().parent / "catalog_identifier_lookup.tsv"

# Catalog names already warned about during this upload.
_unmapped_names = set()


def lookup_file():
    return str(LOOKUP_FILE)


@cache
def load_mapping():
    """{catalog name: identifier}."""
    with open(LOOKUP_FILE, newline="", encoding="utf-8") as f:
        mapping = {row["Catalog Name"].strip(): row["Identifier"].strip() for row in csv.DictReader(f, delimiter="\t")}
    logger.info("Loaded %s catalog identifiers", len(mapping))
    return mapping


def reset_caches():
    load_mapping.cache_clear()
    _unmapped_names.clear()


def add_catalog_identifiers(docs, source):
    """Set `identifier` on each `includedInDataCatalog` entry whose name is in the lookup."""
    mapping = load_mapping()
    for doc in docs:
        for catalog in dict_entries(doc, "includedInDataCatalog"):
            name = catalog.get("name")
            if not isinstance(name, str):
                continue
            identifier = mapping.get(name.strip())
            if identifier:
                catalog["identifier"] = identifier
            elif name not in _unmapped_names:
                _unmapped_names.add(name)
                logger.warning(
                    "No catalog identifier for %r (source %s): add a row to %s", name, source, LOOKUP_FILE.name
                )
        yield doc
