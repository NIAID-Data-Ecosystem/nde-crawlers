"""DDE-specific curation of terms that already carry an ontology URL.

Data Discovery Engine submitters pick terms from a set of ontologies, so a term
that arrives with a URL from one of those ontologies can be curated directly
instead of being guessed at by the `terms` stage. The ontology-per-property
table lives in a Google Sheet.

Used only by the `dde` source.
"""

import csv
import datetime
from functools import cache
from io import StringIO

import requests
from config import logger

from .common import as_list
from .terms import get_species_details, query_condition

PROPERTIES_CSV_URL = "https://docs.google.com/spreadsheets/d/107WVX39r_a6xBGZ_gCku0LBNRWWBmk9x7Dg53wj1SiI/export?format=csv"

HEALTH_CONDITION_TERM_SETS = ("MeSH", "DOID", "NCIT", "MONDO")

_species_cache = {}
_health_condition_cache = {}


def _curated_by():
    return {
        "curatedBy": {
            "name": "Data Discovery Engine",
            "url": "https://discovery.biothings.io/",
            "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
        }
    }


def process_csv_data(csv_content):
    """Parse the DDE sheet into {property: {ontology base url: term set}}."""
    properties = {}
    current_property = None

    for row in csv.DictReader(StringIO(csv_content)):
        if ".inDefinedTermSet" in row["Old button text"]:
            # Start of a new property section
            current_property = row["Old button text"].split(".")[0]
            properties[current_property] = {}
        elif row["base url (for mapping purposes)"] and current_property:
            base_url = row["base url (for mapping purposes)"]
            if "http" not in base_url:
                continue
            properties[current_property][base_url] = row["InDefinedTermSet"]

    return properties


@cache
def load_properties():
    """Fetch (once per upload) the ontology-per-property table."""
    properties = process_csv_data(requests.get(PROPERTIES_CSV_URL).text)
    logger.info("Loaded DDE term sets for %s properties", len(properties))
    return properties


def species_func(species, term_set):
    """Curate a species term from its UniProt URL."""
    if term_set != "UniProt":
        return species

    key = species.get("name")
    if key in _species_cache:
        return _species_cache[key]

    identifier = species.get("url").split("_")[-1]
    try:
        species.update(get_species_details(species.get("name"), identifier))
        species.update(_curated_by())
        _species_cache[key] = species
    except Exception as e:
        logger.error("Error fetching species details: %s", e)
    return species


def health_condition_func(health_condition, term_set):
    """Curate a health condition term against its ontology."""
    key = (health_condition.get("name"), term_set)
    if key in _health_condition_cache:
        return _health_condition_cache[key]

    if term_set in HEALTH_CONDITION_TERM_SETS:
        if result := query_condition(health_condition.get("name")):
            health_condition.update(result)
            health_condition.update(_curated_by())
            _health_condition_cache[key] = health_condition
    else:
        health_condition["inDefinedTermSet"] = "Other"
        _health_condition_cache[key] = health_condition
    return health_condition


def get_term_set_helper(url, property_dict):
    """Return the term set whose base URL matches `url`."""
    for prop, term_set in property_dict.items():
        if prop in url:
            return term_set
    return None


def de_duplicate_dicts(dict_list):
    """Deduplicate terms by identifier, keeping the ones that have none."""
    unique_dicts = []
    seen = set()
    for obj in dict_list:
        identifier = obj.get("identifier")
        if identifier is None:
            unique_dicts.append(obj)
        elif identifier not in seen:
            seen.add(identifier)
            unique_dicts.append(obj)
    return unique_dicts


def get_in_defined_term_set(doc, properties_dict):
    """Curate every term in `doc` whose URL belongs to a known ontology."""
    logger.debug("Processing doc: %s", doc.get("_id"))
    nde_properties = {
        "species": species_func,
        "infectiousAgent": species_func,
        "healthCondition": health_condition_func,
        # TODO add other properties
    }

    # Species and infectiousAgent are collected together so a term can move
    # between them based on how it classifies.
    organisms = {"species": [], "infectiousAgent": []}

    for nde_property, curate in nde_properties.items():
        doc_property = doc.get(nde_property)
        if not doc_property:
            continue
        is_organism = nde_property in organisms

        for item in as_list(doc_property):
            url = item.get("url")
            term_set = get_term_set_helper(url, properties_dict.get(nde_property, {})) if url else None
            if not term_set:
                # Without a known URL, leave the term for the pubtator lookup to curate.
                if is_organism:
                    organisms[nde_property].append(item)
                continue

            result = curate(item, term_set)
            if is_organism:
                target = "infectiousAgent" if result.get("classification") == "infectiousAgent" else "species"
                organisms[target].append(result)
            elif isinstance(doc_property, list):
                doc_property[doc_property.index(item)] = result
            else:
                doc[nde_property] = result

    for field, terms in organisms.items():
        if terms:
            doc[field] = de_duplicate_dicts(terms)

    return doc


def handle_dde_docs(docs):
    """Curate DDE-submitted terms for every record in one batch."""
    properties = load_properties()
    for doc in docs:
        yield get_in_defined_term_set(doc, properties)
