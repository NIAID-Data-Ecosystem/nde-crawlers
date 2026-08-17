"""Curated measurementTechnique mappings.

Runs for records with a `measurementTechnique` when the source has a mapping
CSV in `/data/nde-hub/standardizers/measurement_technique_lookup/`. Each row
maps a repository technique to an ontology term, optionally routing it to a
different field (`Field` column, e.g. `keywords`).
"""

import csv
from functools import cache
from urllib.parse import quote

import requests
from config import logger
from rdflib import Graph, URIRef
from rdflib.namespace import RDFS, SKOS

from .common import dict_entries

LOOKUP_DIR = "/data/nde-hub/standardizers/measurement_technique_lookup"


def lookup_file(source):
    return f"{LOOKUP_DIR}/{source}.csv"


def get_identifier(url):
    """Extract the ontology term identifier from its URL."""
    parts = url.split("_")
    return parts[-1] if parts[-1] else "0000"


@cache
def fetch_term_name_from_url(url):
    """Fetch an ontology term's label, preferring OLS then the raw RDF. None if unavailable.

    Unbounded: the URLs come from the curated technique CSVs, so the key space is
    bounded by those files, and the values are short labels.
    """
    try:
        # Ontobee URLs carry the real IRI as a query parameter.
        if "ontobee.org" in url:
            if "?iri=" in url:
                return _ols_label(url.split("?iri=")[-1])
            return None

        if label := _ols_label(url):
            return label

        response = requests.get(
            url,
            timeout=10,
            headers={"Accept": "application/rdf+xml, text/turtle, application/ld+json"},
        )
        response.raise_for_status()

        graph = Graph()
        for rdf_format in ("xml", "turtle"):
            try:
                graph.parse(data=response.text, format=rdf_format)
                break
            except Exception:
                continue
        else:
            return None

        term_uri = URIRef(url)
        for predicate in (RDFS.label, SKOS.prefLabel):
            for _, _, label in graph.triples((term_uri, predicate, None)):
                if label and str(label).strip():
                    return str(label)
        return None
    except Exception as e:
        logger.warning("Failed to fetch term from %s: %s", url, e)
        return None


def _ols_label(iri):
    """Look a term's label up in the EBI Ontology Lookup Service."""
    response = requests.get(f"https://www.ebi.ac.uk/ols4/api/terms?iri={quote(iri, safe='')}", timeout=10)
    if response.status_code != 200:
        return None
    terms = response.json().get("_embedded", {}).get("terms") or []
    return terms[0].get("label") if terms else None


@cache
def load_mapping(source):
    """Load a source's technique mapping: {repository technique: [target terms]}.

    Term names come from the ontology URL when it resolves, falling back to the
    manually mapped term and then to the repository's own wording. Cached, since
    resolving the URLs costs a request each.
    """
    mapping = {}
    with open(lookup_file(source), "r", newline="", encoding="utf-8") as csvfile:
        for row in csv.DictReader(csvfile):
            repo_technique = row["Repository Technique"].strip()
            manually_mapped = row["Manually Mapped Term"].strip() if row["Manually Mapped Term"] else ""
            url = row["URL"].strip()

            entry = {
                "@type": "DefinedTerm",
                "name": fetch_term_name_from_url(url) or manually_mapped or repo_technique,
                "inDefinedTermSet": row["Ontology"].strip(),
                "url": url,
                "identifier": get_identifier(url),
                "isCurated": True,
                "field": row.get("Field", "measurementTechnique").strip() or "measurementTechnique",
            }
            mapping.setdefault(repo_technique, []).append(entry)
    logger.info("Loaded %s measurementTechnique mappings for %s", len(mapping), source)
    return mapping


def append_to_field(doc, field, new_entry):
    """Append `new_entry` to `doc[field]`, making it a list if it isn't one."""
    if field in doc:
        if isinstance(doc[field], list):
            doc[field].append(new_entry)
        else:
            doc[field] = [doc[field], new_entry]
    else:
        doc[field] = [new_entry]


def _apply_mapping(doc, mapping):
    new_mt = []
    for item in dict_entries(doc, "measurementTechnique"):
        original_name = item.get("name")
        if original_name not in mapping:
            new_mt.append(item)
            continue
        for map_record in mapping[original_name]:
            new_entry = dict(map_record)
            target_field = new_entry.pop("field", "measurementTechnique")
            new_entry["originalName"] = original_name
            if target_field == "measurementTechnique":
                new_mt.append(new_entry)
            elif target_field == "keywords":
                append_to_field(doc, target_field, new_entry["name"])
            else:
                append_to_field(doc, target_field, new_entry)

    if new_mt:
        doc["measurementTechnique"] = new_mt
    else:
        doc.pop("measurementTechnique", None)
    return doc


def process_measurement_technique(docs, source):
    """Standardize the measurementTechnique of every record in one batch."""
    mapping = load_mapping(source)
    for doc in docs:
        yield _apply_mapping(doc, mapping)
