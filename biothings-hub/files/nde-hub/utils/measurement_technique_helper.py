import csv
import os
from functools import lru_cache

import orjson
import requests
from rdflib import Graph


def get_identifier(url):
    """
    Extracts the identifier from the given URL.
    """
    parts = url.split("_")
    return parts[-1] if parts[-1] else "0000"


@lru_cache(maxsize=500)
def fetch_term_name_from_url(url):
    """
    Fetches the proper term name from an ontology URL.
    Caches results to minimize API calls.
    Returns None if the fetch fails.
    """
    try:
        from urllib.parse import quote

        from rdflib import URIRef
        from rdflib.namespace import RDFS, SKOS

        # Special handling for Ontobee URLs (extract IRI and use OLS)
        if 'ontobee.org' in url:
            # Extract the actual ontology IRI from the URL
            if '?iri=' in url:
                actual_iri = url.split('?iri=')[-1]
                # Try OLS API for Ontobee URLs
                ols_url = f"https://www.ebi.ac.uk/ols4/api/terms?iri={quote(actual_iri, safe='')}"
                response = requests.get(ols_url, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    if data.get('_embedded') and data['_embedded'].get('terms'):
                        terms = data['_embedded']['terms']
                        if terms and len(terms) > 0:
                            label = terms[0].get('label')
                            if label:
                                return label
            return None

        # Try OLS API first for all URLs (it's the most reliable and comprehensive)
        ols_url = f"https://www.ebi.ac.uk/ols4/api/terms?iri={quote(url, safe='')}"
        response = requests.get(ols_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data.get('_embedded') and data['_embedded'].get('terms'):
                terms = data['_embedded']['terms']
                if terms and len(terms) > 0:
                    label = terms[0].get('label')
                    if label:
                        return label

        # Fallback: Try to fetch RDF/OWL content directly
        response = requests.get(url, timeout=10, headers={'Accept': 'application/rdf+xml, text/turtle, application/ld+json'})
        response.raise_for_status()

        # Parse the RDF/OWL content
        g = Graph()
        try:
            g.parse(data=response.text, format='xml')
        except Exception:
            # Try turtle format if XML fails
            try:
                g.parse(data=response.text, format='turtle')
            except Exception:
                return None

        # Create a URIRef for the term we're looking for
        term_uri = URIRef(url)

        # Try to find the label for this specific term (rdfs:label or skos:prefLabel)
        for s, p, o in g.triples((term_uri, RDFS.label, None)):
            if o and str(o).strip():  # Make sure we have a non-empty value
                return str(o)

        for s, p, o in g.triples((term_uri, SKOS.prefLabel, None)):
            if o and str(o).strip():  # Make sure we have a non-empty value
                return str(o)

        return None
    except Exception as e:
        print(f"Warning: Failed to fetch term from {url}: {e}")
        return None


def load_mapping(name):
    """
    Loads mappings from the CSV and groups multiple mappings (if any) for the same repository technique.
    If the "Field" column is missing or empty, defaults to "measurementTechnique".
    If "Manually Mapped Term" is empty, uses the "Repository Technique" for the name.
    Fetches the actual term name from the URL if possible, falling back to the manually mapped term.
    """
    csv_file = f".ark_data/measurement_technique_lookup/{name}.csv"
    mapping = {}
    with open(csv_file, "r", newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            repo_technique = row["Repository Technique"].strip()
            manually_mapped = row["Manually Mapped Term"].strip() if row["Manually Mapped Term"] else ""
            url = row["URL"].strip()

            # Try to fetch the term name from the URL
            term_name = fetch_term_name_from_url(url)

            # Fallback to manually mapped term if URL fetch failed
            if not term_name:
                term_name = manually_mapped if manually_mapped else repo_technique

            entry = {
                "@type": "DefinedTerm",
                "name": term_name,
                "inDefinedTermSet": row["Ontology"].strip(),
                "url": url,
                "identifier": get_identifier(url),
                "isCurated": True,
                "field": row.get("Field", "measurementTechnique").strip() or "measurementTechnique",
            }
            if repo_technique in mapping:
                mapping[repo_technique].append(entry)
            else:
                mapping[repo_technique] = [entry]
    return mapping


def append_to_field(doc, field, new_entry):
    """
    Appends new_entry to doc[field]. If the field exists and is not a list,
    it converts it into a list first.
    """
    if field in doc:
        if isinstance(doc[field], list):
            doc[field].append(new_entry)
        else:
            doc[field] = [doc[field], new_entry]
    else:
        doc[field] = [new_entry]


def process_measurement_technique(data, name):
    """
    Updates the measurementTechnique field in each document based on the provided mapping.
    Each document's measurementTechnique property is standardized as a list.
    If a repository technique maps to more than one target field, each mapping entry is
    added to the appropriate field.
    If measurementTechnique ends up as an empty list, it is removed from the document.
    For the 'keywords' field, only the mapping name string is appended.
    """
    mapping = load_mapping(name)

    # Load documents from file or list.
    if isinstance(data, str):
        with open(os.path.join(data, "data.ndjson"), "rb") as f:
            doc_list = [orjson.loads(line) for line in f]
    else:
        doc_list = list(data)

    for doc in doc_list:
        new_mt = []  # This will store standardized measurementTechnique entries
        if "measurementTechnique" in doc:
            mt = doc["measurementTechnique"]
            # Normalize mt to a list for easier processing.
            if isinstance(mt, dict):
                mt = [mt]
            elif not isinstance(mt, list):
                mt = []

            for item in mt:
                original_name = item.get("name")
                if original_name in mapping:
                    # Process all mapping entries for this repository technique.
                    for map_record in mapping[original_name]:
                        target_field = map_record.get("field", "measurementTechnique")
                        new_entry = dict(map_record)
                        new_entry["originalName"] = original_name
                        if target_field == "keywords":
                            # Append only the name string for keywords.
                            append_to_field(doc, target_field, new_entry["name"])
                        elif target_field == "measurementTechnique":
                            new_mt.append(new_entry)
                        else:
                            append_to_field(doc, target_field, new_entry)
                else:
                    new_mt.append(item)

        # Overwrite or remove measurementTechnique based on whether new_mt is empty.
        if new_mt:
            doc["measurementTechnique"] = new_mt
        else:
            doc.pop("measurementTechnique", None)

    return doc_list
