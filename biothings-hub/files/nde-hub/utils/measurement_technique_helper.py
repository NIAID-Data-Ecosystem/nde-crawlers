import csv
import os

import orjson


def get_identifier(url):
    """
    Extracts the identifier from the given URL.
    """
    parts = url.split("_")
    return parts[-1] if parts[-1] else "0000"


def load_mapping(name):
    """
    Loads mappings from the CSV and groups multiple mappings (if any) for the same repository technique.
    If the "Field" column is missing or empty, defaults to "measurementTechnique".
    """
    csv_file = f"/data/nde-hub/standardizers/measurement_technique_lookup/{name}.csv"
    mapping = {}
    with open(csv_file, "r", newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            repo_technique = row["Repository Technique"]
            entry = {
                "name": row["Manually Mapped Term"],
                "inDefinedTermSet": row["Ontology"],
                "url": row["URL"],
                "identifier": get_identifier(row["URL"]),
                "isCurated": True,
                "field": row.get("Field", "measurementTechnique") or "measurementTechnique",
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
                        if target_field == "measurementTechnique":
                            new_mt.append(new_entry)
                        else:
                            # Append the new entry to the corresponding target field.
                            append_to_field(doc, target_field, new_entry)
                else:
                    # If no mapping found, keep the original item in measurementTechnique.
                    new_mt.append(item)

        # Overwrite or remove measurementTechnique based on whether new_mt is empty.
        if new_mt:
            doc["measurementTechnique"] = new_mt
        else:
            doc.pop("measurementTechnique", None)

    return doc_list
