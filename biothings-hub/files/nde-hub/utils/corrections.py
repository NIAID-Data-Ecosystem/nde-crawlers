import json
import logging
import os

import orjson
import requests

logging.basicConfig(level=logging.INFO)


def load_documents(data):
    """
    :param data: a list of documents or a path to a directory containing data.ndjson file
    :return: a list of documents
    """
    doc_list = []
    count = 0

    if isinstance(data, str):
        ndjson_file = os.path.join(data, "data.ndjson")
        with open(ndjson_file, "rb") as f:
            for line in f:
                doc = orjson.loads(line)
                doc_list.append(doc)
                count += 1
                if count % 1000 == 0:
                    logging.info(f"Processed {count} documents")
    else:
        doc_list = list(data)

    return doc_list


def get_github_file_content(owner, repo, path):
    url = f"https://api.github.com/repos/{owner}/{repo}/contents/{path}"
    headers = {"Accept": "application/vnd.github.v3.raw"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Error fetching {path}: {response.status_code} {response.text}")


def get_record_ids(records_content):
    ids = []
    for line in records_content.strip().splitlines():
        line = line.strip()
        if line:
            if "id=" in line:
                id_value = line.split("id=")[-1]
                ids.append(id_value.lower())
    return ids


def update_source_organization(record_metadata, correction_organizations):
    existing_orgs = record_metadata.get("sourceOrganization", [])

    # Ensure existing_orgs is a list
    if not isinstance(existing_orgs, list):
        existing_orgs = [existing_orgs]

    # Build a set of existing organization identifiers
    existing_org_identifiers = set()
    for org in existing_orgs:
        identifier = org.get("name") or org.get("url")
        if identifier:
            existing_org_identifiers.add(identifier.lower())

    # Append new organizations that are not already in existing_orgs
    for new_org in correction_organizations:
        identifier = new_org.get("name") or new_org.get("url")
        if identifier and identifier.lower() not in existing_org_identifiers:
            existing_orgs.append(new_org)
            existing_org_identifiers.add(identifier.lower())

    record_metadata["sourceOrganization"] = existing_orgs
    return record_metadata


def update_documents_with_corrections(documents, correction_name):
    owner = "NIAID-Data-Ecosystem"
    repo = "nde-metadata-corrections"
    corrections_dir = "collections_corrections"

    correction_file_path = f"{corrections_dir}/{correction_name}_correction.json"
    records_file_path = f"{corrections_dir}/{correction_name}_records.txt"

    try:
        # Fetch correction content
        correction_content = get_github_file_content(owner, repo, correction_file_path)
        correction_json = json.loads(correction_content)
        correction_organizations = correction_json["source_organization"]

        # Fetch records content
        records_content = get_github_file_content(owner, repo, records_file_path)
        record_ids = set(get_record_ids(records_content))
    except Exception as e:
        logging.error(f"Error fetching correction files: {e}")
        return documents

    # Update documents
    updated_documents = []
    for doc in documents:
        doc_id = doc.get("_id")
        if doc_id:
            doc_id = doc_id.lower()
        if doc_id and doc_id in record_ids:
            updated_doc = update_source_organization(doc, correction_organizations)
            updated_documents.append(updated_doc)
            logging.info(f"Updated document with ID: {doc_id}")
        else:
            updated_documents.append(doc)

    return updated_documents


def corrections(data, correction_name):
    # Load documents from NDJSON file or list
    documents = load_documents(data)

    # Update documents with corrections
    updated_documents = update_documents_with_corrections(documents, correction_name)

    return updated_documents
