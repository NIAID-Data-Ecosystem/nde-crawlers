import json
import logging
import os

import orjson
import requests

logging.basicConfig(level=logging.INFO)

def list_github_files(owner, repo, directory):
    """
    List files in a GitHub repository directory using the GitHub contents API.
    """
    url = f"https://api.github.com/repos/{owner}/{repo}/contents/{directory}"
    headers = {"Accept": "application/vnd.github.v3+json"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = response.json()
        return [item["name"] for item in data if item["type"] == "file"]
    else:
        raise Exception(f"Error listing directory {directory}: {response.status_code} {response.text}")

def get_github_file_content(owner, repo, path):
    """
    Fetch file content from GitHub.
    """
    url = f"https://api.github.com/repos/{owner}/{repo}/contents/{path}"
    headers = {"Accept": "application/vnd.github.v3.raw"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Error fetching {path}: {response.status_code} {response.text}")

def fetch_correction_files(correction_name):
    """
    Try to fetch the correction JSON and records file from production first.
    If not found, then try staging.
    Returns:
      - correction_json: Parsed correction JSON.
      - approved: Boolean flag.
      - records_file_path: The full path to the records file to fetch.
    """
    owner = "NIAID-Data-Ecosystem"
    repo = "nde-metadata-corrections"

    # Build paths for production and staging folders
    prod_corrections_dir = "collections_corrections_production"
    staging_corrections_dir = "collections_corrections_staging"

    prod_correction_file = f"{prod_corrections_dir}/{correction_name}_correction.json"
    prod_records_file = f"{prod_corrections_dir}/{correction_name}_records.txt"

    staging_correction_file = f"{staging_corrections_dir}/{correction_name}_correction.json"
    staging_records_file = f"{staging_corrections_dir}/{correction_name}_records.txt"

    try:
        # Try production folder first.
        correction_content = get_github_file_content(owner, repo, prod_correction_file)
        correction_json = json.loads(correction_content)
        # If the JSON contains an "approved" field, use it; otherwise assume production is approved.
        approved = correction_json.get("approved", True)
        records_file_path = prod_records_file
    except Exception as prod_error:
        logging.info(f"Production file for '{correction_name}' not found or error encountered: {prod_error}. Trying staging folder.")
        correction_content = get_github_file_content(owner, repo, staging_correction_file)
        correction_json = json.loads(correction_content)
        approved = correction_json.get("approved", False)
        records_file_path = staging_records_file

    return correction_json, approved, records_file_path

def get_record_ids(records_content):
    """
    Extract record ids from the records.txt file content.
    Assumes each non-empty line contains an URL with an "id=" query parameter.
    """
    ids = []
    for line in records_content.strip().splitlines():
        line = line.strip()
        if line:
            if "id=" in line:
                id_value = line.split("id=")[-1]
                ids.append(id_value.lower())
    return ids

def update_source_organization(record_metadata, correction_organizations, approved=True):
    existing_orgs = record_metadata.get("sourceOrganization", [])
    if not isinstance(existing_orgs, list):
        existing_orgs = [existing_orgs]

    # Gather identifiers already in the document (using name or url)
    existing_org_identifiers = {
        (org.get("name") or org.get("url")).lower()
        for org in existing_orgs if (org.get("name") or org.get("url"))
    }

    # Add each new correction organization, attaching the approval flag to the object
    for new_org in correction_organizations:
        identifier = new_org.get("name") or new_org.get("url")
        if identifier and identifier.lower() not in existing_org_identifiers:
            new_org["correctionApproved"] = approved
            existing_orgs.append(new_org)
            existing_org_identifiers.add(identifier.lower())

    # Save back into the document; note that we no longer set a top-level flag.
    record_metadata["sourceOrganization"] = existing_orgs
    return record_metadata

def load_documents(data):
    """
    Load documents from a list or a path to a directory containing a data.ndjson file.
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

def update_documents_with_corrections(documents):
    """
    For every correction available in our corrections repo, update any documents whose
    _id appears in that correction's records file.
    """
    owner = "NIAID-Data-Ecosystem"
    repo = "nde-metadata-corrections"
    prod_dir = "collections_corrections_production"
    staging_dir = "collections_corrections_staging"

    corrections_dict = {}

    # List production corrections (files ending with "_records.txt")
    try:
        production_files = list_github_files(owner, repo, prod_dir)
        production_records = [fname for fname in production_files if fname.endswith("_records.txt")]
        for fname in production_records:
            correction_name = fname[:-len("_records.txt")]
            corrections_dict[correction_name] = None  # placeholder for later details
    except Exception as e:
        logging.error(f"Error listing production corrections: {e}")

    # List staging corrections (if any) that are not already in production.
    try:
        staging_files = list_github_files(owner, repo, staging_dir)
        staging_records = [fname for fname in staging_files if fname.endswith("_records.txt")]
        for fname in staging_records:
            correction_name = fname[:-len("_records.txt")]
            if correction_name not in corrections_dict:
                corrections_dict[correction_name] = None
    except Exception as e:
        logging.error(f"Error listing staging corrections: {e}")

    logging.info(f"Found correction names: {list(corrections_dict.keys())}")

    # For each correction, fetch the correction JSON and records file, and parse the record IDs.
    for correction_name in list(corrections_dict.keys()):
        try:
            correction_json, approved, records_file_path = fetch_correction_files(correction_name)
            correction_organizations = correction_json["sourceOrganization"]
            records_content = get_github_file_content(owner, repo, records_file_path)
            record_ids = set(get_record_ids(records_content))
            corrections_dict[correction_name] = (record_ids, correction_organizations, approved)
            logging.info(f"Loaded correction '{correction_name}' with {len(record_ids)} record IDs (approved: {approved})")
        except Exception as e:
            logging.error(f"Error fetching correction files for '{correction_name}': {e}")
            corrections_dict.pop(correction_name, None)

    # Update each document: for every correction, if the document _id is in the records list, update it.
    for doc in documents:
        doc_id = doc.get("_id", "").lower()
        for correction_name, correction_data in corrections_dict.items():
            if correction_data is None:
                continue
            record_ids, correction_organizations, approved = correction_data
            if doc_id and doc_id in record_ids:
                doc = update_source_organization(doc, correction_organizations, approved=approved)
                logging.info(f"Updated document {doc_id} with correction '{correction_name}' (approved: {approved})")
    return documents

def corrections(data):
    """
    Load documents and then update them by applying every correction found in our corrections repo.
    """
    documents = load_documents(data)
    updated_documents = update_documents_with_corrections(documents)
    return updated_documents
