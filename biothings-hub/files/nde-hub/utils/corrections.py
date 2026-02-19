"""
Corrections module for applying sourceOrganization metadata from the
nde-metadata-corrections GitHub repository.

Documents are matched to corrections using two strategies:
  1. ID-based:      Match document _id against pre-curated _records.txt lists.
  2. Funding-based: Match document funding.identifier against grant patterns
                    stored in correction JSONs (fundingIdentifiers field).

Strategy 2 ensures that NEW records are corrected immediately at upload time
without waiting for the records list to be regenerated externally.

Corrections data is cached at the module level so it is fetched from GitHub
once per build (or per CACHE_TTL window) and reused across all sources.
"""

import json
import logging
import math
import os
import threading
import time

import orjson
import requests
from config import token

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
OWNER = "NIAID-Data-Ecosystem"
REPO = "nde-metadata-corrections"
PROD_DIR = "collections_corrections_production"
STAGING_DIR = "collections_corrections_staging"
CACHE_TTL = 3600  # seconds

# ---------------------------------------------------------------------------
# Module-level cache (thread-safe)
# ---------------------------------------------------------------------------
_cache_lock = threading.Lock()
_corrections_cache = None
_corrections_cache_time = 0.0


# ---------------------------------------------------------------------------
# GitHub helpers
# ---------------------------------------------------------------------------
def get_auth_headers(accept_header):
    """
    Return headers for GitHub API requests, including the Authorization header
    if a token is available in the environment variable GITHUB_TOKEN.
    """
    headers = {"Accept": accept_header}
    if token:
        headers["Authorization"] = f"token {token}"
    return headers


def list_github_files(owner, repo, directory):
    """
    List files in a GitHub repository directory using the GitHub contents API.
    """
    url = f"https://api.github.com/repos/{owner}/{repo}/contents/{directory}"
    headers = get_auth_headers("application/vnd.github.v3+json")
    response = requests.get(url, headers=headers, timeout=30)
    if response.status_code == 200:
        data = response.json()
        return [item["name"] for item in data if item["type"] == "file"]
    else:
        raise Exception(f"Error listing directory {directory}: {response.status_code} {response.text}")


def get_github_file_content(owner, repo, path):
    """
    Fetch raw file content from GitHub.
    """
    url = f"https://api.github.com/repos/{owner}/{repo}/contents/{path}"
    headers = get_auth_headers("application/vnd.github.v3.raw")
    response = requests.get(url, headers=headers, timeout=30)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Error fetching {path}: {response.status_code} {response.text}")


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------
def get_record_ids(records_content):
    """
    Extract record ids from the records.txt file content.
    Assumes each non-empty line contains a URL with an "id=" query parameter.
    """
    ids = []
    for line in records_content.strip().splitlines():
        line = line.strip()
        if line and "id=" in line:
            id_value = line.split("id=")[-1]
            ids.append(id_value.lower())
    return ids


def sanitize_org(org):
    """
    Replace NaN values in the 'url' field with an empty string.
    """
    if "url" in org:
        url_val = org["url"]
        if isinstance(url_val, float) and math.isnan(url_val):
            org["url"] = ""
    return org


# ---------------------------------------------------------------------------
# Source-organization merging
# ---------------------------------------------------------------------------
def update_source_organization(record_metadata, correction_organizations, approved=True):
    """
    Merge correction organizations into the document's sourceOrganization,
    sanitizing any NaN values in the 'url' fields.

    Deduplication is based on the lowercase name (or url) of each organization,
    so calling this multiple times with the same correction is safe/idempotent.
    """
    existing_orgs = record_metadata.get("sourceOrganization", [])

    if not isinstance(existing_orgs, list):
        existing_orgs = [existing_orgs]

    # Sanitize existing organizations.
    existing_orgs = [sanitize_org(org) for org in existing_orgs]

    # Gather identifiers already present (name or url, lowercased).
    existing_org_identifiers = {
        (org.get("name") or org.get("url")).lower()
        for org in existing_orgs
        if (org.get("name") or org.get("url"))
    }

    # Sanitize and merge new correction organizations.
    for new_org in correction_organizations:
        new_org = sanitize_org(new_org)
        identifier = new_org.get("name") or new_org.get("url")
        if identifier and identifier.lower() not in existing_org_identifiers:
            new_org["correctionApproved"] = approved
            existing_orgs.append(new_org)
            existing_org_identifiers.add(identifier.lower())

    record_metadata["sourceOrganization"] = existing_orgs
    return record_metadata


# ---------------------------------------------------------------------------
# Correction file fetching
# ---------------------------------------------------------------------------
def fetch_correction_files(correction_name):
    """
    Try to fetch the correction JSON and records file from production first.
    If not found, then try staging.

    Returns:
      - correction_json:   Parsed correction JSON dict.
      - approved:          Boolean flag (True for production, False for staging by default).
      - records_file_path: The full GitHub path to the records file.
    """
    prod_correction_file = f"{PROD_DIR}/{correction_name}_correction.json"
    prod_records_file = f"{PROD_DIR}/{correction_name}_records.txt"

    staging_correction_file = f"{STAGING_DIR}/{correction_name}_correction.json"
    staging_records_file = f"{STAGING_DIR}/{correction_name}_records.txt"

    try:
        # Try production folder first.
        correction_content = get_github_file_content(OWNER, REPO, prod_correction_file)
        correction_json = json.loads(correction_content)
        approved = correction_json.get("approved", True)
        records_file_path = prod_records_file
    except Exception as prod_error:
        logger.info(
            f"Production file for '{correction_name}' not found ({prod_error}). "
            "Trying staging folder."
        )
        correction_content = get_github_file_content(OWNER, REPO, staging_correction_file)
        correction_json = json.loads(correction_content)
        approved = correction_json.get("approved", False)
        records_file_path = staging_records_file

    return correction_json, approved, records_file_path


# ---------------------------------------------------------------------------
# Corrections index — built once, cached for CACHE_TTL
# ---------------------------------------------------------------------------
def _build_corrections_index():
    """
    Fetch all correction data from GitHub and build an efficient lookup index.

    Returns dict::

        {
            "by_id": {
                "<record_id>": [
                    {"name": str, "organizations": list, "approved": bool}, ...
                ]
            },
            "by_funding": [
                {
                    "name": str,
                    "patterns": ["AI123456", ...],   # uppercased
                    "organizations": list,
                    "approved": bool
                }, ...
            ]
        }
    """
    index = {"by_id": {}, "by_funding": []}
    correction_names = {}  # name -> source_folder ("production" | "staging")

    # --- Discover correction names from production and staging ---
    try:
        prod_files = list_github_files(OWNER, REPO, PROD_DIR)
        for fname in prod_files:
            if fname.endswith("_records.txt"):
                name = fname[: -len("_records.txt")]
                correction_names[name] = "production"
    except Exception as e:
        logger.error(f"Error listing production corrections: {e}")

    try:
        staging_files = list_github_files(OWNER, REPO, STAGING_DIR)
        for fname in staging_files:
            if fname.endswith("_records.txt"):
                name = fname[: -len("_records.txt")]
                if name not in correction_names:
                    correction_names[name] = "staging"
    except Exception as e:
        logger.error(f"Error listing staging corrections: {e}")

    logger.info(f"Discovered {len(correction_names)} correction(s) to load")

    # --- Load each correction ---
    for correction_name in list(correction_names.keys()):
        try:
            correction_json, approved, records_file_path = fetch_correction_files(correction_name)
            correction_orgs = correction_json.get("sourceOrganization", [])

            if not correction_orgs:
                logger.warning(f"Correction '{correction_name}' has no sourceOrganization - skipping")
                continue

            # -- ID-based index --
            try:
                records_content = get_github_file_content(OWNER, REPO, records_file_path)
                record_ids = get_record_ids(records_content)
            except Exception as rec_err:
                logger.warning(f"Could not load records for '{correction_name}': {rec_err}")
                record_ids = []

            for rid in record_ids:
                index["by_id"].setdefault(rid, []).append(
                    {
                        "name": correction_name,
                        "organizations": correction_orgs,
                        "approved": approved,
                    }
                )

            # -- Funding-based index --
            funding_ids = correction_json.get("fundingIdentifiers", [])
            if funding_ids:
                # Normalize patterns to uppercase for case-insensitive matching.
                index["by_funding"].append(
                    {
                        "name": correction_name,
                        "patterns": [fid.upper() for fid in funding_ids if fid],
                        "organizations": correction_orgs,
                        "approved": approved,
                    }
                )

            logger.info(
                f"Loaded correction '{correction_name}': "
                f"{len(record_ids)} record IDs, "
                f"{len(funding_ids)} funding patterns "
                f"(approved={approved})"
            )
        except Exception as e:
            logger.error(f"Error loading correction '{correction_name}': {e}")

    id_count = sum(len(v) for v in index["by_id"].values())
    logger.info(
        f"Corrections index built: {id_count} ID mappings, "
        f"{len(index['by_funding'])} funding-based entries"
    )
    return index


def get_corrections_index():
    """
    Return the cached corrections index, rebuilding it if stale or missing.
    Thread-safe; falls back to a stale cache on network errors.
    """
    global _corrections_cache, _corrections_cache_time

    with _cache_lock:
        now = time.time()
        if _corrections_cache is not None and (now - _corrections_cache_time) < CACHE_TTL:
            return _corrections_cache
        stale_cache = _corrections_cache

    # Build outside the lock (network I/O can be slow).
    try:
        index = _build_corrections_index()
    except Exception as e:
        logger.error(f"Failed to build corrections index: {e}")
        if stale_cache is not None:
            logger.warning("Falling back to stale corrections cache")
            return stale_cache
        # Return empty index so documents can still flow through.
        return {"by_id": {}, "by_funding": []}

    with _cache_lock:
        _corrections_cache = index
        _corrections_cache_time = time.time()

    return index


# ---------------------------------------------------------------------------
# Funding-based matching
# ---------------------------------------------------------------------------
def _match_funding(doc_funding, funding_patterns):
    """
    Return True if any of the document's funding identifiers contain
    any of the correction's funding patterns as a substring.

    Args:
        doc_funding:      The document's ``funding`` field (dict, list, or None).
        funding_patterns: List of **uppercased** grant patterns to match against.
    """
    if not doc_funding or not funding_patterns:
        return False

    if isinstance(doc_funding, dict):
        doc_funding = [doc_funding]

    if not isinstance(doc_funding, list):
        return False

    for fund_entry in doc_funding:
        if not isinstance(fund_entry, dict):
            continue
        identifier = fund_entry.get("identifier")
        if not identifier or not isinstance(identifier, str):
            continue
        identifier_upper = identifier.upper()
        for pattern in funding_patterns:
            if pattern in identifier_upper:
                return True
    return False


# ---------------------------------------------------------------------------
# Per-document correction — the main public API
# ---------------------------------------------------------------------------
def apply_corrections(doc):
    """
    Apply all matching corrections to a single document.

    Matching strategies (applied in order):
      1. **ID-based** - doc ``_id`` appears in a correction's ``_records.txt``.
         This is the backward-compatible path that uses the pre-curated lists.
      2. **Funding-based** - doc ``funding.identifier`` contains a grant pattern
         from a correction's ``fundingIdentifiers`` list.
         This catches newly-added records immediately, without waiting for the
         external workflow to regenerate records lists.

    Deduplication ensures that the same correction is never applied twice
    (even if both strategies match), and ``update_source_organization``
    prevents duplicate organizations from being added.

    Args:
        doc: A document dict.

    Returns:
        The same document dict, potentially with sourceOrganization updated.
    """
    try:
        index = get_corrections_index()
    except Exception as e:
        logger.error(f"Failed to get corrections index - skipping corrections: {e}")
        return doc

    doc_id = doc.get("_id", "").lower()
    applied_corrections = set()

    # --- Strategy 1: ID-based matching ---
    if doc_id and doc_id in index["by_id"]:
        for correction in index["by_id"][doc_id]:
            cname = correction["name"]
            if cname not in applied_corrections:
                doc = update_source_organization(doc, correction["organizations"], correction["approved"])
                applied_corrections.add(cname)
                logger.debug(
                    f"Applied correction '{cname}' to '{doc_id}' "
                    f"(matched by ID, approved={correction['approved']})"
                )

    # --- Strategy 2: Funding-based matching ---
    doc_funding = doc.get("funding")
    if doc_funding:
        for correction in index["by_funding"]:
            cname = correction["name"]
            if cname not in applied_corrections:
                if _match_funding(doc_funding, correction["patterns"]):
                    doc = update_source_organization(doc, correction["organizations"], correction["approved"])
                    applied_corrections.add(cname)
                    logger.debug(
                        f"Applied correction '{cname}' to '{doc_id}' "
                        f"(matched by funding, approved={correction['approved']})"
                    )

    return doc


# ---------------------------------------------------------------------------
# Document loading helper
# ---------------------------------------------------------------------------
def load_documents(data):
    """
    Yield documents from *data*, which is either:
      - a directory path (str) containing a ``data.ndjson`` file, or
      - an iterable of document dicts.

    Unlike the legacy version, this is a generator - it does NOT load
    everything into memory at once.
    """
    count = 0

    if isinstance(data, str):
        ndjson_file = os.path.join(data, "data.ndjson")
        with open(ndjson_file, "rb") as f:
            for line in f:
                doc = orjson.loads(line)
                yield doc
                count += 1
                if count % 10000 == 0:
                    logging.info(f"Loaded {count} documents")
    else:
        yield from data


# ---------------------------------------------------------------------------
# Legacy / backward-compatible entry point
# ---------------------------------------------------------------------------
def corrections(data):
    """
    Load documents and apply corrections.

    This is the **legacy** entry point kept for backward compatibility.
    Prefer calling ``apply_corrections(doc)`` per-document in the upload
    wrapper instead.  If both are used, deduplication in
    ``update_source_organization`` prevents duplicate organizations.
    """
    for doc in load_documents(data):
        yield apply_corrections(doc)
