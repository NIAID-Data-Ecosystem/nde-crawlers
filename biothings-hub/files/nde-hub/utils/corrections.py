import fcntl
import json
import math
import os
import tempfile
import time

import config
import requests
from config import logger, token

from .common import as_list

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
OWNER = "NIAID-Data-Ecosystem"
REPO = "nde-metadata-corrections"
PROD_DIR = "collections_corrections_production"
STAGING_DIR = "collections_corrections_staging"

# Each correction is `<name>_records.txt` plus `<name>_correction.json`.
RECORDS_SUFFIX = "_records.txt"

# Every uploader runs in its own process. Keep a short-lived shared snapshot so
# one process fetches correction definitions from GitHub and the others reuse
# exactly the same index for the upload wave.
_CORRECTIONS_CACHE_VERSION = 1
_CORRECTIONS_CACHE_FILENAME = "nde-corrections-index-v1.json"
_CORRECTIONS_CACHE_TTL_SECONDS = 24 * 60 * 60

# ---------------------------------------------------------------------------
# Module-level cache, loaded once per uploader process
# ---------------------------------------------------------------------------
_corrections_cache = None


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
    # Sanitize existing organizations.
    existing_orgs = [sanitize_org(org) for org in as_list(record_metadata.get("sourceOrganization"))]

    # Gather identifiers already present (name or url, lowercased).
    existing_org_identifiers = {
        (org.get("name") or org.get("url")).lower() for org in existing_orgs if (org.get("name") or org.get("url"))
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
    prod_records_file = f"{PROD_DIR}/{correction_name}{RECORDS_SUFFIX}"

    staging_correction_file = f"{STAGING_DIR}/{correction_name}_correction.json"
    staging_records_file = f"{STAGING_DIR}/{correction_name}{RECORDS_SUFFIX}"

    try:
        # Try production folder first.
        correction_content = get_github_file_content(OWNER, REPO, prod_correction_file)
        correction_json = json.loads(correction_content)
        approved = correction_json.get("approved", True)
        records_file_path = prod_records_file
    except Exception as prod_error:
        logger.debug(f"Production file for '{correction_name}' not found ({prod_error}). " "Trying staging folder.")
        correction_content = get_github_file_content(OWNER, REPO, staging_correction_file)
        correction_json = json.loads(correction_content)
        approved = correction_json.get("approved", False)
        records_file_path = staging_records_file

    return correction_json, approved, records_file_path


# ---------------------------------------------------------------------------
# Corrections index — built once per uploader process
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

    # --- Discover correction names from production, then staging ---
    # Each correction is a pair of files sharing a base name; the records file is
    # what we list on, so every correction is discovered exactly once. Production
    # is scanned first and setdefault keeps it winning over a staging namesake.
    for directory, folder in ((PROD_DIR, "production"), (STAGING_DIR, "staging")):
        try:
            for fname in list_github_files(OWNER, REPO, directory):
                if fname.endswith(RECORDS_SUFFIX):
                    correction_names.setdefault(fname.removesuffix(RECORDS_SUFFIX), folder)
        except Exception as e:
            logger.error("Error listing %s corrections: %s", folder, e)

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

            logger.debug(
                f"Loaded correction '{correction_name}': "
                f"{len(record_ids)} record IDs, "
                f"{len(funding_ids)} funding patterns "
                f"(approved={approved})"
            )
        except Exception as e:
            logger.error(f"Error loading correction '{correction_name}': {e}")

    id_count = sum(len(v) for v in index["by_id"].values())
    logger.info(
        f"Corrections index built: {id_count} ID mappings, " f"{len(index['by_funding'])} funding-based entries"
    )
    return index


def _shared_cache_path():
    cache_folder = os.environ.get("CORRECTIONS_CACHE_FOLDER") or getattr(config, "CACHE_FOLDER", None)
    if not cache_folder:
        return None
    return os.path.join(cache_folder, _CORRECTIONS_CACHE_FILENAME)


def _valid_corrections_index(index):
    return (
        isinstance(index, dict) and isinstance(index.get("by_id"), dict) and isinstance(index.get("by_funding"), list)
    )


def _read_shared_cache(path, allow_stale=False):
    with open(path, encoding="utf-8") as cache_file:
        payload = json.load(cache_file)

    if not isinstance(payload, dict):
        return None
    if payload.get("version") != _CORRECTIONS_CACHE_VERSION:
        return None
    index = payload.get("index")
    if not _valid_corrections_index(index):
        return None

    created_at = float(payload["created_at"])
    age = max(0.0, time.time() - created_at)
    if not allow_stale and age > _CORRECTIONS_CACHE_TTL_SECONDS:
        return None
    return index, age


def _write_shared_cache(path, index):
    cache_folder = os.path.dirname(path)
    fd, temporary_path = tempfile.mkstemp(prefix=".nde-corrections-", suffix=".json", dir=cache_folder)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as cache_file:
            json.dump(
                {
                    "version": _CORRECTIONS_CACHE_VERSION,
                    "created_at": time.time(),
                    "index": index,
                },
                cache_file,
                separators=(",", ":"),
            )
        os.replace(temporary_path, path)
    except Exception:
        try:
            os.unlink(temporary_path)
        except FileNotFoundError:
            pass
        raise


def _load_or_build_corrections_index():
    """Load a fresh shared index, or build it once under a process lock."""
    path = _shared_cache_path()
    if path is None:
        return _build_corrections_index()

    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        lock_file = open(f"{path}.lock", "a+", encoding="utf-8")
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
    except OSError as e:
        if "lock_file" in locals():
            lock_file.close()
        logger.warning("Corrections shared cache is unavailable (%s); building a process-local index", e)
        return _build_corrections_index()

    try:
        try:
            cached = _read_shared_cache(path)
        except (FileNotFoundError, json.JSONDecodeError, KeyError, TypeError, ValueError, OSError):
            cached = None
        if cached is not None:
            index, age = cached
            logger.info("Loaded corrections index from shared cache (%.0fs old)", age)
            return index

        try:
            stale = _read_shared_cache(path, allow_stale=True)
        except (FileNotFoundError, json.JSONDecodeError, KeyError, TypeError, ValueError, OSError):
            stale = None

        logger.info("Corrections shared cache is missing or stale; refreshing it from GitHub")
        index = _build_corrections_index()
        if not index["by_id"] and not index["by_funding"]:
            if stale is not None:
                logger.warning("Correction refresh returned an empty index; using the stale shared snapshot")
                return stale[0]
            return index

        try:
            _write_shared_cache(path, index)
        except OSError as e:
            logger.warning("Could not write corrections shared cache: %s", e)
        return index
    finally:
        lock_file.close()


def get_corrections_index():
    """Return the corrections index, loading it once per process.

    The on-disk snapshot is shared by uploader processes for 24 hours. Until a
    refresh succeeds an empty index (or a stale shared index) is returned,
    letting documents flow through rather than failing the upload.
    """
    global _corrections_cache

    if _corrections_cache is None:
        try:
            _corrections_cache = _load_or_build_corrections_index()
        except Exception as e:
            logger.error("Failed to build corrections index: %s", e)
            _corrections_cache = {"by_id": {}, "by_funding": []}

    return _corrections_cache


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

    for fund_entry in as_list(doc_funding):
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
                    f"Applied correction '{cname}' to '{doc_id}' " f"(matched by ID, approved={correction['approved']})"
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
