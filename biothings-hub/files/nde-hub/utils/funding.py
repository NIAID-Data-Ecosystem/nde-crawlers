"""Funding standardization.

Runs whenever a batch of records carries `funding`:

  * `funding.identifier` is looked up in the funding cache and, on a hit, the
    whole entry is replaced with the curated NIH RePORTER grant,
  * otherwise each `funding.funder` name is standardized against CrossRef.

Both lookups are cached in SQLite, so a repeat run costs no HTTP requests. The
cache is populated offline by `update_funding` -- upload never calls the
RePORTER API itself.
"""

import json
import re

import requests
from config import logger

from .cache import SqliteCache
from .common import as_list

DB_PATH = "/data/nde-hub/standardizers/funding_lookup/funding_lookup.db"

# normalize=False: both tables are keyed by the caller's own spelling -- grants by
# `_funding_key` (which strips every space, not just the ends) and funders by the
# raw name the source supplied.
FUNDING_LOOKUP = SqliteCache(DB_PATH, "funding_lookup", "funding_id", "funding", memoize=False, normalize=False)
FUNDERS = SqliteCache(DB_PATH, "funder_cache", "funder_name", "funder_data", normalize=False)


def reset_caches():
    """Drop the process-wide caches so a new upload sees fresh lookup data."""
    FUNDING_LOOKUP.reset()
    FUNDERS.reset()


def _funding_key(identifier):
    return identifier.replace(" ", "").lower()


def standardize_funder(funder):
    """Standardize a funder name against CrossRef, caching hits and misses."""
    if cached := FUNDERS.get(funder):
        # Older rows predate @type always being set.
        cached.setdefault("@type", "Organization")
        return cached

    funder_name = re.sub(r"\([^)]*\)", "", funder).strip().replace("&", "and")
    funder_acronym_list = re.findall(r"\(([^)]*)\)", funder)
    funder_acronym = funder_acronym_list[0] if funder_acronym_list else ""
    url = f"https://api.crossref.org/funders?query={funder_name}"

    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
    except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
        logger.error("Error fetching funder from CrossRef API: %s", e)
        return {"name": funder, "@type": "Organization"}

    funder_dict = {}
    if "message" not in data or "items" not in data["message"]:
        logger.debug("No message in response for %s, %s", funder_name, url)
    else:
        for item in data["message"]["items"]:
            alt_names = [alt_name.lower() for alt_name in item.get("alt-names", [])]
            if item["name"].lower() == funder_name.lower() or (funder_acronym and funder_acronym.lower() in alt_names):
                funder_dict = {
                    "@type": "Organization",
                    "name": item["name"],
                    "alternateName": item.get("alt-names", []),
                    "identifier": item["id"],
                }
                break

    if not funder_dict:
        logger.debug("NO FUNDING INFORMATION FOUND FOR %s, %s", funder_name, url)
        funder_dict = {"name": funder, "@type": "Organization"}

    FUNDERS.put(funder, funder_dict)
    return funder_dict


def _standardize_funders(entry):
    """Standardize `entry["funder"]`, which may be a single funder or a list."""
    funder = entry.get("funder")
    if isinstance(funder, dict):
        if name := funder.get("name"):
            entry["funder"] = standardize_funder(name)
    elif isinstance(funder, list):
        for i, funder_dict in enumerate(funder):
            if isinstance(funder_dict, dict) and (name := funder_dict.get("name")):
                funder[i] = standardize_funder(name)


def standardize_funding(docs):
    """Standardize the funding of one batch of records."""
    docs = list(docs)

    funding_ids = {
        _funding_key(entry["identifier"])
        for doc in docs
        for entry in as_list(doc.get("funding"))
        if isinstance(entry, dict) and entry.get("identifier")
    }
    funding_cache = FUNDING_LOOKUP.get_many(funding_ids)
    logger.debug("Found %s of %s funding records in the cache", len(funding_cache), len(funding_ids))

    for doc in docs:
        funding = doc.get("funding")
        if not funding:
            yield doc
            continue

        entries = as_list(funding)
        for i, entry in enumerate(entries):
            if not isinstance(entry, dict):
                continue
            if identifier := entry.get("identifier"):
                # A curated grant replaces the whole entry; a miss leaves it untouched.
                if cached := funding_cache.get(_funding_key(identifier)):
                    entries[i] = cached
                else:
                    logger.debug("Not in cache: %s, skipping API lookup", _funding_key(identifier))
                continue
            try:
                _standardize_funders(entry)
            except Exception as e:
                logger.error("ERROR standardizing funder for %s, skipping: %s", doc.get("_id", "unknown id"), e)

        doc["funding"] = entries if isinstance(funding, list) else entries[0]
        yield doc


# ---------------------------------------------------------------------------
# Offline cache population (NIH RePORTER)
# ---------------------------------------------------------------------------
def update_sqlite_db(funding_id, new_funding):
    """Store a curated grant in the funding cache."""
    logger.debug("Updating funding information for %s in the funding cache.", funding_id)
    FUNDING_LOOKUP.put(_funding_key(funding_id), new_funding)


def is_valid_nih_funding_id(funding_id):
    """Quick plausibility check that filters European grants, prose and markdown junk."""
    stripped = funding_id.strip()

    # Too short or too long to be a real project number
    if len(stripped) < 5 or len(stripped) > 50:
        return False
    # Non-ASCII means European/foreign grant IDs
    if not stripped.isascii():
        return False
    # Must contain at least one digit
    if not re.search(r"\d", stripped):
        return False
    # Starts with junk characters
    if stripped[0] in "\"'#_":
        return False
    # Too many words -- it's a description, not an ID
    if len(stripped.split()) > 4:
        return False
    # After the same cleanup update_funding does, should be alphanumeric
    return stripped.replace("NIH", "").replace("-", "").replace(" ", "").isalnum()


def update_funding(funding_id):
    """Fetch the parent project for `funding_id` from NIH RePORTER."""
    if not is_valid_nih_funding_id(funding_id):
        logger.debug("INVALID FUNDING ID FOR %s", funding_id)
        return None

    funding_id = funding_id.replace("NIH", "")
    # Remove dashes between alphanumeric segments (e.g. U01-HL -> U01HL)
    # but keep the suffix dash before digits (e.g. -01, -01A1)
    funding_id = re.sub(r"-(?=[A-Za-z])", "", funding_id)

    data = []
    offset = 0
    count = 0
    while True:
        count += 1
        params = {
            "criteria": {"project_nums": [funding_id]},
            "include_fields": [
                "ProjectTitle",
                "ApplId",
                "SubprojectId",
                "FiscalYear",
                "Organization",
                "ProjectNum",
                "OrgCountry",
                "ProjectNumSplit",
                "PrincipalInvestigators",
                "AllText",
                "FullStudySection",
                "ProjectStartDate",
                "ProjectEndDate",
                "FullFoa",
                "ProgramOfficers",
                "AwardAmount",
                "AgencyIcFundings",
                "AgencyIcAdmin",
            ],
            "offset": offset,
            "limit": 500,
            "sort_field": "fiscal_year",
            "sort_order": "asc",
        }
        response = requests.post("https://api.reporter.nih.gov/v2/projects/search", json=params)
        try:
            funding_data = response.json()
        except Exception as e:
            logger.error("ERROR for %s request, skipping: %s", count, e)
            continue
        if funding_data == ["Invalid project number"]:
            logger.debug("INVALID PROJECT NUMBER FOR %s", funding_id)
            return None
        if funding_data["meta"]["total"] == 0:
            logger.debug("NO RESULTS FOUND FOR %s", funding_id)
            return None
        data.extend(funding_data["results"])
        offset += 500
        if offset >= funding_data["meta"]["total"]:
            break

    logger.debug("FOUND %s RESULTS FOR %s", len(data), funding_id)
    year = data[0]["fiscal_year"]
    parent = [obj for obj in data if obj["fiscal_year"] == year and obj["subproject_id"] is None]

    if not parent:
        logger.debug("NO PARENT PROJECT FOUND FOR %s", funding_id)
        return None
    if len(parent) == 1:
        logger.debug("FOUND PARENT PROJECT FOR %s", funding_id)
        return build_funding_dict(parent[0])

    logger.debug("MULTIPLE PARENT PROJECTS FOUND FOR %s", funding_id)
    parent = sorted(
        (i for i in parent if i["award_amount"] is not None),
        key=lambda x: x["award_amount"],
        reverse=True,
    )
    if not parent:
        logger.debug("NO PARENT PROJECTS WITH AWARD AMOUNT FOUND FOR %s", funding_id)
        return None

    largest_amount_parents = [i for i in parent if i["award_amount"] == parent[0]["award_amount"]]
    if len(largest_amount_parents) == 1:
        logger.debug("USING LARGEST AWARD AMOUNT PARENT FOR %s", funding_id)
        return build_funding_dict(largest_amount_parents[0])

    logger.debug("MULTIPLE PARENTS WITH LARGEST AMOUNT FOUND FOR %s", funding_id)
    for large_parent in largest_amount_parents:
        split = large_parent.get("project_num_split")
        if split is not None and split["appl_type_code"] == "1":
            logger.debug("USING PARENT WITH APPLICATION TYPE 1 FOR %s", funding_id)
            return build_funding_dict(large_parent)
    logger.debug("NO PARENT WITH APPLICATION TYPE 1 FOUND FOR %s, USING FIRST INDEX", funding_id)
    return build_funding_dict(largest_amount_parents[0])


def build_funding_dict(funding_info):
    """Map a RePORTER project onto our MonetaryGrant shape."""
    funding_dict = {"@type": "MonetaryGrant"}
    if appl_id := funding_info.get("appl_id"):
        funding_dict["url"] = f"https://reporter.nih.gov/project-details/{appl_id}"
    if project_num := funding_info.get("project_num"):
        funding_dict["identifier"] = project_num
    if project_title := funding_info.get("project_title"):
        funding_dict["name"] = project_title

    funders = []
    for ic_funding in funding_info.get("agency_ic_fundings") or []:
        standardized_funder_dict = standardize_funder(ic_funding["name"])
        if standardized_funder_dict["name"] == funding_info.get("agency_ic_admin", {}).get("name"):
            if employees := _build_employees(funding_info.get("program_officers")):
                standardized_funder_dict["employee"] = employees
        funders.append(standardized_funder_dict)

    if project_start_date := funding_info.get("project_start_date"):
        funding_dict["startDate"] = project_start_date.split("T")[0]
    if project_end_date := funding_info.get("project_end_date"):
        funding_dict["endDate"] = project_end_date.split("T")[0]
    if full_foa := funding_info.get("full_foa"):
        funding_dict["isBasedOn"] = {"identifier": full_foa}
    if funders:
        funding_dict["funder"] = funders
    return funding_dict


def _build_employees(program_officers):
    employees = []
    for officer in program_officers or []:
        employee = {}
        first_name = officer.get("first_name")
        last_name = officer.get("last_name")
        if first_name:
            employee["givenName"] = first_name
        if last_name:
            employee["familyName"] = last_name
        if first_name and last_name:
            employee["name"] = f"{first_name} {last_name}"
        if (email := officer.get("email")) and re.match(r"[^@]+@[^@]+\.[^@]+", email):
            employee["email"] = email
        employees.append(employee)
    return employees
