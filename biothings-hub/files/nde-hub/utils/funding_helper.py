import json
import os
import re
import sqlite3

import orjson
import requests
from config import logger

from .utils import retry

DB_PATH = "/data/nde-hub/standardizers/funding_lookup/funding_lookup.db"


def create_sqlite_db(conn):
    """
    Create the SQLite database tables if they don't exist.
    :param conn: SQLite connection object
    :return: None
    """
    c = conn.cursor()

    # Create the funding_lookup table
    c.execute(
        """CREATE TABLE IF NOT EXISTS funding_lookup (
                 funding_id TEXT PRIMARY KEY,
                 funding TEXT
           )"""
    )

    # Create the funder_cache table
    c.execute(
        """CREATE TABLE IF NOT EXISTS funder_cache (
                 funder_name TEXT PRIMARY KEY,
                 funder_data TEXT
           )"""
    )

    conn.commit()


def update_sqlite_db(conn, funding_id, new_funding):
    """
    Update the SQLite database with the new funding information.
    :param conn: SQLite connection object
    :param funding_id: the funding id
    :param new_funding: the new funding information
    :return: None
    """
    logger.info(f"Updating funding information for {funding_id} in SQLite database.")
    funding_id = funding_id.replace(" ", "").lower()

    c = conn.cursor()
    c.execute(
        "INSERT OR REPLACE INTO funding_lookup (funding_id, funding) VALUES (?, ?)",
        (funding_id, orjson.dumps(new_funding).decode("utf-8")),
    )
    conn.commit()

    logger.info(f"Successfully updated funding information for {funding_id}.")


def batch_sqlite_lookup(conn, funding_ids):
    """
    Perform a batch lookup of funding IDs.
    :param conn: SQLite connection object
    :param funding_ids: a set of funding IDs
    :return: a dictionary mapping funding IDs to funding data
    """
    logger.info(f"Performing batch lookup for {len(funding_ids)} funding IDs.")
    placeholders = ",".join("?" for _ in funding_ids)
    c = conn.cursor()
    c.execute(f"SELECT funding_id, funding FROM funding_lookup WHERE funding_id IN ({placeholders})", list(funding_ids))
    result = c.fetchall()
    funding_cache = {row[0]: orjson.loads(row[1]) for row in result}
    logger.info(f"Found {len(funding_cache)} funding records in the database.")
    return funding_cache


def standardize_funding(data):
    """
    Standardize funding information for all documents in a data source.
    :param data: a list of documents or a path to a data.ndjson file
    :return: a list of documents with standardized funding information
    """
    count = 0

    # Open a persistent SQLite connection
    conn = sqlite3.connect(DB_PATH)

    try:
        # Ensure the database exists
        create_sqlite_db(conn)

        # Check if data is a file path to an NDJSON file
        if isinstance(data, str):
            doc_list = []
            with open(os.path.join(data, "data.ndjson"), "rb") as f:
                for line in f:
                    doc = orjson.loads(line)
                    doc_list.append(doc)
                    count += 1
                    if count % 1000 == 0:
                        logger.info(f"Processed {count} documents")
        else:
            doc_list = list(data)

        # First we collect all unique funding IDs
        funding_ids = set()
        for doc in doc_list:
            funding_field = doc.get("funding", {})
            if isinstance(funding_field, list):
                for funding_dict in funding_field:
                    if "identifier" in funding_dict:
                        funding_id = funding_dict["identifier"].replace(" ", "").lower()
                        funding_ids.add(funding_id)
            elif isinstance(funding_field, dict):
                if "identifier" in funding_field:
                    funding_id = funding_field["identifier"].replace(" ", "").lower()
                    funding_ids.add(funding_id)
        logger.info(f"Unique funding IDs collected: {len(funding_ids)}")

        # Batch query the database
        funding_cache_dict = batch_sqlite_lookup(conn, funding_ids)

        # Standardize funding information
        for doc in doc_list:
            count += 1
            if count % 1000 == 0:
                logger.info(f"Processed {count} documents")

            if isinstance(doc.get("funding", {}), list):
                for i, funding_dict in enumerate(doc["funding"]):
                    if "identifier" in funding_dict:
                        funding_id = funding_dict["identifier"].replace(" ", "").lower()
                        funding_cache = funding_cache_dict.get(funding_id)
                        if funding_cache:
                            doc["funding"][i] = funding_cache
                        else:
                            logger.info(f"Not in cache: {funding_id}, skipping...")
                            continue
                            try:
                                new_funding = update_funding(funding_id)
                            except Exception as e:
                                logger.error(f"ERROR for {doc['_id']} request, skipping...")
                                logger.error(e)
                                continue
                            if new_funding:
                                update_sqlite_db(conn, funding_id, new_funding)
                                doc["funding"][i] = new_funding
                    elif isinstance(funding_dict.get("funder", {}), dict):
                        if funder_name := funding_dict.get("funder", {}).get("name"):
                            try:
                                doc["funding"][i]["funder"] = standardize_funder(funder_name)
                            except Exception as e:
                                logger.error(f"ERROR for {doc['_id']} request, skipping...")
                                logger.error(e)
                                continue
                    elif isinstance(funding_dict.get("funder", {}), list):
                        for j, funder_dict in enumerate(funding_dict["funder"]):
                            if funder_name := funder_dict.get("name"):
                                try:
                                    doc["funding"][i]["funder"][j] = standardize_funder(funder_name)
                                except Exception as e:
                                    logger.error(f"ERROR for {doc['_id']} request, skipping...")
                                    logger.error(e)
                                    continue

            elif funding_id := doc.get("funding", {}).get("identifier"):
                funding_id = funding_id.replace(" ", "").lower()
                funding_cache = funding_cache_dict.get(funding_id)
                if funding_cache:
                    doc["funding"] = funding_cache
                else:
                    logger.info(f"Not in cache: {funding_id}, skipping...")
                    continue
                    try:
                        new_funding = update_funding(funding_id)
                    except Exception as e:
                        logger.error(f"ERROR for {doc['_id']} request, skipping...")
                        logger.error(e)
                        continue
                    if new_funding:
                        update_sqlite_db(conn, funding_id, new_funding)
                        doc["funding"] = new_funding
            elif isinstance(doc.get("funding", {}).get("funder", {}), dict):
                if funder_name := doc.get("funding", {}).get("funder", {}).get("name"):
                    try:
                        doc["funding"]["funder"] = standardize_funder(funder_name)
                    except Exception as e:
                        logger.error(f"ERROR for {doc['_id']} request, skipping...")
                        logger.error(e)
                        continue
            elif isinstance(doc.get("funding", {}).get("funder", {}), list):
                for i, funder_dict in enumerate(doc["funding"]["funder"]):
                    if funder_name := funder_dict.get("name"):
                        try:
                            doc["funding"]["funder"][i] = standardize_funder(funder_name)
                        except Exception as e:
                            logger.error(f"ERROR for {doc['_id']} request, skipping...")
                            logger.error(e)
                            continue

        return doc_list
    finally:
        conn.close()


@retry(3, 5)
def update_funding(funding_id):
    """
    Update funding information for a given funding id.
    :param funding_id: the funding id
    :return: the funding information
    """
    count = 0

    # TODO - remove this if statement when funding ids are fixed
    if len(funding_id) < 5:
        logger.info(f"INVALID FUNDING ID FOR {funding_id}")
        return None
    if "NIH" in funding_id:
        funding_id = funding_id.replace("NIH", "")

    data = []
    offset = 0
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
            logger.error(f"ERROR for {count} request, skipping...")
            logger.error(e)
            continue
        if funding_data == ["Invalid project number"]:
            logger.info(f"INVALID PROJECT NUMBER FOR {funding_id}")
            return None
        elif funding_data["meta"]["total"] == 0:
            logger.info(f"NO RESULTS FOUND FOR {funding_id}")
            return None
        else:
            data.extend(funding_data["results"])
            offset += 500
            if offset >= funding_data["meta"]["total"]:
                break

    logger.info(f"FOUND {len(data)} RESULTS FOR {funding_id}")

    parent = []
    year = data[0]["fiscal_year"]
    for obj in data:
        if year != obj["fiscal_year"]:
            break
        if obj["subproject_id"] == None:
            parent.append(obj)
    if len(parent) == 0:
        logger.info(f"NO PARENT PROJECT FOUND FOR {funding_id}")
    elif len(parent) == 1:
        logger.info(f"FOUND PARENT PROJECT FOR {funding_id}")
        return build_funding_dict(parent[0])
    elif len(parent) > 1:
        logger.info(f"MULTIPLE PARENT PROJECTS FOUND FOR {funding_id}")
        parent = [i for i in parent if i["award_amount"] != None]
        parent.sort(key=lambda x: x["award_amount"], reverse=True)
        if len(parent) == 0:
            logger.info(f"NO PARENT PROJECTS WITH AWARD AMOUNT FOUND FOR {funding_id}")
            logger.info(f"USING APPLICATION TYPE 1 FOR {funding_id}")
            parent = [i for i in parent if i["project_num_split"]["appl_type_code"] == "1"]
            return build_funding_dict(parent[0])
        largest_amount = parent[0]["award_amount"]
        # get all parents with the largest amount
        largest_amount_parents = [i for i in parent if i["award_amount"] == largest_amount]
        # if there is only one parent with the largest amount, use it
        if len(largest_amount_parents) == 1:
            logger.info(f"USING LARGEST AWARD AMOUNT PARENT FOR {funding_id}")
            return build_funding_dict(largest_amount_parents[0])
        # if there are more than one parent with the largest amount, use the one application_type is 1
        elif len(largest_amount_parents) > 1:
            logger.info(f"MULTIPLE PARENTS WITH LARGEST AMOUNT FOUND FOR {funding_id}")
            for large_parent in largest_amount_parents:
                if (
                    large_parent["project_num_split"] is not None
                    and large_parent["project_num_split"]["appl_type_code"] == "1"
                ):
                    logger.info(f"USING PARENT WITH APPLICATION TYPE 1 FOR {funding_id}")
                    return build_funding_dict(large_parent)
            logger.info(f"NO PARENT WITH APPLICATION TYPE 1 FOUND FOR {funding_id}, USING FIRST INDEX")
            return build_funding_dict(largest_amount_parents[0])


@retry(3, 5)
def standardize_funder(funder):
    """
    Standardize the funder dictionary using a SQLite cache.
    :param funder: the funder name
    :return: the funder dictionary
    """
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    cursor.execute("SELECT funder_data FROM funder_cache WHERE funder_name = ?", (funder,))
    result = cursor.fetchone()

    if result:
        logger.info(f"FOUND FUNDING INFORMATION FOR {funder} in SQLite cache")
        return json.loads(result[0])

    funder_dict = {}
    funder_name = re.sub(r"\([^)]*\)", "", funder).strip()
    funder_name = funder_name.replace("&", "and")
    funder_acronym_list = re.findall(r"\(([^)]*)\)", funder)
    funder_acronym = funder_acronym_list[0] if funder_acronym_list else ""

    url = f"https://api.crossref.org/funders?query={funder_name}"
    response = requests.get(url)
    data = response.json()

    if "message" not in data or "items" not in data["message"]:
        logger.info(f"No message in response for {funder_name}, {url}")
        funder_dict = {"name": funder, "@type": "Organization"}
    else:
        for item in data["message"]["items"]:
            if item["name"].lower() == funder_name.lower() or funder_acronym.lower() in [
                alt_name.lower() for alt_name in item.get("alt-names", [])
            ]:
                funder_dict["name"] = item["name"]
                funder_dict["alternateName"] = item.get("alt-names", [])
                funder_dict["identifier"] = item["id"]
                break

    if not funder_dict:
        logger.info(f"NO FUNDING INFORMATION FOUND FOR {funder_name}, {url}")
        funder_dict = {"name": funder, "@type": "Organization"}

    cursor.execute(
        "INSERT OR REPLACE INTO funder_cache (funder_name, funder_data) VALUES (?, ?)",
        (funder, json.dumps(funder_dict)),
    )
    conn.commit()

    logger.info(f"Standardized and saved funding information for {funder_name}")
    return funder_dict


def build_funding_dict(funding_info):
    """
    Build a correctly mapped funding dictionary from the funding information.
    :param funding_info: the funding information
    :return: the funding dictionary
    """

    funding_dict = {}
    funders = []

    if appl_id := funding_info.get("appl_id"):
        funding_dict["url"] = f"https://reporter.nih.gov/project-details/{appl_id}"
    if project_num := funding_info.get("project_num"):
        funding_dict["identifier"] = project_num
    if project_title := funding_info.get("project_title"):
        funding_dict["name"] = project_title
    if agency_ic_fundings := funding_info.get("agency_ic_fundings"):
        for ic_funding in agency_ic_fundings:
            standardized_funder_dict = standardize_funder(ic_funding["name"])
            if standardized_funder_dict["name"] == funding_info.get("agency_ic_admin").get("name"):
                if program_officers := funding_info.get("program_officers"):
                    employees = []
                    for officer in program_officers:
                        employee = {}
                        if first_name := officer.get("first_name"):
                            employee["givenName"] = first_name
                        if last_name := officer.get("last_name"):
                            employee["familyName"] = last_name
                        if first_name and last_name:
                            employee["name"] = f"{first_name} {last_name}"
                        if email := officer.get("email"):
                            if re.match(r"[^@]+@[^@]+\.[^@]+", email):
                                employee["email"] = email
                        employees.append(employee)
                    if len(employees) > 0:
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
