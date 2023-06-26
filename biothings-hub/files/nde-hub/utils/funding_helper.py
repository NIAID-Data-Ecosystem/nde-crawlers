import os
import re
import sqlite3

import orjson
import requests
from config import logger

DB_PATH = "data/nde-hub/standardizers/funding_lookup/funding_lookup.db"


def create_sqlite_db(DB_PATH):
    """
    Create the sqlite database.
    :param DB_PATH: the path to the sqlite database
    :return: None
    """
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute(
        """CREATE TABLE funding_lookup
                 (funding_id text, funding text)"""
    )
    conn.commit()
    conn.close()


def update_sqlite_db(funding_id, new_funding):
    """
    Update the sqlite database with the new funding information.
    :param funding_id: the funding id
    :param new_funding: the new funding information
    :return: None
    """
    logger.info(f"Updating funding information for {funding_id} in sqlite database.")
    funding_id = funding_id.replace(" ", "").lower()

    if not os.path.exists(DB_PATH):
        logger.info("Creating sqlite database.")
        create_sqlite_db(DB_PATH)

    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()

    c.execute("INSERT INTO funding_lookup VALUES (?, ?)", (funding_id, orjson.dumps(new_funding).decode("utf-8")))

    conn.commit()
    conn.close()

    logger.info(f"Successfully updated funding information for {funding_id}.")


def sqlite_lookup(funding_id):
    """
    Look up the funding information for a given funding id.
    :param funding_id: the funding id
    :return: the funding information
    """
    logger.info(f"Looking up funding information for {funding_id} in sqlite database.")
    funding_id = funding_id.replace(" ", "").lower()

    if not os.path.exists(DB_PATH):
        logger.info("Creating sqlite database.")
        create_sqlite_db(DB_PATH)

    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()

    c.execute("SELECT funding FROM funding_lookup WHERE funding_id=?", (funding_id,))
    result = c.fetchone()

    conn.close()

    if result:
        logger.info(f"Successfully found funding information for {funding_id}.")
        return orjson.loads(result[0])
    else:
        logger.info(f"No funding information found for {funding_id}.")
        return None



def standardize_funding(data):
    """
    Standardize funding information for all documents in a data source.
    :param data: a list of documents or a path to a data.ndjson file
    :return: a list of documents with standardized funding information
    """
    count = 0

    if isinstance(data, str):
        with open(os.path.join(data, "data.ndjson"), "rb") as f:

            doc_list = []

            for line in f:
                doc = orjson.loads(line)
                count += 1
                if count % 1000 == 0:
                    logger.info(f"Processed {count} documents")

                if isinstance(doc.get("funding", {}), list):
                    for i, funding_dict in enumerate(doc["funding"]):
                        if "identifier" in funding_dict:
                            funding_id = funding_dict["identifier"]
                            if sqlite_lookup(funding_id):
                                doc["funding"][i] = sqlite_lookup(funding_id)
                            else:
                                new_funding = update_funding(funding_id)
                                if new_funding:
                                    update_sqlite_db(funding_id, new_funding)
                                    doc["funding"][i] = new_funding
                                else:
                                    update_sqlite_db(funding_id, None)

                elif funding_id := doc.get("funding", {}).get("identifier"):
                    if sqlite_lookup(funding_id):
                        doc["funding"] = sqlite_lookup(funding_id)
                    else:
                        new_funding = update_funding(funding_id)
                        if new_funding:
                            update_sqlite_db(funding_id, new_funding)
                            doc["funding"] = new_funding
                        else:
                            update_sqlite_db(funding_id, None)

                doc_list.append(doc)

            logger.info(f"Finished processing {count} documents")

            return doc_list

    else:
        doc_list = list(data)

        for doc in doc_list:
            count += 1
            if count % 1000 == 0:
                logger.info(f"Processed {count} documents")

            if isinstance(doc.get("funding", {}), list):
                for i, funding_dict in enumerate(doc["funding"]):
                    if "identifier" in funding_dict:
                        funding_id = funding_dict["identifier"]
                        if sqlite_lookup(funding_id):
                            doc["funding"][i] = sqlite_lookup(funding_id)
                        else:
                            new_funding = update_funding(funding_id)
                            if new_funding:
                                update_sqlite_db(funding_id, new_funding)
                                doc["funding"][i] = new_funding
                            else:
                                update_sqlite_db(funding_id, None)

            elif funding_id := doc.get("funding", {}).get("identifier"):
                if sqlite_lookup(funding_id):
                    doc["funding"] = sqlite_lookup(funding_id)
                else:
                    new_funding = update_funding(funding_id)
                    if new_funding:
                        update_sqlite_db(funding_id, new_funding)
                        doc["funding"] = new_funding
                    else:
                        update_sqlite_db(funding_id, None)

        return doc_list


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
                "AbstractText",
                "Terms",
                "FullFoa",
                "ProgramOfficers",
                "AwardAmount",
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
            logger.error(f"ERROR for #{count} request, skipping...")
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


def build_funding_dict(funding_info):
    """
    Build a correctly mapped funding dictionary from the funding information.
    :param funding_info: the funding information
    :return: the funding dictionary
    """

    ic_lookup = {
        "CC": "Clinical Center",
        "RG": "Center for Scientific Review",
        "CIT": "Center for Information Technology",
        "TW": "John E. Fogarty International Center",
        "TR": "National Center for Advancing Translational Sciences (NCATS)",
        "AT": "National Center for Complementary and Integrative Health",
        "CA": "National Cancer Institute",
        "RR": "National Center for Research Resources (dissolved 12/2011)",
        "EY": "National Eye Institute",
        "HG": "National Human Genome Research Institute",
        "HL": "National Heart, Lung, and Blood Institute",
        "AG": "National Institute on Aging",
        "AA": "National Institute on Alcohol Abuse and Alcoholism",
        "AI": "National Institute of Allergy and Infectious Diseases",
        "AR": "National Institute of Arthritis and Musculoskeletal and Skin Diseases",
        "EB": "National Institute of Biomedical Imaging and Bioengineering",
        "HD": "Eunice Kennedy Shriver National Institute of Child Health and Human Development",
        "DA": "National Institute on Drug Abuse",
        "DC": "National Institute on Deafness and Other Communication Disorders",
        "DE": "National Institute of Dental and Craniofacial Research",
        "DK": "National Institute of Diabetes and Digestive and Kidney Diseases",
        "ES": "National Institute of Environmental Health Sciences",
        "GM": "National Institute of General Medical Sciences",
        "MH": "National Institute of Mental Health",
        "MD": "National Institute on Minority Health and Health Disparities",
        "NS": "National Institute of Neurological Disorders and Stroke",
        "NR": "National Institute of Nursing Research",
        "LM": "National Library of Medicine",
        "OD": "Office of the Director",
    }

    funding_dict = {}
    funder_dict = {}

    if appl_id := funding_info.get("appl_id"):
        funding_dict["url"] = f"https://reporter.nih.gov/project-details/{appl_id}"
    if project_num := funding_info.get("project_num"):
        funding_dict["identifier"] = project_num
    if project_title := funding_info.get("project_title"):
        funding_dict["name"] = project_title
    if abstract_text := funding_info.get("abstract_text"):
        funding_dict["description"] = abstract_text
    if project_num_split := funding_info.get("project_num_split"):
        if ic_code := project_num_split.get("ic_code"):
            try:
                funder_dict["name"] = ic_lookup[ic_code]
            except KeyError:
                logger.info(f"IC CODE {ic_code} NOT FOUND IN LOOKUP")
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
            funder_dict["employee"] = employees
    if project_start_date := funding_info.get("project_start_date"):
        funding_dict["startDate"] = project_start_date.split("T")[0]
    if project_end_date := funding_info.get("project_end_date"):
        funding_dict["endDate"] = project_end_date.split("T")[0]
    if terms := funding_info.get("terms"):
        terms = terms.replace("<", "").replace(">", " ").split()
        funding_dict["keywords"] = terms

    if full_foa := funding_info.get("full_foa"):
        funding_dict["isBasedOn"] = {"identifier": full_foa}

    if bool(funder_dict):
        funding_dict["funder"] = funder_dict
    return funding_dict
