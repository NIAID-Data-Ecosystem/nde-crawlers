# Helper file to batch call pmids to get citations and funding
# Helpful links to documentation of Biopython Package for writing this file
# https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec162
# https://biopython.org/docs/1.76/api/Bio.Entrez.html
# https://www.nlm.nih.gov/bsd/mms/medlineelements.html
# https://dataguide.nlm.nih.gov/eutilities/utilities.html#efetch
# https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/

import csv
import functools
import gzip
import json
import os
import os.path
import sqlite3
import time
import urllib.error
from copy import copy
from datetime import datetime
from itertools import islice
from typing import Dict, Iterable, Optional

import orjson
import requests
from Bio import Entrez, Medline
from config import GEO_API_KEY, GEO_EMAIL, logger

from .funding_helper import standardize_funder
from .pubtator import DB_PATH, get_species_details, query_condition
from .utils import retry

PMID_DB_PATH = "/data/nde-hub/standardizers/pmid_lookup/pmid_lookup.db"


def remove_first_by_name(lst, target):
    target_lower = target.lower()
    for i, item in enumerate(lst):
        if item.get("name", "").lower() == target_lower:
            lst.pop(i)
            break


def update_species_and_disease(rec, pmids):
    species = get_data_for_pmids(pmids, "species")
    if species:
        update_record_species(rec, species)
    diseases = get_data_for_pmids(pmids, "disease")
    if diseases:
        update_record_disease(rec, diseases)
    return rec


def create_db():
    os.makedirs(os.path.dirname(PMID_DB_PATH))
    # Create or open a SQLite database
    conn = sqlite3.connect(PMID_DB_PATH)
    cur = conn.cursor()

    # Create table with a unique constraint on pmid and entity_id
    cur.execute(
        """CREATE TABLE IF NOT EXISTS species_data
                   (pmid TEXT, entity_id TEXT, names TEXT,
                    PRIMARY KEY (pmid, entity_id))"""
    )
    cur.execute(
        """CREATE TABLE IF NOT EXISTS disease_data
                   (pmid TEXT, entity_id TEXT, names TEXT,
                    PRIMARY KEY (pmid, entity_id))"""
    )

    conn.commit()
    conn.close()


def stream_and_store(filename, entity_type):
    conn = sqlite3.connect(PMID_DB_PATH)
    cur = conn.cursor()
    file_path = "/data/nde-hub/standardizers/pmid_lookup/"
    full_path = os.path.join(file_path, filename)

    with gzip.open(full_path, "rt") as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            pmid = row[0]
            entity_id = row[2]
            names = "|".join(row[3].split("|"))

            # Upsert query: update if exists, else insert
            cur.execute(
                f"""
                INSERT INTO {entity_type}_data (pmid, entity_id, names)
                VALUES (?, ?, ?)
                ON CONFLICT(pmid, entity_id) DO UPDATE SET names=excluded.names;
            """,
                (pmid, entity_id, names),
            )

    conn.commit()
    conn.close()


def get_data_for_pmids(pmids, entity_type):
    conn = sqlite3.connect(PMID_DB_PATH)
    cur = conn.cursor()
    placeholder = "?"  # SQLite placeholder
    placeholders = ", ".join(placeholder for unused in pmids)
    query = f"""SELECT entity_id, names FROM {entity_type}_data WHERE pmid IN ({placeholders})"""
    cur.execute(query, pmids)
    rows = cur.fetchall()

    data = {}
    for entity_id, names in rows:
        if entity_id not in data:
            data[entity_id] = []
        data[entity_id].extend(names.split("|"))

    conn.close()
    return data


def pubtator_lookup(name, table):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute(f"SELECT * FROM {table} WHERE original_name=?", (name.lower().strip(),))
    result = c.fetchone()
    conn.close()
    if result:
        logger.info(f"Skipping {name}, already in database")
        return json.loads(result[1])
    return None


def pubtator_add(name, table, standard_dict):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute(f"INSERT INTO {table} VALUES (?, ?)", (name.lower().strip(), standard_dict))
    conn.commit()
    conn.close()


@retry(7, 5)
def get_disease_details(identifier, original_name):
    """
    Retrieves disease details from nih for a given MeSH id.

    Parameters:
    - original_name (str): The original disease name.
    - identifier (str): The UniProt identifier for the disease.

    Returns:
    - standard_dict (dict): The standardized disease dictionary.
    """
    identifier = identifier.split(":")[-1]
    lookup_result = pubtator_lookup(original_name, "health_conditions")
    if lookup_result:
        logger.info(f"Skipping {original_name}, already in database")
        lookup_result["fromPMID"] = True
        lookup_result["isCurated"] = False
        if "curatedBy" in lookup_result:
            lookup_result.pop("curatedBy")
        else:
            logger.info(f"originalName not found in {lookup_result}")
        if "originalName" in lookup_result:
            lookup_result.pop("originalName")
        else:
            logger.info(f"originalName not found in {lookup_result}")
        return lookup_result
    logger.info(f"Converting {original_name} from MeSH {identifier} to standard format")
    non_mesh_result = query_condition(original_name, identifier)
    if non_mesh_result:
        logger.info(f"Converted {original_name} from MeSH {identifier} to standard format")
        pubtator_add(original_name, "health_conditions", json.dumps(non_mesh_result))
        non_mesh_result["fromPMID"] = True
        non_mesh_result["isCurated"] = False
        non_mesh_result.pop("curatedBy")
        if "originalName" in non_mesh_result:
            non_mesh_result.pop("originalName")
        else:
            logger.info(f"originalName not found in {non_mesh_result}")
        return non_mesh_result

    logger.info(f"Fetching details for {original_name} with ID {identifier}")
    # Fetch details from the MeSH API
    disease_info = requests.get(f"https://id.nlm.nih.gov/mesh/{identifier}.json")
    disease_info.raise_for_status()
    disease_info = disease_info.json()
    standard_dict = {
        "@type": "DefinedTerm",
        "identifier": identifier,
        "inDefinedTermSet": "MeSH",
        "url": f"https://id.nlm.nih.gov/mesh/{identifier}.html",
        "isCurated": False,
    }
    if terms := disease_info.get("terms"):
        alternative_names = []
        for term in terms:
            if term["preferred"] == True:
                standard_dict["name"] = term["label"]
            else:
                alternative_names.append(term["label"])
        if alternative_names:
            standard_dict["alternateName"] = alternative_names
    if label := disease_info.get("label"):
        standard_dict["name"] = label["@value"]
    if "name" not in standard_dict:
        raise Exception(f"No name found for {identifier}")

    pubtator_add(original_name, "health_conditions", json.dumps(standard_dict))
    standard_dict["fromPMID"] = True
    return standard_dict


def download_file(url, local_filename):
    file_path = "/data/nde-hub/standardizers/pmid_lookup/"
    full_path = os.path.join(file_path, local_filename)
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(full_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename


def file_needs_update(url, local_filename):
    file_path = "/data/nde-hub/standardizers/pmid_lookup/"
    full_path = os.path.join(file_path, local_filename)
    response = requests.head(url)
    if response.status_code == 200:
        remote_last_modified = response.headers.get("Last-Modified")
        if remote_last_modified:
            remote_last_modified = datetime.strptime(remote_last_modified, "%a, %d %b %Y %H:%M:%S GMT")
            if os.path.exists(full_path):
                local_last_modified = datetime.utcfromtimestamp(os.path.getmtime(full_path))
                return remote_last_modified > local_last_modified
            else:
                return True  # File does not exist locally, needs download
        else:
            print("Cannot determine remote file's last modified time.")
            return False
    else:
        print(f"Error {response.status_code} - Cannot check remote file's last modified time.")
        return False


def unzip_and_parse(filename, entity_type):
    data = {}
    with gzip.open(filename, "rt") as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            pmid = row[0]
            entity_id = row[2]
            names = row[3].split("|")
            # Initialize the PMID entry if not present
            if pmid not in data:
                data[pmid] = {"species": {}, "disease": {}}
            # Add the entity data to the PMID entry
            data[pmid][entity_type][entity_id] = names
    return data


def update_record_disease(rec, disease_data):
    if isinstance(rec.get("healthCondition"), dict):
        rec["healthCondition"] = [rec["healthCondition"]]

    existing_diseases = {disease["name"].lower(): i for i, disease in enumerate(rec.get("healthCondition", []))}

    abstract = rec.get("abstract", "")
    description = rec.get("description", "")
    title = rec.get("name", "")

    incorrect_terms = [
        "novel tumor",
        "hypervirulent covs",
        "covr-covs",
        "cancer stemness affords novel cancer",
        "2 tumors",
        "novel disease",
        "novel icos deficiency",
        "kucap-2 tumors",
        "hif-1/hif-2",
        "ncp",
    ]

    for mesh_id, diseases in disease_data.items():
        if "MESH" not in mesh_id:
            logger.warning(f"Invalid MeSH ID {mesh_id}")
            continue
        update_made = False
        for disease in diseases:
            name = disease.strip()
            if name.lower() in incorrect_terms and mesh_id == "MESH:C000657245":
                logger.warning(f"Incorrect Covid-19 mapping found for {name}")
                continue
            if mesh_id.strip() == "":
                logger.warning(f"Empty mesh_id found for {name}")
                continue
            if name.strip() == "":
                logger.warning(f"Empty name found for {mesh_id}")
                continue
            if mesh_id == "MESH:C000657245":
                mesh_id = "MESH:D000086382"

            if name.lower() in abstract.lower() or name.lower() in description.lower() or name.lower() in title.lower():
                logger.info(f"Found {name} in abstract, description, or title")
                if name.isupper():
                    logger.info(f"Possible Acronym: {name} in record: {rec['_id']}")
                logger.info(f"Adding {name} to record {rec['_id']}")
                try:
                    standardized_dict = get_disease_details(mesh_id, name)
                except Exception as e:
                    logger.warning(f"Could not get details for {name} with ID {mesh_id}: {e}")
                    logger.warning(f"URL: https://id.nlm.nih.gov/mesh/{mesh_id}.json")
                    continue

                if any(d.get("name", "").lower() == name.lower() for d in rec.get("healthCondition", [])):
                    remove_first_by_name(rec["healthCondition"], name)

                rec["healthCondition"] = rec.get("healthCondition", []) + [standardized_dict]
                logger.info(f"Added disease {name} to record {rec['_id']}")
                update_made = True
                break
            else:
                logger.info(f"{name} is not in abstract, description, or title. Not adding to record {rec['_id']}")
        if update_made:
            continue


def update_record_species(rec, species_data):
    """
    Updates the species in a given record based on abstract, description, and title.
    """
    if isinstance(rec.get("species"), dict):
        rec["species"] = [rec["species"]]
    if isinstance(rec.get("infectiousAgent"), dict):
        rec["infectiousAgent"] = [rec["infectiousAgent"]]

    existing_species = {spec["name"].lower(): i for i, spec in enumerate(rec.get("species", []))}
    existing_infectious_agents = {spec["name"].lower(): i for i, spec in enumerate(rec.get("infectiousAgent", []))}

    abstract = rec.get("abstract", "")
    description = rec.get("description", "")
    title = rec.get("name", "")

    blacklist = ["PERCH", "D-FISH"]

    for taxonomy_id, species_names in species_data.items():
        update_made = False
        for name in species_names:
            name = name.strip()
            if taxonomy_id.strip() == "":
                logger.warning(f"Empty mesh_id found for {taxonomy_id}")
                continue
            if name.strip() == "":
                logger.warning(f"Empty name found for {name}")
                continue
            if name.lower() in abstract.lower() or name.lower() in description.lower() or name.lower() in title.lower():
                logger.info(f"Found {name} in abstract, description, or title")
                if name.isupper():
                    logger.info(f"Possible Acronym: {name} in record: {rec['_id']}")
                if name in blacklist:
                    logger.info(f"Blacklisted: {name} in record: {rec['_id']}, skipping")
                    continue
                logger.info(f"Adding {name} to record {rec['_id']}")
                lookup_result = pubtator_lookup(name, "species")
                if lookup_result:
                    standardized_dict = lookup_result
                    standardized_dict["fromPMID"] = True
                    standardized_dict["isCurated"] = False
                    standardized_dict.pop("curatedBy")
                    standardized_dict.pop("originalName")
                else:
                    try:
                        standardized_dict = get_species_details(name, taxonomy_id)
                        pubtator_add(name, "species", json.dumps(standardized_dict))
                        standardized_dict.pop("originalName")
                        standardized_dict.pop("curatedBy")
                        standardized_dict["fromPMID"] = True
                        standardized_dict["isCurated"] = False
                    except Exception as e:
                        logger.warning(f"Could not get details for {name} with ID {taxonomy_id}: {e}")
                        logger.warning(f"URL: https://rest.uniprot.org/taxonomy/{taxonomy_id}")
                        continue

                if any(spec.get("name", "").lower() == name.lower() for spec in rec.get("species", [])):
                    remove_first_by_name(rec["species"], name)
                elif any(spec.get("name", "").lower() == name.lower() for spec in rec.get("infectiousAgent", [])):
                    remove_first_by_name(rec["infectiousAgent"], name)

                if "classification" not in standardized_dict:
                    logger.warning(f"Could not classify {name} with ID {taxonomy_id}")
                    rec["species"] = rec.get("species", []) + [standardized_dict]
                    logger.info(f"Added species {name} to record {rec['_id']}")
                elif standardized_dict["classification"] == "host":
                    rec["species"] = rec.get("species", []) + [standardized_dict]
                    logger.info(f"Added species {name} to record {rec['_id']}")
                elif standardized_dict["classification"] == "infectiousAgent":
                    rec["infectiousAgent"] = rec.get("infectiousAgent", []) + [standardized_dict]
                    logger.info(f"Added infectious agent {name} to record {rec['_id']}")
                update_made = True
                break
            else:
                logger.info(f"{name} is not in abstract, description, or title. Not adding to record {rec['_id']}")
        if update_made:
            continue


@retry(7, 5)
def _convert_pmc(pmc_list, pmc_dict):
    base_url = (
        "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=my_email@example.com&format=json&"
    )
    # convert all unique pmcs into a comma separated string
    request_url = "ids=" + ",".join(list(set(pmc_list)))
    # make the request
    request_url = base_url + request_url
    request = requests.get(request_url).json()
    # key pmcid value pmid
    for record in request.get("records"):
        pmc_dict[record.get("pmcid")] = record.get("pmid")
    # reset the list
    pmc_list.clear()
    time.sleep(0.5)


@retry(7, 5)
def _convert_doi(doi, doi_dict):
    api_key = GEO_API_KEY
    email = GEO_EMAIL
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    if doi not in doi_dict:
        term = f"{doi}[DOI]"
        handle = Entrez.esearch(db="pubmed", term=term, retmode="json")
        data = json.loads(handle.read())
        handle.close()
        pmids = data.get("esearchresult", {}).get("idlist", [])
        doi_dict[doi] = pmids[0] if pmids else None
        if api_key:
            time.sleep(0.05)
        else:
            time.sleep(0.35)


def _get_pub_date(date: str):
    """helper method to solve the problem transforming dates such as "2000 Spring" into date().isoformat dates

    https://www.nlm.nih.gov/bsd/mms/medlineelements.html#dp
    Returns:
        An isoformatted date depending on context:

        Seasons use metrological start
        Winter: December 1
        Spring: March 1
        Summer: June 1
        Fall: September 1

        Dates with Y/M/D-D take only the beginning day
        Dates with only Y/M or Y/M-M take the first day of that month
        Dates with only Y or Y-Y take first day of that year
        TODO: Not important since only one instance so far but fix edge case "2016 11-12"
    """

    months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    seasons = {"spring": " mar 1", "summer": " jun 1", "fall": " sep 1", "winter": " dec 1"}

    s_date = date.lower().split()
    date_len = len(s_date)
    # if length is 1 can either be year or year-year
    if date_len == 1:
        return datetime.strptime(s_date[0].split("-")[0], "%Y").date().isoformat()
    # if length is 2 can either be year season or year month or year month-month
    elif date_len == 2:
        if s_date[1][:3] in months:
            return datetime.strptime(s_date[0] + " " + s_date[1][:3], "%Y %b").date().isoformat()
        elif season := seasons.get(s_date[1]):
            return datetime.strptime(s_date[0] + season, "%Y %b %d").date().isoformat()
        else:
            logger.warning("Need to update isoformat transformation: %s", date)
            return None
    # if length is 3 should be year month day or year month day-day or year month-month day
    elif date_len == 3:
        year = s_date[0]
        og_month = s_date[1].split("-")[0]
        day = s_date[2].split("-")[0]
        month = next((month for month in months if month in og_month), None)
        if month:
            return datetime.strptime(year + " " + month + " " + day, "%Y %b %d").date().isoformat()
        else:
            logger.warning("Need to update isoformat transformation: %s", date)
    # exception case there are quite a few entries with this case "2020 Jan - Feb"
    elif date_len == 4:
        if s_date[1] in months and s_date[3] in months and s_date[2] == "-":
            return datetime.strptime(s_date[0] + " " + s_date[1], "%Y %b").date().isoformat()
        else:
            logger.warning("Need to update isoformat transformation %s", date)
    else:
        logger.warning("Need to update isoformat transformation: %s", date)
        return None


@retry(7, 30)
def batch_get_pmid_eutils(pmids: Iterable[str], email: str, api_key: Optional[str] = None) -> Dict:
    """Use pmid to retrieve both citation and funding info in batch
    :param pmids: A list of PubMed PMIDs
    :param api_key: API Key from NCBI to access E-utilities
    :return: A dictionary containing the pmids which hold citations and funding.
    """
    # probably dont need this line. Using python package, should work both ways.
    # if pmids is str:
    #     warnings.warn(f"Got str:{pmids} as parameter, expecting an Iterable of str", RuntimeWarning)

    # set up Entrez variables. Email is required.
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    ct_fd = {}

    # api query to parse citations
    try:
        handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
    except urllib.error.HTTPError as err:
        logger.error("This is the length of the pmids %s", len(pmids))
        logger.error("The list of pmids %s", pmids)
        logger.error("HTTP url: %s", err.url)
        raise err

    records = Medline.parse(handle)
    for record in records:
        citation = {}
        # rename keys
        if name := record.get("TI"):
            citation["name"] = name
        if pmid := record.get("PMID"):
            citation["pmid"] = pmid
            citation["identifier"] = "PMID:" + pmid
            citation["url"] = "https://pubmed.ncbi.nlm.nih.gov/" + pmid + "/"
        if journal_name := record.get("JT"):
            citation["journalName"] = journal_name
        if date_published := record.get("DP"):
            if date := _get_pub_date(date_published):
                citation["datePublished"] = date

        # make an empty list if there is some kind of author
        if record.get("AU") or record.get("CN"):
            citation["author"] = []
        if authors := record.get("AU"):
            for author in authors:
                citation["author"].append({"@type": "Person", "name": author})
        if corp_authors := record.get("CN"):
            for corp_author in corp_authors:
                citation["author"].append({"@type": "Organization", "name": corp_author})
        if citation:
            citation["fromPMID"] = True
        # put citation in dictionary
        ct_fd[pmid] = {"citation": citation}

    # throttle request rates, NCBI says up to 10 requests per second with API Key, 3/s without.
    if api_key:
        time.sleep(0.1)
    else:
        time.sleep(0.35)

    # get the funding using xml file because of problems parsing the medline file
    # https://www.nlm.nih.gov/bsd/mms/medlineelements.html#gr
    handle = Entrez.efetch(db="pubmed", id=pmids, retmode="xml")

    # Have to use Entrez.read() instead of Entrez.parse(). The problem is discussed here: https://github.com/biopython/biopython/issues/1027
    # This can get an incompleteread error we rerun this at the top layer
    records = Entrez.read(handle)
    records = records["PubmedArticle"]

    funding = []
    for record in records:
        if grants := record["MedlineCitation"]["Article"].get("GrantList"):
            for grant in grants:
                fund = {}
                if grant_id := grant.get("GrantID"):
                    grant_id = str(grant_id)
                    fund["identifier"] = grant_id
                if agency := grant.get("Agency"):
                    agency = str(agency)
                    fund_dict = standardize_funder(agency)
                    if fund_dict:
                        fund["funder"] = fund_dict
                    else:
                        fund["funder"] = {"@type": "Organization", "name": agency}
                if agency or grant_id:
                    fund["fromPMID"] = True
                funding.append(fund)
        if pmid := record["MedlineCitation"].get("PMID"):
            if funding:
                ct_fd[pmid]["funding"] = funding
            funding = []

    return ct_fd


def load_pmid_ctfd(data):
    """Takes 1000 documents at a time and batch queries all of the pmids in the documents to improve runtime.
    If there are any pmcids, convert all of them into pmids before running the batch query.
    Loads the citation and funding into the documents.
      Returns: A generator with the completed documents
    """

    # a temporary solution to make bigger batch api call instead of multiple smaller calls in crawler to improve runtime
    # TODO: figure out how to make a batch api call in crawler perferrably

    api_key = GEO_API_KEY
    email = GEO_EMAIL

    # if no access to config file comment out above and enter your own email
    # email = myemail@gmail.com

    # URLs for the files
    species_url = "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/species2pubtatorcentral.gz"
    disease_url = "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/disease2pubtatorcentral.gz"

    # Filenames for the local copies of the files
    species_filename = "species2pubtatorcentral.gz"
    disease_filename = "disease2pubtatorcentral.gz"

    if not os.path.exists(PMID_DB_PATH):
        create_db()

    # Check if files need to be updated and download if necessary
    if file_needs_update(species_url, species_filename):
        logger.info(f"Downloading species data from {species_url}...")
        download_file(species_url, species_filename)
        logger.info("Unzipping and storing species data...")
        stream_and_store(species_filename, "species")

    if file_needs_update(disease_url, disease_filename):
        logger.info(f"Downloading disease data from {disease_url}...")
        download_file(disease_url, disease_filename)
        logger.info("Unzipping and storing disease data...")
        stream_and_store(disease_filename, "disease")

    # process the file
    def process_data(d, is_file: bool = False):
        count = 0
        while True:
            # dict to convert pmcs to pmids
            pmc_pmid = {}
            # pmc list for batch query
            pmc_list = []
            # pmid list for batch query
            pmid_list = []
            # docs to yield for each batch query
            doc_list = []
            # dict to convert dois to pmids
            doi_pmid = {}

            # to make batch api query take the next 1000 docs and collect all the pmids
            next_n_lines = list(islice(d, 1000))
            if not next_n_lines:
                break
            for line in next_n_lines:
                count += 1
                if count % 1000 == 0:
                    logger.info("Processed %s documents", count)
                if is_file:
                    doc = orjson.loads(line)
                else:
                    doc = line
                doc_list.append(doc)
                if pmcs := doc.get("pmcs"):
                    pmcs = [pmc.strip() for pmc in pmcs.split(",")]
                    # limit reached need to make request to pmc convertor
                    if (len(pmc_list) + len(pmcs)) >= 200:
                        _convert_pmc(pmc_list, pmc_pmid)
                    pmc_list += pmcs

                if citations := doc.get("citation"):
                    if not isinstance(citations, list):
                        citations = [citations]
                    for citation in citations:
                        if doi := citation.get("doi"):
                            _convert_doi(doi, doi_pmid)

            # make request to pmc convertor
            if pmc_list:
                _convert_pmc(pmc_list, pmc_pmid)

            # loop through doc_list and convert
            for loc, doc in enumerate(doc_list):
                if pmcs := doc.pop("pmcs", None):
                    pmcs = [pmc.strip() for pmc in pmcs.split(",")]
                    for pmc in pmcs:
                        if pmid := pmc_pmid.get(pmc):
                            doc["pmids"] = doc.get("pmids") + "," + pmid if doc.get("pmids") else pmid
                        else:
                            logger.info("There is an issue with this PMCID. PMCID: %s, rec_id: %s", pmc, doc["_id"])
                if citations := doc.get("citation", None):
                    if not isinstance(citations, list):
                        citations = [citations]
                    for citation in citations:
                        if doi := citation.get("doi"):
                            if pmid := doi_pmid.get(doi):
                                doc["pmids"] = doc.get("pmids") + "," + pmid if doc.get("pmids") else pmid

                    doc_list[loc] = doc
                if pmids := doc.get("pmids"):
                    pmid_list += [pmid.strip() for pmid in pmids.split(",")]

            # check if there are any pmids before requesting
            if pmid_list:
                # same thing as list(set(pmid_list))
                pmid_list = [*set(pmid_list)]
                # batch request retry up to 3 times
                eutils_info = batch_get_pmid_eutils(pmid_list, email, api_key)
                # throttle request rates, NCBI says up to 10 requests per second with API Key, 3/s without.
                if api_key:
                    time.sleep(0.1)
                else:
                    time.sleep(0.35)

            # add in the citation and funding to each doc in doc_list and yield
            for rec in doc_list:
                if pmids := rec.pop("pmids", None):
                    pmids = [pmid.strip().lstrip("0") for pmid in pmids.split(",")]
                    pmids = list(set(pmids))
                    species = get_data_for_pmids(pmids, "species")
                    if species:
                        update_record_species(rec, species)
                    diseases = get_data_for_pmids(pmids, "disease")
                    if diseases:
                        update_record_disease(rec, diseases)
                    # Sequential processing for the rest
                    for pmid in pmids:
                        if not eutils_info.get(pmid):
                            logger.info("There is an issue with this pmid. PMID: %s, rec_id: %s", pmid, rec["_id"])
                        # this fixes the error where there is no pmid
                        # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41964
                        if eutils_info.get(pmid):
                            if citation := eutils_info[pmid].get("citation"):
                                if rec_citation := rec.get("citation"):
                                    # if the user originally had a citation field that is not a list change citation to list
                                    if not isinstance(rec_citation, list):
                                        rec["citation"] = [rec_citation]
                                    rec["citation"].append(citation)
                                else:
                                    rec["citation"] = [citation]
                            if funding := eutils_info[pmid].get("funding"):
                                if rec_funding := rec.get("funding"):
                                    # if the user originally had a funding field that is not a list change funding to list
                                    if not isinstance(rec_funding, list):
                                        rec["funding"] = [rec_funding]
                                    rec["funding"] += funding
                                else:
                                    rec["funding"] = copy(funding)
                yield rec

    if isinstance(data, str):
        with open(os.path.join(data, "data.ndjson"), "rb") as f:
            yield from process_data(f, True)
    else:
        yield from process_data(data, False)  # yields each document individually


def load_pmid_ctfd_wrapper(func):
    """Wrapper function that takes in a generator and yields a generator that adds citations and funding to the documents using pmids.
    Takes 1000 documents at a time and batch queries all of the pmids in the documents to improve runtime.
    If there are any pmcids, convert all of them into pmids before running the batch query.
    Loads the citation and funding into the documents.
      Returns: A generator with the completed documents
    """

    # a temporary solution to make bigger batch api call instead of multiple smaller calls in crawler to improve runtime
    # TODO: figure out how to make a batch api call in crawler perferrably

    api_key = GEO_API_KEY
    email = GEO_EMAIL

    # if no access to config file comment out above and enter your own email
    # email = myemail@gmail.com
    # URLs for the files
    species_url = "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/species2pubtatorcentral.gz"
    disease_url = "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/disease2pubtatorcentral.gz"

    # Filenames for the local copies of the files
    species_filename = "species2pubtatorcentral.gz"
    disease_filename = "disease2pubtatorcentral.gz"

    if not os.path.exists(PMID_DB_PATH):
        create_db()

    # Check if files need to be updated and download if necessary
    if file_needs_update(species_url, species_filename):
        logger.info(f"Downloading species data from {species_url}...")
        download_file(species_url, species_filename)
        logger.info("Unzipping and storing species data...")
        stream_and_store(species_filename, "species")

    if file_needs_update(disease_url, disease_filename):
        logger.info(f"Downloading disease data from {disease_url}...")
        download_file(disease_url, disease_filename)
        logger.info("Unzipping and storing disease data...")
        stream_and_store(disease_filename, "disease")

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        count = 0
        while True:
            # dict to convert pmcs to pmids
            pmc_pmid = {}
            # pmc list for batch query
            pmc_list = []
            # pmid list for batch query
            pmid_list = []
            # docs to yield for each batch query
            doc_list = []
            # dict to convert dois to pmids
            doi_pmid = {}

            # to make batch api query take the next 1000 docs and collect all the pmids
            next_n_lines = list(islice(data, 1000))
            if not next_n_lines:
                break
            for line in next_n_lines:
                count += 1
                if count % 1000 == 0:
                    logger.info("Processed %s documents", count)
                doc = orjson.loads(line)
                doc_list.append(doc)
                if pmcs := doc.get("pmcs"):
                    pmcs = [pmc.strip() for pmc in pmcs.split(",")]
                    # limit reached need to make request to pmc convertor
                    if (len(pmc_list) + len(pmcs)) >= 200:
                        _convert_pmc(pmc_list, pmc_pmid)
                    pmc_list += pmcs
                if citations := doc.get("citation"):
                    if not isinstance(citations, list):
                        citations = [citations]
                    for citation in citations:
                        if doi := citation.get("doi"):
                            _convert_doi(doi, doi_pmid)

            # make request to pmc convertor
            if pmc_list:
                _convert_pmc(pmc_list, pmc_pmid)

            # loop through doc_list and convert
            for loc, doc in enumerate(doc_list):
                if pmcs := doc.pop("pmcs", None):
                    pmcs = [pmc.strip() for pmc in pmcs.split(",")]
                    for pmc in pmcs:
                        if pmid := pmc_pmid.get(pmc):
                            doc["pmids"] = doc.get("pmids") + "," + pmid if doc.get("pmids") else pmid
                        else:
                            logger.info("There is an issue with this PMCID. PMCID: %s, rec_id: %s", pmc, doc["_id"])
                if citations := doc.get("citation", None):
                    if not isinstance(citations, list):
                        citations = [citations]
                    for citation in citations:
                        if doi := citation.get("doi"):
                            if pmid := doi_pmid.get(doi):
                                doc["pmids"] = doc.get("pmids") + "," + pmid if doc.get("pmids") else pmid

                    doc_list[loc] = doc
                if pmids := doc.get("pmids"):
                    pmid_list += [pmid.strip() for pmid in pmids.split(",")]

            # check if there are any pmids before requesting
            if pmid_list:
                # same thing as list(set(pmid_list))
                pmid_list = [*set(pmid_list)]
                # batch request retry up to 3 times
                eutils_info = batch_get_pmid_eutils(pmid_list, email, api_key)
                # throttle request rates, NCBI says up to 10 requests per second with API Key, 3/s without.
                if api_key:
                    time.sleep(0.1)
                else:
                    time.sleep(0.35)

            # add in the citation and funding to each doc in doc_list and yield
            for rec in doc_list:
                if pmids := rec.pop("pmids", None):
                    pmids = [pmid.strip() for pmid in pmids.split(",")]
                    # fixes issue where pmid numbers under 10 is read as 04 instead of 4
                    pmids = [pmid.lstrip("0") for pmid in pmids]
                    pmids = [*set(pmids)]
                    species = get_data_for_pmids(pmids, "species")
                    if species:
                        update_record_species(rec, species)
                    diseases = get_data_for_pmids(pmids, "disease")
                    if diseases:
                        update_record_disease(rec, diseases)
                    if not eutils_info.get(pmid):
                        logger.info("There is an issue with this pmid. PMID: %s, rec_id: %s", pmid, rec["_id"])
                    # this fixes the error where there is no pmid
                    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41964
                    if eutils_info.get(pmid):
                        if citation := eutils_info[pmid].get("citation"):
                            if rec_citation := rec.get("citation"):
                                # if the user originally had a citation field that is not a list change citation to list
                                if not isinstance(rec_citation, list):
                                    rec["citation"] = [rec_citation]
                                rec["citation"].append(citation)
                            else:
                                rec["citation"] = [citation]
                        if funding := eutils_info[pmid].get("funding"):
                            if rec_funding := rec.get("funding"):
                                # if the user originally had a funding field that is not a list change funding to list
                                if not isinstance(rec_funding, list):
                                    rec["funding"] = [rec_funding]
                                rec["funding"] += funding
                            else:
                                rec["funding"] = copy(funding)
            yield rec

    return wrapper


def standardize_fields(docs):
    api_key = GEO_API_KEY
    email = GEO_EMAIL

    fields = ["isBasedOn", "isBasisFor", "citedBy", "isPartOf", "hasPart"]
    while True:
        # pmid list for batch query
        pmid_list = []
        # docs to yield for each batch query
        doc_list = []

        next_n_docs = list(islice(docs, 1000))
        if not next_n_docs:
            break
        for count, doc in enumerate(next_n_docs, start=1):
            if count % 1000 == 0:
                logger.info("Processed %s documents", count)

            doc_list.append(doc)

            for field in fields:
                if objects := doc.get(field):
                    if not isinstance(objects, list):
                        objects = [objects]
                    for object in objects:
                        if object.get("@type") == "ScholarlyArticle":
                            if pmid := object.get("pmid"):
                                doc["pmids"] = doc.get("pmids") + "," + str(pmid) if doc.get("pmids") else str(pmid)

            if pmids := doc.pop("pmids", None):
                pmid_list += [pmid.strip() for pmid in pmids.split(",")]

        # check if there are any pmids before requesting
        if pmid_list:
            # same thing as list(set(pmid_list))
            pmid_list = [*set(pmid_list)]
            # batch request retry up to 3 times
            eutils_info = batch_get_pmid_eutils(pmid_list, email, api_key)
            # throttle request rates, NCBI says up to 10 requests per second with API Key, 3/s without.
            if api_key:
                time.sleep(0.1)
            else:
                time.sleep(0.35)

        for doc in doc_list:
            for field in fields:
                if objects := doc.get(field):
                    if not isinstance(objects, list):
                        objects = [objects]
                    field_pmids = [
                        str(object.get("pmid"))
                        for object in objects
                        if object.get("pmid") and object.get("@type") == "ScholarlyArticle"
                    ]
                    doc[field] = [
                        object
                        for object in objects
                        if not (object.get("pmid") and object.get("@type") == "ScholarlyArticle")
                    ]
                    if field_pmids:
                        for pmid in field_pmids:
                            if eutils_info.get(pmid):
                                if citation := eutils_info[pmid].get("citation"):
                                    citation["type"] = "ScholarlyArticle"
                                    doc[field].append(citation)
            yield doc
