# Helper file to batch call pmids to get citations and funding
# Helpful links to documentation of Biopython Package for writing this file
# https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec162
# https://biopython.org/docs/1.76/api/Bio.Entrez.html
# https://www.nlm.nih.gov/bsd/mms/medlineelements.html
# https://dataguide.nlm.nih.gov/eutilities/utilities.html#efetch
# https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/

import functools
import os
import os.path
import time
from copy import copy
from datetime import datetime
from ftplib import FTP
from itertools import islice
from typing import Dict, Iterable, Optional

import orjson
import requests
from Bio import Entrez, Medline
from config import GEO_API_KEY, GEO_EMAIL, logger

from .funding_helper import standardize_funder
from .pubtator import classify_as_host_or_agent
from .utils import retry

SPECIES_CACHE = {}
DISEASE_CACHE = {}


@retry(3, 5)
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
    if identifier in DISEASE_CACHE:
        logger.info(f"Fetching details from cache for {original_name}")
        return DISEASE_CACHE[identifier]
    logger.info(f"Getting details for {original_name}")

    # Fetch details from the MeSH API
    disease_info = requests.get(f"https://id.nlm.nih.gov/mesh/{identifier}.json")
    disease_info.raise_for_status()
    disease_info = disease_info.json()
    standard_dict = {
        "identifier": identifier,
        "inDefinedTermSet": "MeSH",
        "url": f"https://id.nlm.nih.gov/mesh/{identifier}.html",
        "fromPMID": True,
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": datetime.now().strftime("%Y-%m-%d"),
        },
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

    DISEASE_CACHE[identifier] = standard_dict
    return standard_dict


@retry(3, 5)
def get_species_from_api(pmids):
    url = f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids={','.join(pmids)}&concepts=species"
    response = requests.get(url)
    response.raise_for_status()
    lines = response.text.split("\n")
    taxonomy_info = []
    for line in lines:
        columns = line.split("\t")
        if len(columns) > 4 and columns[4] == "Disease":
            taxonomy_info.append({"identifier": columns[5], "name": columns[3].lower()})
    return taxonomy_info


@retry(3, 5)
def get_species_details(identifier, original_name):
    """
    Retrieves species details from UniProt for a given species name.

    Parameters:
    - original_name (str): The original species name.
    - identifier (str): The UniProt identifier for the species.

    Returns:
    - standard_dict (dict): The standardized species dictionary.
    """
    identifier = identifier.split("*")[-1]
    if identifier in SPECIES_CACHE:
        logger.info(f"Fetching details from cache for {original_name}")
        return SPECIES_CACHE[identifier]

    logger.info(f"Getting details for {original_name}")
    # Fetch details from the UniProt API

    species_info = requests.get(f"https://rest.uniprot.org/taxonomy/{identifier}")
    species_info.raise_for_status()
    species_info = species_info.json()
    standard_dict = {
        "identifier": identifier,
        "inDefinedTermSet": "UniProt",
        "url": f"https://www.uniprot.org/taxonomy/{identifier}",
        "fromPMID": True,
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": datetime.now().strftime("%Y-%m-%d"),
        },
    }
    if scientific_name := species_info.get("scientificName"):
        standard_dict["name"] = scientific_name
    else:
        standard_dict["name"] = original_name

    alternative_names = []
    if common_name := species_info.get("commonName"):
        standard_dict["commonName"] = common_name
        alternative_names.append(common_name)

        standard_dict["displayName"] = f"{common_name} | {scientific_name}"

    if other_names := species_info.get("otherNames"):
        alternative_names.extend(other_names)

    if alternative_names:
        standard_dict["alternateName"] = alternative_names

    if lineage := species_info.get("lineage"):
        standard_dict["classification"] = classify_as_host_or_agent(lineage)
    else:
        logger.warning(f"No lineage found for {original_name} {identifier}")
        standard_dict["classification"] = "host"

    SPECIES_CACHE[identifier] = standard_dict
    return standard_dict


@retry(3, 5)
def get_disease_from_api(pmids):
    url = f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids={','.join(pmids)}&concepts=disease"
    response = requests.get(url)
    response.raise_for_status()
    lines = response.text.split("\n")
    mesh_info = []
    for line in lines:
        columns = line.split("\t")
        if len(columns) > 4 and columns[4] == "Disease":
            if columns[5] == "":
                logger.warning(f"No MeSH ID found for {columns[3]}")
                logger.info(f"url: {url}")
                continue
            mesh_info.append({"identifier": columns[5], "name": columns[3].lower()})
    return mesh_info


def update_record_disease(rec, disease_data):
    if isinstance(rec.get("healthCondition"), dict):
        rec["healthCondition"] = [rec["healthCondition"]]

    existing_diseases = {spec["name"].lower(): i for i, spec in enumerate(rec.get("healthCondition", []))}

    # Fetch the abstract, description, and title for current record
    abstract = rec.get("abstract", "")
    description = rec.get("description", "")
    title = rec.get("name", "")

    added_mesh_ids = set()

    # Iterate over new disease data and update the record
    for disease_dict in disease_data:
        # for identifier, name in disease_data.items():
        identifier = disease_dict["identifier"]
        name = disease_dict["name"]
        # Skip if mesh ID has already been added
        if identifier in added_mesh_ids:
            logger.info(f"Skipping {identifier} because it has already been added")
            continue

        if name.lower() in abstract.lower() or name.lower() in description.lower() or name.lower() in title.lower():
            logger.info(f"Found {name} in abstract, description, or title")
            if name.isupper():
                logger.info(f"Possible Acronym: {name} in record: {rec['_id']}")
            logger.info(f"Adding {name} to record {rec['_id']}")
            try:
                standardized_dict = get_disease_details(identifier, name)
            except Exception as e:
                logger.warning(f"Could not get details for {name} with ID {identifier}: {e}")
                logger.warning(f"URL: https://id.nlm.nih.gov/mesh/{identifier}.json")
                continue

            if name.lower() in existing_diseases:
                del rec["healthCondition"][existing_diseases[name.lower()]]

            rec["healthCondition"] = rec.get("healthCondition", []) + [standardized_dict]
            logger.info(f"Added disease {name} to record {rec['_id']}")

            added_mesh_ids.add(identifier)

        else:
            logger.info(f"{name} is not in abstract, description, or title. Not adding to record {rec['_id']}")


def update_record_species(rec, species_data):
    """
    Updates the species in a given record based on abstract, description, and title.

    Parameters:
    - rec (dict): The record to be updated.
    - species_data (dict): Dictionary containing species data, with species id as key and lowercase names separated by a "|" as the value.
    """

    # Convert species and infectiousAgent to lists if they are not already
    if isinstance(rec.get("species"), dict):
        rec["species"] = [rec["species"]]
    if isinstance(rec.get("infectiousAgent"), dict):
        rec["infectiousAgent"] = [rec["infectiousAgent"]]

    existing_species = {spec["name"].lower(): i for i, spec in enumerate(rec.get("species", []))}
    existing_infectious_agents = {spec["name"].lower(): i for i, spec in enumerate(rec.get("infectiousAgent", []))}

    # Fetch the abstract, description, and title for current record
    abstract = rec.get("abstract", "")
    description = rec.get("description", "")
    title = rec.get("name", "")

    added_taxonomy_ids = set()

    blacklist = ["PERCH", "D-FISH"]

    # Iterate over new species data and update the record
    for species_dict in species_data:
        identifier = species_dict["identifier"]
        name = species_dict["name"]
        # Skip if taxonomy ID has already been added
        if identifier in added_taxonomy_ids:
            logger.info(f"Skipping {identifier} because it has already been added")
            continue

        if name.lower() in abstract.lower() or name.lower() in description.lower() or name.lower() in title.lower():
            logger.info(f"Found {name} in abstract, description, or title")
            if name.isupper():
                logger.info(f"Possible Acronym: {name} in record: {rec['_id']}")
            if name in blacklist:
                logger.info(f"Blacklisted: {name} in record: {rec['_id']}, skipping")
                continue
            logger.info(f"Adding {name} to record {rec['_id']}")

            try:
                standardized_dict = get_species_details(identifier, name)
            except Exception as e:
                logger.warning(f"Could not get details for {name} with ID {identifier}: {e}")
                logger.warning(f"URL: https://rest.uniprot.org/taxonomy/{identifier}")
                continue

            if name.lower() in existing_species:
                del rec["species"][existing_species[name.lower()]]
            elif name.lower() in existing_infectious_agents:
                del rec["infectiousAgent"][existing_infectious_agents[name.lower()]]

            if "classification" not in standardized_dict:
                logger.warning(f"Could not classify {name} with ID {identifier}")
                rec["species"] = rec.get("species", []) + [standardized_dict]
                logger.info(f"Added species {name} to record {rec['_id']}")

            elif standardized_dict["classification"] == "host":
                rec["species"] = rec.get("species", []) + [standardized_dict]
                logger.info(f"Added species {name} to record {rec['_id']}")

            elif standardized_dict["classification"] == "infectiousAgent":
                rec["infectiousAgent"] = rec.get("infectiousAgent", []) + [standardized_dict]
                logger.info(f"Added infectious agent {name} to record {rec['_id']}")

            added_taxonomy_ids.add(identifier)

        else:
            logger.info(f"{name} is not in abstract, description, or title. Not adding to record {rec['_id']}")


# retry 3 times sleep 5 seconds between each retry
@retry(3, 5)
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
        return (
            datetime.strptime(s_date[0] + " " + s_date[1].split("-")[0] + " " + s_date[2].split("-")[0], "%Y %b %d")
            .date()
            .isoformat()
        )
    # exception case there are quite a few entries with this case "2020 Jan - Feb"
    elif date_len == 4:
        if s_date[1] in months and s_date[3] in months and s_date[2] == "-":
            return datetime.strptime(s_date[0] + " " + s_date[1], "%Y %b").date().isoformat()
        else:
            logger.warning("Need to update isoformat transformation %s", date)
    else:
        logger.warning("Need to update isoformat transformation: %s", date)
        return None


# retry 3 times sleep 5 seconds between each retry
@retry(3, 5)
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
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")

    records = Medline.parse(handle)
    # This can get an incompleteread error we rerun this at the top layer
    records = list(records)
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
                    fund["identifier"] = grant_id
                if agency := grant.get("Agency"):
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


def load_pmid_ctfd(data_folder):
    """Takes 1000 documents at a time and batch queries all of the pmids in the documents to improve runtime.
    If there are any pmcids, convert all of them into pmids before running the batch query.
    Loads the citation and funding into the documents along with uploading the date field.
      Returns: A generator with the completed documents
    """

    # a temporary solution to make bigger batch api call instead of multiple smaller calls in crawler to improve runtime
    # TODO: figure out how to make a batch api call in crawler perferrably

    api_key = GEO_API_KEY
    email = GEO_EMAIL

    # if no access to config file comment out above and enter your own email
    # email = myemail@gmail.com

    with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
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

            # to make batch api query take the next 1000 docs and collect all the pmids
            next_n_lines = list(islice(f, 1000))
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

                    # get species tax ids from pubtator
                    try:
                        species_data = get_species_from_api(pmids)
                    except Exception as e:
                        logger.warning(f"Could not get species data for {pmids}: {e}")
                        species_data = []
                    # get disease mesh ids from pubtator
                    try:
                        disease_data = get_disease_from_api(pmids)
                    except Exception as e:
                        logger.warning(f"Could not get disease data for {pmids}: {e}")
                        disease_data = []
                    if species_data:
                        # Update the species field in the record
                        update_record_species(rec, species_data)
                    if disease_data:
                        # Update the healthCondition field in the record
                        update_record_disease(rec, disease_data)

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


def load_pmid_ctfd_wrapper(func):
    """Wrapper function that takes in a generator and yields a generator that adds citations and funding to the documents using pmids.
    Takes 1000 documents at a time and batch queries all of the pmids in the documents to improve runtime.
    If there are any pmcids, convert all of them into pmids before running the batch query.
    Loads the citation and funding into the documents along with uploading the date field.
      Returns: A generator with the completed documents
    """

    # a temporary solution to make bigger batch api call instead of multiple smaller calls in crawler to improve runtime
    # TODO: figure out how to make a batch api call in crawler perferrably

    api_key = GEO_API_KEY
    email = GEO_EMAIL

    # if no access to config file comment out above and enter your own email
    # email = myemail@gmail.com

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

                    # get species tax ids from pubtator
                    try:
                        species_data = get_species_from_api(pmids)
                    except Exception as e:
                        logger.warning(f"Could not get species data for {pmids}: {e}")
                        species_data = []
                    # get disease mesh ids from pubtator
                    try:
                        disease_data = get_disease_from_api(pmids)
                    except Exception as e:
                        logger.warning(f"Could not get disease data for {pmids}: {e}")
                        disease_data = []
                    if species_data:
                        # Update the species field in the record
                        update_record_species(rec, species_data)
                    if disease_data:
                        # Update the healthCondition field in the record
                        update_record_disease(rec, disease_data)

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

    return wrapper
