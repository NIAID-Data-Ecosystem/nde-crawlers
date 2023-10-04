# Helper file to batch call pmids to get citations and funding
# Helpful links to documentation of Biopython Package for writing this file
# https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec162
# https://biopython.org/docs/1.76/api/Bio.Entrez.html
# https://www.nlm.nih.gov/bsd/mms/medlineelements.html
# https://dataguide.nlm.nih.gov/eutilities/utilities.html#efetch
# https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/
import os
import time
from copy import copy
from datetime import datetime
from itertools import islice
from typing import Dict, Iterable, Optional

import orjson
import requests
from Bio import Entrez, Medline
from config import GEO_API_KEY, GEO_EMAIL, logger

from .utils import retry

SPECIES_CACHE = {}


def classify_as_host_or_agent(lineage):
    # Extracting scientific names for easy processing
    scientific_names = [item["scientificName"] for item in lineage]

    # Check for host species conditions
    if "Deuterostomia" in scientific_names:
        new_classification = "host"
    elif "Embryophyta" in scientific_names and not any(
        parasite in scientific_names for parasite in ["Arceuthobium", "Cuscuta", "Orobanche", "Striga", "Phoradendron"]
    ):
        new_classification = "host"
    elif "Arthropoda" in scientific_names:
        if "Acari" in scientific_names:
            if "Ixodida" in scientific_names:
                new_classification = "host"
            else:
                new_classification = "infectiousAgent"
        else:
            new_classification = "host"
    else:
        # If not falling under the above host conditions, classify as infectiousAgent
        new_classification = "infectiousAgent"

    return new_classification


def get_species_details(original_name, identifier):
    if identifier in SPECIES_CACHE:
        logger.info(f"Fetching details from cache for {original_name}")
        return SPECIES_CACHE[identifier]

    logger.info(f"Getting details for {original_name}")
    # Fetch details from the UniProt API
    identifier = identifier.split("*")[-1]

    species_info = requests.get(f"https://rest.uniprot.org/taxonomy/{identifier}")
    species_info.raise_for_status()
    species_info = species_info.json()
    standard_dict = {
        "identifier": identifier,
        "inDefinedTermSet": "UniProt",
        "url": f"https://www.uniprot.org/taxonomy/{identifier}",
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
        standard_dict["classification"] = classify_as_host_or_agent(lineage, identifier)
    else:
        logger.warning(f"No lineage found for {original_name} {identifier}")
        standard_dict["classification"] = "host"

    SPECIES_CACHE[identifier] = standard_dict
    return standard_dict


@retry(3, 5)
def get_species_from_pubtator(pmids):
    """
    Retrieves species data from PubTator for a given list of PMIDs.

    Parameters:
    - pmids (list): List of PMIDs to retrieve species data for.
    """
    url = f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids={','.join(pmids)}&concepts=species"
    response = requests.get(url)
    lines = response.text.split("\n")
    species_info = {}
    for line in lines:
        parts = line.split("\t")
        if len(parts) >= 6 and parts[4] == "Species":
            species_name = parts[3].lower()
            species_id = parts[5]
            if species_id not in species_info:
                species_info[species_id] = species_name
                logger.info(f"Species found: {species_name} with ID: {species_id}")

    return species_info


def update_record_species(rec, species_data):
    """
    Updates the species in a given record based on abstract, description, and title.

    Parameters:
    - rec (dict): The record to be updated.
    - species_data (dict): Dictionary containing species data, with species name as key and ID as value.
    """

    # Convert species and infectiousAgent to lists if they are not already
    if isinstance(rec.get("species"), dict):
        rec["species"] = [rec["species"]]
    if isinstance(rec.get("infectiousAgent"), dict):
        rec["infectiousAgent"] = [rec["infectiousAgent"]]

    existing_species = {spec["name"]: spec for spec in rec.get("species", [])}
    existing_infectious_agents = {spec["name"]: spec for spec in rec.get("infectiousAgent", [])}

    # Fetch the abstract, description, and title for current record
    abstract = rec.get("abstract", "")
    description = rec.get("description", "")
    title = rec.get("name", "")

    # Iterate over new species data and update the record
    for id, name in species_data.items():
        if name.upper() in abstract or name.upper() in description or name.upper() in title:
            # avoid acronyms ex. PERCH (Pneumonia Etiology Research for Child Health) not a species
            continue
        # If species name is in abstract, description or title and not in existing species, add it to the record
        if (
            (name in abstract.lower() or name in description.lower() or name in title.lower())
            and name not in existing_species
            and name not in existing_infectious_agents
        ):
            standardized_dict = get_species_details(name, id)

            if standardized_dict["classification"] == "host":
                rec["species"] = rec.get("species", []) + [standardized_dict]
                logger.info(f"Added species {name} to record {rec['_id']}")

            elif standardized_dict["classification"] == "infectiousAgent":
                rec["infectiousAgent"] = rec.get("infectiousAgent", []) + [standardized_dict]
                logger.info(f"Added infectious agent {name} to record {rec['_id']}")

            else:
                logger.warning(f"Unknown classification: {standardized_dict['classification']}")


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
                    fund["funder"] = {"@type": "Organization", "name": agency}
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
                    # get species names from pubtator
                    species_data = get_species_from_pubtator(pmids)
                    if species_data:
                        # Update the species field in the record
                        update_record_species(rec, species_data)
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
