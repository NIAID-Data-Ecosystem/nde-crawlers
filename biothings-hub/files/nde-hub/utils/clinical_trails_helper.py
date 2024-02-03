import functools
import time
from itertools import islice
from typing import Dict, Generator, List

import requests
from config import logger

from .utils import retry


def get_creators(row: Dict) -> List:
    """
    Takes a single study from the clinical trials API and returns a list of creators.
    """
    creators = []
    responsible_party = row.get("sponsorCollaboratorsModule", {}).get("responsibleParty")
    if responsible_party and responsible_party.get("type"):
        if responsible_party["type"].casefold() != "sponsor":
            obj = {}
            obj["@type"] = "Person"
            obj["name"] = row["sponsorCollaboratorsModule"]["responsibleParty"]["investigatorFullName"]
            obj["affiliation"] = [
                {
                    "@type": "Organization",
                    "name": row["sponsorCollaboratorsModule"]["responsibleParty"]["investigatorAffiliation"],
                }
            ]
            obj["title"] = row["sponsorCollaboratorsModule"]["responsibleParty"]["investigatorTitle"]
            obj["role"] = row["sponsorCollaboratorsModule"]["responsibleParty"]["type"]
            creators.append(obj)
    if "contactsLocationsModule" in row.keys():
        if row.get("contactsLocationsModule"):
            if "centralContacts" in row["contactsLocationsModule"].keys():
                contacts = row["contactsLocationsModule"]["centralContacts"]
                for contact in contacts:
                    obj = {}
                    obj["@type"] = "Person"
                    obj["name"] = contact["name"]
                    obj["role"] = contact["role"]
                    creators.append(obj)
            if "overallOfficials" in row["contactsLocationsModule"].keys():
                contacts = row["contactsLocationsModule"]["overallOfficials"]
                for contact in contacts:
                    obj = {}
                    obj["@type"] = "Person"
                    obj["name"] = contact["name"]
                    # Covers Exception id NCT01693562
                    if role := contact.get("role"):
                        obj["role"] = role
                    if "affiliation" in contact.keys():
                        obj["affiliation"] = [{"@type": "Organization", "name": contact["affiliation"]}]
                    creators.append(obj)

    return creators


@retry(3, 5)
def batch_ct(ct_ids: List[str]) -> Dict[str, List[Dict]]:
    """
    Takes a list of clinicaltrials.gov ids and returns a dictionary of creators.
    """

    creators = {}
    string_ct_ids = ",".join(ct_ids)
    query = {"query.id": string_ct_ids, "countTotal": "true", "pageSize": 1000}
    url = "https://clinicaltrials.gov/api/v2/studies"
    request = requests.get(url, params=query)
    clinical_trials = request.json()

    assert clinical_trials.get("totalCount") == len(
        ct_ids
    ), f"Clinical Trials response does not match ids given to request. Response length: {clinical_trials.get('totalCount')}. URL: {request.url}, {query}"

    for study in clinical_trials.get("studies"):
        study = study.get("protocolSection")
        _id = study["identificationModule"].get("nctId")

        try:
            creators[_id] = get_creators(study)
        except Exception as e:
            logger.error("This id: %s cannot be converted into an creator", _id)
            raise e

    return creators


def load_ct_wrapper(func: Generator) -> Generator:
    """
    Takes a generator function and yields a generator function that adds creators from clinicaltrials.gov.
    Each record needs to have a sdPublisher field with a name field that contains "clinicaltrials.gov" and an identifier field that contains the clinicaltrials.gov id.
    TODO BATCH QUERY WILL FAIL IF THERE IS EVEN 1 INCORRECT ID.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        count = 0
        data = func(*args, **kwargs)
        while True:
            # take a chunk/slice of records
            doc_list = list(islice(data, 300))
            if not doc_list:
                break

            # get only the clinicaltrail.gov ids
            ct_ids = [
                sd_publisher.get("identifier")
                for doc in doc_list
                if "sdPublisher" in doc
                for sd_publisher in doc.get("sdPublisher")
                if sd_publisher.get("name") and "clinicaltrials.gov" in sd_publisher.get("name").casefold()
            ]
            # remove None values and whitespace
            ct_ids = [_id.strip() for _id in ct_ids if _id]
            # remove duplicates
            ct_ids = [*set(ct_ids)]

            if ct_ids:
                # make request
                creators_lookup = batch_ct(ct_ids)
                time.sleep(0.2)

            for doc in doc_list:
                creators = []
                if sd_publishers := doc.get("sdPublisher"):
                    for sd_publisher in sd_publishers:
                        if sd_publisher.get("name") and "clinicaltrials.gov" in sd_publisher.get("name").casefold():
                            creator = creators_lookup.get(sd_publisher.get("identifier"))
                            creators += creator
                if creators:
                    if doc_creator := doc.get("creator"):
                        if isinstance(doc_creator, list):
                            doc["creator"] += creators
                        else:
                            doc["creator"] = list(doc_creator) + creators
                    else:
                        doc["creator"] = creators

                count += 1
                if count % 1000 == 0:
                    logger.info("Processed %s documents", count)
                yield doc
        logger.info("Finished processing %s documents", count)

    return wrapper
