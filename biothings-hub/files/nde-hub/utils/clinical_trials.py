"""Creators pulled from clinicaltrials.gov.

Each record needs an `sdPublisher` whose name contains "clinicaltrials.gov" and
whose identifier is the trial's NCT id. Used only by the `vivli` source, as a
`load_data` decorator underneath `@nde_upload_wrapper`.
"""

import functools
import time
from itertools import islice
from typing import Dict, Generator, List

import requests
from config import logger

from .common import retry

STUDIES_URL = "https://clinicaltrials.gov/api/v2/studies"
_BATCH_SIZE = 300


def get_creators(row: Dict) -> List:
    """Build the creator list for one clinical trial's protocol section."""
    creators = []

    responsible_party = row.get("sponsorCollaboratorsModule", {}).get("responsibleParty")
    if responsible_party and responsible_party.get("type") and responsible_party["type"].casefold() != "sponsor":
        creators.append(
            {
                "@type": "Person",
                "name": responsible_party["investigatorFullName"],
                "affiliation": [{"@type": "Organization", "name": responsible_party["investigatorAffiliation"]}],
                "title": responsible_party["investigatorTitle"],
                "role": responsible_party["type"],
            }
        )

    contacts_module = row.get("contactsLocationsModule") or {}
    for contact in contacts_module.get("centralContacts", []):
        creators.append({"@type": "Person", "name": contact["name"], "role": contact["role"]})

    for contact in contacts_module.get("overallOfficials", []):
        creator = {"@type": "Person", "name": contact["name"]}
        # Covers Exception id NCT01693562
        if role := contact.get("role"):
            creator["role"] = role
        if affiliation := contact.get("affiliation"):
            creator["affiliation"] = [{"@type": "Organization", "name": affiliation}]
        creators.append(creator)

    return creators


@retry(3, 5)
def batch_ct(ct_ids: List[str]) -> Dict[str, List[Dict]]:
    """Fetch the creators for a batch of clinicaltrials.gov ids."""
    query = {"query.id": ",".join(ct_ids), "countTotal": "true", "pageSize": 1000}
    request = requests.get(STUDIES_URL, params=query)
    clinical_trials = request.json()

    assert clinical_trials.get("totalCount") == len(ct_ids), (
        "Clinical Trials response does not match ids given to request. "
        f"Response length: {clinical_trials.get('totalCount')}. URL: {request.url}, {query}"
    )

    creators = {}
    for study in clinical_trials.get("studies"):
        study = study.get("protocolSection")
        _id = study["identificationModule"].get("nctId")
        try:
            creators[_id] = get_creators(study)
        except Exception as e:
            logger.error("This id: %s cannot be converted into an creator", _id)
            raise e
    return creators


def _trial_ids(doc):
    for sd_publisher in doc.get("sdPublisher") or []:
        name = sd_publisher.get("name")
        if name and "clinicaltrials.gov" in name.casefold() and sd_publisher.get("identifier"):
            yield sd_publisher["identifier"].strip()


def load_ct_wrapper(func: Generator) -> Generator:
    """Add clinicaltrials.gov creators to the records a `load_data` yields.

    TODO BATCH QUERY WILL FAIL IF THERE IS EVEN 1 INCORRECT ID.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        count = 0
        data = func(*args, **kwargs)
        while True:
            doc_list = list(islice(data, _BATCH_SIZE))
            if not doc_list:
                break

            ct_ids = sorted({ct_id for doc in doc_list for ct_id in _trial_ids(doc)})
            creators_lookup = batch_ct(ct_ids) if ct_ids else {}
            if ct_ids:
                time.sleep(0.2)

            for doc in doc_list:
                creators = [creator for ct_id in _trial_ids(doc) for creator in creators_lookup.get(ct_id, [])]
                if creators:
                    doc_creator = doc.get("creator")
                    doc["creator"] = (list(doc_creator) if doc_creator else []) + creators

                count += 1
                if count % 1000 == 0:
                    logger.info("Processed %s documents", count)
                yield doc
        logger.info("Finished processing %s documents", count)

    return wrapper
