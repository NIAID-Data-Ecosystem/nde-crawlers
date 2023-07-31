import datetime
import functools
import time
import traceback
from typing import Dict, Generator, Iterable

from config import logger


def retry(retry_num, retry_sleep_sec):
    """
    retry help decorator.
    :param retry_num: the retry num; retry sleep sec
    :return: decorator
    """

    def decorator(func):
        """decorator"""

        # preserve information about the original function, or the func name will be "wrapper" not "func"
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """wrapper"""
            for attempt in range(retry_num):
                try:
                    return func(*args, **kwargs)  # should return the raw function's return value
                except Exception as err:
                    logger.error(err)
                    logger.error(traceback.format_exc())
                    time.sleep(retry_sleep_sec)
                logger.info("Retrying failed func %s. Trying attempt %s of %s.", func, attempt + 1, retry_num)
            logger.error("func %s retry failed", func)
            raise Exception("Exceed max retry num: {} failed".format(retry_num))

        return wrapper

    return decorator


def check_schema(doc: Dict) -> Dict:
    """
    check doc before inserting into MongoDb
    """

    assert isinstance(doc, dict), "doc is not a dict"
    assert doc.get("_id"), "_id is None"
    assert doc.get("@type"), "@type is None"
    assert doc.get("includedInDataCatalog"), "includedInDataCatalog is None"
    if coa := doc.get("conditionsOfAccess"):
        enum = ["Open", "Restricted", "Closed", "Embargoed"]
        assert coa in enum, "%s is not a valid conditionsOfAccess. Allowed conditionsOfAccess: %s" % (coa, enum)
    return doc


def add_date(doc: Dict) -> Dict:
    """
    The date field is the latest date from the following fields:
    date, dateCreated, dateModified, datePublished
    :param func: a generator function that yields documents
    """

    dates = []
    if doc.get("date"):
        dates.append(doc.get("date"))
    if doc.get("dateCreated"):
        dates.append(doc.get("dateCreated"))
    if doc.get("dateModified"):
        dates.append(doc.get("dateModified"))
    if doc.get("datePublished"):
        dates.append(doc.get("datePublished"))
    if dates:
        dates.sort()
        date = datetime.datetime.fromisoformat(dates[-1]).date().isoformat()
        doc["date"] = date

    return doc


def merge_duplicates(doc: Dict) -> Dict:
    if _id := doc.get("doi"):
        if isinstance(_id, list):
            if len(_id) == 1:
                assert isinstance(_id[0], str), "Doi is not a string %s" % _id
                doc["_id"] = _id[0]
        else:
            assert isinstance(_id, str), "Doi is not a string %s" % _id
            if _id.casefold() == "none":
                doc.pop("doi")
            else:
                doc["_id"] = _id

    return doc

MAPPING_SCORES = {
    "abstract": 0.4,
    "alternateName": 0.4,
    "author": {
        "familyName": 0.2,
        "givenName": 0.2,
        "name": 0.4,
    },
    "citation": {
        "author": {
            "familyName": 0.1,
            "givenName": 0.1,
            "name": 0.2,
        },
        "doi": 0.2,
        "name": 0.1,
        "pmid": 0.1,
        "url": 0.1,
    },
    "citedBy": {
        "doi": 0.1,
        "name": 0.1,
        "pmid": 0.1,
        "url": 0.1,
    },
    "contentUrl": 0.2,
    "dateCreated": 0.2,
    "dateModified": 0.2,
    "datePublished": 0.2,
    "distribution": {
        "contentUrl": 0.1,
    },
    "doi": 0.3,
    "funding": {
        "url": 0.3,
        "name": 0.1,
        "funder": {
            "name": 0.1,
            "url": 0.3,
        },
    },
    "healthCondition": {
        "name": 0.3,
        "isCurated": 0.3,
    },
    "infectiousAgent": {
        "name": 0.3,
        "isCurated": 0.3,
    },
    "isBasedOn": {
        "doi": 0.3,
        "name": 0.1,
        "pmid": 0.3,
        "url": 0.1,
    },
    "keywords": 0.1,
    "measurementTechnique": {
        "name": 0.1,
    },
    "sdPublisher": {
        "name": 0.3,
        "url": 0.3,
    },
    "species": {
        "name": 0.3,
        "isCurated": 0.3,
    },
}
def calculate_score(data, mapping):
    score = 0
    for key, value in data.items():
        if key == "description":
            score += 0.5 * min(1, len(value) / 500)  # Normalized based on an arbitrary max length of 500
        elif key in mapping:
            if isinstance(value, dict):
                score += calculate_score(value, mapping[key])
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        score += calculate_score(item, mapping[key])
            else:
                score += mapping[key]
    return score

def add_metadata_score(document, mapping):
    document["metadata_score"] = calculate_score(document, mapping)
    return document


def nde_upload_wrapper(func: Iterable[Dict]) -> Generator[dict, dict, Generator]:
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        gen = func(*args, **kwargs)
        for doc in gen:
            # dictionaries are mutable so we dont need to reassign
            add_date(doc)
            # TODO FOR TESTING
            # merge_duplicates(doc)
            add_metadata_score(doc, MAPPING_SCORES)
            # This will always be last
            check_schema(doc)
            yield doc

    return wrapper


def zenodo_upload_wrapper(func: Iterable[Dict]) -> Generator[dict, dict, Generator]:
    """REMOVE THIS WHEN WE CAN GET A SUCCESSFUL RUN OF ZENODO"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        gen = func(*args, **kwargs)
        for doc in gen:
            add_date(doc)
            # TODO FOR TESTING
            # merge_duplicates(doc)
            yield doc

    return wrapper
