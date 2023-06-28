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


# def check_schema(func: Iterable[Dict]) -> Generator[dict, dict, Generator]:
#     """
#     check doc before inserting into MongoDb
#     :param func: a generator function that yields documents
#     """

#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         gen = func(*args, **kwargs)
#         for doc in gen:
#             assert isinstance(doc, dict), "doc is not a dict"
#             assert doc.get("_id"), "_id is None"
#             assert doc.get("@type"), "@type is None"
#             assert doc.get("includedInDataCatalog"), "includedInDataCatalog is None"
#             if coa := doc.get("conditionsOfAccess"):
#                 enum = ["Open", "Restricted", "Closed", "Embargoed"]
#                 assert coa in enum, "%s is not a valid conditionsOfAccess. Allowed conditionsOfAccess: %s" % (coa, enum)
#             yield doc

#     return wrapper


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
                doc["_id"] = _id[0]
        else:
            assert isinstance(_id, str), "Doi is not a string %s" % _id
            doc["_id"] = _id

    return doc


def nde_upload_wrapper(func: Iterable[Dict]) -> Generator[dict,dict,Generator]:
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        gen = func(*args, **kwargs)
        for doc in gen:
            add_date(doc)
            merge_duplicates(doc)
            check_schema(doc)
            yield doc

    return wrapper
