import datetime
import functools
import re
import time
import traceback
from typing import Dict, Generator, Iterable

import bson
from config import logger
from lxml import html
from scores import RECOMMENDED_AUGMENTED_FIELDS, RECOMMENDED_FIELDS, REQUIRED_AUGMENTED_FIELDS, REQUIRED_FIELDS


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
    _id = doc.get("_id")
    assert _id, "_id is None"
    unsafe_chars_pattern = re.compile(r"[^a-zA-Z0-9_-]")
    # assert not re.search(unsafe_chars_pattern, doc.get("_id")), "_id: %s contains unsafe URL characters" % doc.get(
    #    "_id"
    # )
    assert doc.get("@type"), "@type is None"
    assert doc.get("includedInDataCatalog"), "includedInDataCatalog is None"
    assert doc.get("version", None) is None, "Remove version field"
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


def is_purely_augmented(field, field_content):
    if field == "includedInDataCatalog":
        return False
    elif isinstance(field_content, (str, bool)):
        return False

    if isinstance(field_content, dict):
        field_content = [field_content]

    if isinstance(field_content, list):
        return all(
            isinstance(item, dict)
            and (item.get("fromPMID", False) or item.get("fromGPT", False) or item.get("fromEXTRACT", False))
            for item in field_content
        )

    return False


def check_augmented_fields(document, augmented_field_list):
    augmented_fields_found = []

    for field in augmented_field_list:
        value = document.get(field, None)

        if isinstance(value, dict) and (
            value.get("fromPMID", False) or value.get("fromGPT", False) or value.get("fromEXTRACT", False)
        ):
            augmented_fields_found.append(field)

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict) and (
                    item.get("fromPMID", False) or item.get("fromGPT", False) or item.get("fromEXTRACT", False)
                ):
                    augmented_fields_found.append(field)
                    break  # Once we find one instance in the list, we can break out

    return augmented_fields_found


def add_metadata_score(document: Dict) -> Dict:
    local_required_fields = list(REQUIRED_FIELDS)
    local_recommended_fields = list(RECOMMENDED_FIELDS)

    # Handle ResourceCatalog types separately to add additional required and recommended fields
    includedInDataCatalog = document.get("includedInDataCatalog")
    if includedInDataCatalog:
        if isinstance(includedInDataCatalog, list):
            for catalog in includedInDataCatalog:
                if catalog.get("@type") == "ResourceCatalog":
                    local_required_fields.append("collectionType")
                    local_recommended_fields.extend(["collectionSize", "hasAPI", "hasDownload"])
                    break
        elif includedInDataCatalog.get("@type") == "ResourceCatalog":
            local_required_fields.append("collectionType")
            local_recommended_fields.extend(["collectionSize", "hasAPI", "hasDownload"])

    required_score = sum(
        1 for field in local_required_fields if field in document and not is_purely_augmented(field, document[field])
    )
    recommended_score = sum(
        1 for field in local_recommended_fields if field in document and not is_purely_augmented(field, document[field])
    )

    required_augmented_fields = check_augmented_fields(document, REQUIRED_AUGMENTED_FIELDS)
    recommended_augmented_fields = check_augmented_fields(document, RECOMMENDED_AUGMENTED_FIELDS)

    total_required = len(local_required_fields)
    total_recommended = len(local_recommended_fields)
    total_required_augmented = len(REQUIRED_AUGMENTED_FIELDS)
    total_recommended_augmented = len(RECOMMENDED_AUGMENTED_FIELDS)

    existing_required_fields = [
        field
        for field in local_required_fields
        if field in document and not is_purely_augmented(field, document[field])
    ]
    existing_recommended_fields = [
        field
        for field in local_recommended_fields
        if field in document and not is_purely_augmented(field, document[field])
    ]

    document.setdefault("_meta", {}).update(
        {
            "required_augmented_fields": required_augmented_fields,
            "recommended_augmented_fields": recommended_augmented_fields,
            "required_fields": existing_required_fields,
            "recommended_fields": existing_recommended_fields,
            "completeness": {
                "total_score": required_score + recommended_score,
                "total_max_score": total_required + total_recommended,
                "required_score": required_score,
                "required_ratio": round(required_score / total_required, 2) if total_required > 0 else 0,
                "required_max_score": total_required,
                "recommended_score": recommended_score,
                "recommended_score_ratio": (
                    round(recommended_score / total_recommended, 2) if total_recommended > 0 else 0
                ),
                "recommended_max_score": total_recommended,
                "augmented_required_ratio": (
                    round(len(required_augmented_fields) / total_required, 2) if total_required > 0 else 0
                ),
                "augmented_recommended_ratio": (
                    round(len(recommended_augmented_fields) / total_recommended, 2) if total_recommended > 0 else 0
                ),
                "total_required_augmented": total_required_augmented,
                "total_recommended_augmented": total_recommended_augmented,
            },
        }
    )

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
            add_metadata_score(doc)

            # Remove HTML tags from description field
            if doc.get("description"):
                doc["description"] = html.fromstring(doc["description"]).text_content()

            # Everything past this point should always be last in order
            check_schema(doc)
            doc["_id"] = doc["_id"].casefold()
            bson_size = len(bson.BSON.encode(doc))
            if bson_size < (16 * 1024 * 1024):  # Check if less than 16 MB
                yield doc
            else:
                # Handle oversized document: log, skip, or process differently
                # For example, log a warning and skip this document
                logger.warning("Document %s exceeds MongoDB's size limit: %s bytes", doc["_id"], bson_size)

    return wrapper
