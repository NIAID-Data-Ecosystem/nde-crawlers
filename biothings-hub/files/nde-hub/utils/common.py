"""Small helpers shared by the upload pipeline stages."""

import functools
import os
import sqlite3
import time
import traceback
from contextlib import contextmanager

import orjson
from config import logger


DESCRIPTION_ENRICHMENT_TYPES = frozenset({"Dataset", "DataCollection", "ResourceCatalog"})


def retry(retry_num, retry_sleep_sec):
    """Retry a function `retry_num` times, sleeping `retry_sleep_sec` between attempts."""

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_exc = None
            for attempt in range(retry_num):
                try:
                    return func(*args, **kwargs)
                except Exception as err:
                    last_exc = err
                    logger.error(err)
                    logger.error(traceback.format_exc())
                    time.sleep(retry_sleep_sec)
                logger.debug("Retrying failed func %s. Trying attempt %s of %s.", func, attempt + 1, retry_num)
            logger.error("func %s retry failed", func)
            raise Exception(
                "Exceed max retry num: {} failed. Last error: {!r}".format(retry_num, last_exc)
            ) from last_exc

        return wrapper

    return decorator


def iter_ndjson(data_folder, filename="data.ndjson"):
    """Yield documents from `<data_folder>/<filename>`."""
    with open(os.path.join(os.fspath(data_folder), filename), "rb") as f:
        for line in f:
            yield orjson.loads(line)


@contextmanager
def sqlite(path, *setup):
    """Open `path`, run any `setup` DDL, commit on success and always close.

    A connection per call rather than a shared one: the species resolvers write
    from a thread pool and SQLite connections cannot cross threads.
    """
    conn = sqlite3.connect(path)
    try:
        with conn:
            for statement in setup:
                conn.execute(statement)
            yield conn
    finally:
        conn.close()


def as_list(value):
    """Normalize a schema field that may be a single object, a list, or missing, to a list."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def supports_description_enrichment(doc):
    """True for resource records and BioSample-flavoured Sample records."""
    record_types = {value for value in as_list(doc.get("@type")) if isinstance(value, str)}
    if record_types & DESCRIPTION_ENRICHMENT_TYPES:
        return True

    additional_types = {value for value in as_list(doc.get("additionalType")) if isinstance(value, str)}
    return "Sample" in record_types and "BioSample" in additional_types


def dict_entries(doc, field):
    """Yield the dict entries of `doc[field]`, tolerating a single dict or a list."""
    for entry in as_list(doc.get(field)):
        if isinstance(entry, dict):
            yield entry
