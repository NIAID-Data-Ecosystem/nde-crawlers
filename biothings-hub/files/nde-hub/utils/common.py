"""Small helpers shared by the upload pipeline stages."""

import functools
import os
import time
import traceback

import orjson
from config import logger


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
                logger.info("Retrying failed func %s. Trying attempt %s of %s.", func, attempt + 1, retry_num)
            logger.error("func %s retry failed", func)
            raise Exception("Exceed max retry num: {} failed. Last error: {!r}".format(retry_num, last_exc)) from last_exc

        return wrapper

    return decorator


def batched(iterable, batch_size):
    """Yield lists of at most `batch_size` items from `iterable`."""
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def iter_ndjson(data_folder, filename="data.ndjson"):
    """Yield documents from `<data_folder>/<filename>`."""
    with open(os.path.join(os.fspath(data_folder), filename), "rb") as f:
        for line in f:
            yield orjson.loads(line)


def as_list(value):
    """Normalize a schema field that may be a single object, a list, or missing, to a list."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def dict_entries(doc, field):
    """Yield the dict entries of `doc[field]`, tolerating a single dict or a list."""
    for entry in as_list(doc.get(field)):
        if isinstance(entry, dict):
            yield entry
