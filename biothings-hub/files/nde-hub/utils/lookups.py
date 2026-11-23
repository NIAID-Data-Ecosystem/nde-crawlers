"""Lookup collections that hub sources upload to the src database.

A source re-upload replaces its collection with two renames, so between them
the collection does not exist. `find_in` reads each chunk between two checks of
the collection's UUID and reads it again when the collection was missing or
replaced.
"""

import time
from itertools import batched

from biothings.utils.mongo import get_src_db
from config import logger
from pymongo.errors import OperationFailure

_CHUNK_SIZE = 5_000
_ATTEMPTS = 5
_RETRY_SECONDS = 2


def _uuid(db, name):
    """The UUID of collection `name`, or None when it does not exist."""
    for spec in db.list_collections(filter={"name": name}):
        return spec.get("info", {}).get("uuid")
    return None


def ready(name, wait=0):
    """True when collection `name` holds a document.

    With `wait`, a missing or empty collection is checked again every second for
    up to `wait` seconds.
    """
    deadline = time.monotonic() + wait
    while True:
        if get_src_db()[name].find_one({}, {"_id": 1}) is not None:
            return True
        if time.monotonic() >= deadline:
            return False
        time.sleep(1)


def find_in(name, field, values, projection=None):
    """The documents of collection `name` whose `field` is one of `values`."""
    db = get_src_db()
    collection = db[name]
    found = []
    for chunk in batched(dict.fromkeys(values), _CHUNK_SIZE):
        for attempt in range(1, _ATTEMPTS + 1):
            before = _uuid(db, name)
            if before is not None:
                try:
                    docs = list(collection.find({field: {"$in": list(chunk)}}, projection))
                except OperationFailure as error:
                    logger.warning("Reading %s failed (%s); retrying", name, error)
                else:
                    if _uuid(db, name) == before:
                        found.extend(docs)
                        break
            if attempt == _ATTEMPTS:
                raise RuntimeError(f"Lookup collection {name} is missing or was replaced during every read")
            time.sleep(_RETRY_SECONDS)
    return found
