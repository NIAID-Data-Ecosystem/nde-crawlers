import datetime
import functools
from typing import Dict, Generator, Iterable

from biothings import config


def add_date(func: Iterable[Dict]) -> Generator[dict, dict, Generator]:
    """Decorator to add a date field to a document.
    The date field is the latest date from the following fields:
    date, dateCreated, dateModified, datePublished
    :param func: a generator function that yields documents
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        gen = func(*args, **kwargs)
        for doc in gen:
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

            yield doc

    return wrapper
