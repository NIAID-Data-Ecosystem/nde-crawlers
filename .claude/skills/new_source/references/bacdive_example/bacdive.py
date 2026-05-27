import csv
import re
import time

import bacdive
import requests
import datetime
import dateutil.parser
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

SPARQL_URL = "https://sparql.dsmz.de/api/bacdive"
STRAIN_QUERY = (
    "SELECT ?s WHERE { ?s a <https://purl.dsmz.de/schema/Strain> }"
)
BATCH_SIZE = 100    # v2 fetch endpoint cap
SLEEP = 0.05        # seconds between fetch calls

def insert_value(d, key, value, extend=False):
    """ Insert a value into a dictionary, handling existing keys by converting to lists or extending strings as needed.
    """

    if key in d and not extend:
        if isinstance(d[key], list):
            if isinstance(value, list):
                for item in value:
                    if item not in d[key]:
                        d[key].append(item)
            elif value not in d[key]:
                d[key].append(value)
        else:
            if isinstance(value, list):
                d[key] = [d[key]] + [v for v in value if v != d[key]]
            elif d[key] != value:
                d[key] = [d[key], value]
    elif d.get(key) and extend:
        d[key] = (d.get(key) + " " + value).strip()
    else:
        d[key] = value

def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(val, ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt

def get_all_bacdive_ids():
    """Return a sorted list of every BacDive strain ID via SPARQL."""
    resp = requests.get(
        SPARQL_URL,
        params={"query": STRAIN_QUERY},
        headers={"Accept": "application/sparql-results+json"},
        timeout=120,
    )
    resp.raise_for_status()
    bindings = resp.json()["results"]["bindings"]
    ids = []
    for b in bindings:
        uri = b["s"]["value"]
        ids.append(int(uri.rsplit("/", 1)[-1]))
    ids.sort()
    return ids


def iter_bacdive_records(ids=None, sleep=SLEEP):
    """Yield (bacdive_id, record_dict) for every ID (all of BacDive by default)."""
    if ids is None:
        ids = get_all_bacdive_ids()
    client = bacdive.BacdiveClient()  # public=True -> no auth
    for i in range(0, len(ids), BATCH_SIZE):
        print(f"Fetching IDs {min(i + BATCH_SIZE, len(ids))} out of {len(ids)}...")
        chunk = ids[i:i + BATCH_SIZE]
        if client.search(id=chunk) == 0:
            continue
        for entry in client.retrieve():
            bid = entry.get("General", {}).get("BacDive-ID")
            yield bid, entry
        time.sleep(sleep)


def parse():
    for bid, record in iter_bacdive_records():
        general = record.get("General") or {}
        natc = record.get("Name and taxonomic classification") or {}
        isolation_root = record.get("Isolation, sampling and environmental information") or {}
        sequence_info = record.get("Sequence information") or {}
        literature = record.get("Literature") or {}
        references = record.get("Reference") or []
        morphology = record.get("Morphology") or {}

        _id = general.get("BacDive-ID")
        if not _id:
            continue
        url = f"https://bacdive.dsmz.de/strain/{_id}"

        output = {
            "@context": "http://schema.org/",
            "@type": "Sample",
            "_id": f"bacdive_{_id}",
            "identifier": str(_id),
            "url": url,
            "distribution": [{"@type": "DataDownload", "contentUrl": url}],
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "BacDive",
                "url": "https://bacdive.dsmz.de/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
        }

