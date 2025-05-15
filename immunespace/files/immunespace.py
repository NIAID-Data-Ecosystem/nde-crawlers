import datetime
import logging

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


def fetch_source_organization():
    """
    Fetch the sourceOrganization from the HIPC_correction.json file.
    """
    url = "https://raw.githubusercontent.com/NIAID-Data-Ecosystem/nde-metadata-corrections/main/collections_corrections_production/HIPC_correction.json"
    logger.info("Fetching sourceOrganization from: %s", url)
    response = requests.get(url)
    response.raise_for_status()  # Raise an error if the request fails
    data = response.json()
    return data.get("sourceOrganization")


def parse():
    # Fetch the sourceOrganization
    source_organization = fetch_source_organization()

    url = "https://immunespace.org/api_kb/get_study_id_dropdown"
    logging.info("Making request: %s", url)
    request = requests.get(url)
    logger.info("Request made. HTTP STATUS: %d", request.status_code)
    records = request.json()
    logger.info("Parsing records...")
    for count, record in enumerate(records, start=1):
        if count % 100 == 0:
            logger.info("Parsed %d records", count)

        output = {
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "ImmuneSpace",
                "url": "https://immunespace.org",
                "versionDate": datetime.date.today().isoformat(),
                "dataset": f"https://immunespace.org/query/study/{record['value']}",
            },
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": record["value"],
            "url": f"https://immunespace.org/query/study/{record['value']}",
            "sourceOrganization": source_organization,  # Add sourceOrganization here
            "species": "Homo sapiens",
        }
        yield output

    logger.info("Finished parsing records")
    logger.info("Total records parsed: %d", count)
