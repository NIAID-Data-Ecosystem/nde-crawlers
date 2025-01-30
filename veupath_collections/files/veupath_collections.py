import datetime
import logging

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


def record_generator():
    sites = [
        "AmoebaDB",
        "CryptoDB",
        "GiardiaDB",
        "HostDB",
        "PlasmoDB",
        "VectorBase",
        "FungiDB",
        "MicrosporidiaDB",
        "ToxoDB",
        "TrichDB",
        "TriTrypDB",
        "PiroplasmaDB",
    ]
    webapps = [
        "amoeba",
        "cryptodb",
        "giardiadb",
        "hostdb",
        "plasmo",
        "vectorbase",
        "fungidb",
        "micro",
        "toxo",
        "trichdb",
        "tritrypdb",
        "piro",
    ]
    for i, site in enumerate(sites):
        url = f"https://{site.casefold()}.org/{webapps[i]}/service/record-types/dataset/searches/AllDatasets/reports/standard"
        logging.info("Making request: %s", url)
        request = requests.get(url)
        logger.info("Request made. HTTP STATUS: %d", request.status_code)
        records = request.json()
        logger.info("Parsing records...")
        for count, record in enumerate(records["records"], start=1):
            if count % 1000 == 0:
                logger.info("Parsed %d records", count)

            url = f"https://{site.casefold()}.org/{webapps[i]}/app/record/dataset/{record['id'][0]['value']}"
            output = {
                "includedInDataCatalog": {
                    "@type": "Dataset",
                    "name": site,
                    "url": f"https://{site.casefold()}.org/{webapps[i]}/app",
                    "versionDate": datetime.date.today().isoformat(),
                    "dataset": url,
                },
                "@context": "http://schema.org/",
                "@type": "Dataset",
                "_id": record["id"][0]["value"],
                "url": url,
            }
            yield output
