import json
import logging
import time
from datetime import datetime, timedelta
from xml.etree import ElementTree

import requests
import xmltodict
from requests.adapters import HTTPAdapter, Retry

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


def clean_metadata(metadata):
    # Used to clean up dictionary keys and convert values to lists
    result = {}
    for k, v in metadata.items():
        if "@" not in k:
            if not isinstance(v, list):
                result[k.split(":")[-1]] = [v]
            else:
                result[k.split(":")[-1]] = v
    return result


def date_range(start_date, end_date, delta):
    current_date = start_date
    while current_date <= end_date:
        yield current_date
        current_date += delta


def oai_helper(last_updated=None):
    start_date = datetime(1975, 1, 1)
    end_date = datetime.now()
    delta = timedelta(days=365)  # One year interval

    for single_date in date_range(start_date, end_date, delta):
        from_date = single_date.strftime("%Y-%m-%dT%H:%M:%SZ")
        until_date = (single_date + delta - timedelta(seconds=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
        logger.info("Fetching records from %s to %s", from_date, until_date)

        if last_updated:
            url = f"https://api.figshare.com/v2/oai?verb=ListRecords&metadataPrefix=uketd_dc&from={last_updated}&set=item_type_3"
            logger.info("Retrieving records from %s", last_updated)
        else:
            url = f"https://api.figshare.com/v2/oai?verb=ListRecords&metadataPrefix=uketd_dc&set=item_type_3&from={from_date}&until={until_date}"

        count = 0

        # Retry logic to handle 500 errors
        s = requests.Session()
        retries = Retry(
            total=60, backoff_factor=60, status_forcelist=[500, 502, 503, 504]
        )  # Retry every minute for 60 minutes
        s.mount("http://", HTTPAdapter(max_retries=retries))

        now = datetime.now()
        response = s.get(url)

        logger.info("Response size: %s", len(response.content))
        logger.info("Response status: %s", response.status_code)

        response_string = response.text
        response_dict = xmltodict.parse(response_string, dict_constructor=dict)

        # Check if records exist
        try:
            records = response_dict["OAI-PMH"]["ListRecords"]["record"]
        except KeyError:
            logger.info("Could not find records in response")
            logger.info(response_dict)
            continue

        if not isinstance(records, list):
            records = [records]

        for record in records:
            # Here we check if the record is a deleted record, if it is we ignore it
            header = record["header"]
            identifier = header["identifier"]
            logger.info("Current Record: %s", identifier)

            if "status" in header and header["status"] == "deleted":
                logger.info("Record %s has been deleted. Not saving record to cache.", header["identifier"])
            else:
                count += 1
                metadata = record["metadata"]["uketd_dc:uketddc"]
                metadata = clean_metadata(metadata)

                xml = ElementTree.fromstring(response_string)
                doc = {"header": header, "metadata": metadata, "xml": ElementTree.tostring(xml, encoding="unicode")}

                yield (identifier, json.dumps(doc))
        end = datetime.now()
        logger.info("Time taken to request and store response: %s", end - now)
        try:
            resumptionToken = response_dict["OAI-PMH"]["ListRecords"]["resumptionToken"]["#text"]
        except KeyError:
            logger.info("Could not find resumptionToken in response")
            logger.info(response_dict)
            continue

        # We will repeat the above since we do not have the resumptionToken defined until AFTER the first request
        while True:
            logger.info("Continuing to next page...")
            try:
                # With the resumptionToken defined we can now start looping through each request
                now = datetime.now()
                response = s.get(f"https://api.figshare.com/v2/oai?verb=ListRecords&resumptionToken={resumptionToken}")

                logger.info("Response size: %s", len(response.content))
                logger.info("Response status: %s", response.status_code)

                response_string = response.text
                response_dict = xmltodict.parse(response_string, dict_constructor=dict)

                try:
                    records = response_dict["OAI-PMH"]["ListRecords"]["record"]
                except KeyError:
                    logger.info("Could not find records in response")
                    logger.info(response_dict)
                    logger.info("Checking if there is another page")
                    try:
                        resumptionToken = response_dict["OAI-PMH"]["ListRecords"]["resumptionToken"]["#text"]
                        logger.info("Found resumptionToken: %s", resumptionToken)
                        continue
                    except KeyError:
                        logger.info("Could not find resumptionToken in response")
                        logger.info(response_dict)
                        break
                if not isinstance(records, list):
                    records = [records]

                for record in records:
                    try:
                        header = record["header"]
                    except (KeyError, TypeError):
                        logger.info("Could not find header in record")
                        logger.info(record)
                        continue
                    identifier = header["identifier"]

                    logger.info("Current Record: %s", identifier)

                    if header and "status" in header and header["status"] == "deleted":
                        logger.info("Record %s has been deleted. Not saving record to cache.", header["identifier"])
                    else:
                        count += 1

                        metadata = record["metadata"]["uketd_dc:uketddc"]
                        metadata = clean_metadata(metadata)

                        xml = ElementTree.fromstring(response_string)
                        doc = {
                            "header": header,
                            "metadata": metadata,
                            "xml": ElementTree.tostring(xml, encoding="unicode"),
                        }

                        yield (identifier, json.dumps(doc))
                end = datetime.now()
                logger.info("Time taken to request and store response: %s", end - now)

                if count % 100 == 0:
                    time.sleep(1)
                    logger.info("Sleeping for 1 second")
                    logger.info("Retrieved %s records", count)

                try:
                    logger.info("Getting resumptionToken...")
                    resumptionToken = response_dict["OAI-PMH"]["ListRecords"]["resumptionToken"]["#text"]
                    logger.info("Found resumptionToken: %s", resumptionToken)
                except KeyError:
                    logger.info("Could not find resumptionToken in response")
                    logger.info(response_dict)
                    break

            except StopIteration:
                logger.info("Finished retrieving %s records", count)
                break
