import time
import json
import requests
import logging
import xmltodict
import requests

from datetime import datetime
from requests.adapters import HTTPAdapter, Retry
from xml.etree import ElementTree

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('nde-logger')


def clean_metadata(metadata):
    # Used to clean up dictionary keys and convert values to lists
    result = {}
    for k, v in metadata.items():
        if '@' not in k:
            if not isinstance(v, list):
                result[k.split(':')[-1]] = [v]
            else:
                result[k.split(':')[-1]] = v
    return result


def oai_helper(last_updated=None):
    # last_updated is a date defined in the update_cache function to retrieve records from a specific date
    # last_updated will default to None when called in load_cache
    if last_updated:
        url = f'https://api.figshare.com/v2/oai?verb=ListRecords&metadataPrefix=uketd_dc&from={last_updated}'
        logger.info('Retrieving records from %s', last_updated)

    url = 'https://api.figshare.com/v2/oai?verb=ListRecords&metadataPrefix=uketd_dc'
    count = 0

    # Retry logic to handle 500 errors

    s = requests.Session()

    retries = Retry(total=10,
                    backoff_factor=0.1,
                    status_forcelist=[500, 502, 503, 504])

    s.mount('http://', HTTPAdapter(max_retries=retries))

    now = datetime.now()
    response = s.get(url)

    logger.info("Response size: %s", len(response.content))
    logger.info("Response status: %s", response.status_code)

    response_string = response.text

    response_dict = xmltodict.parse(
        response_string, dict_constructor=dict)

    # Check if records exist
    try:
        records = response_dict['OAI-PMH']['ListRecords']['record']
    except KeyError:
        logger.info("Could not find records in response")
        logger.info(response_dict)
        return

    for record in records:
        # Here we check if the record is a deleted record, if it is we ignore it
        header = record['header']
        identifier = header['identifier']
        logger.info('Current Record: %s', identifier)

        if 'status' in header and header['status'] == 'deleted':
            logger.info('Record %s has been deleted. Not saving record to cache.',
                        header['identifier'])
        else:
            count += 1
            metadata = record['metadata']['uketd_dc:uketddc']
            metadata = clean_metadata(metadata)

            xml = ElementTree.fromstring(response_string)

            doc = {'header': header,
                   'metadata': metadata,
                   'xml': ElementTree.tostring(xml, encoding='unicode')
                   }

            yield (identifier, json.dumps(doc))
    end = datetime.now()
    logger.info(
        'Time taken to request and store response: %s', end - now)
    try:
        resumptionToken = response_dict['OAI-PMH']['ListRecords']['resumptionToken']['#text']
    except KeyError:
        logger.info("Could not find resumptionToken in response")
        logger.info(response_dict)
        return

    # We will repeat the above since we do not have the resumptionToken defined until AFTER the first request

    while True:
        logger.info("Continuing to next page...")
        try:
            # With the resumptionToken defined we can now start looping through each request
            now = datetime.now()
            response = s.get(
                f'https://api.figshare.com/v2/oai?verb=ListRecords&resumptionToken={resumptionToken}')

            logger.info("Response size: %s", len(response.content))
            logger.info("Response status: %s", response.status_code)

            response_string = response.text

            response_dict = xmltodict.parse(
                response_string, dict_constructor=dict)

            try:
                records = response_dict['OAI-PMH']['ListRecords']['record']
            except KeyError:
                logger.info("Could not find records in response")
                logger.info(response_dict)
                raise StopIteration

            for record in records:
                header = record['header']
                identifier = header['identifier']

                logger.info('Current Record: %s', identifier)

                if 'status' in header and header['status'] == 'deleted':
                    logger.info('Record %s has been deleted. Not saving record to cache.',
                                header['identifier'])
                else:
                    count += 1

                    metadata = record['metadata']['uketd_dc:uketddc']
                    metadata = clean_metadata(metadata)

                    xml = ElementTree.fromstring(response_string)

                    doc = {'header': header,
                           'metadata': metadata,
                           'xml': ElementTree.tostring(xml, encoding='unicode')
                           }

                    yield (identifier, json.dumps(doc))
            end = datetime.now()
            logger.info(
                'Time taken to request and store response: %s', end - now)

            if count % 100 == 0:
                time.sleep(1)
                logger.info('Sleeping for 1 second')
                logger.info("Retrieved %s records", count)

            # test small number of records
            # if count % 1000 == 0:
            #     raise StopIteration

            try:
                logger.info('Getting resumptionToken...')
                resumptionToken = response_dict['OAI-PMH']['ListRecords']['resumptionToken']['#text']
                logger.info('Gottem: %s', resumptionToken)
            except KeyError:
                logger.info("Could not find resumptionToken in response")
                logger.info(response_dict)
                raise StopIteration

        except StopIteration:
            logger.info('Finished retrieving %s records', count)
            break
