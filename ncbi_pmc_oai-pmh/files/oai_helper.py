import json
import requests
import logging
import requests

from requests.adapters import HTTPAdapter, Retry
from xml.etree import ElementTree

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('nde-logger')


def get_header(record):
    # Used to create a dictionary of the header information
    header = record.find(
        './/{http://www.openarchives.org/OAI/2.0/}header')
    header_dict = {}
    for el in header:
        header_dict[el.tag.split('}')[-1]] = el.text
    return header_dict


def get_metadata(record):
    # Used to create a dictionary of the metadata information
    metadata_xml = record.find(
        './/{http://www.openarchives.org/OAI/2.0/}metadata')
    metadata_dict = {}
    for el in metadata_xml.iter():
        key = el.tag.split('}')[-1]
        if key not in metadata_dict:
            metadata_dict[key] = [x.strip()
                                  for x in el.itertext() if x.strip() != '']
        else:
            text_list = [x.strip()
                         for x in el.itertext() if x.strip() != '']
            for text in text_list:
                metadata_dict[key].append(text)
    metadata_dict.pop('metadata')
    metadata_dict.pop('article')
    return metadata_dict


def oai_helper(last_updated=None):
    # last_updated is a date defined in the update_cache function to retrieve records from a specific date
    # last_updated will default to None when called in load_cache
    if last_updated:
        url = f'https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=ListRecords&metadataPrefix=pmc&from={last_updated}'
        logger.info('Retrieving records from %s', last_updated)

    url = 'https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=ListRecords&metadataPrefix=pmc'
    count = 0

    s = requests.Session()

    # Retry logic to handle 500 errors
    retries = Retry(total=5,
                    backoff_factor=0.1,
                    status_forcelist=[500, 502, 503, 504])

    s.mount('http://', HTTPAdapter(max_retries=retries))

    response = s.get(url)

    logger.info("Response size: %s", len(response.content))
    logger.info("Response status: %s", response.status_code)

    response_string = response.text

    # Check if records exist
    xml = ElementTree.fromstring(response_string)
    records = xml.findall(
        './/{http://www.openarchives.org/OAI/2.0/}record')

    if records == []:
        logger.info("Could not find records in response")
        logger.info(response_string)
        return

    for record in records:
        # Here we check if the record is a deleted record, if it is we ignore it
        header_dict = get_header(record)
        if 'status' in header_dict and header_dict['status'] == 'deleted':
            logger.info('Record %s has been deleted. Not saving record to cache.',
                        header_dict['identifier'])
        else:
            count += 1

            metadata_dict = get_metadata(record)

            record_xml = ElementTree.tostring(record, encoding='unicode')
            identifier = header_dict['identifier']
            doc = {'header': header_dict, 'metadata': metadata_dict,
                   'xml': record_xml}

            logger.info('Current Record: %s', identifier)
            yield (identifier, json.dumps(doc))

    resumptionToken = xml.find(
        './/{http://www.openarchives.org/OAI/2.0/}resumptionToken').text
    if resumptionToken is None:
        logger.info("No resumptionToken found")
        logger.info(response_string)
        return

    # We will repeat the above since we do not have the resumptionToken defined until AFTER the first request

    while True:
        try:
            s = requests.Session()

            retries = Retry(total=5,
                            backoff_factor=0.1,
                            status_forcelist=[500, 502, 503, 504])

            s.mount('http://', HTTPAdapter(max_retries=retries))

            # With the resumptionToken defined we can now start looping through each request
            response = s.get(
                f'https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=ListRecords&resumptionToken={resumptionToken}')

            logger.info("Response size: %s", len(response.content))
            logger.info("Response status: %s", response.status_code)

            response_string = response.text

            xml = ElementTree.fromstring(response_string)
            records = xml.findall(
                './/{http://www.openarchives.org/OAI/2.0/}record')

            if records == []:
                logger.info("Could not find records in response")
                logger.info(response_string)
                raise StopIteration

            for record in records:
                header_dict = get_header(record)
                if 'status' in header_dict and header_dict['status'] == 'deleted':
                    logger.info('Record %s has been deleted. Not saving record to cache.',
                                header_dict['identifier'])
                else:
                    count += 1
                    metadata_dict = get_metadata(record)

                    record_xml = ElementTree.tostring(
                        record, encoding='unicode')
                    identifier = header_dict['identifier']
                    doc = {'header': header_dict, 'metadata': metadata_dict,
                           'xml': record_xml}

                    logger.info('Current Record: %s', identifier)

                    if count % 1000 == 0:
                        raise StopIteration

                    yield (identifier, json.dumps(doc))

            resumptionToken = xml.find(
                './/{http://www.openarchives.org/OAI/2.0/}resumptionToken').text

        except StopIteration:
            logger.info('Finished retrieving %s records', count)
            break
