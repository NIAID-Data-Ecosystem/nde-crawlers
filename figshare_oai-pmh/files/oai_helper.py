import time
import json
import requests
import logging
import xmltodict
import requests

from requests.adapters import HTTPAdapter, Retry
from xml.etree import ElementTree

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('nde-logger')


def oai_helper():

    count = 0

    response = requests.get(
        f"https://api.figshare.com/v2/oai?verb=ListRecords&metadataPrefix=uketd_dc")

    logger.info("Response size: %s", len(response.content))
    logger.info("Response status: %s", response.status_code)

    response_string = response.text

    response_dict = xmltodict.parse(response_string, dict_constructor=dict)

    try:
        records = response_dict['OAI-PMH']['ListRecords']['record']
    except KeyError:
        logger.info("Could not find records in response")
        logger.info(response_dict)

    for record in records:
        count += 1

        metadata = record['metadata']['uketd_dc:uketddc']
        metadata = {k: v for k, v in metadata.items() if '@' not in k}
        metadata = {k.split(':')[-1]: v for k, v in metadata.items()}
        metadata = {k: [v] if not isinstance(
            v, list) else v for k, v in metadata.items()}

        xml = ElementTree.fromstring(response_string)

        doc = {'header': record['header'],
               'metadata': metadata,
               'xml': ElementTree.tostring(xml, encoding='unicode')
               }

        identifier = record['header']['identifier']

        yield (identifier, json.dumps(doc))

    try:
        resumptionToken = response_dict['OAI-PMH']['ListRecords']['resumptionToken']['#text']
    except KeyError:
        logger.info("Could not find resumptionToken in response")
        logger.info(response_dict)

    while True:
        try:
            s = requests.Session()

            retries = Retry(total=5,
                            backoff_factor=0.1,
                            status_forcelist=[500, 502, 503, 504])

            s.mount('http://', HTTPAdapter(max_retries=retries))

            response = s.get(
                f"https://api.figshare.com/v2/oai?verb=ListRecords&resumptionToken={resumptionToken}")

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
                break

            for record in records:
                count += 1

                metadata = record['metadata']['uketd_dc:uketddc']
                metadata = {k: v for k, v in metadata.items() if '@' not in k}
                metadata = {k.split(':')[-1]: v for k, v in metadata.items()}
                metadata = {k: [v] if not isinstance(
                    v, list) else v for k, v in metadata.items()}

                xml = ElementTree.fromstring(response_string)

                doc = {'header': record['header'],
                       'metadata': metadata,
                       'xml': ElementTree.tostring(xml, encoding='unicode')
                       }

                identifier = record['header']['identifier']

                yield (identifier, json.dumps(doc))

            if count % 100 == 0:
                time.sleep(1)
                logger.info("Retrieved %s records", count)

            if count % 1000 == 0:
                raise StopIteration

            try:
                resumptionToken = response_dict['OAI-PMH']['ListRecords']['resumptionToken']['#text']
            except KeyError:
                logger.info("Could not find resumptionToken in response")
                logger.info(response_dict)
                break

        except StopIteration:
            logger.info('Finished retrieving %s records', count)
            break

    # by set
    # sickle = Sickle('http://www.duo.uio.no/oai/request')
    # sets = sickle.ListSets()  # gets all sets
    # for recs in sets:
    #     for rec in recs:
    #         if rec[0] == 'setSpec':
    #             try:
    #                 print rec[1][0], self.spec_list[rec[1][0]]
    #                 records = sickle.ListRecords(
    #                     metadataPrefix='xoai', set=rec[1][0], ignore_deleted=True)
    #                 self.write_file_and_metadata()
    #             except Exception as e:
    #                 # simple exception handling if not possible to retrieve record
    #                 print('Exception: {}'.format(e))
