import time
import logging
from datetime import datetime

from sickle import Sickle
from pprint import pprint


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')

properties = ["title", "creator", "subject", "description", "date", "publisher", "type", "identifier", "language", "relation",
              "rights", "license", "abstract", "isReferencedBy", "issued", "institution", "department", "sponsor", "grantnumber", "authoridentifier", "embargodate", "qualificationname", "qualificationlevel", "advisor"]
missing = []
# connect to the website
sickle = Sickle('https://api.figshare.com/v2/oai', max_retries=3)
records = sickle.ListRecords(
    metadataPrefix='uketd_dc', ignore_deleted=True)
# record = sickle.GetRecord(
#     identifier='oai:figshare.com:article/139141', metadataPrefix='uketd_dc')
# pprint(record.metadata)

count = 0
while True:
    try:
        # get the next item
        record = records.next()
        metadata = record.metadata
        for key in metadata:
            if key not in properties:
                missing.append(key)
        count += 1
        if count % 10 == 0:
            time.sleep(1)
            logger.info("Parsed %s records", count)
            if len(missing):
                logger.info(f'Missing {missing}')

        output = {
            "includedInDataCatalog": {"name": "Figshare"},
            'versionDate': datetime.today().isoformat()
        }
        if title := metadata.get('title'):
            output['name'] = title[0]
        if creators := metadata.get('creator'):
            creator_list = []
            for creator in creators:
                creator_list.append({"name": creator})
            output['author'] = creator_list
        if subject := metadata.get('subject'):
            output['keywords'] = subject
        if description := metadata.get('description'):
            output['description'] = description[0]
        if date := metadata.get('date'):
            output['dateModified'] = datetime.strptime(
                date[0], '%Y-%m-%dT%H:%M:%S%z').strftime('%Y-%m-%d')
        if publisher := metadata.get('publisher'):
            output['sdPublisher'] = publisher[0]
        if type := metadata.get('type'):
            output['@type'] = type[0]
        if identifier := metadata.get('identifier'):
            output['identifier'] = identifier[0]
        if distribution := metadata.get('identifier'):
            if len(distribution) > 1:
                output['distribution'] = {'url': distribution[1]}
        if language := metadata.get('language'):
            output['language'] = language[0]
        if relation := metadata.get('relation'):
            output['url'] = relation[0]
            output['_id'] = 'Figshare_' + relation[0].split('/')[-1]
        if license := metadata.get('license'):
            output['license'] = license[0]
        # if reference := metadata.get('isReferencedBy'):
        #     output['url'] = reference[0]
        if issued := metadata.get('issued'):
            output['datePublished'] = issued[0]
        if sponsor := metadata.get('sponsor'):
            output['funding'] = {'funder': {'name': sponsor[0]}}
        if grantnumber := metadata.get('grantnumber'):
            output['funding'] = {'funder': {'identifier': grantnumber[0]}}
    except StopIteration:
        logger.info("Finished Parsing. Total Records: %s", count)
        # if StopIteration is raised, break from loop
        break
