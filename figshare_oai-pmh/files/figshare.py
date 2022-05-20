import json
import time
import logging
from datetime import datetime
from sickle import Sickle


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')

# during testing used this list to check if any properties were NOT in the list, if not then print the property
# properties = ["title", "creator", "subject", "description", "date", "publisher", "type", "identifier", "language", "relation",
#               "rights", "license", "abstract", "isReferencedBy", "issued", "institution", "department", "sponsor", "grantnumber", "authoridentifier", "embargodate", "qualificationname", "qualificationlevel", "advisor"]
# used to add missing properties and print later
# missing = []

# TODO
# known covid categories according to docs https://covid19.figshare.com/f/faqs to siphon off later
# covid_keywords = ['covid-19', '2019-ncov', 'coronavirus', 'sars-cov-2', '2019ncov',
#                   'covid19', 'corona virus', 'sarscov2', 'covid2019', 'covid_19', 'sars-cov2', 'covid 19']


def parse():
    # connect to the website
    logger.info("Parsing records")
    sickle = Sickle('https://api.figshare.com/v2/oai', max_retries=3)
    records = sickle.ListRecords(
        metadataPrefix='uketd_dc', ignore_deleted=True)

    # used to test single record
    # record = sickle.GetRecord(
    #     identifier='oai:figshare.com:article/5849037', metadataPrefix='uketd_dc')

    count = 0
    while True:
        try:
            # get the next item
            record = records.next()
            metadata = record.metadata

            # testing for missing properties
            # for key in metadata:
            #     if key not in properties:
            #         missing.append(key)

            count += 1
            if count % 10 == 0:
                # figshare requires us to parse 10 records a second for the oai-pmh
                time.sleep(1)
            if count % 1000 == 0:
                logger.info("Parsed %s records", count)
                # logging missing properties
                # if len(missing):
                #     logger.info(f'Missing {missing}')

            output = {
                "includedInDataCatalog": {"name": "Figshare", 'versionDate': datetime.today().isoformat(), 'url': "https://figshare.com"},
            }
            if title := metadata.get('title'):
                output['name'] = title[0]

            if creators := metadata.get('creator'):
                creator_list = []
                for creator in creators:
                    creator_list.append({"name": creator})
                output['author'] = creator_list

            # TODO
            # checking if covid related article for outbreak api
                #     for keyword in subject:
                #         if keyword.lower() in covid_keywords:
                #             output['outbreakapi'] = True
                #         else:
                #             output['outbreakapi'] = False
            if subject := metadata.get('subject'):
                output['keywords'] = subject

            if description := metadata.get('description'):
                output['description'] = description[0]
            if description == None:
                if abstract := metadata.get('abstract'):
                    output['description'] = abstract

            if date := metadata.get('date'):
                output['dateModified'] = datetime.strptime(
                    date[0], '%Y-%m-%dT%H:%M:%S%z').strftime('%Y-%m-%d')

            # sdPublsher = publisher > instituion/department > figshare
            if publisher := metadata.get('publisher'):
                output['sdPublisher'] = {'name': publisher[0]}
            if publisher == None:
                institution = metadata.get('institution')
                department = metadata.get('department')
                if institution and department:
                    if institution[0] != department[0]:
                        output['sdPublisher'] = institution[0] + \
                            '/' + department[0]
                    else:
                        output['sdPublisher'] = institution[0]
                else:
                    output['sdPublisher'] = 'figshare'

            if type := metadata.get('type'):
                if len(type) > 1:
                    if type[0] == type[1]:
                        if type[0] == 'Software':
                            output['@type'] = 'ComputationalTool'
                        else:
                            output['@type'] = type[0]
                    else:
                        output['@type'] = 'Collection'
                        output['hasPart'] = {'@type': type}
                else:
                    output['@type'] = 'Collection'
                    output['hasPart'] = {'@type': type[0]}

            if identifier := metadata.get('identifier'):
                for el in identifier:
                    # [None, 'https://ndownloader.figshare.com/files/35101450']
                    if el is not None:
                        if 'ndownloader' in el:
                            output['distribution'] = {'url': el}
                        elif '10.' in el:
                            output['doi'] = el
            if identifier == None:
                if reference := metadata.get('isReferencedBy'):
                    output['doi'] = reference[0]

            if language := metadata.get('language'):
                output['language'] = language[0]
            if relation := metadata.get('relation'):
                output['url'] = relation[0]
                output['_id'] = 'Figshare_' + relation[0].split('/')[-1]
                output['identifier'] = relation[0].split('/')[-1]
            if license := metadata.get('license'):
                output['license'] = license[0]
            if issued := metadata.get('issued'):
                output['datePublished'] = issued[0]
            if sponsor := metadata.get('sponsor'):
                grantnumber = metadata.get('grantnumber')
                if grantnumber:
                    output['funding'] = {'funder': {
                        'name': sponsor[0], 'identifier': grantnumber[0]}}
                else:
                    output['funding'] = {'funder': {'name': sponsor[0]}}

            yield output

        except StopIteration:
            logger.info("Finished Parsing. Total Records: %s", count)
            # if StopIteration is raised, break from loop
            break
