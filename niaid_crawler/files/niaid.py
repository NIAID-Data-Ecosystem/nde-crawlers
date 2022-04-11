import json
import requests
import logging
import validators

# from pprint import pprint


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def parse():
    # schema
    # https://docs.google.com/spreadsheets/d/1XlpvoeSWSiqfw1pPAHHFTjr9Wuowj3TA-eCsyGEj0ng/edit#gid=0

    url = "https://accessclinicaldata.niaid.nih.gov/guppy/graphql"
    query = """query ($filter: JSON) {
        clinical_trials (filter: $filter, first: 10000, accessibility: accessible) {
            title,cmc_unique_id,brief_summary,data_availability_date,most_recent_update,data_available,creator,nct_number,condition,clinical_trial_website,publications,data_available_for_request
        }
    }"""

    r = requests.post(url, json={'query': query})

    # parse json string and convert to dictionary
    json_obj = json.loads(r.text)

    # grab list of trials from dictionary
    trials = json_obj['data']['clinical_trials']

    count = 0

    for trial in trials:
        try:
            trial['name'] = trial.pop('title')
            trial['identifier'] = trial.pop('cmc_unique_id')
            trial['description'] = trial.pop('brief_summary')
            trial['datePublished'] = trial.pop('data_availability_date')
            trial['dateModified'] = trial.pop('most_recent_update')
            trial['additionalType'] = trial.pop('data_available')
            trial['funding'] = [{'funder': {'name': trial.pop('creator')}}]
            trial['nctid'] = trial.pop('nct_number')
            trial['infectiousDisease'] = trial.pop('condition')
            trial['mainEntityOfPage'] = trial.pop('clinical_trial_website')

            # check if url is valid then take pubmed id
            citation_URL = trial.pop('publications')
            if citation_URL is not None and validators.url(citation_URL):
                if 'pubmed' in citation_URL:
                    trial['pmids'] = citation_URL.split('/')[-2]
                # To convert doi id to pubmed id, use requests library to search doi id on pubmed search engine, catch redirect and take the pubmed id from url
                elif 'doi' in citation_URL:
                    doi_id = citation_URL.split('/')[-1]
                    r = requests.get("https://pubmed.ncbi.nlm.nih.gov/?term=" + doi_id)
                    if 'pubmed' in r.url:
                        trial['pmids'] = r.url.split('/')[-2]
                    else:
                        trial['citation'] = None
                else:
                    trial['citation'] = [{'url':citation_URL}]
            else:
                trial['citation'] = None

            if trial['data_available_for_request'] == 'true':
                trial['data_available_for_request'] = "Restricted"
            if trial['data_available_for_request'] == 'false':
                trial['data_available_for_request'] = "Closed"

            trial['conditionsOfAccess'] = trial.pop(
                'data_available_for_request')

            # de-duplication of identifier
            if trial['identifier'] != trial['nctid']:
                trial['identifier'] = [trial['identifier'], trial['nctid']]
            else:
                trial['identifier'] = [trial['nctid']]

            # unique _id appending identifier
            trial['_id'] = "accessclinicaldata_" + trial['identifier'][0]
            trial['includedInDataCatalog'] = {
                "name": "AccessClinicalData@NIAID"
            }
            trial['@type'] = "Dataset"
            trial['url'] = "https://accessclinicaldata.niaid.nih.gov/study-viewer/clinical_trials/" + \
                trial['identifier'][0]

            yield trial

            count += 1
            logger.info("Parsed %s records", count)
        except StopIteration:
            logger.info("Finished Parsing. Total Records: %s", count)
            # if StopIteration is raised, break from loop
            break

# parse()
# for i in parse():
#     pprint(i)
