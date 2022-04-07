import json
import requests
import logging

from pprint import pprint


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def parse():
    # schema
    # https://docs.google.com/spreadsheets/d/1XlpvoeSWSiqfw1pPAHHFTjr9Wuowj3TA-eCsyGEj0ng/edit#gid=0

    url = "https://accessclinicaldata.niaid.nih.gov/guppy/graphql"
    query = """query ($filter: JSON) {
        clinical_trials (filter: $filter, first: 10000, accessibility: accessible) {
            title,cmc_unique_id,brief_summary,data_availability_date,most_recent_update,data_available,creator,nct_number,condition,clinical_trial_website,publications,,data_available_for_request
        }
    }"""

    r = requests.post(url, json={'query': query})

    # parse json string and convert to dictionary
    json_obj = json.loads(r.text)

    # grab list of trials from dictionary
    trials = json_obj['data']['clinical_trials']

    count = 0

    for trial in trials:

        count += 1
        logger.info("Parsed %s records", count)

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

        # TODO make citiation an object using helper function
        trial['citation'] = trial.pop('publications')

        if trial['data_available_for_request'] == 'true':
            trial['data_available_for_request'] = "restricted access"
        if trial['data_available_for_request'] == 'false':
            trial['data_available_for_request'] = "closed access"

        trial['conditionsOfAccess'] = trial.pop('data_available_for_request')

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


for i in parse():
    print(i)
