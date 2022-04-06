import json
import requests
import logging

from pprint import pprint


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def parse():
    # dictionary to convert metadata property to schema property
    # https://docs.google.com/spreadsheets/d/1XlpvoeSWSiqfw1pPAHHFTjr9Wuowj3TA-eCsyGEj0ng/edit#gid=0

    url = "https://accessclinicaldata.niaid.nih.gov/guppy/graphql"
    query = """query ($filter: JSON) {
        clinical_trials (filter: $filter, first: 10000, accessibility: accessible) {
            title,cmc_unique_id,brief_summary,data_availability_date,most_recent_update,data_available,creator,nct_number,condition,clinical_trial_website,publications,,data_available_for_request
        }
    }"""

    r = requests.post(url, json={'query': query})
    json_obj = json.loads(r.text)
    trials = json_obj['data']['clinical_trials']

    for trial in trials:
        trial['name'] = trial.pop('title')
        trial['identifier'] = trial.pop('cmc_unique_id')
        trial['description'] = trial.pop('brief_summary')
        trial['datePublished'] = trial.pop('data_availability_date')
        trial['dateModified'] = trial.pop('most_recent_update')

        if trial['dateModified'] == None:
            trial['dateModified'] = "null"

        trial['additionalType'] = trial.pop('data_available')
        trial['funding'] = [{'funder': {'name': trial.pop('creator')}}]
        trial['nctid'] = trial.pop('nct_number')
        trial['infectiousDisease'] = trial.pop('condition')
        trial['mainEntityOfPage'] = trial.pop('clinical_trial_website')
        trial['citation'] = trial.pop('publications')

        if trial['data_available_for_request'] == None:
            trial['data_available_for_request'] = "null"
        if trial['data_available_for_request'] == 'true':
            trial['data_available_for_request'] = "restricted access"
        if trial['data_available_for_request'] == 'false':
            trial['data_available_for_request'] = "closed access"

        trial['conditionsOfAccess'] = trial.pop('data_available_for_request')

        if trial['identifier'] != trial['nctid']:
            trial['identifier'] = [trial['identifier'], trial['nctid']]
        else:
            trial['identifier'] = [trial['nctid']]

        trial['_id'] = "accessclinicaldata_" + trial['identifier'][0]
        trial['includedInDataCatalog'] = {
            "name": "AccessClinicalData@NIAID"
        }
        trial['@type'] = "Dataset"
        trial['url'] = "https://accessclinicaldata.niaid.nih.gov/study-viewer/clinical_trials/" + \
            trial['identifier'][0]
    pprint(trials)
    # return trial


parse()
