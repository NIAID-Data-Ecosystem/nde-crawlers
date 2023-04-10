import orjson
import os
import requests
import time
import json
import sqlite3
import logging
import datetime

from .date import add_date
from Bio import Entrez
from biothings.utils.dataload import tab2dict
from config import GEO_EMAIL, logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')

DB_PATH = '/data/nde-hub/standardizers/pubtator_lookup/pubtator_lookup.db'


MANUAL_HEALTH_CONDITIONS = [
    {
        "identifier": "D000086382",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D000086382",
        "originalName": "covid-19",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02"
        },
        "name": "COVID-19",
        "alternateName": [
            "2019 Novel Coronavirus Disease",
            "2019 Novel Coronavirus Infection",
            "2019-nCoV Disease",
            "2019-nCoV Infection",
            "COVID-19 Pandemic",
            "COVID-19 Pandemics",
            "COVID-19 Virus Disease",
            "COVID-19 Virus Infection",
            "COVID19",
            "Coronavirus Disease 2019",
            "Coronavirus Disease-19",
            "SARS Coronavirus 2 Infection",
            "SARS-CoV-2 Infection",
            "Severe Acute Respiratory Syndrome Coronavirus 2 Infection"
        ]
    },
    {
        "identifier": "D000086382",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D000086382",
        "originalName": "sars-cov-2",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02"
        },
        "name": "COVID-19",
        "alternateName": [
            "2019 Novel Coronavirus Disease",
            "2019 Novel Coronavirus Infection",
            "2019-nCoV Disease",
            "2019-nCoV Infection",
            "COVID-19 Pandemic",
            "COVID-19 Pandemics",
            "COVID-19 Virus Disease",
            "COVID-19 Virus Infection",
            "COVID19",
            "Coronavirus Disease 2019",
            "Coronavirus Disease-19",
            "SARS Coronavirus 2 Infection",
            "SARS-CoV-2 Infection",
            "Severe Acute Respiratory Syndrome Coronavirus 2 Infection"
        ]
    },
    {
        "identifier": "D006262",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D006262",
        "originalName": "healthy",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02"
        },
        "name": "Health",
        "alternateName": [
            "Individual Health",
            "Normalcy",
            "Normality"
        ]
    },
    {
        "identifier": "D000086002",
        "inDefinedTermSet": "MESH",
        "url": "http://id.nlm.nih.gov/mesh/D000086002",
        "originalName": "malignant mesothelioma",
        "isCurated": True,
        "curatedBy": {
            "name": "PubTator",
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
            "dateModified": "2023-04-02"
        },
        "name": "Mesothelioma, Malignant",
        "alternateName": [
            "Malignant Mesothelioma",
            "Malignant Pleural Mesothelioma",
            "Mesothelioma, Malignant Pleural"
        ]
    }
]

@add_date
def standardize_data(data):
    logger.info('Standardizing data...')
    if isinstance(data, str):
        logger.info('Reading data from file...')
        with open(os.path.join(data, 'data.ndjson'), 'rb') as f:
            health_conditions_list = []
            species_list = []
            infectious_agents_list = []
            doc_list = []
            count = 0
            for line in f:
                count += 1
                if count % 1000 == 0:
                    logger.info(f'Processed {count} lines')
                doc = orjson.loads(line)
                doc_list.append(doc)

                if isinstance(doc.get('healthCondition', {}), list):
                    for health_condition in doc['healthCondition']:
                        health_conditions_list.append(
                            health_condition['name'])
                elif health_condition_name := doc.get('healthCondition', {}).get('name'):
                    health_conditions_list.append(health_condition_name)

                if isinstance(doc.get('species', {}), list):
                    for species in doc['species']:
                        species_list.append(species['name'])
                elif species_name := doc.get('species', {}).get('name'):
                    species_list.append(species_name)

                if isinstance(doc.get('infectiousAgent', {}), list):
                    for infectious_agent in doc['infectiousAgent']:
                        infectious_agents_list.append(
                            infectious_agent['name'])
                elif infectious_agent_name := doc.get('infectiousAgent', {}).get('name'):
                    infectious_agents_list.append(infectious_agent_name)

            if len(health_conditions_list):
                health_conditions_list = list(
                    dict.fromkeys([x.lower().strip() for x in health_conditions_list]))
            if len(species_list):
                species_list = list(
                    dict.fromkeys([x.lower().strip() for x in species_list]))
            if len(infectious_agents_list):
                infectious_agents_list = list(
                    dict.fromkeys([x.lower().strip() for x in infectious_agents_list]))

            logger.info(
                f'Found {len(health_conditions_list)} health conditions, {len(species_list)} species, {len(infectious_agents_list)} infectious agents')
            return update_lookup_dict(health_conditions_list, species_list,
                                      infectious_agents_list, doc_list)

    else:
        logger.info('Reading data from list...')
        health_conditions_list = []
        species_list = []
        infectious_agents_list = []
        doc_list = list(data)
        count = 0
        for doc in doc_list:
            count += 1
            if count % 1000 == 0:
                logger.info(f'Processed {count} records')
            if isinstance(doc.get('healthCondition', {}), list):
                for health_condition in doc['healthCondition']:
                    health_conditions_list.append(health_condition['name'])
            elif health_condition_name := doc.get('healthCondition', {}).get('name'):
                health_conditions_list.append(health_condition_name)

            if isinstance(doc.get('species', {}), list):
                for species in doc['species']:
                    species_list.append(species['name'])
            elif species_name := doc.get('species', {}).get('name'):
                species_list.append(species_name)

            if isinstance(doc.get('infectiousAgent', {}), list):
                for infectious_agent in doc['infectiousAgent']:
                    infectious_agents_list.append(infectious_agent['name'])
            elif infectious_agent_name := doc.get('infectiousAgent', {}).get('name'):
                infectious_agents_list.append(infectious_agent_name)

        if len(health_conditions_list):
            health_conditions_list = list(
                dict.fromkeys([x.lower().strip() for x in health_conditions_list]))
        if len(species_list):
            species_list = list(
                dict.fromkeys([x.lower().strip() for x in species_list]))
        if len(infectious_agents_list):
            infectious_agents_list = list(
                dict.fromkeys([x.lower().strip() for x in infectious_agents_list]))

        logger.info(
            f'Found {len(health_conditions_list)} health conditions, {len(species_list)} species, {len(infectious_agents_list)} infectious agents')
        return update_lookup_dict(health_conditions_list, species_list,
                                  infectious_agents_list, doc_list)


def get_details(original_name, type, identifier):

    logger.info(f'Getting details for {original_name}')

    if type == 'Disease':
        if 'covid-19' in original_name.lower() or 'sars-cov-2' in original_name.lower():
            info = requests.get(
                'https://id.nlm.nih.gov/mesh/lookup/details?descriptor=D000086382')
            identifier = 'MESH:D000086382'
        elif 'malignant mesothelioma' in original_name.lower():
            info = requests.get(
                'https://id.nlm.nih.gov/mesh/lookup/details?descriptor=D000086002')
            identifier = 'MESH:D000086002'
        elif original_name.lower() == 'healthy':
            info = requests.get(
                'https://id.nlm.nih.gov/mesh/lookup/details?descriptor=D006262')
            identifier = 'MESH:D006262'
        elif 'asthma' in original_name.lower():
            info = requests.get(
                'https://id.nlm.nih.gov/mesh/lookup/details?descriptor=D004802')
            identifier = 'MESH:D004802'
        else:
            info = requests.get(
                f'https://id.nlm.nih.gov/mesh/lookup/details?descriptor={identifier.split(":")[1]}')
        try:
            info.raise_for_status()
            info = info.json()

            standard_dict = {
                'identifier': identifier.split(":")[1],
                'inDefinedTermSet': 'MESH',
                'url': info['descriptor'],
                'originalName': original_name,
                'isCurated': True,
                'curatedBy': {
                    'name': 'PubTator',
                    'url': 'https://www.ncbi.nlm.nih.gov/research/pubtator/api.html',
                    'dateModified': datetime.datetime.now().strftime('%Y-%m-%d')
                }
            }
            if len(info['terms']) > 0:
                standard_dict['name'] = info['terms'][0]['label']
            else:
                logger.info(f'No offical name for {original_name}')
                standard_dict['name'] = original_name

            if len(info['terms']) > 1:
                standard_dict['alternateName'] = [
                    x['label'] for x in info['terms'][1:]]

            return standard_dict
        except requests.exceptions.HTTPError as e:
            logger.info(
                f'No information found for {original_name}, {identifier}')
            return None

    elif type == 'Species':
        identifier = identifier.split('*')[-1]

        try:
            species_info = requests.get(
                f"https://rest.uniprot.org/taxonomy/{identifier}")
            species_info.raise_for_status()
            species_info = species_info.json()
            standard_dict = {
                'identifier': identifier,
                'inDefinedTermSet': 'UniProt',
                'url': f'https://www.uniprot.org/taxonomy/{identifier}',
                'originalName': original_name,
                'isCurated': True,
                'curatedBy': {
                    'name': 'PubTator',
                    'url': 'https://www.ncbi.nlm.nih.gov/research/pubtator/api.html',
                    'dateModified': datetime.datetime.now().strftime('%Y-%m-%d')

                }
            }
            if scientific_name := species_info.get('scientificName'):
                standard_dict['name'] = scientific_name
            else:
                standard_dict['name'] = original_name

            alternative_names = []
            if common_name := species_info.get('commonName'):
                standard_dict['commonName'] = common_name
                alternative_names.append(common_name)

            if other_names := species_info.get('otherNames'):
                alternative_names.extend(other_names)

            if alternative_names:
                standard_dict['alternateName'] = alternative_names

            return standard_dict

        except requests.exceptions.HTTPError as e:
            logger.info(
                f'No Uniprot information found for {original_name}, {identifier}. Trying NCBI...')

        Entrez.email = GEO_EMAIL

        while True:
            try:
                handle = Entrez.efetch(
                    db='taxonomy', id=identifier, retmode='xml', max_tries=10)
                record = Entrez.read(handle, validate=False)
                handle.close()
                break
            except Exception as e:
                logger.info(f'Error: {e}')
                logger.info(f'Retrying in 5 seconds...')
                time.sleep(5)
                get_details(original_name, type, identifier)

        standard_dict = {
            'identifier': identifier,
            'inDefinedTermSet': 'NCBI Taxonomy',
            'url': f'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={identifier}',
            'originalName': original_name,
            'isCurated': True,
            'curatedBy': {
                'name': 'PubTator',
                'url': 'https://www.ncbi.nlm.nih.gov/research/pubtator/api.html',
                'dateModified': datetime.datetime.now().strftime('%Y-%m-%d')
            }
        }
        if scientific_name := record[0].get('ScientificName'):
            standard_dict['name'] = scientific_name
        else:
            standard_dict['name'] = original_name

        alternative_names = []
        if common_name := record[0].get('OtherNames', {}).get('GenbankCommonName'):
            standard_dict['commonName'] = common_name
            alternative_names.append(common_name)

        if other_names := record[0].get('OtherNames', {}).get('Name'):
            for name_obj in other_names:
                if name_obj['ClassCDE'] == 'authority':
                    if name_obj['DispName'] not in alternative_names:
                        alternative_names.append(name_obj['DispName'])
        if alternative_names:
            standard_dict['alternateName'] = alternative_names

        return standard_dict


def get_new_health_conditions(health_conditions, sql_data):
    chunks = [health_conditions[x:x+1000]
              for x in range(0, len(health_conditions), 1000)]
    logger.info(f'Splitting into {len(chunks)} chunks')
    chunk_count = 0
    for chunk in chunks:
        chunk_count += 1
        logger.info(f'Processing chunk {chunk_count}')
        data = json.dumps('.    '.join(chunk))
        submit_response = requests.post(
            'https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/submit/disease', data=data)
        logger.info(f'Waiting for response, {submit_response.text}')
        timeout = 0
        retries = 0
        while True:
            time.sleep(10)
            timeout += 10
            if timeout > 100:
                retries += 1
                if retries > 3:
                    raise Exception('Attempted 3 times, giving up')
                logger.info('Timeout, retrying...')
                try:
                    os.remove(f'{submit_response.text}_response.csv')
                except FileNotFoundError:
                    logger.info('Issue removing file: %s_response.csv',
                                submit_response.text)
                get_new_health_conditions(health_conditions, sql_data)
            retrieve_response = requests.get(
                f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/retrieve/{submit_response.text}")
            if retrieve_response.status_code == 200:
                logger.info('Got response')
                break

        result = retrieve_response.text.split(
            '00000|a|-NoAbstract-')[1].strip()

        if len(result.split('\n')) == 1:
            result = result + '\n' + result

        with open(f'{submit_response.text}_response.csv', 'a') as f:
            for line in result.split('\n'):
                if 'MESH' in line:
                    f.write(line + '\n')
                else:
                    logger.info(f'Removed {line}')
        with open(f'{submit_response.text}_response.csv', 'r') as f:
            lines = f.readlines()

        with open(f'{submit_response.text}_response.csv', 'w') as f:
            for line in lines:
                if line.strip():
                    f.write(line)

            f.seek(0, 0)
            f.write('\n')

    if os.stat(f'{submit_response.text}_response.csv').st_size != 0:
        result = tab2dict(f'{submit_response.text}_response.csv', cols=[3, 4, 5],
                          key=0, sep='\t', alwayslist=True)

        remove_dupes = {}

        for key, value in result.items():
            if key.strip().lower() not in remove_dupes:
                remove_dupes[key.strip().lower()] = value

        for original_name in remove_dupes:

            if original_name not in [x[0].lower().strip() for x in sql_data]:

                result_dict = get_details(original_name, remove_dupes[original_name]
                                          [0][0], remove_dupes[original_name][0][1])
                conn = sqlite3.connect(DB_PATH)
                c = conn.cursor()
                c.execute("INSERT INTO health_conditions VALUES (?, ?)",
                          (original_name.lower(), json.dumps(result_dict)))
                logger.info(f'Added {original_name}')
                conn.commit()
                conn.close()

        difference = list(set([x.lower().strip()
                               for x in health_conditions]) - set(remove_dupes.keys()))
        logger.info(f'New health conditions: {difference}')
        try:
            os.remove(f'{submit_response.text}_response.csv')
        except FileNotFoundError:
            logger.info('Issue removing file: %s_response.csv',
                        submit_response.text)
        not_found = []
        for submitted_health_condition in health_conditions:
            if submitted_health_condition.lower() not in remove_dupes.keys():
                not_found.append(submitted_health_condition)

        return (difference, not_found)
    else:
        try:
            os.remove(f'{submit_response.text}_response.csv')
        except FileNotFoundError:
            logger.info('Issue removing file: %s_response.csv',
                        submit_response.text)

        logger.info("No new health conditions found")
        return None


def get_new_species(species, sql_data):
    chunks = [species[x:x+1000]
              for x in range(0, len(species), 1000)]
    logger.info(f'Splitting into {len(chunks)} chunks')
    chunk_count = 0
    for chunk in chunks:
        chunk_count += 1
        logger.info(f'Processing chunk {chunk_count}')
        data = json.dumps('.    '.join(chunk))
        submit_response = requests.post(
            'https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/submit/species', data=data)
        logger.info(f'Waiting for response, {submit_response.text}')
        timeout = 0
        retries = 0
        while True:
            time.sleep(10)
            timeout += 10
            if timeout > 100:
                retries += 1
                if retries > 3:
                    raise Exception('Attempted 3 times, giving up')
                logger.info('Timeout, retrying...')
                try:
                    os.remove(f'{submit_response.text}_response.csv')
                except FileNotFoundError:
                    logger.info('Issue removing file: %s_response.csv',
                                submit_response.text)
                get_new_species(species, sql_data)
            retrieve_response = requests.get(
                f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/annotations/annotate/retrieve/{submit_response.text}")
            if retrieve_response.status_code == 200:
                logger.info('Got response')
                break

        result = retrieve_response.text.split(
            '00000|a|-NoAbstract-')[1].strip()

        if len(result.split('\n')) == 1:
            result = result + '\n' + result

        with open(f'{submit_response.text}_response.csv', 'a') as f:
            for line in result.split('\n'):
                if 'Species' in line:
                    f.write(line + '\n')
                else:
                    logger.info(f'Removed {line}')

        with open(f'{submit_response.text}_response.csv', 'r') as f:
            lines = f.readlines()

        with open(f'{submit_response.text}_response.csv', 'w') as f:
            for line in lines:
                if line.strip():
                    f.write(line)

            f.seek(0, 0)
            f.write('\n')

    if os.stat(f'{submit_response.text}_response.csv').st_size != 0:
        result = tab2dict(f'{submit_response.text}_response.csv', cols=[3, 4, 5],
                          key=0, sep='\t', alwayslist=True)

        remove_dupes = {}

        for key, value in result.items():
            if key.strip().lower() not in remove_dupes:
                remove_dupes[key.strip().lower()] = value
        for original_name in remove_dupes:
            if original_name not in [x[0].lower().strip() for x in sql_data]:

                result_dict = get_details(original_name, remove_dupes[original_name]
                                          [0][0], remove_dupes[original_name][0][1])
                conn = sqlite3.connect(DB_PATH)
                c = conn.cursor()
                c.execute("INSERT INTO species VALUES (?, ?)",
                          (original_name.lower(), json.dumps(result_dict)))
                logger.info(f'Added {original_name}')
                conn.commit()
                conn.close()

        difference = list(set([x.lower().strip()
                               for x in species]) - set(remove_dupes.keys()))
        try:
            os.remove(f'{submit_response.text}_response.csv')
        except FileNotFoundError:
            logger.info('Issue removing file: %s_response.csv',
                        submit_response.text)
        logger.info(f'New species: {difference}')

        not_found = []
        for submitted_species in species:
            if submitted_species.lower() not in remove_dupes.keys():
                not_found.append(submitted_species)

        return (difference, not_found)
    else:
        try:
            os.remove(f'{submit_response.text}_response.csv')
        except FileNotFoundError:
            logger.info('Issue removing file: %s_response.csv',
                        submit_response.text)
        logger.info("No new species found")
        return None


def transform(doc_list):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT * FROM health_conditions")
    health_conditions_data = c.fetchall()
    c.execute("SELECT * FROM species")
    species_data = c.fetchall()
    conn.close()

    for doc in doc_list:
        new_health_conditions_list = []
        if isinstance(doc.get('healthCondition', {}), list):
            health_conditions_count = 0
            for original_health_condition_obj in doc['healthCondition']:
                health_conditions_count += 1
                if original_health_condition_name := original_health_condition_obj.get('name'):
                    for health_condition in health_conditions_data:
                        if original_health_condition_name.lower().strip() == health_condition[0].lower().strip():
                            logger.info(
                                f'Found {original_health_condition_name} in lookup dictionary')
                            if health_condition[1] is not None:
                                new_health_conditions_list.append(
                                    json.loads(health_condition[1]))
                                break
                        if health_condition[1] is not None and 'name' in health_condition[1]:
                            scientific_name = json.loads(
                                health_condition[1])['name']
                            if original_health_condition_name.lower().strip() == scientific_name.lower().strip():
                                logger.info(
                                    f'Found {original_health_condition_name} in lookup dictionary through the scientific name of {health_condition[0]}')
                                new_health_conditions_list.append(
                                    json.loads(health_condition[1]))
                                break
                        if health_condition[1] is not None and 'alternateName' in health_condition[1]:
                            found_alternate_name = False
                            for alternate_name in json.loads(health_condition[1])['alternateName']:
                                if original_health_condition_name.lower().strip() == alternate_name.lower().strip():
                                    logger.info(
                                        f'Found {original_health_condition_name} in lookup dictionary through an alternate name of {health_condition[0]}')
                                    new_health_conditions_list.append(
                                        json.loads(health_condition[1]))
                                    found_alternate_name = True
                                    break
                            if found_alternate_name:
                                break

                    if len(new_health_conditions_list) < health_conditions_count:
                        logger.info(
                            f'No information found for {original_health_condition_name}')
                        new_health_conditions_list.append(
                            original_health_condition_obj)
        elif original_health_condition_name := doc.get('healthCondition', {}).get('name'):
            for health_condition in health_conditions_data:
                if original_health_condition_name.lower().strip() == health_condition[0].lower().strip():
                    logger.info(
                        f'Found {original_health_condition_name} in lookup dictionary')
                    if health_condition[1] is not None:
                        new_health_conditions_list.append(
                            json.loads(health_condition[1]))
                        break
                if health_condition[1] is not None and 'name' in health_condition[1]:
                    scientific_name = json.loads(
                        health_condition[1])['name']
                    if original_health_condition_name.lower().strip() == scientific_name.lower().strip():
                        logger.info(
                            f'Found {original_health_condition_name} in lookup dictionary through the scientific name of {health_condition[0]}')
                        new_health_conditions_list.append(
                            json.loads(health_condition[1]))
                        break
                if health_condition[1] is not None and 'alternateName' in health_condition[1]:
                    found_alternate_name = False
                    for alternate_name in json.loads(health_condition[1])['alternateName']:
                        if original_health_condition_name.lower().strip() == alternate_name.lower().strip():
                            logger.info(
                                f'Found {original_health_condition_name} in lookup dictionary through an alternate name of {health_condition[0]}')
                            new_health_conditions_list.append(
                                json.loads(health_condition[1]))
                            found_alternate_name = True
                            break
                    if found_alternate_name:
                        break
            if not len(new_health_conditions_list):
                logger.info(
                    f'No information found for {original_health_condition_name}')
                new_health_conditions_list.append(doc['healthCondition'])

        new_species_list = []
        if isinstance(doc.get('species', {}), list):
            species_count = 0
            for original_species_obj in doc['species']:
                species_count += 1
                if original_species_name := original_species_obj.get('name'):
                    for species in species_data:
                        if original_species_name.lower().strip() == species[0].lower().strip():
                            if species[1] is not None:
                                logger.info(
                                    f'Found {original_species_name} in lookup dictionary')
                                new_species_list.append(
                                    json.loads(species[1]))
                                break
                        if species[1] is not None and 'name' in species[1]:
                            scientific_name = json.loads(
                                species[1])['name']
                            if original_species_name.lower().strip() == scientific_name.lower().strip():
                                logger.info(
                                    f'Found {original_species_name} in lookup dictionary through the scientific name of {species[0]}')
                                new_species_list.append(
                                    json.loads(species[1]))
                                break
                        if species[1] is not None and 'alternateName' in species[1]:
                            found_alternate_name = False
                            for alternate_name in json.loads(species[1])['alternateName']:
                                if original_species_name.lower().strip() == alternate_name.lower().strip():
                                    logger.info(
                                        f'Found {original_species_name} in lookup dictionary through an alternate name of {species[0]}')
                                    new_species_list.append(
                                        json.loads(species[1]))
                                    found_alternate_name = True
                                    break
                            if found_alternate_name:
                                break
                    if len(new_species_list) < species_count:
                        logger.info(
                            f'No information found for {original_species_name}')
                        new_species_list.append(
                            original_species_obj)

        elif original_species_name := doc.get('species', {}).get('name'):
            for species in species_data:
                if original_species_name.lower().strip() == species[0].lower().strip():
                    logger.info(
                        f'Found {original_species_name} in lookup dictionary')
                    if species[1] is not None:
                        new_species_list.append(
                            json.loads(species[1]))
                        break
                if species[1] is not None and 'name' in species[1]:
                    scientific_name = json.loads(
                        species[1])['name']
                    if original_species_name.lower().strip() == scientific_name.lower().strip():
                        logger.info(
                            f'Found {original_species_name} in lookup dictionary through the scientific name of {species[0]}')
                        new_species_list.append(
                            json.loads(species[1]))
                        break
                if species[1] is not None and 'alternateName' in species[1]:
                    found_alternate_name = False
                    for alternate_name in json.loads(species[1])['alternateName']:
                        if original_species_name.lower().strip() == alternate_name.lower().strip():
                            logger.info(
                                f'Found {original_species_name} in lookup dictionary through an alternate name of {species[0]}')
                            new_species_list.append(
                                json.loads(species[1]))
                            found_alternate_name = True
                            break
                    if found_alternate_name:
                        break
            if not len(new_species_list):
                logger.info(
                    f'No information found for {original_species_name}')
                new_species_list.append(doc['species'])

        new_infectious_agent_list = []
        if isinstance(doc.get('infectiousAgent', {}), list):
            species_count = 0
            for original_species_obj in doc['infectiousAgent']:
                species_count += 1
                if original_species_name := original_species_obj.get('name'):
                    for species in species_data:
                        if original_species_name.lower().strip() == species[0].lower().strip():
                            logger.info(
                                f'Found {original_species_name} in lookup dictionary')
                            if species[1] is not None:
                                new_infectious_agent_list.append(
                                    json.loads(species[1]))
                                break
                        if species[1] is not None and 'name' in species[1]:
                            scientific_name = json.loads(
                                species[1])['name']
                            if original_species_name.lower().strip() == scientific_name.lower().strip():
                                logger.info(
                                    f'Found {original_species_name} in lookup dictionary through the scientific name of {species[0]}')
                                new_infectious_agent_list.append(
                                    json.loads(species[1]))
                                break
                        if species[1] is not None and 'alternateName' in species[1]:
                            for alternate_name in json.loads(species[1])['alternateName']:
                                found_alternate_name = False
                                if original_species_name.lower().strip() == alternate_name.lower().strip():
                                    logger.info(
                                        f'Found {original_species_name} in lookup dictionary through an alternate name of {species[0]}')
                                    new_infectious_agent_list.append(
                                        json.loads(species[1]))
                                    found_alternate_name = True
                                    break
                            if found_alternate_name:
                                break
                    if len(new_infectious_agent_list) < species_count:
                        logger.info(
                            f'No information found for {original_species_name}')
                        new_infectious_agent_list.append(
                            original_species_obj)

        elif original_species_name := doc.get('infectiousAgent', {}).get('name'):
            for species in species_data:
                if original_species_name.lower().strip() == species[0].lower().strip():
                    logger.info(
                        f'Found {original_species_name} in lookup dictionary')
                    if species[1] is not None:
                        new_infectious_agent_list.append(
                            json.loads(species[1]))
                        break
                if species[1] is not None and 'name' in species[1]:
                    scientific_name = json.loads(
                        species[1])['name']
                    if original_species_name.lower().strip() == scientific_name.lower().strip():
                        logger.info(
                            f'Found {original_species_name} in lookup dictionary through the scientific name of {species[0]}')
                        new_infectious_agent_list.append(
                            json.loads(species[1]))
                        break
                if species[1] is not None and 'alternateName' in species[1]:
                    found_alternate_name = False
                    for alternate_name in json.loads(species[1])['alternateName']:
                        if original_species_name.lower().strip() == alternate_name.lower().strip():
                            logger.info(
                                f'Found {original_species_name} in lookup dictionary through an alternate name of {species[0]}')
                            new_infectious_agent_list.append(
                                json.loads(species[1]))
                            found_alternate_name = True
                            break
                    if found_alternate_name:
                        break
            if not len(new_infectious_agent_list):
                logger.info(
                    f'No information found for {original_species_name}')
                new_infectious_agent_list.append(doc['infectiousAgent'])

        if new_health_conditions_list:
            doc['healthCondition'] = new_health_conditions_list
        if new_species_list:
            doc['species'] = new_species_list
        if new_infectious_agent_list:
            doc['infectiousAgent'] = new_infectious_agent_list
        yield doc


def update_lookup_dict(health_conditions_list, species_list,
                       infectious_agents_list, doc_list):

    if not os.path.exists(DB_PATH) or os.stat(DB_PATH).st_size == 0:
        logger.info('No lookup dictionary found, creating new one')
        os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute('''CREATE TABLE health_conditions
                        (original_name text, standard_dict text)''')
        for health_condition in MANUAL_HEALTH_CONDITIONS:
            c.execute(
                f"INSERT INTO health_conditions VALUES ('{health_condition['originalName']}', '{json.dumps(health_condition)}')")
            conn.commit()

        c.execute('''CREATE TABLE species
                        (original_name text, standard_dict text)''')

    if len(health_conditions_list):
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("SELECT * FROM health_conditions")
        data = c.fetchall()
        conn.close()
        no_matches_health_conditions = []
        if len(health_conditions_list):
            for original_name in health_conditions_list:
                found_match = lookup(original_name, data)
                if found_match is False and original_name not in no_matches_health_conditions:
                    logger.info(
                        f'No {original_name} in lookup dictionary, saving for batch...')
                    no_matches_health_conditions.append(original_name)
        if no_matches_health_conditions:
            logger.info('Getting new health conditions')
            no_results = get_new_health_conditions(
                no_matches_health_conditions, data)
            if no_results == None:
                logger.info(
                    f'No results for {no_matches_health_conditions}, adding to lookup dictionary with no detailed information')
                conn = sqlite3.connect(DB_PATH)
                c = conn.cursor()
                for item in no_matches_health_conditions:
                    c.execute("INSERT INTO health_conditions VALUES (?, ?)",
                              (item.lower(), None))
                conn.commit()
                conn.close()
            elif len(no_results[1]):
                for item in no_results[1]:
                    logger.info(
                        f'No results for {item}, adding to lookup dictionary with no detailed information')
                    conn = sqlite3.connect(DB_PATH)
                    c = conn.cursor()
                    c.execute("INSERT INTO health_conditions VALUES (?, ?)",
                              (item.lower(), None))
                    conn.commit()
                    conn.close()

    if len(species_list) or len(infectious_agents_list):
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("SELECT * FROM species")
        data = c.fetchall()
        conn.close()
        no_matches_species = []
        if len(species_list):
            for original_name in species_list:
                found_match = lookup(original_name, data)
                if found_match is False and original_name not in no_matches_species:
                    logger.info(
                        f'No {original_name} in lookup dictionary, saving for batch...')
                    no_matches_species.append(original_name)
        if len(infectious_agents_list):
            for original_name in infectious_agents_list:
                found_match = lookup(original_name, data)
                if found_match is False and original_name not in no_matches_species:
                    logger.info(
                        f'No {original_name} in lookup dictionary, saving for batch...')
                    no_matches_species.append(original_name)
        if no_matches_species:
            logger.info('Getting new species')
            no_results = get_new_species(
                no_matches_species, data)
            if no_results == None:
                logger.info(
                    f'No results for {no_matches_species}, adding to lookup dictionary with no detailed information')
                conn = sqlite3.connect(DB_PATH)
                c = conn.cursor()
                for item in no_matches_species:
                    c.execute("INSERT INTO species VALUES (?, ?)",
                              (item.lower(), None))
                conn.commit()
                conn.close()
            elif len(no_results[1]):
                for item in no_results[1]:
                    logger.info(
                        f'No results for {item}, adding to lookup dictionary with no detailed information')
                    conn = sqlite3.connect(DB_PATH)
                    c = conn.cursor()
                    c.execute("INSERT INTO species VALUES (?, ?)",
                              (item.lower(), None))
                    conn.commit()
                    conn.close()

    logger.info('Finished updating lookup dictionary')
    logger.info('Standardizing document')
    transformed_docs = transform(doc_list)
    return transformed_docs


def lookup(original_name, data):
    for item in data:
        if original_name.lower().strip() == item[0].lower().strip():
            if item[1] is not None:
                logger.info(f'Found {original_name} in lookup dictionary')
                return True
            else:
                logger.info(
                    f'{original_name} exists in lookup dictionary, but no detailed information')
                return True
        elif item[1] is not None and 'alternateName' in item[1]:
            for alternate_name in json.loads(item[1])['alternateName']:
                if original_name.lower().strip() == alternate_name.lower().strip():
                    logger.info(
                        f'Found {original_name} in alternate names of {item[0]}')
                    return True
    return False
