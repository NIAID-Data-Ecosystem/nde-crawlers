import re
import datetime
import logging
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


COOKIE = 'user="2|1:0|10:1665435969|4:user|32:ImR5bGFud2VsemVsQGdtYWlsLmNvbSI=|a81f6d43704c37a1422da8da927b936505257850f1b5708364a98d163d4f3ad3"'


def retrieve_study_metdata():
    logger.info('Retrieving study metadata from Qiita')

    url = 'https://qiita.ucsd.edu/study/list_studies/?&user=dylanwelzel@gmail.com&visibility=public&sEcho=344&query=&_=1665437236850'
    headers = {'Cookie': COOKIE}
    r = requests.get(url, headers=headers)
    return r.json()['aaData']


def retrieve_study_samples(study_id):
    logger.info('Retrieving sample information from each individual study')

    url = f'https://qiita.ucsd.edu/study/description/sample_template/columns/?study_id={study_id}'
    headers = {'Cookie': COOKIE}
    r = requests.get(url, headers=headers)

    result = {}

    for column_name in r.json()['values']:
        url = f'https://qiita.ucsd.edu/study/description/sample_template/columns/?study_id={study_id}&column={column_name}'
        headers = {'Cookie': COOKIE}
        r = requests.get(url, headers=headers)

        result[column_name] = {}

        for value in r.json()['values']:
            if value in result[column_name]:
                result[column_name][value] += 1
            else:
                result[column_name][value] = 1

    return result


def parse():
    # Parse the metadata
    studies = retrieve_study_metdata()

    for study in studies:
        study_id = study['study_id']
        # study_samples = retrieve_study_samples(study_id)

        output = {'includedInDataCatalog': {
            '@type': 'Dataset',
            'name': 'Qiita',
            'url': 'https://qiita.ucsd.edu/',
            'versionDate': datetime.date.today().isoformat()
        },
            '@type': 'Dataset',
        }

        if study_abstract := study.get('study_abstract'):
            output['description'] = study_abstract

        if study_id := study.get('study_id'):
            output['url'] = f'https://qiita.ucsd.edu/study/description/{study_id}'

        if study_alias := study.get('study_alias'):
            output['_id'] = study_alias

        if study_title := study.get('study_title'):
            output['name'] = study_title

        if study_tags := study.get('study_tags'):
            output['keywords'] = study_tags

        # if owner := study.get('owner'):
        #     print(owner)

        if pi := study.get('pi'):
            name = re.findall(r'>(.*?)<', pi)
            output['author'] = {'name': name[0]}

        if pubs := study.get('pubs'):
            publications = re.findall(r'>(.*?)<', pubs)
            doi = []
            pmid = []
            for publication in publications:
                if '10.' in publication:
                    if 'http' in publication:
                        doi.append(publication.replace(
                            'https://doi.org/', ''))
                    else:
                        doi.append(publication)
                if len(publication) == 8:
                    pmid.append(publication)
            output['pmids'] = pmid
            output['doi'] = doi

        if ebi_study_accession := study.get('ebi_study_accession'):
            output['mainEntityOfPage'] = f'https://www.ebi.ac.uk/ena/browser/view/{ebi_study_accession}'

        yield output
