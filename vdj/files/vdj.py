import re
import datetime
import logging
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def retrieve_study_metadata():
    logger.info('Retrieving study metadata from VDJ')
    urls = ['https://vdj-staging.tacc.utexas.edu/airr/v1/repertoire', 'https://ipa3.ireceptor.org/airr/v1/repertoire',
            'https://vdjserver.org/airr/v1/repertoire', 'https://ipa4.ireceptor.org/airr/v1/repertoire', 'https://ipa1.ireceptor.org/airr/v1/repertoire', 'https://ipa2.ireceptor.org/airr/v1/repertoire', 'https://scireptor.dkfz.de/airr/v1/repertoire']

    studies = []
    for url in urls:
        r = requests.post(url, data='{}')
        data = r.json()['Repertoire']
        for dict in data:
            id = dict['study']['study_id']
            studies.append({id: dict})
    # print(metadata)
    # print(len(metadata))
    return studies


def parse():
    studies = retrieve_study_metadata()
    count = 0
    logger.info('Parsing study metadata')
    for study in studies:
        count += 1
        logger.info(f'Parsing study {count} of {len(studies)}')
        if count % 10 == 0:
            logger.info('Parsed %s studies', count)

        output = {'includedInDataCatalog': {
            '@type': 'Dataset',
            'name': 'VDJ',
            'url': 'https://vdj-staging.tacc.utexas.edu/community/',
            'versionDate': datetime.date.today().isoformat()
        },
            '@type': 'Dataset',
        }

        if study_abstract := study.get('study_description'):
            output['description'] = study_abstract

        if study_id := study.get('uuid'):
            # output['url'] = f'https://qiita.ucsd.edu/study/description/{study_id}'
            output['_id'] = f'qiita_{study_id}'

        if study_title := study.get('study_title'):
            output['name'] = study_title

        if study_tags := study.get('keywords_study'):
            output['keywords'] = study_tags

        # if owner := study.get('owner'):
        #     print(owner)

        if author := study.get('lab_name'):
            output['author'] = {'name': author}

        if pub_ids := study.get('pub_ids'):
            output['pmids'] = [pub_ids]

        # if ebi_study_accession := study.get('ebi_study_accession'):
        #     output['mainEntityOfPage'] = f'https://www.ebi.ac.uk/ena/browser/view/{ebi_study_accession}'

        yield output
