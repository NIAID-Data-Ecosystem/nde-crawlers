import time
import os
import logging
import pandas as pd
import json
import datetime
import wget

from kingfisher import annotate
from sql_database import NDEDatabase
from concurrent.futures import ThreadPoolExecutor

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


class NCBI_SRA(NDEDatabase):
    # override variables
    SQL_DB = "ncbi_sra.db"
    EXPIRE = datetime.timedelta(days=90)
    NO_CACHE = True

    # Used for testing small chunks of data
    DATA_LIMIT = None

    # API
    def query_sra(self, study_info):

        study_acc = study_info[0]
        if study_acc == '-':
            return (study_acc, json.dumps(None))
        try:
            annotate(run_identifiers=None, run_identifiers_file=None,
                     bioproject_accession=study_acc, output_file=f'{study_acc}.json', output_format='json', all_columns=True)
            with open(f'{study_acc}.json') as f:
                data = json.load(f)
                ftp_info = {}
                if study_info[0] != '-':
                    ftp_info['ftp_acc'] = study_info[0]
                if study_info[1] != '-':
                    ftp_info['ftp_type'] = study_info[1]
                if study_info[3] != '-':
                    ftp_info['ftp_updated'] = study_info[3]
                if study_info[4] != '-':
                    ftp_info['ftp_published'] = study_info[4]
                if study_info[5] != '-':
                    ftp_info['ftp_experiment'] = study_info[5]
                if study_info[6] != '-':
                    ftp_info['ftp_sample'] = study_info[6]
                if study_info[7] != '-':
                    ftp_info['ftp_bioProject'] = study_info[7]
                if study_info[8] != '-':
                    ftp_info['ftp_replacedBy'] = study_info[8]
                data.append(ftp_info)
                os.remove(f'{study_acc}.json')
                return (study_acc, json.dumps(data))
        except KeyError as e:
            logger.error(f'KeyError for {study_acc}: {e}')
            return (study_acc, json.dumps(None))
        except Exception as e:
            if "HTTP Failure" in str(e):
                logger.error(f'HTTP Failure for {study_acc}: {e}')
                return (study_acc, json.dumps(None))
            else:
                logger.error(f'Unknown Error for {study_acc}: {e}')
                return (study_acc, json.dumps(None))




    def load_cache(self):
        logger.info('Starting FTP Download')

        fileloc = 'https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab'
        wget.download(fileloc, out='SRA_Accessions.tab')

        logger.info('FTP Download Complete')

        logger.info('Retrieving Studies from SRA_Accessions.tab')

        df = pd.read_csv(r"SRA_Accessions.tab", sep="\t",
                         usecols=['Accession', 'Type', 'Status', 'Updated', 'Published', 'Experiment', 'Sample', 'BioProject', 'ReplacedBy'])

        only_live = df[df['Status'] == 'live']
        filtered = only_live[only_live['Type'] == 'STUDY']
        accession_list = filtered[['Accession', 'Type', 'Status', 'Updated', 'Published',
                                   'Experiment', 'Sample', 'BioProject', 'ReplacedBy']].values.tolist()

        # Used for testing small chunks of data
        if self.DATA_LIMIT:
            accession_list = accession_list[:self.DATA_LIMIT]

        logger.info('Total Studies Found: {}'.format(len(accession_list)))
        count = 0

        logger.info('Retrieving Individual Study Metadata from API')

        start = time.time()
        with ThreadPoolExecutor(max_workers=3) as pool:
            data = pool.map(self.query_sra, accession_list)
            for item in data:
                yield item
                count += 1
                if count % 1000 == 0:
                    end = time.time()
                    logger.info('Retrieved {} Studies in {} seconds'.format(
                        count, end - start))
                    start = time.time()

        logger.info('Removing SRA_Accessions.tab')
        os.remove("SRA_Accessions.tab")
        logger.info('Removed SRA_Accessions.tab')

    def parse(self, studies):
        logger.info('Parsing Individual Study Metadata')

        count = 0

        for study in studies:

            count += 1
            if count % 1000 == 0:
                logger.info('{} Studies Parsed'.format(count))

            metadata = json.loads(study[1])

            if metadata is None:
                logger.info(f'No Metadata for {study[0]}')
                continue
            if len(metadata) == 0:
                logger.info(f'No Metadata for {study[0]}')
                continue

            # get last item in list, which is the ftp info
            ftp_info = metadata[-1]
            output = {
                '@context': "https://schema.org/",
                'includedInDataCatalog': {
                    '@type': 'Dataset',
                    'name': 'NCBI SRA',
                    'url': 'https://www.ncbi.nlm.nih.gov/sra/',
                    'versionDate': datetime.date.today().isoformat()
                },
                '@type': 'Dataset',
            }
            # Top Level Study Data
            if accession := ftp_info.get('ftp_acc'):
                output['_id'] = 'NCBI_SRA_' + accession
                output['url'] = 'https://www.ncbi.nlm.nih.gov/sra/' + accession
            if updated := ftp_info.get('ftp_updated'):
                output['dateModified'] = datetime.datetime.strptime(
                    updated, '%Y-%m-%dT%H:%M:%SZ').strftime('%Y-%m-%d')
            if published := ftp_info.get('ftp_published'):
                output['datePublished'] = datetime.datetime.strptime(
                    published, '%Y-%m-%dT%H:%M:%SZ').strftime('%Y-%m-%d')
            if replaced_by := ftp_info.get('ftp_replacedBy'):
                output['sameAs'] = replaced_by

            for run in study:

                if study_title := run.get('study_title'):
                    output['name'] = study_title[0]
                if study_abstract := metadata.get('study_abstract'):
                    output['description'] = study_abstract[0]
                if contact_name := metadata.get('contact_name'):
                    output['author'] = {
                        'name': contact_name[0],
                    }

                # distribution
                distribution_list = []
                if gcp_urls := metadata.get('GCP_url'):
                    if gcp_free_egress := metadata.get('GCP_free_egress'):
                        for url, string in zip(gcp_urls, gcp_free_egress):
                            distribution_dict = {}
                            if string is not None:
                                output['isAccessibleForFree'] = True
                            if url is not None:
                                distribution_dict['contentUrl'] = url
                            if bool(distribution_dict) and distribution_dict not in distribution_list:
                                distribution_list.append(distribution_dict)
                if aws_urls := metadata.get('AWS_url'):
                    if aws_free_egress := metadata.get('AWS_free_egress'):
                        for url, string in zip(aws_urls, aws_free_egress):
                            distribution_dict = {}
                            if string is not None:
                                output['isAccessibleForFree'] = True
                            if url is not None:
                                distribution_dict['contentUrl'] = url
                            if bool(distribution_dict) and distribution_dict not in distribution_list:
                                distribution_list.append(distribution_dict)
                if len(distribution_list):
                    output['distribution'] = distribution_list

                # species
                species_list = []
                organism_taxids = metadata.get('organism_taxid')
                organism_names = metadata.get('organism_name')
                if organism_taxids and organism_names:
                    for taxid, name in zip(organism_taxids, organism_names):
                        species_dict = {}
                        species_dict['name'] = name
                        species_dict['identifier'] = taxid
                        species_dict['additionalType'] = {
                            'name': 'Species', 'url': 'http://purl.obolibrary.org/obo/NCIT_C45293'}
                        if species_dict not in species_list:
                            species_list.append(species_dict)
                elif organism_taxids:
                    for taxid in organism_taxids:
                        species_dict = {}
                        species_dict['identifier'] = taxid
                        species_dict['additionalType'] = {
                            'name': 'Species', 'url': 'http://purl.obolibrary.org/obo/NCIT_C45293'}
                        if species_dict not in species_list:
                            species_list.append(species_dict)
                elif organism_names:
                    for name in organism_names:
                        species_dict = {}
                        species_dict['name'] = name
                        species_dict['additionalType'] = {
                            'name': 'Species', 'url': 'http://purl.obolibrary.org/obo/NCIT_C45293'}
                        if species_dict not in species_list:
                            species_list.append(species_dict)
                output['species'] = species_list

                # isBasedOn
                is_based_on = []
                # runs
                if run_accessions := metadata.get('run_accession'):
                    for run_accession in run_accessions:
                        run_dict = {}
                        run_dict['identifier'] = run_accession
                        run_dict['additionalType'] = {
                            'name': 'Run', 'url': 'http://purl.obolibrary.org/obo/NCIT_C47911'}
                        run_dict_filtered = {
                            k: v for k, v in run_dict.items() if v is not None}

                        if run_dict_filtered not in is_based_on:
                            is_based_on.append(run_dict_filtered)
                # bioproject
                if bio_project := metadata.get('BioProject'):
                    bio_project_dict = {}
                    bio_project_dict['identifier'] = bio_project
                    bio_project_dict['additionalType'] = {
                        'name': 'BioProject', 'url': 'http://purl.obolibrary.org/obo/NCIT_C45293'}
                    bio_project_dict_filtered = {
                        k: v for k, v in bio_project_dict.items() if v is not None}

                    if bio_project_dict_filtered not in is_based_on:
                        is_based_on.append(bio_project_dict_filtered)

                # experiments
                if experiment_accessions := metadata.get('experiment_accession'):
                    if experiment_title := metadata.get('experiment_title'):
                        if experiment_desc := metadata.get('experiment_desc'):
                            for experiment_accession, experiment_title, experiment_desc in zip(experiment_accessions, experiment_title, experiment_desc):
                                experiment_dict = {}
                                experiment_dict['identifier'] = experiment_accession
                                experiment_dict['name'] = experiment_title
                                experiment_dict['description'] = experiment_desc
                                experiment_dict['additionalType'] = {
                                    'name': 'Experiment', 'url': 'http://purl.obolibrary.org/obo/NCIT_C42790'}
                                experiment_dict_filtered = {
                                    k: v for k, v in experiment_dict.items() if v is not None}

                                if experiment_dict_filtered not in is_based_on:
                                    is_based_on.append(
                                        experiment_dict_filtered)
                # samples
                if sample_accessions := metadata.get('sample_accession'):
                    if sample_title := metadata.get('sample_title'):
                        if sample_comment := metadata.get('sample comment'):
                            for sample_accession, sample_title, sample_comment in zip(sample_accessions, sample_title, sample_comment):
                                sample_dict = {}
                                sample_dict['identifier'] = sample_accession
                                sample_dict['name'] = sample_title
                                sample_dict['description'] = sample_comment
                                sample_dict['additionalType'] = {
                                    'name': 'Sample', 'url': 'http://purl.obolibrary.org/obo/NCIT_C70699'}
                                sample_dict_filtered = {
                                    k: v for k, v in sample_dict.items() if v is not None}

                                if sample_dict_filtered not in is_based_on:
                                    is_based_on.append(sample_dict_filtered)
                # instruments
                if instruments := metadata.get('instrument'):
                    if instrument_models := metadata.get('instrument_model'):
                        if instrument_model_descriptions := metadata.get('instrument_model_desc'):
                            for instrument, instrument_model, instrument_model_desc in zip(instruments, instrument_models, instrument_model_descriptions):
                                instrument_dict = {}
                                instrument_dict['name'] = instrument
                                if instrument != instrument_model:
                                    instrument_dict['identifier'] = instrument_model
                                instrument_dict['description'] = instrument_model_desc
                                instrument_dict['additionalType'] = {
                                    'name': 'Instrument', 'url': 'http://purl.obolibrary.org/obo/NCIT_C16742'}
                                instrument_dict_filtered = {
                                    k: v for k, v in instrument_dict.items() if v is not None}

                                if instrument_dict_filtered not in is_based_on:
                                    is_based_on.append(
                                        instrument_dict_filtered)
                # cells
                if cell_lines := metadata.get('cell line'):
                    if cell_strain := metadata.get('strain'):
                        if cell_type := metadata.get('cell type'):
                            if cell_line_name := metadata.get('cell line name'):
                                for cell_line, cell_strain, cell_type, cell_line_name in zip(cell_lines, cell_strain, cell_type, cell_line_name):
                                    cell_dict = {}
                                    cell_dict['name'] = cell_line
                                    cell_dict['name'] = cell_line_name
                                    cell_dict['identifier'] = cell_strain
                                    cell_dict['additionalType'] = {
                                        'name': 'Cell', 'url': 'http://purl.obolibrary.org/obo/NCIT_C12508'}
                                cell_dict_filtered = {
                                    k: v for k, v in cell_dict.items() if v is not None}

                                if cell_dict_filtered not in is_based_on:
                                    is_based_on.append(cell_dict_filtered)
                # hapmap
                if hapmap_ids := metadata.get('HapMap sample ID'):
                    if cell_lines := metadata.get('Cell line'):
                        if sexes := metadata.get('sex'):
                            for hapmap_id, cell_line, sex in zip(hapmap_ids, cell_lines, sexes):
                                hapmap_dict = {}
                                hapmap_dict['identifier'] = hapmap_id
                                hapmap_dict['name'] = cell_line
                                hapmap_dict['gender'] = sex
                                hapmap_dict['additionalType'] = {
                                    'name': 'HapMap', 'url': 'http://purl.obolibrary.org/obo/NCIT_C70979'}
                                hapmap_dict_filtered = {
                                    k: v for k, v in hapmap_dict.items() if v is not None}

                                if hapmap_dict_filtered not in is_based_on:
                                    is_based_on.append(hapmap_dict_filtered)

                if len(is_based_on):
                    output['isBasedOn'] = is_based_on

                yield output

        logger.info(f'Finished Parsing {count} Studies')

    def update_cache(self):
        start = time.time()
        last_updated = self.retreive_last_updated()
        # logger.info('Starting FTP Download')

        fileloc = 'https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab'
        wget.download(fileloc, out='SRA_Accessions.tab')

        # logger.info('FTP Download Complete')

        logger.info('Retrieving Studies from SRA_Accessions.tab')

        df = pd.read_csv(r"SRA_Accessions.tab", sep="\t",
                         usecols=['Accession', 'Type', 'Status', 'Updated', 'Published', 'Experiment', 'Sample', 'BioProject', 'ReplacedBy'])
        only_live = df[df['Status'] == 'live']
        filtered = only_live[only_live['Type'] == 'STUDY']
        new_studies = filtered[filtered['Updated'] > last_updated]
        accession_list = new_studies[['Accession', 'Type', 'Status', 'Updated', 'Published',
                                      'Experiment', 'Sample', 'BioProject', 'ReplacedBy']].values.tolist()

        # Used for testing small chunks of data
        if self.DATA_LIMIT:
            accession_list = accession_list[:self.DATA_LIMIT]

        logger.info('Total Studies Found: {}'.format(len(accession_list)))
        count = 0

        logger.info('Retrieving Individual Study Metadata from API')

        with ThreadPoolExecutor(max_workers=3) as pool:
            data = pool.map(self.query_sra, accession_list)
            for item in data:
                yield item
                count += 1
                if count % 1000 == 0:
                    end = time.time()
                    logger.info('Retrieved {} Studies in {} seconds'.format(
                        count, end - start))
                    start = time.time()

        logger.info('Removing SRA_Accessions.tab')
        os.remove("SRA_Accessions.tab")
        logger.info('Removed SRA_Accessions.tab')
        self.insert_last_updated()
