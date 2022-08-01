import datetime
import requests
import numpy as np
import json
import ftplib
import pandas as pd
from pysradb.sraweb import SRAweb
import logging
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def parse():
    logger.info('Starting FTP Download')
    ftp = ftplib.FTP(
        'ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('/sra/reports/Metadata')
    entries = list(ftp.mlsd())
    ftp_accessions = [
        entry for entry in entries if 'SRA_Accessions.tab' in entry[0]]
    ftp_accessions.sort(key=lambda entry: entry[1]['modify'], reverse=True)
    filename = ftp_accessions[0][0]
    with open(filename, "wb") as file:
        # use FTP's RETR command to download the file
        ftp.retrbinary(f"RETR {filename}", file.write)
    ftp.quit()
    logger.info('FTP Download Complete')

    logger.info('Retrieving Studies')
    db = SRAweb()
    df = pd.read_csv(r"SRA_Accessions.tab", sep="\t",
                     usecols=['Accession', 'Type', 'Status', 'Updated', 'Published', 'Experiment', 'Sample', 'BioProject', 'ReplacedBy'])
    only_live = df[df['Status'] == 'live']
    filtered = only_live[only_live['Type'] == 'STUDY']
    accession_list = filtered[['Accession', 'Type', 'Status', 'Updated', 'Published',
                               'Experiment', 'Sample', 'BioProject', 'ReplacedBy']].values.tolist()
    logger.info('Total Studies Found: {}'.format(len(accession_list)))
    count = 0

    # One at a time
    logger.info('Retrieving Individual Study Metadata')
    metadata_list = []
    for x in accession_list:
        try:
            meta_df = db.sra_metadata(x[0], detailed=True)
            meta_df = meta_df.replace({np.nan: None, '-': None})
            meta_dict = meta_df.to_dict(orient='list')
            meta_dict['Accession'] = x[0]
            drs_url = f'https://locate.be-md.ncbi.nlm.nih.gov/idx/v1/{x[0]}?submitted=true&etl=false'
            drs_response = requests.get(drs_url)
            drs_json = json.loads(drs_response.text)
            if 'response' in drs_json:
                drs_id = drs_json['response'][x[0]]['drs']
                meta_dict['DRS'] = drs_id
            if x[1] != '-':
                meta_dict['Type'] = x[1]
            if x[3] != '-':
                meta_dict['Updated'] = x[3]
            if x[4] != '-':
                meta_dict['Published'] = x[4]
            if x[5] != '-':
                meta_dict['Experiment'] = x[5]
            if x[6] != '-':
                meta_dict['Sample'] = x[6]
            if x[7] != '-':
                meta_dict['BioProject'] = x[7]
            if x[8] != '-':
                meta_dict['ReplacedBy'] = x[8]
            metadata_list.append(meta_dict)
            count += 1
            if count % 1000 == 0:
                break
            if count % 100 == 0:
                logger.info('{} Studies Retrieved'.format(count))
        except KeyError as e:
            print(x[0])
            print(e)
            continue

    logger.info('Parsing Individual Study Metadata')
    count = 0
    for metadata in metadata_list:
        count += 1
        if count % 100 == 0:
            logger.info('{} Studies Parsed'.format(count))
        if test := metadata.get('ena_fastq_ftp'):
            if test[0]:
                logger.info(metadata)
        if test2 := metadata.get('ena_fastq_http'):
            if test2[0]:
                logger.info(metadata)
        output = {
            '@context': "https://schema.org/",
            'includedInDataCatalog': {
                '@type': 'Dataset',
                'name': 'NCBI SRA',
                'url': 'https://www.ncbi.nlm.nih.gov/sra/',
                'versionDate': datetime.date.today().isoformat()
            }
        }
        # top level
        if accession := metadata.get('Accession'):
            output['_id'] = 'NCBI_SRA_' + accession
        if updated := metadata.get('Updated'):
            output['dateModified'] = datetime.datetime.strptime(
                updated, '%Y-%m-%dT%H:%M:%SZ').strftime('%Y-%m-%d')
        if published := metadata.get('Published'):
            output['datePublished'] = datetime.datetime.strptime(
                published, '%Y-%m-%dT%H:%M:%SZ').strftime('%Y-%m-%d')
        if recieved := metadata.get('Recieved'):
            output['dateCreated'] = datetime.datetime.strptime(
                recieved, '%Y-%m-%dT%H:%M:%SZ').strftime('%Y-%m-%d')
        if visibility := metadata.get('Visibility'):
            output['conditionsOfAccess'] = visibility
        if replaced_by := metadata.get('ReplacedBy'):
            output['sameAs'] = replaced_by
        if study_title := metadata.get('study_title'):
            output['name'] = study_title[0]
        if study_abstract := metadata.get('study_abstract'):
            output['description'] = study_abstract[0]
        if contact_name := metadata.get('contact_name'):
            output['author'] = {
                'name': contact_name[0],
            }
        # distribution
        distribution_list = []
        if sra_urls := metadata.get('sra_url'):
            for url in sra_urls:
                distribution_dict = {}
                distribution_dict['url'] = url
                if distribution_dict not in distribution_list:
                    distribution_list.append(distribution_dict)
        if gcp_urls := metadata.get('GCP_url'):
            for url in gcp_urls:
                distribution_dict = {}
                distribution_dict['url'] = url
                if distribution_dict not in distribution_list:
                    distribution_list.append(distribution_dict)
        if aws_urls := metadata.get('AWS_url'):
            for url in aws_urls:
                distribution_dict = {}
                distribution_dict['contentUrl'] = url
                if distribution_dict not in distribution_list:
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

        # isBasedOn
        is_based_on = []
        # runs
        if run_accessions := metadata.get('run_accession'):
            for run_accession in run_accessions:
                run_dict = {}
                run_dict['identifier'] = run_accession
                run_dict['additionalType'] = {
                    'name': 'Run', 'url': 'http://purl.obolibrary.org/obo/NCIT_C47911'}
                is_based_on.append(run_dict)
        # bioproject
        if bio_project := metadata.get('BioProject'):
            is_based_on.append(
                {'identifier': bio_project, 'additionalType': {'name': 'BioProject Accession Number', 'url': 'http://purl.obolibrary.org/obo/NCIT_C175890'}})

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
                        if experiment_dict not in is_based_on:
                            is_based_on.append(experiment_dict)
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
                        if sample_dict not in is_based_on:
                            is_based_on.append(sample_dict)
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
                        if instrument_dict not in is_based_on:
                            is_based_on.append(instrument_dict)
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
                            if cell_dict not in is_based_on:
                                is_based_on.append(cell_dict)
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
                        if hapmap_dict not in is_based_on:
                            is_based_on.append(hapmap_dict)

        if len(is_based_on):
            output['isBasedOn'] = is_based_on

        conditions_of_access = ''
        if gcp_free_egress := metadata.get('GCP_free_egress'):
            if conditions_of_access == '':
                conditions_of_access = gcp_free_egress
            else:
                conditions_of_access += f', {gcp_free_egress}'
        if gcp_access_type := metadata.get('GCP_access_type'):
            if conditions_of_access == '':
                conditions_of_access = gcp_access_type
            else:
                conditions_of_access += f', {gcp_access_type}'
        if aws_free_egress := metadata.get('AWS_free_egress'):
            if conditions_of_access == '':
                conditions_of_access = aws_free_egress
            else:
                conditions_of_access += f', {aws_free_egress}'
        if aws_access_type := metadata.get('AWS_access_type'):
            if conditions_of_access == '':
                conditions_of_access = aws_access_type
            else:
                conditions_of_access += f', {aws_access_type}'
        if conditions_of_access != '':
            output['conditionsOfAccess'] = conditions_of_access
        yield output

    logger.info('Removing SRA_Accessions.tab')
    os.remove("SRA_Accessions.tab")
    logger.info('Complete')
