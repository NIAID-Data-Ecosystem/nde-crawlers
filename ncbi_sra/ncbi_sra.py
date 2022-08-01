import datetime
import requests
import numpy as np
import time
import json
import ftplib
import pandas as pd
import tarfile
from pysradb.sraweb import SRAweb
from pprint import pprint
from pysradb.sradb import SRAdb
import sqlite3


# url = ' https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$query&usehistory=y'
# r = requests.get(url)
# WITH LOCAL DB
# conn = sqlite3.connect("SRAmetadb.sqlite")
# cur = conn.cursor()
# cur.execute('SELECT study_accession FROM sra')
# a_list = cur.fetchall()
# accession_list = []
# for x in a_list:
#     if not x[0].isdigit() and not x[0] == 'SAMEAPERIM' and not x[0] == '-':
#         accession_list.append(x[0])
#         # print(x[0])
# db = SRAdb("SRAmetadb.sqlite")
# # test = db.sra_metadata(accession_list[0:100], detailed=True)
# count = 0
# for x in accession_list[0:10000]:
#     try:
#         count += 1
#         db.sra_metadata(x, detailed=True)
#     except ValueError:
#         print(x)

# WITH WEB
# db = SRAweb()
# df = db.sra_metadata("DRP002325", detailed=True)
# df.to_csv('out.csv', sep='\t', encoding='utf-8')

# df = db.search_sra(search_str='"test"')
# print(df.to_json())
# print(df)
# test = db.sra_metadata(accession_list[0:100], detailed=True)


# ftp = ftplib.FTP(
#     'ftp.ncbi.nlm.nih.gov')
# ftp.login()
# ftp.cwd('/sra/reports/Metadata')
# entries = list(ftp.mlsd())
# metadata = [entry for entry in entries if 'SRA_Accessions.tab' in entry[0]]
# metadata.sort(key=lambda entry: entry[1]['modify'], reverse=True)
# filename = metadata[0][0]
# with open(filename, "wb") as file:
#     # use FTP's RETR command to download the file
#     ftp.retrbinary(f"RETR {filename}", file.write)
# ftp.quit()

# file = tarfile.open(filename)
# file.extract('sample.txt', './Destination_FolderName')
# file.close()
# DRR000533
# DRR000588

db = SRAweb()
# df = db.sra_metadata('DRR000588', detailed=True)
# print(df.to_dict())
print('Reading SRA_Accessions')
df = pd.read_csv(r"SRA_Accessions.tab", sep="\t",
                 usecols=['Accession', 'Type', 'Status', 'Updated', 'Published', 'Experiment', 'Sample', 'BioProject', 'ReplacedBy'])
only_live = df[df['Status'] == 'live']
filtered = only_live[only_live['Type'] == 'STUDY']
accession_list = filtered[['Accession', 'Type', 'Status', 'Updated', 'Published',
                           'Experiment', 'Sample', 'BioProject', 'ReplacedBy']].values.tolist()
print(len(accession_list))
print('Getting Metadata')
count = 0

# Multiple at a time
# df = db.sra_metadata(accession_list[0:100], detailed=True)
# for i in range(100, len(accession_list), 100):
#     count += 100
#     try:
#         start = time.time()
#         to_merge = db.sra_metadata(accession_list[i:i+100], detailed=True)
#         df = pd.concat((df, to_merge), axis=1)
#         print(f'Time: {time.time() - start}')
#         print(count)
#         if count % 100 == 0:
#             break
#     except KeyError:
#         print(accession_list[i:i+100])
#         continue

# One at a time
metadata_list = []
for x in accession_list:
    try:
        start = time.time()
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
        print(f'Time: {time.time() - start}')
        count += 1
        if count % 1000 == 0:
            break
        if count % 100 == 0:
            print(count)
    except KeyError as e:
        print(x[0])
        print(e)
        continue

for metadata in metadata_list:
    if test := metadata.get('ena_fastq_ftp'):
        if test[0]:
            pprint(metadata)
    if test2 := metadata.get('ena_fastq_http'):
        if test2[0]:
            pprint(metadata)
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
            species_dict['additionalType'] = 'http://purl.obolibrary.org/obo/NCIT_C45293'
            if species_dict not in species_list:
                species_list.append(species_dict)
    elif organism_taxids:
        for taxid in organism_taxids:
            species_dict = {}
            species_dict['identifier'] = taxid
            species_dict['additionalType'] = 'http://purl.obolibrary.org/obo/NCIT_C45293'
            if species_dict not in species_list:
                species_list.append(species_dict)
    elif organism_names:
        for name in organism_names:
            species_dict = {}
            species_dict['name'] = name
            species_dict['additionalType'] = 'http://purl.obolibrary.org/obo/NCIT_C45293'
            if species_dict not in species_list:
                species_list.append(species_dict)

    # isBasedOn
    is_based_on = []
    # runs
    if run_accessions := metadata.get('run_accession'):
        for run_accession in run_accessions:
            run_dict = {}
            run_dict['identifier'] = run_accession
            run_dict['additionalType'] = 'http://purl.obolibrary.org/obo/NCIT_C47911'
            is_based_on.append(run_dict)
    # bioproject
    if bio_project := metadata.get('BioProject'):
        is_based_on.append(
            {'identifier': bio_project, 'additionalType': 'http://purl.obolibrary.org/obo/NCIT_C175890'})

    # experiments
    if experiment_accessions := metadata.get('experiment_accession'):
        if experiment_title := metadata.get('experiment_title'):
            if experiment_desc := metadata.get('experiment_desc'):
                for experiment_accession, experiment_title, experiment_desc in zip(experiment_accessions, experiment_title, experiment_desc):
                    experiment_dict = {}
                    experiment_dict['identifier'] = experiment_accession
                    experiment_dict['name'] = experiment_title
                    experiment_dict['description'] = experiment_desc
                    experiment_dict['additionalType'] = 'http://purl.obolibrary.org/obo/NCIT_C42790'
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
                    sample_dict['additionalType'] = 'http://purl.obolibrary.org/obo/NCIT_C70699'
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
                    instrument_dict['additionalType'] = 'http://purl.obolibrary.org/obo/NCIT_C16742'
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
                        cell_dict['additionalType'] = 'http://purl.obolibrary.org/obo/NCIT_C12508'
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

    # pprint(output)

    # each run
    # pprint(metadata)

    # with open('output3.txt', 'w') as f:
    #     for record in metadata_list:
    #         f.write(json.dumps(record) + '\n')

    # df = df.fillna(np.nan).replace([np.nan], [None])
    # with open('output.txt', 'w') as f:
    #     for row in df.to_dict('records'):
    #         f.write(json.dumps(row) + '\n')
    # df = df.fillna(np.nan).replace([np.nan], [None])
    # df_dict = df.to_dict('records')
    # with open('output2.txt', 'w') as f:
    #     f.write(json.dumps(df_dict) + '\n')
    # print(metadata_list)
    # One at a time
    # start = time.time()
    # df = db.sra_metadata(accession_list[0:100], detailed=True)
    # df2 = db.sra_metadata(accession_list[100:200], detailed=True)
    # merged = pd.concat([df, df2])
    # print(f'Time: {time.time() - start}')
    # # metadata_list.append(df.to_json(orient='records', lines=True))
    # d = merged.to_dict(orient='records')
    # print(d)
    # print(d[0])
    # with open('sra_metadata.json', 'w') as f:
    #     json.dump(metadata_list, f)
    # for x in accession_list:
    #     try:
    #         count += 1
    #         print(count)
    #         if count % 10 == 0:
    #             break
    #     except ValueError:
    #         print(x)
    #         continue
    # print(metadata_list)
    # try:
    #     count += 1
    #     db.sra_metadata(x, detailed=True)
    # except ValueError:
    #     print(x)
print('finished')
# for accession_id in accession_list:
#     count += 1
#     df = db.sra_metadata(srp=accession_id, detailed=True,
#                          expand_sample_attributes=True, sample_attribute=True)
#     if count % 10 == 0:
#         print(count)
#     # print(df.to_json())
# print("Parsing Metadata")
# df = db.sra_metadata(srp=accession_list, detailed=True,
#                      expand_sample_attributes=True, sample_attribute=True)

# df = db.sra_metadata(srp="SRX15798918", detailed=True,
#                      expand_sample_attributes=True, sample_attribute=True)
# df = db.sra_metadata(accession_list[0:100], detailed=True)
# print(df.to_json())
# print(df.to_json())
# test = db.sra_metadata(accession_list[0:100], detailed=True)
