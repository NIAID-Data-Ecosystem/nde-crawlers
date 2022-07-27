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
# df = db.sra_metadata(srp="DRP000196", detailed=True)
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
filtered = only_live[only_live['Type'] == 'STUDY'].replace({'-': None})
accession_list = filtered[['Accession', 'Type', 'Status', 'Updated', 'Published',
                           'Experiment', 'Sample', 'BioProject', 'ReplacedBy']].values.tolist()
print(len(accession_list))
print('Getting Metadata')
count = 100

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
        meta_dict = meta_df.replace({'<NA>': None}).to_dict()
        meta_dict['Accession'] = x[0]
        drs_url = f'https://locate.be-md.ncbi.nlm.nih.gov/idx/v1/{x[0]}?submitted=true&etl=false'
        drs_response = requests.get(drs_url)
        drs_json = json.loads(drs_response.text)
        drs_id = drs_json['response'][x[0]]['drs']
        meta_dict['DRS'] = drs_id
        meta_dict['Type'] = x[1]
        meta_dict['Updated'] = x[3]
        meta_dict['Published'] = x[4]
        meta_dict['Experiment'] = x[5]
        meta_dict['Sample'] = x[6]
        meta_dict['BioProject'] = x[7]
        meta_dict['ReplacedBy'] = x[8]
        metadata_list.append(meta_dict)
        print(f'Time: {time.time() - start}')
        count += 1
        if count % 105 == 0:
            break
    except KeyError:
        print(x)
        continue

for metadata in metadata_list:
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
        output['dateModified'] = updated
    if published := metadata.get('Published'):
        output['datePublished'] = published
    if recieved := metadata.get('Recieved'):
        output['dateCreated'] = recieved
    if visibility := metadata.get('Visibility'):
        output['conditionsOfAccess'] = visibility
    if bio_project := metadata.get('BioProject'):
        output['isBasedOn'] = {
            'identifier': bio_project
        }
    if replaced_by := metadata.get('ReplacedBy'):
        output['sameAs'] = replaced_by
    pprint(output)

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
