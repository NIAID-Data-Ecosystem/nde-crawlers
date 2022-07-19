from datetime import datetime
import json
import unicodedata
import xmltodict
from sickle import Sickle
from lxml import etree
from pprint import pprint

sickle = Sickle('https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi',
                max_retries=10, default_retry_after=20)

# records = sickle.ListRecords(
#     metadataPrefix='pmc', ignore_deleted=True)

record = sickle.GetRecord(
    identifier='oai:pubmedcentral.nih.gov:8524328', metadataPrefix='pmc')

xml_string = etree.tostring(record.xml)
root = etree.fromstring(xml_string)
metadata = record.metadata
header = record.header

output = {"includedInDataCatalog":
          {"name": "NCBI PMC",
           'versionDate': datetime.today().isoformat(),
           'url': "https://www.ncbi.nlm.nih.gov/pmc/"},
          'dateModified': datetime.strptime(record.header.datestamp, '%Y-%m-%d').isoformat(),
          "@type": "Dataset"
          }

if journal_title := metadata.get('journal-title'):
    output['journalName'] = journal_title[0]

if issn := metadata.get('issn'):
    output['issueNumber'] = issn[0]

if license := metadata.get('license_ref'):
    output['license'] = license[0]

if keywords := metadata.get('kwd'):
    output['keywords'] = keywords

funding_group = root.findall(
    ".//{https://jats.nlm.nih.gov/ns/archiving/1.3/}funding-source")
funder_list = []
for funder in funding_group:
    if funder.text is not None:
        funder_list.append(
            {'funder': {'name': funder.text}})
    else:
        name = funder.find(
            './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}institution')
        # institution_id = funder.find(
        #     './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}institution-id')
        if name is not None:
            funder_list.append(
                {'funder': {'name': name.text}})
if len(funder_list):
    output['funding'] = funder_list

publisher_sec = root.find(
    ".//{https://jats.nlm.nih.gov/ns/archiving/1.3/}publisher").getchildren()
for publisher in publisher_sec:
    if publisher.tag == '{https://jats.nlm.nih.gov/ns/archiving/1.3/}publisher-name':
        output['sdPublisher'] = {'name': publisher.text}

if title := metadata.get('article-title'):
    output['name'] = title

article_meta_sec = root.find(
    ".//{https://jats.nlm.nih.gov/ns/archiving/1.3/}article-meta").getchildren()
for article_meta in article_meta_sec:
    if article_meta.tag == '{https://jats.nlm.nih.gov/ns/archiving/1.3/}title-group':
        for child in article_meta.getchildren():
            if child.tag == '{https://jats.nlm.nih.gov/ns/archiving/1.3/}article-title':
                output['name'] = child.text

if identifiers := metadata.get('article-id'):
    for identifier in identifiers:
        if identifier.startswith('PMC'):
            output['identifier'] = identifier
            output['url'] = f'https://www.ncbi.nlm.nih.gov/pmc/articles/{identifier}'
            output['_id'] = 'NCBI_PMC_' + identifier
        if identifier.startswith('10.'):
            output['doi'] = identifier

pmid = root.find(
    ".//{https://jats.nlm.nih.gov/ns/archiving/1.3/}article-id[@pub-id-type='pmid']").text
if pmid is not None:
    output['pmids'] = pmid

contrib_list = root.findall(
    './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}contrib')
author_list = []
for contrib in contrib_list:
    author_dict = {}
    surname = contrib.find(
        './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}surname')
    if surname is not None:
        author_dict['familyName'] = surname.text
    given_name = contrib.find(
        './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}given-names')
    if given_name is not None:
        author_dict['givenName'] = given_name.text
    if bool(author_dict):
        author_list.append(author_dict)
if len(author_list):
    output['author'] = author_list

# Supplemental Data
supplemental_data_arr = root.findall(
    './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}supplementary-material')
for supplemental_data in supplemental_data_arr:
    caption = supplemental_data.find(
        '{https://jats.nlm.nih.gov/ns/archiving/1.3/}caption')
    if caption is not None:
        file_name = caption.find(
            '{https://jats.nlm.nih.gov/ns/archiving/1.3/}p').text
        file_name = unicodedata.normalize('NFKD', file_name)
    media = supplemental_data.find(
        '{https://jats.nlm.nih.gov/ns/archiving/1.3/}media')
    if media is not None:
        file_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{output['identifier']}/bin/"+media.get(
            '{http://www.w3.org/1999/xlink}href')
    if 'file_name' in locals() and 'file_url' in locals():
        test = {'file_name': file_name, 'file_url': file_url}

# Data Availability
data_availability_sec = root.find(
    './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}sec[@sec-type="data-availability"]')
if data_availability_sec is not None:
    data_availability_titles = data_availability_sec.findall(
        '{https://jats.nlm.nih.gov/ns/archiving/1.3/}p')
    if data_availability_titles is not None:
        for title in data_availability_titles:
            data_availability_dict = {}
            content_list = []
            data_availability_text = title.text
            if data_availability_text is not None:
                data_availability_dict['title'] = data_availability_text

            data_availability_content_titles = title.findall(
                './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}p')
            for title in data_availability_content_titles:
                data_availability_content_link = title.find(
                    './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}ext-link')
                content_list.append({
                    'link_title': title.text.replace(' (', ''), 'link': data_availability_content_link.text})
            data_availability_dict['content'] = content_list

# abstract
abstract_sec = root.find(
    './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}abstract')
if abstract_sec is not None:
    description_string = ''
    abstract_title = abstract_sec.findall(
        './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}title')
    abstract_text = abstract_sec.findall(
        './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}p')
    for title, text in zip(abstract_title, abstract_text):
        if title is not None:
            if description_string == '':
                description_string = title.text
            else:
                description_string += '\n' + title.text
        if text is not None:
            if description_string == '':
                description_string = text.text
            else:
                description_string += '\n' + text.text
    if description_string != '':
        output['description'] = description_string
print(output['description'])
