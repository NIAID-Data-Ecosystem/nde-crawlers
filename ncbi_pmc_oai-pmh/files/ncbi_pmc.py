import datetime
import json
import unicodedata
import logging

from sickle import Sickle
from lxml import etree
from sql_database import NDEDatabase


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('nde-logger')


class NCBI_PMC(NDEDatabase):
    # override variables
    SQL_DB = "ncbi_pmc.db"
    EXPIRE = datetime.timedelta(days=90)

    sickle = Sickle('https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi',
                    max_retries=10, default_retry_after=20)

    def load_cache(self):
        """Retrives the raw data using a sickle request and formats so dump can store it into the cache table
        Returns:
            A tuple (_id, data formatted as a json string)
        """

        records = self.sickle.ListRecords(
            metadataPrefix='pmc', ignore_deleted=True)
        count = 0
        while True:
            try:
                # get the next item
                record = records.next()
                count += 1
                if count % 100 == 0:
                    logger.info("Loading cache. Loaded %s records", count)

                # in each doc we want record.identifier and record stored
                doc = {'header': dict(record.header), 'metadata': record.metadata,
                       'xml': etree.tostring(record.xml, encoding='unicode')}

                yield (record.header.identifier, json.dumps(doc))

            except StopIteration:
                logger.info("Finished Loading. Total Records: %s", count)
                # if StopIteration is raised, break from loop
                break

    def parse(self, records):
        for record in records:
            data = json.loads(record[1])
            root = etree.fromstring(data['xml'])
            metadata = data['metadata']
            header = data['header']

            output = {"includedInDataCatalog":
                      {"name": "NCBI PMC",
                       'versionDate': datetime.date.today().strftime('%Y-%m-%d'),
                       'url': "https://www.ncbi.nlm.nih.gov/pmc/"},
                      'dateModified': datetime.datetime.strptime(header['datestamp'], '%Y-%m-%d').strftime('%Y-%m-%d'),
                      "@type": "Dataset"
                      }

            citation_dict = {}
            if identifiers := metadata.get('article-id'):
                for identifier in identifiers:
                    if identifier.startswith('PMC'):
                        output['identifier'] = identifier
                        output['url'] = f'https://www.ncbi.nlm.nih.gov/pmc/articles/{identifier}'
                        output['_id'] = 'NCBI_PMC_' + identifier
                    if identifier.startswith('10.'):
                        citation_dict['doi'] = identifier
            notes = root.find(
                './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}notes')
            if notes is not None:
                volume = notes.find(
                    './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}volume')
                if volume is not None:
                    citation_dict['volume'] = volume.text

            # if bool(citation_dict):
            #     output['citation'] = [citation_dict]

            # Supplemental Data
            supplemental_data_arr = root.findall(
                './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}supplementary-material')
            distribuiton_list = []
            for supplemental_data in supplemental_data_arr:
                distribution_obj = {}
                caption = supplemental_data.find(
                    '{https://jats.nlm.nih.gov/ns/archiving/1.3/}caption')
                if caption is not None:
                    file_name = caption.find(
                        '{https://jats.nlm.nih.gov/ns/archiving/1.3/}title')
                    if file_name is not None:
                        if file_name.text is not None:
                            name = unicodedata.normalize(
                                'NFKD', file_name.text)
                            distribution_obj['name'] = name
                    description = caption.find(
                        '{https://jats.nlm.nih.gov/ns/archiving/1.3/}p')
                    if description is not None:
                        distribution_obj['description'] = description.text
                media = supplemental_data.find(
                    '{https://jats.nlm.nih.gov/ns/archiving/1.3/}media')
                if media is not None:
                    file_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{output['identifier']}/bin/"+media.get(
                        '{http://www.w3.org/1999/xlink}href')
                    distribution_obj['url'] = file_url
                if bool(distribution_obj):
                    distribuiton_list.append(distribution_obj)
            if len(distribuiton_list):
                output['distribution'] = distribuiton_list

            if journal_title := metadata.get('journal-title'):
                citation_dict['journalName'] = journal_title[0]

            pub_date = root.find(
                './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}pub-date[@pub-type="epub"]')
            if pub_date is not None:
                date_string = ''
                year = pub_date.find(
                    '{https://jats.nlm.nih.gov/ns/archiving/1.3/}year')
                month = pub_date.find(
                    '{https://jats.nlm.nih.gov/ns/archiving/1.3/}month')
                day = pub_date.find(
                    '{https://jats.nlm.nih.gov/ns/archiving/1.3/}day')
                if year is not None:
                    date_string = year.text
                if month is not None:
                    if len(month.text) == 1:
                        date_string += f'-0{month.text}'
                    else:
                        date_string += '-' + month.text
                if day is not None:
                    if len(day.text) == 1:
                        date_string += f'-0{day.text}'
                    else:
                        date_string += '-' + day.text
                if date_string != '':
                    citation_dict['datePublished'] = date_string

            # if issn := metadata.get('issn'):
            #     output['issueNumber'] = issn[0]

            # if license := metadata.get('license_ref'):
            #     output['license'] = license[0]

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
                ".//{https://jats.nlm.nih.gov/ns/archiving/1.3/}publisher")
            if publisher_sec is not None:
                publisher_list = publisher_sec.getchildren()
                for publisher in publisher_list:
                    if publisher.tag == '{https://jats.nlm.nih.gov/ns/archiving/1.3/}publisher-name':
                        output['sdPublisher'] = {'name': publisher.text}

            # if title := metadata.get('article-title'):
            #     if title[0] is not None:
            #         output['name'] = 'Supplementary materials in ' + title[0]

            article_meta_sec = root.find(
                ".//{https://jats.nlm.nih.gov/ns/archiving/1.3/}article-meta")
            if article_meta_sec is not None:
                if 'volume' not in citation_dict:
                    volume = article_meta_sec.find(
                        './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}volume')
                    if volume is not None:
                        citation_dict['volume'] = volume.text

                article_meta_list = article_meta_sec.getchildren()
                for article_meta in article_meta_list:
                    if article_meta.tag == '{https://jats.nlm.nih.gov/ns/archiving/1.3/}title-group':
                        for child in article_meta.getchildren():
                            if child.tag == '{https://jats.nlm.nih.gov/ns/archiving/1.3/}article-title':
                                if child.text is not None:
                                    output['name'] = 'Supplementary materials in ' + \
                                        '"'+child.text+'"'
                                    citation_dict['name'] = child.text

            pmid = root.find(
                ".//{https://jats.nlm.nih.gov/ns/archiving/1.3/}article-id[@pub-id-type='pmid']")
            if pmid is not None:
                citation_dict['pmid'] = pmid.text
                citation_dict['identifier'] = 'PMID:' + pmid.text
                citation_dict['url'] = f'https://pubmed.ncbi.nlm.nih.gov/{pmid.text}'

            contrib_list = root.findall(
                './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}contrib')
            author_list = []
            for contrib in contrib_list:
                author_dict = {}
                author_id = contrib.find(
                    '{https://jats.nlm.nih.gov/ns/archiving/1.3/}contrib-id')
                if author_id is not None:
                    author_dict['identifier'] = author_id.text
                surname = contrib.find(
                    './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}surname')
                if surname is not None:
                    author_dict['familyName'] = surname.text
                given_name = contrib.find(
                    './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}given-names')
                if given_name is not None:
                    author_dict['givenName'] = given_name.text
                if bool(author_dict):
                    author_dict['@type'] = 'Person'
                    author_list.append(author_dict)
            if len(author_list):
                output['author'] = author_list
                citation_dict['author'] = author_list

            # Data Availability
            # data_availability_sec = root.find(
            #     './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}sec[@sec-type="data-availability"]')
            # if data_availability_sec is not None:
            #     data_availability_titles = data_availability_sec.findall(
            #         '{https://jats.nlm.nih.gov/ns/archiving/1.3/}p')
            #     if data_availability_titles is not None:
            #         for title in data_availability_titles:
            #             data_availability_dict = {}
            #             content_list = []
            #             data_availability_text = title.text
            #             if data_availability_text is not None:
            #                 data_availability_dict['title'] = data_availability_text

            #             data_availability_content_titles = title.findall(
            #                 './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}p')
            #             for title in data_availability_content_titles:
            #                 data_availability_content_link = title.find(
            #                     './/{https://jats.nlm.nih.gov/ns/archiving/1.3/}ext-link')
            #                 content_list.append({
            #                     'link_title': title.text.replace(' (', ''), 'link': data_availability_content_link.text})
            #             data_availability_dict['content'] = content_list

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
                    if title.text is not None:
                        if description_string == '':
                            description_string = title.text
                        else:
                            description_string += '\n' + title.text
                    if text.text is not None:
                        if description_string == '':
                            description_string = text.text
                        else:
                            description_string += '\n' + text.text
                if description_string != '':
                    output['description'] = description_string
            if bool(citation_dict):
                output['citation'] = citation_dict

            if 'distribution' in output:
                yield output

    def update_cache(self):
        """If cache is not expired get the new records to add to the cache since last_updated"""

        last_updated = self.retreive_last_updated()

        # get all the records since last_updated to add into current cache
        logger.info("Updating cache from %s", last_updated)
        records = self.sickle.ListRecords(**{
            'metadataPrefix': 'pmc',
            'ignore_deleted': True,
            'from': last_updated
        }
        )

        # Very similar to load_cache()
        count = 0
        while True:
            try:
                # get the next item
                record = records.next()
                count += 1
                if count % 100 == 0:
                    logger.info(
                        "Updating cache. %s new updated records", count)

                # in each doc we want record.identifier and record stored
                doc = {'header': dict(record.header), 'metadata': record.metadata,
                       'xml': etree.tostring(record.xml, encoding='unicode')}

                yield (record.header.identifier, json.dumps(doc))

            except StopIteration:
                logger.info(
                    "Finished updating cache. Total new records: %s", count)
                # if StopIteration is raised, break from loop
                break

        self.insert_last_updated()
