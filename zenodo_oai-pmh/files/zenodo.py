# to query one document: https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:5680920&metadataPrefix=oai_datacite
# to query all documents: https://zenodo.org/oai2d?verb=ListRecords&metadataPrefix=oai_datacite
# import re
import time
import logging
import datetime
import shelve
import sqlite3
import json

from sql_database import NDEDatabase
from sickle import Sickle
from pprint import pprint
from xml.etree import ElementTree

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


class Zenodo(NDEDatabase):
    # override variables
    DBM_NAME = "zenodo.db"
    EXPIRE = datetime.timedelta(days=90)

    # connect to the website
    sickle = Sickle('https://zenodo.org/oai2d', max_retries=3)

    def new_cache(self):
        """Creates a new tables: metadata and cache. Upserts two entries: date_created, date_updated in metadata table"""
        con = sqlite3.connect(self.path + '/' + self.DBM_NAME)
        c = con.cursor()

        c.execute("""CREATE TABLE IF NOT EXISTS metadata (
                name text NOT NULL PRIMARY KEY,
                date text NOT NULL
                )""")

        today = datetime.date.today().isoformat()

        # used for testing
        # today = datetime.date(2022, 6, 1).isoformat()

        # upserting in sqlite https://www.sqlite.org/lang_UPSERT.html
        # https://stackoverflow.com/questions/62274285/sqlite3-programmingerror-incorrect-number-of-bindings-supplied-the-current-stat
        c.execute("""INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """, ('date_created', today))
        c.execute("""INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """, ('date_updated', today))
        con.commit()

        c.execute("DROP TABLE IF EXISTS cache")
        c.execute("""CREATE TABLE cache (
                      _id text NOT NULL PRIMARY KEY,
                      data text NOT NULL
                     )""")
        con.commit()

        con.close()

    def dump(self):
        """Connects to sickle and stores raw data into the cache table, only runs with cache is expired"""
        con = sqlite3.connect(self.path + '/' + self.DBM_NAME)
        c = con.cursor()

        records = self.sickle.ListRecords(metadataPrefix='oai_datacite', ignore_deleted=True)

        # used to test individual records
        # record = self.sickle.GetRecord(identifier='oai:zenodo.org:1202262', metadataPrefix='oai_datacite')
        # pprint(record.metadata)
        # pprint(type(record))
        # print(type(record.metadata))
        # print(type(record.header))
        # print(vars(record.header))
        count = 0
        while True:
            try:
                # get the next item
                record = records.next()
                count += 1
                if count % 100 == 0:
                    time.sleep(.5)
                    logger.info("Dumped %s records", count)

                # in each doc we want record.identifier and record stored
                doc = {'header': dict(record.header), 'metadata': record.metadata,
                       'xml': ElementTree.tostring(record.xml, encoding='unicode')}

                c.execute("INSERT INTO cache VALUES(?, ?)", (record.header.identifier, json.dumps(doc)))
                con.commit()

            except StopIteration:
                logger.info("Finished Dumping. Total Records: %s", count)
                con.close()
                # if StopIteration is raised, break from loop
                break

    def parse(self):
        """Transforms/pipeline data to the nde schema before writing the information into the ndson file"""
        # dictionary to convert general type to @type
        # https://docs.google.com/spreadsheets/d/1DOwMjvFL3CGPkdoaFCveKNb_pPVpboEgUjRTT_3eAW0/edit#gid=0
        gen_type = {'annotationcollection': 'Collection', 'article': 'ScholarlyArticle', 'audiovisual': 'MediaObject',
                    'book': 'Book', 'bookchapter': 'Chapter', 'collection': 'Collection',
                    'conferencepaper': 'ScholarlyArticle',
                    'datamanagementplan': 'CreativeWork', 'dataset': 'Dataset', 'deliverable': 'CreativeWork',
                    'diagram': 'ImageObject', 'drawing': 'Drawing', 'figure': 'ImageObject', 'image': 'ImageObject',
                    'interactiveresource': 'CreativeWork', 'journalarticle': 'ScholarlyArticle',
                    'lesson': 'LearningResource',
                    'other': 'CreativeWork', 'outputmanagementplan': 'CreativeWork', 'patent': 'CreativeWork',
                    'photo': 'Photograph',
                    'physicalobject': 'Thing', 'plot': 'ImageObject', 'poster': 'Poster',
                    'preprint': 'ScholarlyArticle',
                    'presentation': 'PresentationDigitalDocument', 'projectdeliverable': 'CreativeWork',
                    'projectmilestone': 'CreativeWork', 'proposal': 'CreativeWork', 'publication': 'ScholarlyArticle',
                    'report': 'Report', 'section': 'CreativeWork', 'software': 'ComputationalTool',
                    'softwaredocumentation': 'TechArticle', 'taxonomictreatment': 'ScholarlyArticle',
                    'technicalnote': 'TechArticle',
                    'thesis': 'ScholarlyArticle', 'video': 'VideoObject', 'workflow': 'CreativeWork', 'workingpaper': 'ScholarlyArticle'}
        # dictionary log the types that cannot be converted
        missing_types = {}

        con = sqlite3.connect(self.path + '/' + self.DBM_NAME)
        c = con.cursor()

        c.execute("SELECT * from cache")
        items = c.fetchall()
        for item in items:
            record = json.loads(item[1])
            identifier = item[0]

            # format the identifier for _id, and identifier
            identifier_split = identifier.rsplit(':', 1)
            identifier = identifier_split[-1]
            # use as much of the metadata variable or header variable to format the transformation
            output = {
                '@context': "https://schema.org/",
                'includedInDataCatalog': {
                    '@type': 'Organization',
                    'name': 'Zenodo',
                    'url': 'https://zenodo.org/',
                    'versionDate': datetime.date.today().isoformat()
                },
                '_id': 'ZENODO_' + identifier,
                'name': record['metadata'].get('title')[0],
                'author': [],
                'identifier': "zenodo." + identifier,
                'dateModified': datetime.datetime.fromisoformat(
                    record['header']['datestamp'][:-1]).astimezone(datetime.timezone.utc)
                    .date().isoformat(),
                'url': "https://zenodo.org/record/" + identifier_split[-1],
            }

            if description := record['metadata'].get('description'):
                output['description'] = description[0]

            if date_published := record['metadata'].get('date'):
                output['datePublished'] = date_published[0]

            if language := record['metadata'].get('language'):
                output['inLanguage'] = {'name': language[0]}

            # zenodo uses different delimiters for keywords. See examples: oai:zenodo.org:1188946, oai:zenodo.org:1204780
            # there may be more delimiters than just '; ' and ', '
            # since zenodo does process the delimiters, for now we will not either example: oai:zenodo.org:1204780
            if keywords := record['metadata'].get('subject'):
                output['keywords'] = []
                for keyword in keywords:
                    # keyword = re.split('; |, ', keyword)
                    # for k in keyword:
                    output['keywords'].append(keyword)

            # used for testing to print out xml tags
            # for element in root.iter():
            #     print("%s - %s" % (element.tag, element.text))

            # use xml to query doi
            root = ElementTree.fromstring(record['xml'])

            doi = root.find(".//{http://datacite.org/schema/kernel-3}identifier[@identifierType='DOI']")
            if doi is not None:
                output['doi'] = doi.text

            # use xml to get the type
            zenodo_type = root.find(".//{http://datacite.org/schema/kernel-3}resourceType[@resourceTypeGeneral]").get('resourceTypeGeneral')
            zenodo_type2 = root.find(".//{http://datacite.org/schema/kernel-3}resourceType[@resourceTypeGeneral]").text

            # format the types to be case insensitive and query the dictionary for the transformation
            zenodo_type = zenodo_type.lower()
            if zenodo_type2:
                zenodo_type2 = zenodo_type2.lower()
                zenodo_type2 = zenodo_type2.replace(" ", "")

            if zenodo_type is not None:
                if zenodo_type in gen_type.keys():
                    output['@type'] = gen_type[zenodo_type]
                elif zenodo_type2 in gen_type.keys():
                    output['@type'] = gen_type[zenodo_type2]
                else:
                    print(identifier_split)
                    missing_types[(zenodo_type, zenodo_type2)] = (zenodo_type, zenodo_type2)

            # xml to find all the authors and format them
            # we cannot use metadata to format due to creator and affiliation being in separate lists
            creators = root.findall(".//{http://datacite.org/schema/kernel-3}creator")
            for creator in creators:
                author = {}
                name = creator.find('./{http://datacite.org/schema/kernel-3}creatorName')
                affiliation = creator.find('./{http://datacite.org/schema/kernel-3}affiliation')
                orcid_id = creator.find("./{http://datacite.org/schema/kernel-3}nameIdentifier[@nameIdentifierScheme='ORCID']")
                if name is not None:
                    author['name'] = name.text
                if affiliation is not None:
                    # elasticsearch cannot index strings greater than 32766. Someone mistakenly put a WHOLE ARTICLE into the affiliation
                    # https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:4675716&metadataPrefix=oai_datacite
                    if len(affiliation.text) < 30000:
                        author['affiliation'] = {'name': affiliation.text}
                if orcid_id is not None:
                    author['identifier'] = orcid_id.get('schemeURI') + orcid_id.text
                if author:
                    output['author'].append(author)

            # use xml to find the license and conditionOfAccess
            rights = root.findall(".//{http://datacite.org/schema/kernel-3}rights")
            for right in rights:
                if "access" in right.text.lower():
                    output['conditionsOfAccess'] = right.text.split(' ', 1)[0]
                else:
                    output['license'] = right.get('rightsURI')

            # use xml to find citedBy field
            related_ids = root.findall(".//{http://datacite.org/schema/kernel-3}relatedIdentifier")
            cited_by = []
            for related_id in related_ids:
                if related_id.get('relationType') == "IsCitedBy" and related_id.get('relatedIdentifierType') == "URL":
                    cited_by.append({'url': related_id.text})
            if cited_by:
                output['citedBy'] = cited_by

            # use xml to find funding
            contributors = root.findall(".//{http://datacite.org/schema/kernel-3}contributor[@contributorType='Funder']")
            funding = []
            for contributor in contributors:
                name = contributor.find("./{http://datacite.org/schema/kernel-3}contributorName")
                if name is not None:
                    funding.append({'funder': {'name': name.text}})
            if funding:
                output['funding'] = funding

            """ TODO try to get the codeRepository field. Not confident how to get it. 
            "if this resource is a 'software' @type and has a relatedIdentifier with the isSupplementTo field".
            Examples: 
            * 	https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:5680920&metadataPrefix=oai_datacite
                    <relatedIdentifier relatedIdentifierType="URL" relationType="IsDerivedFrom" >https://gitlab.com/astron-idg/idg-fpga/</relatedIdentifier>
            * 	https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:1135290&metadataPrefix=oai_datacite
                    <relatedIdentifier relatedIdentifierType="URL" relationType="IsSupplementTo" >https://github.com/ljcohen/planets/tree/v0.1</relatedIdentifier>
            """

            yield output

        # output the missing transformations for @type
        if len(missing_types.keys()) > 0:
            logger.warning("Missing type transformation: {}".format(str(missing_types.keys())))

    def update_cache(self):
        """If cache is not expired get the new records to add to the cache since last_updated"""
        con = sqlite3.connect(self.path + '/' + self.DBM_NAME)
        c = con.cursor()

        c.execute("SELECT date from metadata WHERE name='date_updated'")
        last_updated = c.fetchall()
        assert len(last_updated) <= 1, "There is more than one last_created."
        last_updated = last_updated[0][0]

        logger.info("Updating cache from %s", last_updated)
        records = self.sickle.ListRecords(**{
                'metadataPrefix': 'oai_datacite',
                'ignore_deleted': True,
                'from': last_updated
            }
        )

        count = 0
        while True:
            try:
                # get the next item
                record = records.next()
                count += 1
                if count % 100 == 0:
                    time.sleep(.5)
                    logger.info("Updating cache. %s new updated records", count)

                # in each doc we want record.identifier and record stored
                doc = {'header': dict(record.header), 'metadata': record.metadata,
                       'xml': ElementTree.tostring(record.xml, encoding='unicode')}

                c.execute("""INSERT INTO cache VALUES(?, ?)
                                ON CONFLICT(_id) DO UPDATE SET data=excluded.data
                          """, (record.header.identifier, json.dumps(doc)))
                con.commit()

            except StopIteration:
                logger.info("Finished updating cache. Total new records: %s", count)
                # if StopIteration is raised, break from loop
                break

        today = datetime.date.today().isoformat()
        c.execute("""INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """, ('date_updated', today))

        logger.info("Cache last_updated on: %s", today)
        con.commit()

        con.close()
