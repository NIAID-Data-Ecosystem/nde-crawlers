import datetime
# from sqlite3 import DataError
import time
import json
import logging
import requests
from html.parser import HTMLParser

from sql_database import NDEDatabase

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('nde-logger')

class Dataverse(NDEDatabase):
    SQL_DB = "dataverse.db"
    EXPIRE = datetime.timedelta(days=90)
    NO_CACHE = True

    DATAVERSE_SERVER = "https://dataverse.harvard.edu/api"
    DATA_URL = f"{DATAVERSE_SERVER}/search?q=data"
    EXPORT_URL = f"{DATAVERSE_SERVER}/datasets/export?exporter=schema.org"


    def scrape_schema_representation(self, url):
        """
        when the schema.org export of the dataset fails
        this will grab it from the url
        by looking for <script type="application/ld+json">
        """
        class SchemaScraper(HTMLParser):
            def __init__(self):
                super().__init__()
                self.readingSchema = False
                self.schema = None

            def handle_starttag(self, tag, attrs):
                if tag == 'script' and 'type' in attrs and attrs.get('type') == "application/ld+json":
                    self.readingSchema = True

            def handle_data(self, data):
                if self.readingSchema:
                    self.schema = data
                    self.readingSchema = False
        try:
            req = requests.get(url)
        except Exception as requestException:
            logger.info(f"failed to get {url} due to {requestException}")
            return False
        if not req.ok:
            logger.info(f"unable to retrieve {url}, status code:{req.status_code}")
            return False
        parser = SchemaScraper()
        parser.feed(req.text)
        if parser.schema:
            return parser.schema
        return False


    def get_schema_document(self, gid, url, verbose=False):
        schema_export_url = f"{self.EXPORT_URL}&persistentId={gid}"
        try:
            req = requests.get(schema_export_url)
        except Exception as requestException:
            logger.info(f"Failed to get {schema_export_url} due to {requestException}")
            return False
        try:
            res = req.json()
        except json.decoder.JSONDecodeError:
            logger.info(f"failed to get {schema_export_url}, {url}, {req.status_code}")
            return False
        if res.get('status') and res.get('status') == 'ERROR':
            #logger.warning(f"document export failed, scraping {url} instead")
            document = self.scrape_schema_representation(url)
            if document:
                logger.info(f"successfully scrapped document, {url}")
                return document
            else:
                logger.info(f"schema.org export failed on {schema_export_url}, {url},  {req.status_code}, {res.get('status')}, {res.get('message')}")
                logger.info("will attempt to yield original metadata for dataset....")
                return False
        else:
            # success, response is the document
            #logger.info(f"schema export successful {gid}")
            return res


    def compile_paginated_data(self, query_endpoint, per_page=1000, verbose=False):
        """ Extract Data By Page
        pages through data, compiling all response['data']['items']
        and returning them.
        per_page max is 1000
        """
        continue_paging = True
        start   = 0
        retries = 0
        #data_pages = []
        
        while continue_paging:
            url = f"{query_endpoint}&per_page={per_page}&start={start}"
            #logger.info(f"querying endpoint {url}")
            try:
                req = requests.get(url)
            except Exception as requestException:
                logger.error(f"Failed to get {url} due to {requestException}")
                if retries > 5:
                    logger.error(f"Failed to retrieve too many times on endpoint {url}")
                retries += 1
                # after 1 retry, limit per-page to 200, after 2, limit to 50
                per_page = 200 if retries == 1 else 50
            try:
                response = req.json()
            except ValueError:
                logger.error(f"Failed to get a JSON response from {url}")
            if response:
                try:
                    total = response.get('data').get('total_count')
                    page_data = response.get('data').get('items')
                    start += per_page
                    for page in page_data:
                        # dataset pages use "global_id" key
                        if "global_id" in page.keys():
                            yield (page['global_id'], page['url'], page)
                        # dataverse pages have "identifier" key
                        elif "identifier" in page.keys():
                            yield (page["identifier"], page['url'], page) # EXTRACT ID HERE 
                    #data_pages.extend(page_data)
                    continue_paging = total and start < total
                except Exception as exception:
                    logger.info("passing datapage because of exception: ", exception)
            else:
                pass


    def load_cache(self):
        """Retrives the raw data using a sickle request and formats so dump can store it into the cache table
        Returns:
            A tuple (_id, data formatted as a json string)
        grabs all datasets and files related to QUERIES both by querying
        and by grabbing everything in related dataverses
        extracts their global_id, which in this case is a DOI
        returns a dictionary mapping global_id -> dataset
        """
        logger.info("Load cache process starting....")
        start_time = time.process_time()
        dataset_ids = [] #set([None])
        #datasets    = []
        handle_ct = 0
        dataset_ct = 0
        dataverse_ct = 0  
        schema_ct = 0
        skipped_data = 0   
        # dataverse type: data --> (identifier, data_dict)
        logger.info("Starting extraction of dataverse type data....")
        for data in self.compile_paginated_data(self.DATA_URL+"&type=dataverse"):
            dataverse_query = f"{self.DATAVERSE_SERVER}/search?q=*&type=dataset&type=file&subtree={data[0]}"
            for dv_data in self.compile_paginated_data(dataverse_query):
                if "https://hdl.handle.net/" in dv_data[1]:
                    schema_record = self.scrape_schema_representation(dv_data[1])
                    handle_ct += 1
                else:
                    schema_record = self.get_schema_document(dv_data[0], dv_data[1])
                if schema_record:
                    dataset_ids.append(dv_data[0])
                    dataverse_ct += 1
                    schema_ct += 1
                    yield (schema_record['@id'], json.dumps(schema_record))
                elif schema_record == False:
                    if "https://hdl.handle.net/" in dv_data[1]:
                        pass
                    elif dv_data[2]:
                        logger.info(f"yielding original dataset for {dv_data[1]}")
                        yield (dv_data[0], json.dumps(dv_data[2]))
                    else:
                        skipped_data += 1
                        logger.info(f"error retrieving schema export for document with id, {dv_data[0]}, skipping dataset") 

        logger.info(f"Extraction of dataverse data complete, {dataverse_ct} dataset schemas were extracted, and {len(set(dataset_ids))} unique ids total.")
        logger.info("starting extraction of dataset type data....")
        # # dataset type: data --> (global_id, url)
        for data in self.compile_paginated_data(self.DATA_URL+"&type=dataset"):
            if data[0] not in dataset_ids:
                dataset_ids.append(data[0])
                if "https://hdl.handle.net/" in data[1]:
                    self.scrape_schema_representation(data[1])
                    handle_ct += 1
                else:
                    schema_record = self.get_schema_document(data[0], data[1])
                    if schema_record:
                        schema_ct += 1
                        dataset_ct += 1
                        #print(json.dumps(schema_record, indent=4))
                        yield (schema_record['@id'], json.dumps(schema_record))
                    elif schema_record == False:
                        if "https://hdl.handle.net/" in data[1]:
                            pass
                        elif data[2]:
                            logger.info(f"yielding original dataset for {data[1]}")
                            yield (data[0], json.dumps(data[2]))
                        else:
                            skipped_data += 1
                            logger.info(f"error retrieving schema export for document with id, {data[0]}, skipping dataset") 
        process_time = time.process_time() - start_time
        logger.info(f"Extraction of dataset data complete, {dataset_ct} dataset schemas were extracted, {len(set(dataset_ids))} unique ids total.")
        logger.info(f"{schema_ct} schemas extracted in {process_time} seconds.")
        logger.info ("LOAD CACHE PROCESS COMPLETE.")


    def parse(self, records):
        start_time = time.process_time()
        count = 0
        logger.info(f"Starting metadata parser...")
        for record in records:
            try:
                dataset=json.loads(record[1])
                # schema.org exported data
                if "@context" in dataset or "@type" in dataset:
                    ...
                    dataset=json.loads(record[1])
                    dataset['url'] = dataset['identifier']
                    dataset['doi'] = dataset.pop('identifier')
                    dataset['_id'] = dataset['doi'].replace("/","_").replace('https:__doi.org','Dataverse')
                    dataset['dateModified'] = datetime.datetime.strptime(dataset['dateModified'] , "%Y-%m-%d").date().isoformat()

                    if 'datePublished' in dataset:
                        dataset['datePublished'] = datetime.datetime.strptime(dataset['datePublished'] , "%Y-%m-%d").date().isoformat()

                    if dataset['author'] is None:
                        dataset.pop('author')
                        if dataset['creator']:
                            dataset['author'] = dataset.pop('creator')
                            for data_dict in dataset["author"]:
                                if "affiliation" in data_dict.keys():
                                    data_dict["affiliation"] = {"name": data_dict.pop("affiliation")}
                        else:
                            dataset.pop('creator')
                    else:
                        for data_dict in dataset["author"]:
                            if "affiliation" in data_dict.keys():
                                data_dict["affiliation"] = {"name": data_dict.pop("affiliation")}
                        dataset.pop('creator')

                    if dataset['publisher'] is None:
                        dataset.pop('publisher')
                        if dataset['provider']:
                            dataset['sdPublisher'] = [dataset.pop('provider')]
                        else:
                            dataset['sdPublisher'] = [dataset.pop('provider')]
                    else:
                        dataset['sdPublisher'] = [dataset.pop("publisher")]
                        dataset.pop('provider')

                    if "funder" in dataset:
                        dataset['funding'] = {"funder": dataset.pop('funder')}

                    if "description" in dataset:
                        dataset["description"] = dataset.pop("description")[0]

                    if "keywords" in dataset:
                        if dataset['keywords']:
                            dataset["keywords"]  = dataset.pop("keywords")[0]

                    if dataset["license"]:
                        if type(dataset["license"]) is str:
                            pass
                        elif type(dataset["license"]) is dict:
                            if "url" in dataset['license'].keys():
                                dataset["license"] = dataset['license']['url']
                            elif "text" in dataset['license'].keys():
                                dataset["license"] = dataset['license']['text']
                            else:
                                dataset.pop('license')

                    if "temporalCoverage" in dataset and dataset["temporalCoverage"]:
                        dataset["temporalCoverage"] = [{
                            "temporalInterval": {
                                "duration": dataset["temporalCoverage"][0]
                            }
                        }]

                    if "spatialCoverage" in dataset:
                        dataset["spatialCoverage"] = [{ "name": dataset.pop("spatialCoverage")[0] }]

                    if "distribution" in dataset:
                        for data_dict in dataset["distribution"]:
                            if "fileFormat" in data_dict.keys():
                                data_dict["encodingFormat"] = data_dict.pop("fileFormat")
                            if "identifier" in data_dict.keys():
                                data_dict["contentUrl"] = data_dict.pop("identifier")
                            if "contentSize" in  data_dict.keys():
                                data_dict.pop("contentSize")

                    pmids = []
                    if "citation" in dataset and dataset['citation']:
                        for data_dict in dataset['citation']:
                            if "@id" in data_dict.keys():
                                pmid = data_dict['@id'].split("/")[-1]
                                pmids.append(pmid)
                        dataset['pmids'] = ','.join(pmids)
                        dataset.pop('citation')

                    if "includedInDataCatalog" in dataset:
                        dataset["includedInDataCatalog"] =  [dataset["includedInDataCatalog"]]

                    count += 1

                    dataset.pop("@id")

                    if "version" in dataset:
                        dataset.pop("version")

                    yield dataset

                else:
                    # here is where the non-exported metadata is translated
                    dataset["context"] = "http://schema.org"
                    dataset['doi'] = dataset.pop('global_id')
                    dataset['_id'] = dataset['doi'].replace("/","_").replace('doi:','Dataverse')
                    #dataset['dateModified'] = datetime.datetime.strptime(dataset['updatedAt'] , "%Y-%m-%d").date().isoformat()
                    # have to strip Z from time to get into correct format 
                    if 'updatedAt' in dataset:
                        dataset['dateModified'] = dataset['updatedAt'].strip("Z")
                        dataset.pop('updatedAt') 
                    if "type" in dataset:
                        dataset["@type"] = dataset.pop('type')
                    if "published_at" in dataset:
                        dataset["datePublished"]= dataset['published_at'].strip("Z")
                        dataset.pop('published_at')
                    if "subjects" in dataset:
                        if "keywords" in dataset:
                            dataset["keywords"] += dataset.pop('subjects')
                        else:
                            dataset["keywords"] = dataset.pop('subjects')    
                    if "createdAt" in dataset:
                        dataset["dateCreated"] = dataset['createdAt'].strip('Z')
                        dataset.pop('createdAt')
                    if "authors" in dataset:
                        dataset["author"]=[{'name': dataset.pop('authors')}]
                    
                    dataset.pop("citationHtml")
                    dataset.pop("identifier_of_dataverse")
                    dataset.pop("name_of_dataverse")
                    dataset.pop("storageIdentifier")
                    dataset.pop("fileCount")
                    dataset.pop("versionId")
                    dataset.pop("versionState")

                    yield dataset

            except Exception as error:
                logger.info(f"skipping record with error - {error}")

        process_time = time.process_time() - start_time
        logger.info(f"Completed parsing individual metadata, {count} records parsed in {process_time:.2f} seconds")
        logger.info("PARSE PROCESS COMPLETE.")
