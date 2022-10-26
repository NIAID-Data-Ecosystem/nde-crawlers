import datetime
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

    def get_all_datasets_from_dataverses(self):
        logger.info("finding all dataverses data available")
        #dataverses = self.find_relevant_dataverses()
        dataverse_endpoint = "https://dataverse.harvard.edu/api/search?q=*&type=dataverse"
        datasets = []
        logger.info("extracting dataverse datasets...")
        for dv_data in self.compile_paginated_data(dataverse_endpoint):
            dv_id = dv_data['identifier']
            dv_id_query = f"{self.DATAVERSE_SERVER}/search?q=*&type=dataset&subtree={dv_id}"
            datasets_and_files = self.compile_paginated_data(dv_id_query)
            datasets.extend(datasets_and_files)
        logger.info(f"{len(datasets)} dataverse datasets extracted")
        return datasets


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
                    print('herexs')

            def handle_data(self, data):
                if self.readingSchema:
                    self.schema = data
                    self.readingSchema = False
        try:
            req = requests.get(url)
        except Exception as requestException:
            logger.error(f"Failed to get {url} due to {requestException}")
            return False
        if not req.ok:
            logger.error(f"failed to get {url}")
            return False
        parser = SchemaScraper()
        parser.feed(req.text)
        if parser.schema:
            return parser.schema
        return False


    def get_schema_document(self, gid, url, verbose=False):
        if "https://hdl.handle.net/" in url:
            return False
        else:
            schema_export_url = f"{self.EXPORT_URL}&persistentId={gid}"
            if verbose == True:
                logger.info(f"exporting schema.org document for {gid}")
            try:
                req = requests.get(schema_export_url)
            except Exception as requestException:
                logger.error(f"Failed to get {schema_export_url} due to {requestException}")
                return False
            try:
                res = req.json()
            except json.decoder.JSONDecodeError:
                return False
            if res.get('status') and res.get('status') == 'ERROR':
                logger.warning(f"document export failed, scraping {url} instead")
                document = self.scrape_schema_representation(url)
                if document:
                    logger.info("successfully scrapped document")
                    return document
                else:
                    return False
            else:
                # success, response is the document
                logger.info(f"schema export successful {gid}")
                return res


    def compile_paginated_data(self, query_endpoint, per_page=50, verbose=False):
        """ Extract Data By Page
        pages through data, compiling all response['data']['items']
        and returning them.
        per_page max is 1000
        """
        continue_paging = True
        start   = 0
        retries = 0
        data_pages = []
        
        while continue_paging:
            url = f"{query_endpoint}&per_page={per_page}&start={start}"
            if verbose == True:
                logger.info(f"querying endpoint {url}")
            try:
                req = requests.get(url)
            except Exception as requestException:
                logger.error(f"Failed to get {url} due to {requestException}")
                if retries > 5:
                    logger.error(f"Failed too many times on endpoint {url}")
                retries += 1
                # after 1 retry, limit per-page to 200, after 2, limit to 50
                per_page = 200 if retries == 1 else 50
            try:
                response = req.json()
            except ValueError:
                logger.error(f"Failed to get a JSON response from {url}") 
            total = response.get('data').get('total_count')
            page_data = response.get('data').get('items')
            start += per_page
            data_pages.extend(page_data)
            continue_paging = total and start < 200 #total
        return data_pages


    def fetch_datasets(self):
        """Retrives the raw data using a sickle request and formats so dump can store it into the cache table
        Returns:
            A tuple (_id, data formatted as a json string)
        grabs all datasets and files related to QUERIES both by querying
        and by grabbing everything in related dataverses
        extracts their global_id, which in this case is a DOI
        returns a dictionary mapping global_id -> dataset
        """
        dataset_ids = set([None])
        datasets    = []

        logger.info("starting extraction of all search data")
        new_data = self.compile_paginated_data(self.DATA_URL)
        unique_data = [dataset for dataset in new_data if dataset.get('global_id') not in dataset_ids]
        datasets.extend(unique_data)
        logger.info(f"{len(unique_data)} unique search dataset records extracted")

        dataset_ids |= set([dataset.get('global_id') for dataset in unique_data])
        data_for_gid = {d.get('global_id'): d for d in datasets}

        logger.info("starting extraction of all dataverse data")
        dataverse_data = self.get_all_datasets_from_dataverses()
        dataverse_data_for_gid = {d.get('global_id'): d for d in dataverse_data}

        total_datasets = {
            **data_for_gid,
            **dataverse_data_for_gid
        }

        try:
            total_datasets.pop('')
        except KeyError:
            pass
        return total_datasets


    def load_cache(self):
        logger.info("Loading cache....")
        start_time = time.process_time()
        datasets = self.fetch_datasets()
        record_ct = 0
        # get the schemas for the data
        for gid, data in datasets.items():
            schema_record = self.get_schema_document(gid, data.get('url'))
            if schema_record:
                record_ct += 1
                yield (schema_record['@id'], json.dumps(schema_record))
        process_time = time.process_time() - start_time
        logger.info(f"load cache successful, {record_ct} datasets loaded in {process_time:.2f} seconds")


    def parse(self, records):
        
        start_time = time.process_time()
        logger.info(f"Starting metadata parser on {len(records)} records...")
        count = 0

        for record in records:
            dataset=json.loads(record[1])

            dataset['dateModified'] = datetime.datetime.strptime(dataset['dateModified'] , "%Y-%m-%d").date().isoformat()
            dataset['datePublished'] = datetime.datetime.strptime(dataset['datePublished'] , "%Y-%m-%d").date().isoformat()

            dataset['url'] = dataset.pop('identifier')
            dataset['doi'] = dataset['url'].replace("/","_").replace('https:__doi.org','Dataverse')

            if dataset['author'] is None:
                dataset.pop('author')
                if dataset['creator']:
                    dataset['author'] = dataset.pop('creator')
                else:
                    dataset.pop('creator')
            else:
                dataset.pop('creator')

            if dataset['publisher'] is None:
                dataset.pop('publisher')
                if dataset['provider']:
                    dataset['publisher'] = dataset.pop('provider')
                else:
                    dataset.pop('provider')
            else:
                dataset.pop('provider')
            
            if "funder" in dataset.keys():
                dataset['funding'] = {"funder": dataset.pop('funder')}

            count += 1

            yield dataset
        process_time = time.process_time() - start_time
        logger.info(f"completed parsing individual metadata, {count} records parsed in {process_time:.2f} seconds")
