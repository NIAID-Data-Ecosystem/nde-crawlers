import datetime
import json
import logging
import time
from html.parser import HTMLParser

import requests
from sql_database import NDEDatabase

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")

MAX_RETRIES = 5
BASE_DELAY = 2  # 2 seconds
MAX_DELAY = 120  # 2 minutes


class Dataverse(NDEDatabase):
    SQL_DB = "dataverse.db"
    EXPIRE = datetime.timedelta(days=90)
    NO_CACHE = True

    DATAVERSE_SERVER = "https://dataverse.harvard.edu/api"
    DATA_URL = f"{DATAVERSE_SERVER}/search?q=data"
    EXPORT_URL = f"{DATAVERSE_SERVER}/datasets/export?exporter=schema.org"

    def scrape_schema_representation(self, url):
        """
        When the schema.org export of the dataset fails,
        this will grab it from the URL by looking for <script type="application/ld+json">
        """

        class SchemaScraper(HTMLParser):
            def __init__(self):
                super().__init__()
                self.readingSchema = False
                self.schema = None

            def handle_starttag(self, tag, attrs):
                attrs = dict(attrs)
                if tag == "script" and attrs.get("type") == "application/ld+json":
                    self.readingSchema = True

            def handle_data(self, data):
                if self.readingSchema:
                    self.schema = data
                    self.readingSchema = False

        retries = 0
        backoff_time = BASE_DELAY
        logger.info(f"Scraping {url} for schema representation")

        while retries <= MAX_RETRIES:
            try:
                req = requests.get(url)
                if not req.ok:
                    logger.info(f"unable to retrieve {url}, status code: {req.status_code}")
                    return False
                parser = SchemaScraper()
                parser.feed(req.text)
                if parser.schema:
                    return parser.schema
                return False
            except requests.RequestException as requestException:
                retries += 1
                if retries > MAX_RETRIES:
                    logger.info(f"Failed to scrape {url} after {MAX_RETRIES} attempts due to {requestException}")
                    return False
                logger.info(f"Scraping failed due to {requestException}, retrying in {backoff_time} seconds...")
                time.sleep(backoff_time)
                backoff_time = min(MAX_DELAY, backoff_time * 2)  # double the wait time, but cap at MAX_DELAY

    def get_schema_document(self, gid, url, verbose=False):
        schema_export_url = f"{self.EXPORT_URL}&persistentId={gid}"

        retries = 0
        backoff_time = BASE_DELAY

        while retries <= MAX_RETRIES:
            try:
                req = requests.get(schema_export_url)
                res = req.json()
                if res.get("status") and res.get("status") == "ERROR":
                    document = self.scrape_schema_representation(url)
                    if document:
                        logger.info(f"successfully scrapped document, {url}")
                        return document
                    else:
                        logger.info(
                            f"schema.org export failed on {schema_export_url}, {url}, {req.status_code}, {res.get('status')}, {res.get('message')}"
                        )
                        logger.info("will attempt to yield original metadata for dataset....")
                        return False
                else:
                    return res
            except (requests.RequestException, json.decoder.JSONDecodeError) as e:
                retries += 1
                if retries > MAX_RETRIES:
                    logger.info(f"Failed to get {schema_export_url} after {MAX_RETRIES} attempts due to {e}")
                    return False
                logger.info(f"Request failed due to {e}, retrying in {backoff_time} seconds...")
                time.sleep(backoff_time)
                backoff_time = min(MAX_DELAY, backoff_time * 2)  # double the wait time, but cap at MAX_DELAY

    def compile_paginated_data(self, query_endpoint, per_page=1000, verbose=False):
        logger.info(f"Compiling paginated data from {query_endpoint}")
        continue_paging = True
        start = 0
        retries = 0

        while continue_paging:
            url = f"{query_endpoint}&per_page={per_page}&start={start}"
            logger.info(f"Requesting {url}")
            try:
                req = requests.get(url)
                req.raise_for_status()  # Raise an exception for HTTP errors
                response = req.json()
                retries = 0  # Reset the retry counter after a successful request
            except (requests.RequestException, ValueError) as e:
                logger.error(f"Error accessing {url}: {str(e)}")
                retries += 1
                if retries > MAX_RETRIES:
                    logger.error(f"Max retries reached for {url}. Skipping.")
                    return
                delay = min(BASE_DELAY * (2**retries), MAX_DELAY)
                logger.info(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            if response:
                try:
                    total = response.get("data").get("total_count")
                    page_data = response.get("data").get("items")
                    logger.info(f"Retrieved {len(page_data)} items")
                    start += per_page
                    for page in page_data:
                        # dataset pages use "global_id" key
                        if "global_id" in page.keys():
                            yield (page["global_id"], page["url"], page)
                        # dataverse pages have "identifier" key
                        elif "identifier" in page.keys():
                            yield (page["identifier"], page["url"], page)  # EXTRACT ID HERE
                    # data_pages.extend(page_data)
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
        dataset_ids = []  # set([None])
        # datasets    = []
        handle_ct = 0
        dataset_ct = 0
        dataverse_ct = 0
        schema_ct = 0
        skipped_data = 0
        # dataverse type: data --> (identifier, data_dict)
        logger.info("Starting extraction of type 'dataverse' data....")
        for data in self.compile_paginated_data(self.DATA_URL + "&type=dataverse"):
            logger.info(f"Extracting dataverse data from {data[1]}")
            dataverse_query = f"{self.DATAVERSE_SERVER}/search?q=*&type=dataset&type=file&subtree={data[0]}"
            for dataverse_data in self.compile_paginated_data(dataverse_query):
                logger.info(f"Extracting dataset data from {dataverse_data[1]}")
                # unique case: if url is handle.net/ data scrape the schema else proceed to export schema document
                if "https://hdl.handle.net/" in dataverse_data[1]:
                    logger.info(f"scraping schema for {dataverse_data[1]}")
                    schema_record = self.scrape_schema_representation(dataverse_data[1])
                    handle_ct += 1
                else:
                    logger.info(f"exporting schema for {dataverse_data[1]}")
                    schema_record = self.get_schema_document(dataverse_data[0], dataverse_data[1])
                # if the schema extraction was successful yield the exported schema document
                # else handle unique cases
                if schema_record:
                    logger.info(f"schema export successful for {dataverse_data[1]}")
                    dataset_ids.append(dataverse_data[0])
                    dataverse_ct += 1
                    schema_ct += 1
                    yield (schema_record["@id"], json.dumps(schema_record))
                    logger.info(f"caching original dataset for {dataverse_data[1]}")
                elif schema_record == False:
                    logger.info(f"schema export failed for {dataverse_data[1]}")
                    if "https://hdl.handle.net/" in dataverse_data[1]:
                        logger.info(f"will attempt to yield original metadata for dataset....")
                        pass
                    elif dataverse_data[2]:
                        logger.info(f"caching original dataset for {dataverse_data[1]}")
                        yield (dataverse_data[0], json.dumps(dataverse_data[2]))
                    else:
                        skipped_data += 1
                        logger.info(
                            f"error retrieving schema export for document with id, {dataverse_data[0]}, skipping dataset"
                        )

        logger.info(
            f"Extraction of dataverse data complete, {dataverse_ct} dataset schemas were extracted, and {len(set(dataset_ids))} unique ids total."
        )
        logger.info("Starting extraction of type 'dataset' data....")
        # # dataset type: data --> (global_id, url)
        for data in self.compile_paginated_data(self.DATA_URL + "&type=dataset"):
            logger.info(f"Extracting dataset data from {data[1]}")
            if data[0] not in dataset_ids:
                dataset_ids.append(data[0])
                logger.info(f"exporting schema for {data[1]}")
                if "https://hdl.handle.net/" in data[1]:
                    logger.info(f"scraping schema for {data[1]}")
                    self.scrape_schema_representation(data[1])
                    handle_ct += 1
                else:
                    schema_record = self.get_schema_document(data[0], data[1])
                    if schema_record:
                        logger.info(f"schema export successful for {data[1]}")
                        schema_ct += 1
                        dataset_ct += 1
                        # print(json.dumps(schema_record, indent=4))
                        yield (schema_record["@id"], json.dumps(schema_record))
                    elif schema_record == False:
                        logger.info(f"schema export failed for {data[1]}")
                        if "https://hdl.handle.net/" in data[1]:
                            logger.info(f"will attempt to yield original metadata for dataset....")
                            pass
                        elif data[2]:
                            logger.info(f"caching original dataset for {data[1]}")
                            yield (data[0], json.dumps(data[2]))
                        else:
                            skipped_data += 1
                            logger.info(
                                f"error retrieving schema export for document with id, {data[0]}, skipping dataset"
                            )
        process_time = time.process_time() - start_time
        logger.info(
            f"Extraction of dataset data complete, {dataset_ct} dataset schemas were extracted, {len(set(dataset_ids))} unique ids total."
        )
        logger.info(f"{schema_ct} schemas extracted in {process_time:.2f} seconds.")
        logger.info("LOAD CACHE PROCESS COMPLETE.")

    def parse(self, records):
        start_time = time.process_time()
        parse_ct = 0
        logger.info("Starting metadata parser...")
        # rec = ('doi:10.18738/T8/YJMLKO', '{"name": "ChIP-seq peak calls for epigenetic marks in GBM tumors", "type": "dataset", "url": "https://doi.org/10.18738/T8/YJMLKO", "global_id": "doi:10.18738/T8/YJMLKO", "description": "MACS2 narrowPeak files from ChIP-seq experiments for 11 primary GBM tumors, each targeting CTCF transcription factor marks and H3K27Ac, H3K27Me3, H3K4Me1, H3K4Me3, H3K9Ac, and H3K9Me3 histone modifications. See Methods section of doi:10.1158/0008-5472.CAN-17-1724 for more information.", "published_at": "2018-11-05T05:17:42Z", "publisher": "Texas Data Repository Harvested Dataverse", "citationHtml": "Battenhouse, Anna; Hall, Amelia Weber, 2018, \\"ChIP-seq peak calls for epigenetic marks in GBM tumors\\", <a href=\\"https://doi.org/10.18738/T8/YJMLKO\\" target=\\"_blank\\">https://doi.org/10.18738/T8/YJMLKO</a>, Texas Data Repository Dataverse", "identifier_of_dataverse": "tdr_harvested", "name_of_dataverse": "Texas Data Repository Harvested Dataverse", "citation": "Battenhouse, Anna; Hall, Amelia Weber, 2018, \\"ChIP-seq peak calls for epigenetic marks in GBM tumors\\", https://doi.org/10.18738/T8/YJMLKO, Texas Data Repository Dataverse", "storageIdentifier": "s3://10.18738/T8/YJMLKO", "keywords": ["Medicine, Health and Life Sciences", "glioblastoma", "bivalent", "enhancer", "epigenetic", "histone modification"], "subjects": [], "fileCount": 84, "versionId": 146549, "versionState": "RELEASED", "createdAt": "2018-11-05T05:17:42Z", "updatedAt": "2018-11-05T05:17:42Z", "contacts": [{"name": "", "affiliation": ""}], "authors": ["Battenhouse, Anna", "Hall, Amelia Weber"]}')
        # records = [rec]
        logger.info(f"Starting parsing of {len(records)} records...")
        for record in records:
            logger.info(f"parsing record {record[0]}")
            try:
                dataset = json.loads(record[1])
                # here is where the exported schema.org data is parsed
                if "@context" in dataset or "@type" in dataset:
                    dataset = json.loads(record[1])
                    dataset["url"] = dataset["identifier"]
                    dataset["doi"] = dataset.pop("identifier")
                    dataset["identifier"] = dataset["doi"].strip("https://doi.org")
                    dataset["_id"] = dataset["doi"].replace("https://doi.org", "Dataverse").replace("/", "_")
                    dataset["dateModified"] = (
                        datetime.datetime.strptime(dataset["dateModified"], "%Y-%m-%d").date().isoformat()
                    )

                    if "datePublished" in dataset:
                        dataset["datePublished"] = (
                            datetime.datetime.strptime(dataset["datePublished"], "%Y-%m-%d").date().isoformat()
                        )

                    if dataset["author"] is None:
                        dataset.pop("author")
                        if dataset["creator"]:
                            dataset["author"] = dataset.pop("creator")
                            for data_dict in dataset["author"]:
                                if "affiliation" in data_dict.keys():
                                    data_dict["affiliation"] = {"name": data_dict.pop("affiliation")}
                        else:
                            dataset.pop("creator")
                    else:
                        for data_dict in dataset["author"]:
                            if "affiliation" in data_dict.keys():
                                data_dict["affiliation"] = {"name": data_dict.pop("affiliation")}
                        dataset.pop("creator")

                    if dataset["publisher"] is None:
                        dataset.pop("publisher")
                        if dataset["provider"]:
                            dataset["sdPublisher"] = [dataset.pop("provider")]
                        else:
                            dataset["sdPublisher"] = [dataset.pop("provider")]
                    else:
                        dataset["sdPublisher"] = [dataset.pop("publisher")]
                        dataset.pop("provider")

                    if "funder" in dataset:
                        dataset["funding"] = {"funder": dataset.pop("funder")}
                        for data_dict in dataset["funding"]["funder"]:
                            data_dict.pop("@type")

                    if "description" in dataset:
                        dataset["description"] = dataset.pop("description")[0]

                    if "keywords" in dataset:
                        if dataset["keywords"]:
                            dataset["keywords"] = dataset.pop("keywords")[0]

                    if dataset["license"]:
                        if type(dataset["license"]) is str:
                            pass
                        elif type(dataset["license"]) is dict:
                            if "url" in dataset["license"].keys():
                                dataset["license"] = dataset["license"]["url"]
                            elif "text" in dataset["license"].keys():
                                dataset["license"] = dataset["license"]["text"]
                            else:
                                dataset.pop("license")

                    if "temporalCoverage" in dataset and dataset["temporalCoverage"]:
                        dataset["temporalCoverage"] = [
                            {"temporalInterval": {"duration": dataset["temporalCoverage"][0]}}
                        ]

                    if "spatialCoverage" in dataset:
                        dataset["spatialCoverage"] = [{"name": dataset.pop("spatialCoverage")[0]}]

                    if "distribution" in dataset:
                        for data_dict in dataset["distribution"]:
                            if "@id" in data_dict:
                                data_dict["@id"] = data_dict["@id"].strip("https://doi.org/")
                            if "fileFormat" in data_dict:
                                data_dict["encodingFormat"] = data_dict.pop("fileFormat")
                            if "identifier" in data_dict:
                                data_dict["contentUrl"] = data_dict.pop("identifier")
                            if "contentSize" in data_dict:
                                data_dict.pop("contentSize")
                            if "description" in data_dict and "https://" in data_dict["description"]:
                                data_dict.pop("description")

                    if "citation" in dataset:
                        cit_list = []
                        for data_dict in dataset["citation"]:
                            if "text" in data_dict:
                                cit_list.append({"citation": data_dict["text"]})
                        if cit_list:
                            dataset["citation"] = cit_list
                        else:
                            dataset.pop("citation")

                    if "includedInDataCatalog" in dataset:
                        dataset["includedInDataCatalog"] = dataset["includedInDataCatalog"]
                        dataset["includedInDataCatalog"]["versionDate"] = datetime.datetime.today().strftime("%Y-%m-%d")

                    dataset.pop("@id")

                    if "version" in dataset:
                        dataset.pop("version")

                    yield dataset
                    parse_ct += 1

                else:
                    # here is where the non-exported metadata is parsed
                    # add schema source info variables
                    dataset["@context"] = "http://schema.org"
                    if "type" in dataset:
                        dataset["@type"] = dataset.pop("type")
                    else:
                        dataset["@type"] = "dataset"
                    dataset["doi"] = dataset.pop("global_id")
                    dataset["identifier"] = dataset["doi"].strip("https://doi.org")
                    dataset["_id"] = dataset["doi"].replace("doi:", "Dataverse_").replace("/", "_")
                    dataset["includedInDataCatalog"] = {
                        "@type": "dataset",
                        "name": "Harvard Dataverse",
                        "url": "https://dataverse.harvard.edu/",
                        "versionDate": datetime.datetime.today().strftime("%Y-%m-%d"),
                    }

                    # have to strip Z from time to get into correct format
                    if "updatedAt" in dataset:
                        dataset["dateModified"] = dataset["updatedAt"].strip("Z")
                        dataset.pop("updatedAt")
                        dataset["dateModified"] = datetime.datetime.strptime(
                            dataset["dateModified"], "%Y-%m-%dT%H:%M:%S"
                        ).strftime("%Y-%m-%d")

                    if "published_at" in dataset:
                        dataset["datePublished"] = dataset["published_at"].strip("Z")
                        dataset.pop("published_at")
                        dataset["datePublished"] = datetime.datetime.strptime(
                            dataset["datePublished"], "%Y-%m-%dT%H:%M:%S"
                        ).strftime("%Y-%m-%d")

                    if "identifier_of_dataverse" in dataset:
                        dataset["sdPublisher"] = {
                            "identifier": dataset.pop("identifier_of_dataverse"),
                            "name": dataset.pop("name_of_dataverse"),
                        }

                    if "subjects" in dataset:
                        if not dataset["subjects"]:
                            pass
                        elif len(dataset["subjects"]) > 1:
                            dataset["topicCategory"] = []
                            for subject in dataset["subjects"]:
                                dataset["topicCategory"].append({"description": subject})
                        else:
                            dataset["topicCategory"] = {"description": dataset["subjects"][0]}
                        dataset.pop("subjects")

                    if "keywords" in dataset:
                        dataset["keywords"] = ",".join(dataset["keywords"])

                    if "createdAt" in dataset:
                        dataset["dateCreated"] = dataset["createdAt"].strip("Z")
                        dataset.pop("createdAt")
                        dataset["dateCreated"] = datetime.datetime.strptime(
                            dataset["dateCreated"], "%Y-%m-%dT%H:%M:%S"
                        ).strftime("%Y-%m-%d")

                    if "authors" in dataset:
                        # print(len(dataset['authors']))
                        if not dataset["authors"]:
                            # print("[case1]")
                            dataset.pop("authors")
                        elif len(dataset["authors"]) > 1:
                            # print("[case2]")
                            dataset["author"] = []
                            for author in dataset["authors"]:
                                author_dict = {"name": author}
                                dataset["author"].append(author_dict)
                            dataset.pop("authors")
                        else:
                            # print("[case3]")
                            dataset["author"] = {"name": dataset.pop("authors")[0]}
                    if "citation" in dataset:
                        dataset["citation"] = {"citation": dataset.pop("citation")}

                    dataset.pop("publisher")
                    dataset.pop("citationHtml")
                    dataset.pop("storageIdentifier")
                    dataset.pop("fileCount")
                    dataset.pop("versionId")
                    dataset.pop("versionState")
                    dataset.pop("contacts")

                    yield dataset
                    parse_ct += 1

            except Exception as error:
                logger.info(f"skipping with error - {error}, record id - {record[0]} \n{record}")

        process_time = time.process_time() - start_time
        logger.info(f"Completed parsing individual metadata, {parse_ct} records parsed in {process_time:.2f} seconds")
        logger.info("Document parsing complete")
