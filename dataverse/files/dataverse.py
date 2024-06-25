import datetime
import json
import logging
import time
from html.parser import HTMLParser

import psutil
import requests
# from sql_database import NDEDatabase

import pprint


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")

MAX_RETRIES = 10
BASE_DELAY = 2  # 2 seconds
MAX_DELAY = 60  # 1 minute

class Dataverse(): #NDEDatabase):
    SQL_DB = "dataverse.db"
    EXPIRE = datetime.timedelta(days=90)

    DATAVERSE_SERVER = "https://dataverse.harvard.edu/api"
    EXPORT_URL = f"{DATAVERSE_SERVER}/datasets/export?exporter=schema.org"

    def run_dataverse_schema_export(self, gid, url, verbose=False):
        """
        Run the dataverse schema.org export on a harvard dataset
        Returns the schema.org json object
        """
        dv_schema_export = f"{self.EXPORT_URL}&persistentId={gid}"
        logger.info(f"Running export using dataverse api - {dv_schema_export}")

        retries = 0
        backoff_time = BASE_DELAY

        while retries <= MAX_RETRIES:
            try:
                req = requests.get(dv_schema_export)
                res = req.json()
                # odd case of a response with a status of ERROR
                if res.get("status") and res.get("status") == "ERROR":
                    document = self.extract_schema_json(url)
                    if document:
                        logger.info(f"successfully scrapped document, {url}")
                        return document
                    else:
                        logger.info(
                            f"schema.org export failed on {dv_schema_export}, {url}, {req.status_code}, {res.get('status')}, {res.get('message')}"
                        )
                        logger.info("will attempt to yield original metadata for dataset....")
                        return False
                else:
                    return res
            except (requests.RequestException, json.decoder.JSONDecodeError) as e:
                retries += 1
                if retries > MAX_RETRIES:
                    logger.info(f"Failed to get {dv_schema_export} after {MAX_RETRIES} attempts due to {e}")
                    return False
                logger.info(f"Request failed due to {e}, retrying in {backoff_time} seconds...")
                time.sleep(backoff_time)
                backoff_time = min(MAX_DELAY, backoff_time * 2)  # double the wait time, but cap at MAX_DELAY

    def compile_paginated_data(self, query_endpoint, per_page=400, start=0):
        logger.info(f"Compiling paginated data from {query_endpoint}")
        continue_paging = True

        while continue_paging:
            pager = per_page
            response = None
            retries = 0
            backoff_time = BASE_DELAY
            while retries <= MAX_RETRIES and not response:
                try:
                    url = f"{query_endpoint}&per_page={pager}&start={start}"
                    logger.info(f"Requesting {url}")
                    req = requests.get(url)
                    req.raise_for_status()  # Raise an exception for HTTP errors
                    response = req.json()
                    retries = 0  # Reset the retry counter after a successful request
                except (requests.RequestException, ValueError) as e:
                    pager = max(1, pager//2)
                    logger.error(f"Error accessing {url}: {str(e)}")
                    retries += 1
                    if retries > MAX_RETRIES:
                        logger.error(f"Max retries reached for {url}. Skipping.")
                        return
                    logger.info(f"Retrying in {backoff_time} seconds...")
                    time.sleep(backoff_time)
                    backoff_time = min(MAX_DELAY, backoff_time * 2)  # double the wait time, but cap at MAX_DELAY
            if response:
                try:
                    total = response.get("data").get("total_count")
                    page_data = response.get("data").get("items")
                    logger.info(f"Retrieved {len(page_data)} items")
                    start += pager
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

    def log_memory_usage(self):
        process = psutil.Process()
        memory_use = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        logging.info(f"Current memory usage: {memory_use:.2f} MB")

    def load_cache(self):
        """
        # For Dev
        # Testing ~ set initial_start_page to view different sources
        # initial_start_page=0 -- for harvard sources
        #   - https://dataverse.harvard.edu/api/search?q=*&type=dataset&per_page=100&start=0
        # initial_start_page=94000 -- for external sources (this may increase as more data is added to dataverse)
        #   - https://dataverse.harvard.edu/api/search?q=*&type=dataset&per_page=5&start=94000
        # initial_start_page=155000 -- for handle urls
        #   - https://dataverse.harvard.edu/api/search?q=*&type=dataset&per_page=5&start=155000
        """
        logger.info("Runing load_cache process....")

        # Initial memory usage log
        self.log_memory_usage()

        # Adjust the start parameter to skip the desired number of records
        initial_page_start = 0
        records_processed = 0
        handle_url_ct=0
        schemas_gathered_ct=0
        sleep_time = 5 #seconds

        query_endpoint = "https://dataverse.harvard.edu/api/search?q=*&type=dataset"

        # Iterate through the paginated data starting from the adjusted initial position
        for global_id, data_url, data_page in self.compile_paginated_data(query_endpoint, per_page=400, start=initial_page_start):
            records_processed += 1
            if "https://hdl.handle.net/" in data_url:
                if "https://hdl.handle.net/1902.4" in data_url:
                    # https://hdl.handle.net/1902.4 : empty/dead-links (404 error) -- maybe handle 404 instead?
                    pass
                elif "https://hdl.handle.net/" in data_url:
                    handle_url_ct += 1
                    if data_page:
                        yield (data_url, json.dumps(data_page))
                        schemas_gathered_ct += 1
                
            elif global_id.startswith("doi:10.7910"):
                # case: harvard data - doi:10.7910
                # run the built-in dataverse schema export on harvard dataverse sources
                schema_record =  self.run_dataverse_schema_export(global_id, data_url)
                if schema_record and isinstance(schema_record, dict) and schema_record.get("@id"):
                    logger.info(f"schema export passed on {data_url}")
                    yield (schema_record["@id"], json.dumps(schema_record))
                    schemas_gathered_ct += 1

            else:
                # case: outside data (non-harvard registered--therefore not standardized)
                if data_page:
                    # yield data available through dataverse api, not extracting from url
                    yield(data_url, json.dumps(data_page))
                    schemas_gathered_ct += 1

            if records_processed % 1000 == 0:
                logger.info(f"Processed {records_processed} datasets, going to sleep for {sleep_time} seconds to manage load...")
                time.sleep(sleep_time)  # Sleep for 5 seconds every 1000 datasets

            # Optional - log memory usage periodically, e.g., every 100 records
            if records_processed % 100 == 0:
                self.log_memory_usage()

        # Final memory usage log
        self.log_memory_usage()
        logger.info(f"Processed {records_processed} datasets and {schemas_gathered_ct} schemas.")

    def parse(self, records):
        start_time = time.process_time()
        parse_ct = 0
        logger.info("Starting metadata parser...")
        # rec = ('doi:10.18738/T8/YJMLKO', '{"name": "ChIP-seq peak calls for epigenetic marks in GBM tumors", "type": "dataset", "url": "https://doi.org/10.18738/T8/YJMLKO", "global_id": "doi:10.18738/T8/YJMLKO", "description": "MACS2 narrowPeak files from ChIP-seq experiments for 11 primary GBM tumors, each targeting CTCF transcription factor marks and H3K27Ac, H3K27Me3, H3K4Me1, H3K4Me3, H3K9Ac, and H3K9Me3 histone modifications. See Methods section of doi:10.1158/0008-5472.CAN-17-1724 for more information.", "published_at": "2018-11-05T05:17:42Z", "publisher": "Texas Data Repository Harvested Dataverse", "citationHtml": "Battenhouse, Anna; Hall, Amelia Weber, 2018, \\"ChIP-seq peak calls for epigenetic marks in GBM tumors\\", <a href=\\"https://doi.org/10.18738/T8/YJMLKO\\" target=\\"_blank\\">https://doi.org/10.18738/T8/YJMLKO</a>, Texas Data Repository Dataverse", "identifier_of_dataverse": "tdr_harvested", "name_of_dataverse": "Texas Data Repository Harvested Dataverse", "citation": "Battenhouse, Anna; Hall, Amelia Weber, 2018, \\"ChIP-seq peak calls for epigenetic marks in GBM tumors\\", https://doi.org/10.18738/T8/YJMLKO, Texas Data Repository Dataverse", "storageIdentifier": "s3://10.18738/T8/YJMLKO", "keywords": ["Medicine, Health and Life Sciences", "glioblastoma", "bivalent", "enhancer", "epigenetic", "histone modification"], "subjects": [], "fileCount": 84, "versionId": 146549, "versionState": "RELEASED", "createdAt": "2018-11-05T05:17:42Z", "updatedAt": "2018-11-05T05:17:42Z", "contacts": [{"name": "", "affiliation": ""}], "authors": ["Battenhouse, Anna", "Hall, Amelia Weber"]}')
        # records = [rec]
        # logger.info(f"Starting parsing of {len(records)} records...")
        for record in records:
            # logger.info(f"parsing record {record[0]}")
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
                                if "affiliation" in data_dict.keys() and isinstance(data_dict["affiliation"], str):
                                    data_dict["affiliation"] = {"name": data_dict.pop("affiliation")}
                        else:
                            dataset.pop("creator")
                    else:
                        for data_dict in dataset["author"]:
                            if "affiliation" in data_dict.keys() and isinstance(data_dict["affiliation"], str):
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
                    else:
                        dataset["includedInDataCatalog"] = {
                            "@type": "dataset",
                            "name": "Harvard Dataverse",
                            "url": "https://dataverse.harvard.edu/",
                            "versionDate": datetime.datetime.today().strftime("%Y-%m-%d"),
                        }

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

                    keys_to_remove = ['dataSources', 'geographicCoverage', 'majorVersion', 'producers',
                                    'publications', 'relatedMaterial', 'topicCategory', 'publisher', 'citationHtml',
                                    'storageIdentifier', 'fileCount', 'versionId', 'versionState', 'contacts']

                    for key in keys_to_remove:
                        print(f"removing key: {key}")
                        dataset.pop(key, None)

                    yield dataset
                    parse_ct += 1

            except Exception as error:
                logger.info(f"skipping with error - {error}, record id - {record[0]}")

        process_time = time.process_time() - start_time
        logger.info(f"Completed parsing individual metadata, {parse_ct} records parsed in {process_time:.2f} seconds")
        logger.info("Document parsing complete")
