# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh

"""
    NCBI Geo Datasource

    Scrape contents in HTML tables
    Not structured metadata

    https://www.ncbi.nlm.nih.gov/robots.txt
    Declared No Robots and crawl delay as of 02/15/2020.
    Defined long rules so the file content is not included here.

"""

import logging
from ftplib import FTP

import scrapy

logger = logging.getLogger("nde-logger")


# may take some time to start up as getting the ids takes a while
class NCBIGeoSpider(scrapy.Spider):

    name = "ncbi_geo"

    custom_settings = {
        "ITEM_PIPELINES": {
            "pipeline.GeoItemProcessorPipeline": 100,
            "ndjson.NDJsonWriterPipeline": 999,
        }
    }

    # gets the list of gse ids to request
    def get_ids(self):
        gse_li = []
        ftp_host = "ftp.ncbi.nlm.nih.gov"
        ftp = FTP(ftp_host)
        ftp.login()
        ftp.cwd("/geo/series/")
        folder_li = ftp.nlst()
        for i in folder_li:
            logger.info("Retriving gse id list from ftp...")
            ftp.cwd(i)
            gse_li = gse_li + ftp.nlst()
            ftp.cwd("..")
        ftp.close()
        logger.info("Completed gse id list.")
        return gse_li

    # THIS USED FOR TESTING/DEBUGGING
    # for small tests can use (start, end) = (1,20)
    # this should be the most recent assession (GSE) link: https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&display=1&zsort=acc
    # def start_requests(self):
    #     start = 1
    #     end = 200
    #     prefix = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE"
    #     for acc_id in range(start, end + 1):
    #         yield scrapy.Request(url=prefix + str(acc_id))

    def start_requests(self):
        # there may be a connection error sometimes. Retry up to 5 times if there is one.
        tries = 6
        for attempt in range(tries):
            try:
                ids = self.get_ids()
            except EOFError as e:
                if attempt < (tries - 1):
                    continue
                else:
                    raise e
            else:
                break
        prefix = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
        for acc_id in ids:
            yield scrapy.Request(url=prefix + str(acc_id))

    def parse(self, response):

        table = response.xpath(
            "/html/body/table/tr/td/table[6]/tr[3]/td[2]" "/table/tr/td/table/tr/td/table[2]/tr/td" "/table[1]/tr"
        )
        data = {}

        for node in table:
            # extract series id
            if node.attrib.get("bgcolor") == "#cccccc":
                data["_id"] = node.xpath(".//strong").attrib.get("id")
            # remove place holder lines
            elif len(node.xpath("./td")) == 2:
                if node.xpath("string(./td[1])").get().strip():
                    # Handle 'Experiment type' specifically (measurementTechnique)
                    if node.xpath("string(./td[1])").get().strip() == "Experiment type":
                        experiment_types = node.xpath("./td[2]//text()").getall()
                        experiment_types = [etype.strip() for etype in experiment_types if etype.strip()]
                        data[key] = experiment_types
                    # extract multi item entry
                    elif node.xpath("./td[2]").attrib.get("onmouseout"):
                        key = node.xpath("./td[1]/text()").get().split()[0]
                        data[key] = node.xpath("./td[2]//a/text()").getall()
                    # extract single item entry
                    else:
                        key = node.xpath("./td[1]/text()").get()
                        data[key] = node.xpath("string(./td[2])").get().strip().replace("\xa0", " ")

        # self.logger.info("Check if cached: " + " ".join(response.flags))

        yield data if data else None
