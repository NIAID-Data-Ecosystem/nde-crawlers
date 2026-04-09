# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh


import json
import logging
import re

import scrapy

logger = logging.getLogger("nde-logger")


# may take some time to start up as getting the ids takes a while
class TychoSpider(scrapy.Spider):
    name = "tycho"

    custom_settings = {
        "ITEM_PIPELINES": {
            "pipeline.TychoItemProcessorPipeline": 100,
            "ndjson.NDJsonWriterPipeline": 999,
        }
    }

    start_urls = ["https://www.tycho.pitt.edu/data/#datasets"]

    def parse(self, response):

        base_url = "https://www.tycho.pitt.edu"
        # Extract dataset links
        covid_table = response.xpath('//*[@id="covid_dataset-table"]/tbody/tr')
        us_table = response.xpath('//*[@id="us_dataset-table"]/tbody/tr')
        dengue_table = response.xpath('//*[@id="dengue_dataset-table"]/tbody/tr')

        for table in [covid_table, us_table, dengue_table]:
            for row in table:
                dataset_link = row.xpath("@data-href").get()
                country = row.xpath("td[1]//text()").get()
                if dataset_link:
                    yield response.follow(
                        base_url + dataset_link,
                        self.parse_dataset,
                        meta={"country": country},
                    )

    def parse_dataset(self, response):
        logger.info(f"Processing dataset: {response.url}")
        country = response.meta.get("country")
        zenodo_link = response.xpath('//a[contains(text(), "Download data and metadata from Zenodo")]/@href').get()

        # Try new format: metadata embedded in <script id="download-dats-json">
        metadata_script = response.xpath('//script[@id="download-dats-json"]/text()').get()
        if metadata_script:
            metadata = json.loads(metadata_script.strip())
        else:
            # Fall back to old format: metadata in onclick attribute
            metadata_download = response.xpath(
                '//div[@class="service"]//a[contains(text(), "Download metadata in DATS JSON format")]/@onclick'
            ).get()
            if not metadata_download:
                logger.warning(f"No metadata download link found for {response.url}, skipping.")
                return
            match = re.search(r"\{.*\}", metadata_download)
            if not match:
                logger.warning(f"Could not parse metadata from onclick for {response.url}, skipping.")
                return
            metadata = match.group(0)
            metadata = metadata.encode("utf-8").decode("unicode_escape")
            metadata = json.loads(metadata)

        logger.info(f"Identifier: {metadata.get('identifier')}")
        logger.info(f"Zenodo link: {zenodo_link}")
        yield {
            "country": country,
            "zenodo_url": zenodo_link,
            "metadata": metadata,
        }
