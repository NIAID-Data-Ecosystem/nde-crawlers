# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh

import json
import logging
from urllib.parse import urljoin

import scrapy

logger = logging.getLogger("nde-logger")


class CeirrSpider(scrapy.Spider):
    name = "ceirr"
    custom_settings = {
        "ITEM_PIPELINES": {
            "pipeline.CeirrItemProcessorPipeline": 100,
            "ndjson.NDJsonWriterPipeline": 999,
        }
    }

    SOURCE_URL = "https://www.ceirr-network.org/resources/reagents"
    DEFAULT_DATA_URL = "https://www.ceirr-network.org/api/resources/reagents"
    PAGE_SIZE = 500
    SOURCE_NAME = "Centers of Excellence for Influenza Research and Response (CEIRR) Resources"

    start_urls = [SOURCE_URL]

    def parse(self, response):
        next_data_text = response.xpath('//script[@id="__NEXT_DATA__"]/text()').get()
        if not next_data_text:
            raise RuntimeError("Could not find CEIRR Next.js page data.")

        page_data = json.loads(next_data_text)
        page_props = page_data.get("props", {}).get("pageProps", {})
        data_url = urljoin(response.url, page_props.get("dataUrl") or self.DEFAULT_DATA_URL)
        reagents = page_props.get("reagents") or []
        abbreviations = page_props.get("abbreviations") or {}

        if not reagents:
            raise RuntimeError("Could not find CEIRR reagent card definitions.")

        logger.info("Found %s CEIRR reagent categories.", len(reagents))
        for reagent in reagents:
            filters = reagent.get("filter") or []
            if not filters:
                logger.warning("Skipping CEIRR reagent category without filters: %s", reagent.get("name"))
                continue

            body = self.build_body(filters)
            yield self.post_reagent_page(
                data_url=data_url,
                body=body,
                reagent=reagent,
                abbreviations=abbreviations,
                source_url=response.url,
            )

    def build_body(self, filters, pagination_key=""):
        return {
            "keyword": "",
            "pagination_key": pagination_key,
            "page_size": self.PAGE_SIZE,
            "filter": filters,
        }

    def post_reagent_page(self, data_url, body, reagent, abbreviations, source_url):
        return scrapy.Request(
            url=data_url,
            method="POST",
            body=json.dumps(body),
            headers={
                "Accept": "application/json",
                "Content-Type": "application/json",
            },
            callback=self.parse_reagent_page,
            meta={
                "data_url": data_url,
                "body": body,
                "reagent": reagent,
                "abbreviations": abbreviations,
                "source_url": source_url,
            },
        )

    def parse_reagent_page(self, response):
        payload = response.json()
        reagent = response.meta["reagent"]
        reagent_name = reagent.get("name") or "Unknown"
        results = payload.get("results") or []
        total_records = payload.get("total_records")

        logger.info("CEIRR %s: received %s of %s records.", reagent_name, len(results), total_records)

        columns = {}
        for column in reagent.get("column") or []:
            field_id = column.get("id")
            display_name = column.get("display_name")
            if field_id and display_name:
                columns.setdefault(field_id, []).append(display_name)

        for record in results:
            record["_ceirr_source_name"] = self.SOURCE_NAME
            record["_ceirr_reagent_category"] = reagent_name
            record["_ceirr_source_url"] = response.meta["source_url"]
            record["_ceirr_data_url"] = response.meta["data_url"]
            record["_ceirr_filter"] = response.meta["body"].get("filter") or []
            record["_ceirr_category_total_records"] = total_records
            record["_ceirr_column_labels"] = columns
            record["_ceirr_abbreviations"] = response.meta["abbreviations"]
            yield record

        pagination_key = payload.get("pagination_key")
        if pagination_key:
            body = dict(response.meta["body"])
            body["pagination_key"] = pagination_key
            yield self.post_reagent_page(
                data_url=response.meta["data_url"],
                body=body,
                reagent=reagent,
                abbreviations=response.meta["abbreviations"],
                source_url=response.meta["source_url"],
            )
