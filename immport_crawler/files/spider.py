# -*- coding: utf-8 -*-
'''
    ImmPort Spider
    Scrape shared studies on ImmPort

    Entry Point: https://www.immport.org/shared/search
    Example Pages:
        Without Pmid: https://www.immport.org/shared/study/SDY1
        With Pmid: https://www.immport.org/shared/study/SDY1025

    No robots.txt detected.

    More on https://github.com/biothings/biothings.crawler/issues/1
'''

import requests
import scrapy
from extruct.jsonld import JsonLdExtractor
from scrapy_selenium import SeleniumRequest
from selenium.webdriver.common.by import By
from selenium.webdriver.support.expected_conditions import presence_of_element_located


class ImmPortSpider(scrapy.Spider):
    """
    Crawl ImmPort with Selenium

    * Queries API for list of items
    * Uses Selenium Chrome Driver to render dynamic contents
    * Wait until page loads the JS so it has the LD-JSON
    * Extract LD-JSON

    """

    name = 'immport'
    custom_settings = {
        'DOWNLOADER_MIDDLEWARES': {
            'scrapy_selenium.SeleniumMiddleware': 800,
        },
        'ITEM_PIPELINES': {
            'pipeline.Immport2OutbreakDatasetPipeline': 300,
            'ndjson.NDJsonWriterPipeline': 999,
        }
    }
    immport_search_payload = {'pageSize': 1000}

    def start_requests(self):
        base_url = "https://www.immport.org/shared/data/query/search?term="
        try:
            j = requests.get(base_url, timeout=5,
                             params=self.immport_search_payload,).json()
            ids = [h['_id'] for h in j['hits']['hits']]
        except requests.exceptions.Timeout:
            # catch the rather harmless exception, and do nothing
            ids = []
            self.logger.warning("Timeout searching ImmPort IDs")

        prefix = "https://www.immport.org/shared/study/"
        for id_ in ids:
            yield SeleniumRequest(
                url=prefix + id_,
                callback=self.parse,
                wait_time=10,
                wait_until=presence_of_element_located(
                    (By.XPATH,
                     '//script[@type="application/ld+json"]')
                )
            )

    def parse(self, response, _id=None):
        """
        Scrapy Spider Request Callback Function

        * Inject an _id field for database pipeline
        * Use response URL as default _id

        """

        jslds = JsonLdExtractor().extract(response.body)

        for jsld in jslds:
            if _id:
                jsld['_id'] = _id
            else:
                jsld['_id'] = response.url
            yield jsld
