# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file 
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh

import scrapy
import time
import json
from scrapy.spiders import SitemapSpider
from extruct import JsonLdExtractor


class DryadSpider(SitemapSpider):

    name = 'dryad'
    custom_settings = {
        'ITEM_PIPELINES': {
            'pipeline.DryadItemProcessorPipeline': 100,
            'ndjson.NDJsonWriterPipeline': 999,
        }
    }

    # parsing huge XMLs is slow
    # expect a very slow start
    sitemap_urls = ['https://datadryad.org/sitemap.xml']
    sitemap_rules = [('/stash/', 'extract_from_jsonld')]

    # for testing purposes limit requests to urls
    # limit = 70  # Limit entries
    # count = 0  # Entries counter
    #
    # def sitemap_filter(self, entries):
    #     for entry in entries:
    #         if self.count >= self.limit:
    #             continue
    #         self.count += 1
    #         yield entry

    def extract_from_jsonld(self, response, **kwargs):
        jslds = JsonLdExtractor().extract(response.body)

        callback_url = "https://datadryad.org/api/v2/datasets/doi%3A" + response.url.split(':')[3].replace('/', '%2F')

        for jsld in jslds:
            out = jsld
            time.sleep(2)
            request = scrapy.Request(callback_url, callback=self.parse_api)
            request.meta['item'] = out
            yield request

    def parse_api(self, response):
        item = response.meta['item']
        item['callback_url'] = response.url
        item.update(json.loads(response.body.decode('utf-8')))
        yield item


