# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file 
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh

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
    # limit = 2000  # Limit entries
    # count = 0  # Entries counter

    # def sitemap_filter(self, entries):
    #     for entry in entries:
    #         if self.count >= self.limit:
    #             continue
    #         self.count += 1
    #         yield entry

    def extract_from_jsonld(self, response, **kwargs):
        jslds = JsonLdExtractor().extract(response.body)

        for jsld in jslds:
            out = jsld
            date_published = response.xpath("//div[@class='o-metadata__group2-item']/text()").extract()[0]
            assert "Publication date" in date_published, f"Did not get publication date: {date_published} ID: {out['@id']}"

            out['datePublished'] = date_published
            yield out



