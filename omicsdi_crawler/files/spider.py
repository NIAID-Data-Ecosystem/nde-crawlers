# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh

from extruct import JsonLdExtractor
from scrapy.spiders import SitemapSpider


class OmicsdiSpider(SitemapSpider):
    name = "omicsdi"
    custom_settings = {
        "ITEM_PIPELINES": {
            "pipeline.OmicsdiItemProcessorPipeline": 100,
            "ndjson.NDJsonWriterPipeline": 999,
        }
    }

    # parsing huge XMLs is slow and there's like six dozens of them
    # expect a very slow start
    sitemap_urls = ["https://www.omicsdi.org/sitemap.xml"]
    sitemap_rules = [("/dataset/", "extract_from_jsonld")]

    def extract_from_jsonld(self, response, **kwargs):
        jslds = JsonLdExtractor().extract(response.body.decode("utf-8"))

        for jsld in jslds:
            out = jsld["mainEntity"]
            out["_id"] = response.url
            out["pmids"] = response.xpath("//*[@id='pubId']/text()").extract()
            yield out
