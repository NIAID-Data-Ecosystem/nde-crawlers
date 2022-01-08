from scrapy.spiders import SitemapSpider
from extruct import JsonLdExtractor


class OmicsdiSpider(SitemapSpider):

    name = 'omicsdi'
    custom_settings = {
        'ITEM_PIPELINES': {
            # TODO: add something here if transformations needed, which is likely
            'ndjson.NDJsonWriterPipeline': 999,
        }
    }

    # parsing huge XMLs is slow and there's like six dozens of them
    # expect a very slow start
    sitemap_urls = ['https://www.omicsdi.org/sitemap.xml']
    sitemap_rules = [('/dataset/', 'extract_from_jsonld')]

    def extract_from_jsonld(self, response, **kwargs):
        jslds = JsonLdExtractor().extract(response.body)

        for jsld in jslds:
            out = jsld['mainEntity']
            out['_id'] = response.url
            yield out
