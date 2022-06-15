import scrapy  
import requests 
import xmltodict 
import logging

logger = logging.getLogger('nde-logger')

class ClinEpiDBSpider(scrapy.Spider):

    name = 'clinepidb'

    sitemap_urls = ['https://clinepidb.org/ce/sitemap.xml']
    #start_urls = ['https://clinepidb.org/ce/app/record/dataset/DS_010e5612b8']

    # get the dataset ids
    def get_ids(self):
        xml_url = "https://clinepidb.org/ce/sitemap-SitemapDatasets.xml"
        id_response = requests.get(xml_url)
        dict_data = xmltodict.parse(id_response.content)
        url_list = [hit['loc'] for hit in dict_data["urlset"]["url"]]
        #id_dict = [hit['loc'].split("/")[-1] for hit in dict_data['urlset']['url']]
        logger.info("[INFO] Completed clinepidb dataset id list")
        return url_list
    
    def start_requests(self):
        urls = self.get_ids()
        for url in urls[:3]:
            yield scrapy.Request(url)

    def parse(self, response): 
        #body = response.xpath('//li[contains(@class,"main-stack") and contains(a/text(), "Interior")]/div//text()').extract()
        pass
    

