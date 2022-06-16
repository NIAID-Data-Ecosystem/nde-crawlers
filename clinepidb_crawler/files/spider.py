import scrapy  
import requests 
import xmltodict 
import logging
import json
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
        test_url = urls[0]
        self.test_id = test_url.split("/")[-1]

        #for url in urls[:3]:
        yield scrapy.Request(test_url, callback=self.parse)

    def parse(self, response): 
        filename = "/Users/nacosta/Documents/local_testing/clinepidb/local_test_clinepidb_testing_"+str(self.test_id)
        #data = json.loads(response.body)
        #body = response.xpath('//li[contains(@class,"main-stack") and contains(a/text(), "Interior")]/div//text()').extract()
        #body = response.xpath('//div[@class="wdk-DataTableCellContent"]').extract()
        with open(filename, 'wb') as f:
            f.write(response.body)

    

