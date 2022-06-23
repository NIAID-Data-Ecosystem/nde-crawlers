import scrapy  
import requests 
import xmltodict 
import logging
from scrapy_selenium import SeleniumRequest
from selenium.webdriver.common.by import By
from selenium.webdriver.support.expected_conditions import presence_of_element_located

logger = logging.getLogger('nde-logger')

class ClinEpiDBSpider(scrapy.Spider):

    name = 'clinepidb'

    sitemap_urls = ['https://clinepidb.org/ce/sitemap.xml']
    #start_urls = ['https://clinepidb.org/ce/app/record/dataset/DS_010e5612b8']

    def start_requests(self):
        xml_url = "https://clinepidb.org/ce/sitemap-SitemapDatasets.xml"
        sitemap_request = requests.get(xml_url)
        dict_data = xmltodict.parse(sitemap_request.content)
        urls = [hit['loc'] for hit in dict_data["urlset"]["url"]]
        test_url = urls[0]
        #id_dict = [hit['loc'].split("/")[-1] for hit in dict_data['urlset']['url']]
        logger.info("[INFO] Completed clinepidb dataset id list")
        test_id = test_url.split("/")[-1]

        for _url in urls[:3]:
            yield SeleniumRequest(
                url=_url
                wait_time=10,
                callback=self.parse,
                wait_until=presence_of_element_located(
                    (By.XPATH,
                    '')
                )
            )
        yield scrapy.Request(test_url, callback=self.parse)

    def parse(self, response): 
        filename = "/Users/nacosta/Documents/local_testing/clinepidb/local_test_clinepidb_testing_"+str(self.test_id)
        #data = json.loads(response.body)
        #body = response.xpath('//li[contains(@class,"main-stack") and contains(a/text(), "Interior")]/div//text()').extract()
        #body = response.xpath('//div[@class="wdk-DataTableCellContent"]').extract()
        #with open(filename, 'wb') as f:
            #f.write(response.body)

