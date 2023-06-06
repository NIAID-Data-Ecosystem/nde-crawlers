# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh

import json
import logging

import scrapy

logger = logging.getLogger("nde-logger")

""" For testing in scrapy shell
import scrapy
import json
url = "http://flowrepository.org/ajax/list_public_ds"
payload = {'pg':'2'}
req = scrapy.Request( url, method='POST', body=json.dumps(payload), headers={'Content-Type':'application/json'}) 
fetch(req)
response.xpath('//td[@class="repid"]/a/@href').extract()
response.xpath('//a[@href="#"]/text()').getall()



scrapy shell "http://flowrepository.org/id/FR-FCM-ZY4P"  
keys = response.xpath('//table[@class="information-table"]//b/text()').extract()
for key in keys:
    " ".join(response.xpath(f"//td[*='{key}']/following-sibling::td[1]/text()").extract_first().split())
    response.xpath(f"//td[*='{key}']/following-sibling::td[1]/descendant-or-self::text()").extract()

"""

# may take some time to start up as getting the ids takes a while
class FlowRepositorySpider(scrapy.Spider):
    name = "flowrepository"

    custom_settings = {
        "ITEM_PIPELINES": {
            "pipeline.FlowRepositoryItemProcessorPipeline": 100,
            "ndjson.NDJsonWriterPipeline": 999,
        }
    }

    # starting page number

    page = 1
    url = "http://flowrepository.org/ajax/list_public_ds"
    allowed_domains = ["flowrepository.org"]

    def start_requests(self):
        """Simulate the AJAX post request to get the table of ids.
        Yields:
            To a function that will parse the urls from the table.
        """

        payload = {"pg": str(self.page)}
        yield scrapy.Request(
            self.url,
            method="POST",
            body=json.dumps(payload),
            headers={"Content-Type": "application/json"},
            callback=self.get_urls,
        )

    def get_urls(self, response):
        id_urls = response.xpath('//td[@class="repid"]/a/@href').extract()
        for id_url in id_urls:
            yield scrapy.Request(id_url, callback=self.parse_urls)
        # Check for the next page button
        next_page = response.xpath('//a[@href="#"]/text()').extract()
        if "Next >>" in next_page:
            self.page += 1
            payload = {"pg": str(self.page)}
            # create the new AJAX post request for the next page
            yield scrapy.Request(
                self.url,
                method="POST",
                body=json.dumps(payload),
                headers={"Content-Type": "application/json"},
                callback=self.get_urls,
            )

    def parse_urls(self, response):
        data = {"url": response.url}
        keys = response.xpath('//table[@class="information-table"]//b/text()').extract()
        for key in keys:
            values = response.xpath(f"//td[*='{key}']/following-sibling::td[1]/descendant-or-self::text()").extract()
            values = [" ".join(value.split()) for value in values if " ".join(value.split())]
            data[key] = values
        yield data
