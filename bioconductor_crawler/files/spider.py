# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh

'''
    Bioconductor Datasource

    Scrape package badge content

    https://bioconductor.org/robots.txt

'''

import scrapy
import logging

from rpy2.robjects.packages import importr

logger = logging.getLogger('nde-logger')


# may take some time to start up as getting the ids takes a while
class BioconductorSpider(scrapy.Spider):

    name = 'bioconductor'

    custom_settings = {
        'ITEM_PIPELINES': {
            'pipeline.BioconductorItemProcessorPipeline': 100,
            'ndjson.NDJsonWriterPipeline': 999,
        }
    }

    def start_requests(self):
        # there may be a connection error sometimes. Retry up to 5 times if there is one.
        biocPkgTools = importr('BiocPkgTools')
        biocPkgList = biocPkgTools.biocPkgList()
        # get all Packages and save to a list
        biocPkgList = biocPkgList.rx2('Package')

        for package in biocPkgList:
            url = f"https://bioconductor.org/packages/release/bioc/html/{package}.html"
            yield scrapy.Request(url)

    def parse(self, response):

        print(response)

        # table = response.xpath(
        #     '/html/body/table/tr/td/table[6]/tr[3]/td[2]'
        #     '/table/tr/td/table/tr/td/table[2]/tr/td'
        #     '/table[1]/tr')
        # data = {}

        # for node in table:
        #     # extract series id
        #     if node.attrib.get('bgcolor') == '#cccccc':
        #         data['_id'] = node.xpath('.//strong').attrib.get('id')
        #     # remove place holder lines
        #     elif len(node.xpath('./td')) == 2:
        #         if node.xpath('string(./td[1])').get().strip():
        #             # extract multi item entry
        #             if node.xpath('./td[2]').attrib.get('onmouseout'):
        #                 key = node.xpath('./td[1]/text()').get().split()[0]
        #                 data[key] = node.xpath('./td[2]//a/text()').getall()
        #             # extract single item entry
        #             else:
        #                 key = node.xpath('./td[1]/text()').get()
        #                 data[key] = node.xpath(
        #                     'string(./td[2])').get().strip().replace('\xa0', ' ')

        # # self.logger.info("Check if cached: " + " ".join(response.flags))

        # yield data if data else None
