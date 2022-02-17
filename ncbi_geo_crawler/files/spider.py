# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file 
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh
'''
    NCBI Geo Datasource

    Scrape contents in HTML tables
    Not structured metadata

    https://www.ncbi.nlm.nih.gov/robots.txt
    Declared No Robots and crawl delay as of 02/15/2020.
    Defined long rules so the file content is not included here.

'''

import scrapy
from ftplib import FTP


class NCBIGeoSpider(scrapy.Spider):

    name = 'ncbi_geo'
    
    custom_settings = {
        'ITEM_PIPELINES': {
            'pipeline.GeoItemProcessorPipeline': 100,
            'ndjson.NDJsonWriterPipeline': 999,
        }
    }

    # gets the list of gse ids to request
    def get_ids():
        gse_li = []
        ftp_host = 'ftp.ncbi.nlm.nih.gov'
        ftp = FTP(ftp_host)
        ftp.login('','')
        ftp.cwd('/geo/series/')
        folder_li = ftp.nlst()
        for i in folder_li:
            ftp.cwd(i)
            gse_li = gse_li + ftp.nlst()
            ftp.cwd('..')
        print(len(gse_li))
        ftp.close()

    # THIS USED FOR TESTING/DEBUGGING
    # for small tests can use (start, end) = (1,20)
    # this should be the most recent assession (GSE) link: https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&display=1&zsort=acc
    def start_requests(self):
        start = 196607
        end = 196616
        prefix = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE"
        for acc_id in range(start, end + 1):
            yield scrapy.Request(url=prefix + str(acc_id))

    def parse(self, response):

        table = response.xpath(
            '/html/body/table/tr/td/table[6]/tr[3]/td[2]'
            '/table/tr/td/table/tr/td/table[2]/tr/td'
            '/table[1]/tr')
        data = {}

        for node in table:
            # extract series id
            if node.attrib.get('bgcolor') == '#cccccc':
                data['_id'] = node.xpath('.//strong').attrib.get('id')
            # remove place holder lines
            elif len(node.xpath('./td')) == 2:
                if node.xpath('string(./td[1])').get().strip():
                    # extract multi item entry
                    if node.xpath('./td[2]').attrib.get('onmouseout'):
                        key = node.xpath('./td[1]/text()').get().split()[0]
                        data[key] = node.xpath('./td[2]//a/text()').getall()
                    # extract single item entry
                    else:
                        key = node.xpath('./td[1]/text()').get()
                        data[key] = node.xpath('string(./td[2])').get().strip().replace('\xa0', ' ')

        yield data if data else None


