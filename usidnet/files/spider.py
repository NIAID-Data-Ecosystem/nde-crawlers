# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh


import logging
import requests
import scrapy

logger = logging.getLogger("nde-logger")


# may take some time to start up as getting the ids takes a while
class USIDNETSpider(scrapy.Spider):

    name = "usidnet_spider"
    complete_data = {}

    custom_settings = {
        "ITEM_PIPELINES": {
            "pipeline.USIDNETItemProcessorPipeline": 100,
            "ndjson.NDJsonWriterPipeline": 999,
        }
    }

    # gets the list of sample ids to request
    def get_ids(self):
        num_found = requests.get("https://www.coriell.org/Search/APIJson?q=*.*&page=1&fq=&pageSize=0&sort=Gene+asc").json().get("response").get("numFound")
        print(f"Total number of entries: {num_found}")
        page = 1
        for page in range(1, (num_found - 1) // 10000 + 2):
            url = f"https://www.coriell.org/Search/APIJson?q=*.*&page={page}&fq=&pageSize=10000&sort=Gene+asc"
            print(f"Fetching page {page} with URL: {url}")
            response = requests.get(url)
            if response.status_code == 200:
                docs = response.json().get("response").get("docs")
                for doc in docs:
                    yield doc.get("CatalogID"), doc
            else:
                print(f"Failed to fetch page {page}. Status code: {response.status_code}")

    def start_requests(self):
        # Scrapy < 2.13 entry point
        prefix = "https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref="
        for catalog_id, doc in self.get_ids():
            yield scrapy.Request(url=prefix + str(catalog_id), callback=self.parse, meta={"doc": doc})

    async def start(self):
        # Scrapy >= 2.13 calls start() instead of start_requests()
        for request in self.start_requests():
            yield request


    def parse(self, response):
        data = response.meta.get("doc", {})
        tables = response.xpath("//*[@role='presentation' and @class='table grid']")
        for table in tables:
            is_publications = bool(table.xpath("./ancestor::div[@id='Publications']"))
            prev_key = None
            rows = table.xpath("./tr")
            for row in rows:
                if is_publications:
                    text = row.xpath("normalize-space(./td)").get()
                    text = text.replace("\xa0", " ").strip()
                    if text:
                        data.setdefault("Publications", []).append(text)
                    continue
                tds = row.xpath("./td")
                if len(tds) == 2:
                    key = row.xpath("normalize-space(./td[1])").get()
                    key = key.replace("\xa0", "").strip()
                    value = row.xpath("normalize-space(./td[2])").get()
                    if key:
                        data[key] = value.replace("\xa0", " ")
                        prev_key = key
                    else:
                        if not isinstance(data.get(prev_key), list):
                            data[prev_key] = [data[prev_key]]
                        data[prev_key].append(value.replace("\xa0", " "))

        # Pricing (Coriell shows tiered prices as <span class="price">$0.00</span>USD).
        # A price of $0.00 across all tiers means the sample is free.
        prices = [p.strip() for p in response.xpath("//span[@class='price']/text()").getall() if p.strip()]
        if prices:
            data["Prices"] = prices

        yield data if data else None


