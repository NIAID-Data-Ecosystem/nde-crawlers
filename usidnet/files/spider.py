# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh


import logging
import re

import requests
import scrapy
from w3lib.html import remove_tags, replace_entities

logger = logging.getLogger("nde-logger")

PAGE_PREFIX = "https://www.coriell.org/0/Sections/Search/"


def _cell_value(td):
    """Text of a value cell; radio/checkbox groups keep only the checked labels, <br> splits lines into a list."""
    if td.xpath(".//input"):
        lines = td.xpath(".//input[@checked]").xpath("string(following-sibling::node()[1])").getall()
    else:
        lines = [replace_entities(remove_tags(part)) for part in re.split(r"<br\s*/?>", td.get())]
    lines = [line for line in (" ".join(line.split()) for line in lines) if line]
    if len(lines) == 1:
        return lines[0]
    return lines or ""


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
        num_found = requests.get("https://www.coriell.org/Search/APIJson?q=*.*&page=1&fq=&pageSize=0&sort=Gene+asc", timeout=60).json().get("response").get("numFound")
        print(f"Total number of entries: {num_found}")
        samples = {}
        page = 1
        for page in range(1, (num_found - 1) // 10000 + 2):
            url = f"https://www.coriell.org/Search/APIJson?q=*.*&page={page}&fq=&pageSize=10000&sort=Gene+asc"
            print(f"Fetching page {page} with URL: {url}")
            response = requests.get(url, timeout=60)
            # an exception here fires spider_error, so NDJsonWriterPipeline discards the run instead of publishing it without this page
            response.raise_for_status()
            for doc in response.json().get("response").get("docs"):
                sample = samples.setdefault(doc.get("CatalogID"), doc)
                if sample is doc:
                    continue
                # the index has one doc per product (e.g. LCL and DNA) of a CatalogID
                products = sample["Product"] if isinstance(sample.get("Product"), list) else [sample.get("Product")]
                if doc.get("Product") and doc["Product"] not in products:
                    sample["Product"] = [p for p in products if p] + [doc["Product"]]
        yield from samples.items()

    async def start(self):
        for catalog_id, doc in self.get_ids():
            yield scrapy.Request(
                url=f"{PAGE_PREFIX}Sample_Detail.aspx?Ref={catalog_id}", callback=self.parse, meta={"doc": doc}
            )


    def parse(self, response):
        data = response.meta.get("doc", {})
        # panels show "Not Found" on Sample_Detail.aspx; their PageName (Panel_Detail.aspx) renders them
        page_name = data.get("PageName")
        not_found = response.xpath("normalize-space(//span[@id='lblRef'])").get() == "Not Found"
        if not_found and page_name and page_name not in response.url:
            yield scrapy.Request(
                url=f"{PAGE_PREFIX}{page_name}?Ref={data.get('CatalogID')}", callback=self.parse, meta={"doc": data}
            )
            return
        tables = response.xpath("//*[@role='presentation' and @class='table grid']")
        for table in tables:
            is_publications = bool(table.xpath("./ancestor::div[@id='Publications']"))
            prev_key = None
            section = None
            rows = table.xpath("./tr")
            for row in rows:
                if is_publications:
                    text = row.xpath("normalize-space(./td)").get()
                    text = text.replace("\xa0", " ").strip()
                    if text:
                        data.setdefault("Publications", []).append(text)
                    continue
                tds = row.xpath("./td")
                if len(tds) != 2:
                    # section header; a keyless row after it does not continue the previous section's key
                    section = " ".join(row.xpath("string(.)").get().split())
                    prev_key = None
                    continue
                # repeats the donor's disease fields for relatives
                if section == "Family History":
                    continue
                key = row.xpath("normalize-space(./td[1])").get()
                key = key.replace("\xa0", "").strip()
                value = _cell_value(tds[1])
                if key:
                    data[key] = value
                    prev_key = key
                elif prev_key and value:
                    values = data[prev_key] if isinstance(data[prev_key], list) else [data[prev_key]]
                    data[prev_key] = values + (value if isinstance(value, list) else [value])

        # Banner beside the catalog ID, e.g. "Fibroblast" or "DNA from LCL"
        if banner := " ".join(response.xpath("string(//span[@id='product-source'])").get().split()):
            data["Banner"] = banner

        # Pricing (Coriell shows tiered prices as <span class="price">$0.00</span>USD).
        # A price of $0.00 across all tiers means the sample is free.
        prices = [p.strip() for p in response.xpath("//span[@class='price']/text()").getall() if p.strip()]
        if prices:
            data["Prices"] = prices

        data["url"] = response.url
        yield data


