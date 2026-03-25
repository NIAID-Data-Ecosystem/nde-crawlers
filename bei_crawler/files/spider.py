# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh


import html
import logging
import os
import re

import scrapy
from scrapy import FormRequest, Selector
from scrapy_selenium import SeleniumRequest

logger = logging.getLogger("nde-logger")


# may take some time to start up as getting the ids takes a while
class BeiSpider(scrapy.Spider):
    name = "bei"
    custom_settings = {
        "ITEM_PIPELINES": {
            "pipeline.BeiItemProcessorPipeline": 100,
            "ndjson.NDJsonWriterPipeline": 999,
        }
    }

    start_urls = ["https://www.beiresources.org/Catalog.aspx?f_instockflag=In+Stock%23%7e%23Temporarily+Out+of+Stock&pagesize=100"]

    def build_form_payload(self, response):

        form = response.xpath("//form[@id='Form']")
        if not form:
            raise RuntimeError("Could not find the BEI search form needed for pagination.")

        payload = {}
        for inp in form.xpath(".//input[@name]"):
            input_type = (inp.xpath("./@type").get() or "").lower()
            if input_type in {"submit", "button", "image", "file"}:
                continue
            name = inp.xpath("./@name").get()
            value = inp.xpath("./@value").get() or ""
            payload[name] = value

        action = form.xpath("./@action").get() or response.url
        action_url = response.urljoin(action)
        return action_url, payload

    def parse(self, response):

        if "error=An+unexpected+error+has+occurred" in response.url:
            count = response.meta.get("count", 2)
            logger.warning(f"Encountered an error page on page {count} of BEI search results.")
            return

        base_url = "https://www.beiresources.org"

        for link in response.xpath('//*[@class="resultsproductname"]/div[1]/a/@href'):
            # logger.info(f"Found dataset link: {link.get()}")
            _id = os.path.splitext(os.path.basename(link.get()))[0]
            yield scrapy.Request(url=base_url + link.get(), callback=self.parse_dataset, meta={"id": _id})


        count = response.meta.get("count", 2)
        logger.info(f"Querying page {count} of BEI search results.")

        action_url, payload = self.build_form_payload(response)
        payload["__EVENTTARGET"] = "dnn$ctr14176$XPNBEI_SearchResults$rptResults"
        payload["__EVENTARGUMENT"] = f"Page${count}"
        logger.info(f"Submitting form to {action_url} with event argument {payload['__EVENTARGUMENT']}")

        yield scrapy.FormRequest(
            url=action_url,
            formdata=payload,
            callback=self.parse,
            meta={"count": count + 1},
        )

    def clean_texts(self, texts):
        out = []
        for t in texts:
            t = html.unescape(t).replace("\xa0", " ")
            t = re.sub(r"\s+", " ", t).strip()
            if not t:
                continue
            if "{" in t and "}" in t:   # drop CSS blocks
                continue
            out.append(t)
        return out

    def populate_output(self, keys, texts, output):
        key_set = set(keys)
        current_key = None

        for text in texts:
            if text in key_set:
                current_key = text
            elif current_key:
                output[current_key] = (output[current_key] + " " + text) if current_key in output else text

    def parse_dataset(self, response):
        table1_keys = response.xpath(
            '//*[@id="dnn_ctr8620_CTiBEI_PopulateTemplates_ctl10_tblItemDetails"]'
            '//*[contains(@style, "float: right")]//text()[normalize-space()]'
        ).getall()
        table1_keys = [k.strip() for k in table1_keys if k.strip()]
        table1 = response.xpath('//*[@id="dnn_ctr8620_CTiBEI_PopulateTemplates_ctl10_tblItemDetails"]//text()[normalize-space()]').getall()
        table2_keys = response.xpath(
            '//*[@id="dnn_ctr8620_CTiBEI_PopulateTemplates_ctl10_tblPSItemDetails"]'
            '//*[contains(@style, "float: right")]//text()[normalize-space()]'
        ).getall()
        table2_keys = [k.strip() for k in table2_keys if k.strip()]
        table2 = response.xpath('//*[@id="dnn_ctr8620_CTiBEI_PopulateTemplates_ctl10_tblPSItemDetails"]//text()[normalize-space()]').getall()
        table1 = self.clean_texts(table1)
        table2 = self.clean_texts(table2)

        output = {}
        output["_id"] = response.meta["id"]
        output["url"] = response.url
        self.populate_output(table1_keys, table1, output)
        self.populate_output(table2_keys, table2, output)
        yield output
