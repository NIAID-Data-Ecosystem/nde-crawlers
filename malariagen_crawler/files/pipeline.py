# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface
# mapping https://docs.google.com/spreadsheets/d/19hQf4sQ6ZLwvt8ADYyvooYrqmQCsxh386kYnycX3xmg/edit#gid=0

import datetime
import logging

from biothings.utils.dataload import dict_sweep

logger = logging.getLogger("nde-logger")


__all__ = [
    "MalariaGenItemProcessorPipeline",
]


class MalariaGenItemProcessorPipeline:
    def process_item(self, item: dict, spider):
        _id = item["url"].rstrip("/").split("/")[-1].lower()
        item["_id"] = "malariagen_" + _id

        item["includedInDataCatalog"] = {
            "name": "MalariaGEN",
            "url": "https://www.malariagen.net/",
            "@type": "Dataset",
            "versionDate": datetime.datetime.today().strftime("%Y-%m-%d"),
            "dataset": item["url"],
        }
        item["@type"] = "Dataset"

        item["healthCondition"] = {
            "@type": "DefinedTerm",
            "identifier": "0005136",
            "inDefinedTermSet": "MONDO",
            "name": "Malaria",
            "url": "http://purl.obolibrary.org/obo/MONDO_0005136",
        }

        item["isAccessibleForFree"] = True

        cleaned_item = dict_sweep(item, vals=[None, "", [], {}])
        return cleaned_item
