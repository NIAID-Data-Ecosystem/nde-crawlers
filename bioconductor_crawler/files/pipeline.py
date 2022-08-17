# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface

import datetime

__all__ = [
    'BioconductorItemProcessorPipeline',
]


def _set_single_element_or_all_in_list(obj, key, value):
    """mutates obj, list of dict or dict"""
    if isinstance(obj, list):
        for element in obj:
            element[key] = value
    else:
        obj[key] = value
    return obj


class BioconductorItemProcessorPipeline:
    def process_item(self, item: dict, spider):
        print(item)
        # return output
