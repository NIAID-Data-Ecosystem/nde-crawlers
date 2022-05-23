# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface

import datetime
import logging

logger = logging.getLogger('nde-logger')


__all__ = [
    'OmicsdiItemProcessorPipeline',
]


class OmicsdiItemProcessorPipeline:
    def process_item(self, item: dict, spider):
        url = item.pop('_id')
        url_split = url.rsplit('/', 1)
        output = {
            "@context": "http://schema.org/",
            "@type": item.pop('@type'),
            "url": url,
            "_id": "OMICSDI_" + url_split[-1],
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "Omics Discovery Index (OmicsDI)",
                "url": "https://www.omicsdi.org/",
                'versionDate': datetime.date.today().isoformat()
            }
        }

        # collect pmids
        if pmids := item.pop('pmids', None):
            for i in range(len(pmids)):
                pmids[i] = pmids[i].split(': ', 1)[-1]
            output['pmids'] = ','.join(pmids)


        # put publisher into journalName since there is no journalName
        # based on https://discovery.biothings.io/view/outbreak/Publication
        if citation := item.pop('citation', None):
            sd_publisher = {}
            if publisher := citation.pop('publisher', None):
                if name := publisher.pop('name', None):
                    sd_publisher['name'] = name
                if ty := publisher.pop('@type', None):
                    sd_publisher['@type'] = ty
            if url := citation.pop('url', None):
                sd_publisher['url'] = url
            output['sdPublisher'] = sd_publisher


        # remove items from dictionary
        if author := item.pop('creator', None):
            output['author'] = author
        if description := item.pop('description', None):
            output['description'] = description
        if distribution := item.pop('distribution', None):
            output['distribution'] = distribution
        if keywords := item.pop('keywords', None):
            output['keywords'] = keywords
        if name := item.pop('name', None):
            output['name'] = name
        if same_as := item.pop('sameAs', None):
            output['sameAs'] = same_as
        if variable_measured := item.pop('variableMeasured', None):
            output['variableMeasured'] = variable_measured

        if item:
            logger.warning("Haven't parsed all keys in omiscdi_crawler: " + "\t".join(item.keys()))
        
        return output
