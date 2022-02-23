# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface

__all__ = [
    'GeoItemProcessorPipeline',
]

def _set_single_element_or_all_in_list(obj, key, value):
    """mutates obj, list of dict or dict"""
    if isinstance(obj, list):
        for element in obj:
            element[key] = value
    else:
        obj[key] = value
    return obj

class GeoItemProcessorPipeline:
    def process_item(self, item: dict, spider):
        _id = item.pop('_id')
        url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + _id
        output = {
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": url,
            "identifier": _id,
            "distribution": {
                "@type": "dataDownload",
                "contentUrl": url
            },
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "NCBI GEO from Metadataplus",
                "url": "https://www.ncbi.nlm.nih.gov/geo/"
            }
        }

        if authors := item.pop('Contributor(s)',None):
            output['author'] = [{
                '@type': 'Person',
                'name': author
            } for author in authors.split(', ')]


        # Using this because the new GSE files use 'Organization name' instead of 'Organization'
        org = dict(filter(lambda value: 'Organization' in value[0], item.items()))
        org = org.keys()
        assert len(org) <= 1, "There is more than one organization key"
        if len(org) == 1:
            if publisher := item.pop(list(org)[0], None):
                output['publisher'] = {
                '@type': 'Organization',
                'name': publisher
                }

        # rename keys
        if name := item.pop('Title', None):
            output['name'] = name
        if organism := item.pop('Organism', None):
            output['organism'] = organism
        if measurement_technique := item.pop('Experiment type', None):
            output['measurementTechnique'] = measurement_technique
        if description := item.pop('Summary', None):
            output['description'] = description
        if date_published := item.pop('Submission date', None):
            output['datePublished'] = date_published
        if date_modified := item.pop('Last update date', None):
            output['dateModified'] = date_modified

        return output
