# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface

import datetime

__all__ = [
    "GeoItemProcessorPipeline",
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
        _id = item.pop("_id")
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + _id
        output = {
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": _id,
            "identifier": _id,
            "url": url,
            "distribution": {"@type": "dataDownload", "contentUrl": url},
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "NCBI GEO",
                "url": "https://www.ncbi.nlm.nih.gov/geo/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
        }

        if authors := item.pop("Contributor(s)", None):
            output["author"] = [{"@type": "Person", "name": author} for author in authors.split(", ")]

        # Using this because the new GSE files use 'Organization name' instead of 'Organization'
        org = dict(filter(lambda value: "Organization" in value[0], item.items()))
        org = org.keys()
        assert len(org) <= 1, "There is more than one organization key"
        if len(org) == 1:
            if publisher := item.pop(list(org)[0], None):
                output["publisher"] = {"@type": "Organization", "name": publisher}

        # rename keys
        if name := item.pop("Title", None):
            output["name"] = name
        if species := item.pop("Organism", None):
            output["species"] = {"name": species}
        if measurement_technique := item.pop("Experiment type", None):
            if isinstance(measurement_technique, list):
                output["measurementTechnique"] = [{"name": technique} for technique in measurement_technique]
            else:
                output["measurementTechnique"] = {"name": measurement_technique}
        description = ""
        if summary := item.pop("Summary", None):
            description += summary
        if overall_design := item.pop("Overall design", None):
            description += "\n" + overall_design
        if description:
            output["description"] = description
        if date_published := item.pop("Submission date", None):
            output["datePublished"] = datetime.datetime.strptime(date_published, "%b %d, %Y").date().isoformat()
        if date_modified := item.pop("Last update date", None):
            output["dateModified"] = datetime.datetime.strptime(date_modified, "%b %d, %Y").date().isoformat()
        # this will be used to call the api to get citation and funding in the uploader
        if pmids := item.pop("Citation(s)", None):
            output["pmids"] = pmids

        return output
