import datetime
import logging

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def get_ids():
    logger.info("Getting datasets")
    url = "https://search.api.hubmapconsortium.org/v3/portal/search"
    payload = {
        "query": {
            "bool": {
                "must": [
                    {
                        "bool": {
                            "must_not": [
                                {"exists": {"field": "next_revision_uuid"}},
                                {"exists": {"field": "sub_status"}},
                            ]
                        }
                    },
                    {
                        "bool": {
                            "must_not": [
                                {"exists": {"field": "next_revision_uuid"}},
                                {"exists": {"field": "sub_status"}},
                            ]
                        }
                    },
                ]
            }
        },
        "post_filter": {"term": {"entity_type.keyword": "Dataset"}},
        "_source": False,
        "size": 10000,
    }
    r = requests.post(url, json=payload).json()
    ids = [x["_id"] for x in r["hits"]["hits"]]
    logger.info("Total datasets: %s", len(ids))
    count = 0
    datasets = []
    for id in ids:
        count += 1
        if count % 50 == 0:
            logger.info("Retrieved %s out of %s datasets", count, len(ids))
        dataset = requests.get(f"https://portal.hubmapconsortium.org/browse/dataset/{id}.json")
        datasets.append(dataset.json())
    logger.info("Retrieved %s out of %s datasets", count, len(ids))
    return datasets


def parse():
    datasets = get_ids()
    count = 0
    for metadata in datasets:
        count += 1
        if count % 50 == 0:
            logger.info("Parsed %s out of %s datasets", count, len(datasets))
        output = {
            "includedInDataCatalog": {
                "@type": "Dataset",
                "name": "HuBMAP",
                "url": "https://hubmapconsortium.org/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Dataset",
        }

        keywords = []
        if anatomy_1 := metadata.get("anatomy_0"):
            for term in anatomy_1:
                keywords.append(term)
        if anatomy_2 := metadata.get("anatomy_1"):
            for term in anatomy_2:
                keywords.append(term)
        if anatomy_3 := metadata.get("anatomy_2"):
            for term in anatomy_3:
                keywords.append(term)
        if display_subtype := metadata.get("display_subtype"):
            keywords.append(display_subtype)
        if len(keywords):
            output["keywords"] = keywords

        author_list = []
        if contacts := metadata.get("contacts"):
            for contact in contacts:
                author = {}
                if affiliation := contact.get("affiliation"):
                    author["affiliation"] = {"name": affiliation}
                if first_name := contact.get("first_name"):
                    author["givenName"] = first_name
                if last_name := contact.get("last_name"):
                    author["familyName"] = last_name
                if middle_name := contact.get("middle_name_or_initial"):
                    author["givenName"] += f" {middle_name}"
                if name := contact.get("name"):
                    author["name"] = name
                if orcid_id := contact.get("orcid_id"):
                    author["identifier"] = orcid_id
                if bool(author):
                    author_list.append(author)

        if contributors := metadata.get("contributors"):
            for contributor in contributors:
                author = {}
                if affiliation := contributor.get("affiliation"):
                    author["affiliation"] = {"name": affiliation}
                if first_name := contributor.get("first_name"):
                    author["givenName"] = first_name
                if last_name := contributor.get("last_name"):
                    author["familyName"] = last_name
                if middle_name := contributor.get("middle_name_or_initial"):
                    author["givenName"] += f" {middle_name}"
                if name := contributor.get("name"):
                    author["name"] = name
                if orcid_id := contributor.get("orcid_id"):
                    author["identifier"] = orcid_id
                if bool(author):
                    author_list.append(author)

        if len(author_list):
            output["author"] = author_list

        if created_timestamp := metadata.get("created_timestamp"):
            output["dateCreated"] = datetime.datetime.utcfromtimestamp(created_timestamp / 1000).strftime("%Y-%m-%d")

        if data_access_level := metadata.get("data_access_level"):
            if data_access_level == "public":
                output["isAccessibleForFree"] = True
            else:
                output["isAccessibleForFree"] = False

        measurement_technique = {}
        if data_types := metadata.get("data_types"):
            measurement_technique["name"] = data_types[0]
        if dataset_info := metadata.get("dataset_info"):
            measurement_technique["description"] = dataset_info
        if "name" in measurement_technique:
            output["measurementTechnique"] = measurement_technique

        if doi_url := metadata.get("doi_url"):
            output["doi"] = doi_url.removeprefix("https://doi.org/")

        if uuid := metadata.get("uuid"):
            url = f"https://portal.hubmapconsortium.org/browse/dataset/{uuid}"
            output["url"] = url
            output["includedInDataCatalog"]["dataset"] = url
            output["_id"] = "HUBMAP_" + uuid

        if files := metadata.get("files"):
            distribution_list = []
            for file in files:
                distribution_dict = {}
                if file_description := file.get("description"):
                    distribution_dict["name"] = file_description
                # if edam_term := file.get('edam_term'):
                #     distribution_dict['encodingFormat'] = edam_term
                if mapped_description := file.get("mapped_description"):
                    if "name" in distribution_dict:
                        distribution_dict["description"] = mapped_description
                    else:
                        distribution_dict["name"] = mapped_description
                if rel_path := file.get("rel_path"):
                    distribution_dict["contentUrl"] = f"https://assets.hubmapconsortium.org/{uuid}/{rel_path}"
                if size := file.get("size"):
                    distribution_dict["contentSize"] = size
                if file_type := file.get("type"):
                    if file_type != "unknown":
                        distribution_dict["encodingFormat"] = file_type
                if bool(distribution_dict):
                    distribution_list.append(distribution_dict)
            if len(distribution_list):
                output["distribution"] = distribution_list

        if hubmap_id := metadata.get("hubmap_id"):
            output["name"] = hubmap_id

        if version := metadata.get("version"):
            output["version"] = version

        if last_modified_timestamp := metadata.get("last_modified_timestamp"):
            output["dateModified"] = datetime.datetime.utcfromtimestamp(last_modified_timestamp / 1000).strftime(
                "%Y-%m-%d"
            )

        if published_timestamp := metadata.get("published_timestamp"):
            output["datePublished"] = datetime.datetime.utcfromtimestamp(published_timestamp / 1000).strftime(
                "%Y-%m-%d"
            )

        if dataset_metadata := metadata.get("metadata"):
            if origin := dataset_metadata.get("dag_provenance_list"):
                output["isBasedOn"] = {"name": origin[0]["origin"]}
            if nested_metadata := dataset_metadata.get("metadata"):
                if protocols_io_doi := nested_metadata.get("protocols_io_doi"):
                    if "isBasedOn" in output:
                        output["isBasedOn"]["doi"] = protocols_io_doi
                    else:
                        output["isBasedOn"] = {"doi": protocols_io_doi}
                if section_prep_protocols_io_doi := nested_metadata.get("section_prep_protocols_io_doi"):
                    if "isBasedOn" in output:
                        output["isBasedOn"]["doi"] = section_prep_protocols_io_doi
                    else:
                        output["isBasedOn"] = {"doi": section_prep_protocols_io_doi}
                # if analyte_class := nested_metadata.get('analyte_class'):
                #     output['variableMeasured'] = analyte_class
                # if assay_category := nested_metadata.get('assay_category'):
                #     output['measurementTechnique']['name'] = assay_category
                # if antibodies_path := nested_metadata.get('antibodies_path'):
                #     print(antibodies_path)
                #     print(output['url'])
                #     if 'distribution' in output:
                #         print(output['distribution'])
        if title := metadata.get("title"):
            output["description"] = title

        yield output

    logger.info("Finished parsing %s datasets", count)
