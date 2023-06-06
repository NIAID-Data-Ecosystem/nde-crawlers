import datetime
import json
import logging
import re
import urllib.parse

import requests
import validators

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


def get_dataset_ids(study_name):
    dataset_ids = []
    page = 1
    if "'" not in study_name:
        # trying to encode single quotes break the url, even after converting to %27
        study_name = urllib.parse.quote(study_name)
    while True:
        facets = [{"name": "study name", "filters": [f"{study_name}"]}]
        url = (
            f"https://dash.nichd.nih.gov/search/dataset?q=&facets={facets}&page={page}&sortBy=title&asc=true&size=1000"
        )
        try:
            data = requests.get(url).json()
        except:
            logger.error(f"Could not get dataset ids for {study_name}")
            break
        if len(data["hits"]["hits"]) == 0:
            break
        for dataset in data["hits"]["hits"]:
            dataset_ids.append(dataset["_id"])
        page += 1
    return dataset_ids


def get_study_ids():
    study_ids = []
    page = 1
    while True:
        url = f"https://dash.nichd.nih.gov/search/study?size=1000&page={page}"
        data = requests.get(url).json()
        if len(data["hits"]["hits"]) == 0:
            break
        for study in data["hits"]["hits"]:
            study_ids.append(study["_id"])
        page += 1
    return study_ids


def get_dataset_info(dataset_id):
    url = f"https://dash.nichd.nih.gov/api/study/datasetOverview/{dataset_id}"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    return json.loads(response.text)


def parse_study_info(study_info):
    study_dict = {}
    citation_dict = {}
    if citation := study_info.get("citation"):
        citation_dict["name"] = citation.strip()

    for obj in study_info["studyInfo"]["title"]:
        if obj["propertyName"] == "Study Name":
            study_dict["name"] = obj["storedValue"]
        if obj["propertyName"] == "Study Abbreviation":
            study_dict["alternateName"] = obj["storedValue"]

    author_list = []
    for obj in study_info["studyInfo"]["main"]:
        if obj["propertyName"] == "DOI":
            citation_dict["doi"] = obj["storedValue"].strip()
        if obj["propertyName"] == "NICHD Division/Branch/Center":
            study_dict["funding"] = {"funder": {"name": obj["storedValue"]}}
        if obj["propertyName"] == "Study Description":
            study_dict["description"] = obj["storedValue"]
        if obj["propertyName"] == "Clinical Research Network Name":
            if obj["storedValue"]:
                author_list.append({"name": obj["storedValue"]})
            if len(obj["storedArray"]) > 0:
                for name in obj["storedArray"]:
                    author_list.append({"name": name})
        if obj["propertyName"] == "Principal Investigator(s)":
            if obj["storedValue"]:
                author_list.append({"name": obj["storedValue"]})
            elif len(obj["storedArray"]) > 0:
                for name in obj["storedArray"]:
                    author_list.append({"name": name})

    for obj in study_info["studyInfo"]["details"]:
        if obj["propertyName"] == "Keywords":
            if obj["storedValue"]:
                study_dict["keywords"] = obj["storedValue"]
            elif len(obj["storedArray"]) > 0:
                study_dict["keywords"] = obj["storedArray"]
        if obj["propertyName"] == "Topic":
            if "keywords" in study_dict:
                if obj["storedValue"]:
                    study_dict["keywords"].append(obj["storedValue"])
                elif len(obj["storedArray"]) > 0:
                    study_dict["keywords"].extend(obj["storedArray"])
            else:
                if obj["storedValue"]:
                    study_dict["keywords"] = obj["storedValue"]
                elif len(obj["storedArray"]) > 0:
                    study_dict["keywords"] = obj["storedArray"]
        if obj["propertyName"] == "Requires IRB approval to obtain data":
            if obj["storedValue"] == "Yes":
                study_dict["conditionsOfAccess"] = "Restricted"
        if obj["propertyName"] == "Additional Approval Entity":
            if obj["storedValue"] == "Yes":
                study_dict["conditionsOfAccess"] = "Restricted"
        if obj["propertyName"] == "Study Type":
            study_dict["measurementTechnique"] = {"name": obj["storedValue"]}
        if obj["propertyName"] == "ClinicalTrials.gov URL":
            study_dict["mainEntityOfPage"] = obj["storedValue"]
            study_dict["nctid"] = obj["storedValue"].split("/")[-1]
        if obj["propertyName"] == "Publication URLs":
            if obj["storedValue"] is not None and validators.url(obj["storedValue"]):
                citation_dict["url"] = obj["storedValue"]
            pmids = []
            pmcids = []
            for url in obj["storedArray"]:
                if url is not None and validators.url(url):
                    if "pubmed" in url:
                        pmids.append(url.split("/")[-1])
                    elif "PMC" in url:
                        pmc_regex = re.compile(r"PMC\d+")
                        pmcids.append(pmc_regex.search(url).group(0))
                    else:
                        citation_dict["url"] = url
            if len(pmcids) > 0:
                url = f'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={",".join(pmcids)}&format=json'
                try:
                    data = requests.get(url).json()
                    pmcids = [obj["pmid"] for obj in data["records"]]
                except:
                    logger.error(f'PMC ID conversion failed for {study_dict["name"]}')
                pmids.extend(pmcids)
            if len(pmids) > 0:
                study_dict["pmids"] = ",".join(pmids)

    if bool(citation_dict):
        study_dict["citation"] = citation_dict

    if "studyTimeline" in study_info["studyInfo"]:
        temporal_coverage = {}
        for obj in study_info["studyInfo"]["studyTimeline"]:
            if obj["propertyName"] == "StudyCollectionStartDate":
                temporal_coverage["startDate"] = obj["storedValue"]
            if obj["propertyName"] == "StudyCollectionEndDate":
                temporal_coverage["endDate"] = obj["storedValue"]

        if bool(temporal_coverage):
            study_dict["temporalCoverage"] = temporal_coverage

    if approval_date := study_info.get("approvalDate"):
        approval_date = approval_date.split("T")[0]
        study_dict["datePublished"] = approval_date

    if update_date := study_info.get("studyUpdateDate"):
        update_date = update_date.split("T")[0]
        study_dict["dateModified"] = update_date

    if study_population := study_info.get("studyPopulation"):
        if population_description := study_population.get("Population Description"):
            if "description" in study_dict:
                study_dict["description"] += f" {population_description}"
            else:
                study_dict["description"] = population_description

    if descriptive_documents := study_info.get("descriptiveDocuments"):
        for obj in descriptive_documents:
            if obj["documentType"] == "Codebook/Variable Dictionary":
                study_dict["hasPart"] = {
                    "additionalType": {"name": obj["documentType"]},
                    "name": obj["filename"],
                    "url": f'https://dash.nichd.nih.gov/download-api/descriptive/file?id={obj["id"]}',
                }
            else:
                study_dict["isBasedOn"] = {
                    "additionalType": {"name": obj["documentType"]},
                    "name": obj["filename"],
                    "url": f'https://dash.nichd.nih.gov/download-api/descriptive/file?id={obj["id"]}',
                }
    return study_dict


def get_study_info(study_id):
    url = f"https://dash.nichd.nih.gov/api/study/studyOverview/{study_id}"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    study_info = json.loads(response.text)
    study_dict = parse_study_info(study_info)

    logger.info(f"Retrieving dataset ids for study id: {study_id}...")
    dataset_ids = get_dataset_ids(study_dict["name"])
    return (study_dict, dataset_ids)


def parse():
    logger.info("Retrieving study ids...")
    study_ids = get_study_ids()
    logger.info(f"Found {len(study_ids)} studies.")

    logger.info("Retrieving individual study info...")
    study_count = 0
    dataset_count = 0
    for study_id in study_ids:
        study_count += 1
        logger.info(f"Parsing study {study_count} of {len(study_ids)}")

        info = get_study_info(study_id)
        parsed_study = info[0]

        if info == None:
            logger.info(f"Failed to get study info for study id:{study_id}, skipping.")
            continue

        dataset_ids = info[1]
        if len(dataset_ids) == 0:
            logger.info(f"No datasets found for study id: {study_id}")
            continue
        else:
            logger.info(f"Found {len(dataset_ids)} datasets for study id: {study_id}")

        logger.info("Retrieving individual dataset info...")
        related_datasets = []
        for dataset_id in dataset_ids:
            dataset_info = get_dataset_info(dataset_id)
            if dataset_info == None:
                logger.info(f"Failed to get dataset info for {dataset_id}, skipping.")
                continue

            output = parsed_study.copy()

            output["isPartOf"] = [
                {
                    "name": output["name"],
                    "identifier": "NICHD_DASH_Study_" + study_id,
                    "url": f"https://dash.nichd.nih.gov/study/{study_id}",
                }
            ]

            output["includedInDataCatalog"] = {
                "@type": "Dataset",
                "name": "NICHD DASH",
                "url": "https://dash.nichd.nih.gov/",
                "versionDate": datetime.date.today().isoformat(),
            }
            output["@type"] = "Dataset"

            output["_id"] = "NICHD_DASH_Dataset_" + dataset_id
            output["url"] = f"https://dash.nichd.nih.gov/dataset/{dataset_id}"
            distribution_dict = {}
            distribution_dict["contentUrl"] = output["url"]

            if dataset_title := dataset_info.get("datasetTitle"):
                output["name"] = f'{dataset_title} in {output["name"]}'

            if dataset_description := dataset_info.get("datasetDescription"):
                output["description"] = f'{dataset_description}\nStudy Description\n{output["description"]}'

            if dataset_format := dataset_info.get("datasetFormat"):
                distribution_dict["encodingFormat"] = dataset_format
            output["distribution"] = distribution_dict

            related_datasets.append(output)

        for dataset in related_datasets:
            dataset["isRelatedTo"] = [
                {
                    "name": x["name"],
                    "identifier": x["_id"],
                    "hasPart": {"identifier": x["isPartOf"][0]["identifier"]},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "NICHD DASH"},
                    "relationship": "Datasets in the same study",
                }
                for x in related_datasets
                if x != dataset
            ]

            dataset_count += 1

            if dataset_count % 100 == 0:
                logger.info(f"Parsed {dataset_count} datasets")

            yield dataset
    logger.info(f"Finished parsing. Total datasets: {dataset_count}")
