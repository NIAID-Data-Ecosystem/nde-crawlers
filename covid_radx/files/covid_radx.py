import datetime
import logging
import time

import dateutil
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def convert_group_string_to_list(input_string):
    # Remove the curly braces
    stripped_string = input_string.strip("{}")
    # Split by comma and strip double quotes from each element
    result_list = [item.strip('"') for item in stripped_string.split(",")]
    return result_list


def extract_ids(urls_string):
    urls = urls_string.split(",")
    pmids = []
    pmcs = []

    for url in urls:
        if "pubmed" in url:
            pmid = url.rstrip("/").split("/")[-1]
            pmids.append(pmid)
        elif "pmc" in url:
            pmc = url.rstrip("/").split("/")[-1]
            pmcs.append(pmc)

    pmids_string = ",".join(pmids)
    pmcs_string = ",".join(pmcs)

    return pmids_string, pmcs_string


def parse(hit):
    output = {
        "_id": f"radx_{hit['_id']}",
        "@type": "Dataset",
        "includedInDataCatalog": {
            "@type": "DataCatalog",
            "name": "RADx Data Hub",
            "url": "https://radxdatahub.nih.gov/",
            "versionDate": datetime.date.today().isoformat(),
        },
        "conditionsOfAccess": "Restricted",
        "license": "https://sharing.nih.gov/accessing-data/accessing-genomic-data/how-to-request-and-access-datasets-from-dbgap",
    }
    # Add the rest of the fields from _source
    hit = hit["_source"]
    if doi := hit.get("study_doi"):
        output["doi"] = doi

    if phs := hit.get("phs"):
        output["identifier"] = phs

    if version := hit.get("@version"):
        output["@version"] = version

    if date_created := hit.get("created_at"):
        output["dateCreated"] = dateutil.parser.parse(date_created, ignoretz=True).date().isoformat()

    if ibo_url := hit.get("ct_url"):
        output["isBasedOn"] = {"url": ibo_url}

    if sd_name := hit.get("dcc"):
        output["sdPublisher"] = {"name": sd_name}

    if description := hit.get("description"):
        output["description"] = description

    if health_conditions := hit.get("disease_specific_related_conditions"):
        health_conditions = health_conditions.split(", ")
        output["healthCondition"] = [{"name": health_condition for health_condition in health_conditions}]
    else:
        output["healthCondition"] = {"name": "COVID-19"}

    funding = {}
    identifier = hit.get("grant_number").strip("{}").split(",")
    funder = hit.get("institutes_supporting_study_array")

    if funder and identifier:
        if len(funder) == len(identifier):
            for i in range(len(funder)):
                funders = []
                funders.append(
                    {
                        "name": funder[i],
                        "identifier": identifier[i],
                    }
                )
        else:
            for f in funder:
                funders.append({"name": f})
    elif funder:
        for f in funder:
            funders = {"name": f}

    if funders:
        funding["funder"] = funders
        output["funding"] = funding

    authors = []
    if names := hit.get("multi_center_sites"):
        names = names.split("; ")
        for name in names:
            authors.append({"name": name, "@type": "Organization"})

    if name := hit.get("pi_name"):
        name = name.split(", ")
        last_name = name[0]
        first_name = name[1]
        authors.append({"familyName": last_name, "givenName": first_name, "@type": "Person"})

    if authors:
        output["author"] = authors

    if date_published := hit.get("release_date"):
        output["datePublished"] = date_published

    if urls := hit.get("publication_url"):
        pmids, pmcs = extract_ids(urls)
        if pmids:
            output["pmids"] = pmids
        if pmcs:
            output["pmcs"] = pmcs

    if date_published := hit.get("release_date"):
        output["datePublished"] = dateutil.parser.parse(date_published).date().isoformat()

    total_measurement_techniques = []
    mts = [hit.get("source_array"), hit.get("types_array"), hit.get("data_general_types_array")]
    for measurement_technique in mts:
        if measurement_technique:
            for mt in measurement_technique:
                total_measurement_techniques.append({"name": mt})
    if total_measurement_techniques:
        output["measurementTechnique"] = total_measurement_techniques

    base_url = "https://radxdatahub.nih.gov/study/"
    if study_id := hit.get("study_id"):
        url = base_url + str(study_id)
        output["url"] = url
        output["includedInDataCatalog"]["archivedAt"] = url

    total_keywords = []
    if keywords := hit.get("population_focus_array"):
        total_keywords.extend(keywords)
    if keywords := hit.get("topics_array"):
        total_keywords.extend(keywords)
    if total_keywords:
        output["keywords"] = total_keywords

    if main_entity_of_page := hit.get("study_website_url"):
        output["mainEntityOfPage"] = main_entity_of_page.strip("{}")

    temporal_coverage = {"@type": "temporalInterval", "temporalType": "study date"}
    if start_date := hit.get("studystartdate"):
        temporal_coverage["startDate"] = dateutil.parser.parse(start_date).date().isoformat()
    if end_date := hit.get("studyenddate"):
        temporal_coverage["endDate"] = dateutil.parser.parse(end_date).date().isoformat()

    if temporal_coverage.get("startDate") or temporal_coverage.get("endDate"):
        output["temporalCoverage"] = {"temporalInterval": temporal_coverage}

    if variable_measured := hit.get("subject_array"):
        vms = []
        for vm in variable_measured:
            vms.append({"name": vm})
        output["variableMeasured"] = vms

    if name := hit.get("title"):
        output["name"] = name

    if date_modified := hit.get("updated_at"):
        output["dateModified"] = dateutil.parser.parse(date_modified, ignoretz=True).date().isoformat()

    usage_info_keys = ["general_research_group", "health_biomed_group"]
    usage_info_list = []
    for key in usage_info_keys:
        if name := hit.get(key):
            name = convert_group_string_to_list(name)
            usage_info = {"name": name}
            usage_info_list.append(usage_info)
    if usage_info_list:
        output["usageInfo"] = usage_info_list

    return output


def make_requests():
    page = 1
    size = 50
    request = requests.get(
        "https://radxdatahub.nih.gov/_next/data/cx3d2E3xBhLBo0cZuZj8o/studyExplorer.json?sort=desc&prop=relevance&page=1&size=50"
    )
    request = request.json()
    total_hits = request["pageProps"]["searchResults"]["hits"]["total"]["value"]

    while (page - 1) * size < total_hits:
        logging.info(f"Processing page {page} of {total_hits // size + 1}")
        request = requests.get(
            f"https://radxdatahub.nih.gov/_next/data/cx3d2E3xBhLBo0cZuZj8o/studyExplorer.json?sort=desc&prop=relevance&page={page}&size={size}"
        )
        page += 1
        request = request.json()
        for hit in request["pageProps"]["searchResults"]["hits"]["hits"]:
            yield parse(hit)

    time.sleep(1)
