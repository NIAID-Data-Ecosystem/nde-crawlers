import datetime
import json
import logging
import re

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def retrieve_study_metadata():
    logger.info("Retrieving study metadata from VDJ")
    study_urls = [
        "https://vdjserver.org/airr/v1/repertoire",
        "https://vdj-staging.tacc.utexas.edu/airr/v1/repertoire",
        "https://ipa3.ireceptor.org/airr/v1/repertoire",
        "https://ipa4.ireceptor.org/airr/v1/repertoire",
        "https://ipa1.ireceptor.org/airr/v1/repertoire",
        "https://ipa2.ireceptor.org/airr/v1/repertoire",
        "https://scireptor.dkfz.de/airr/v1/repertoire",
    ]

    response_data = {}
    for url in study_urls:
        logger.info("Retrieving studies from %s", url)
        try:
            r = requests.post(url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logger.error("Error retrieving studies from %s", url)
            logger.error(e)
            continue
        except requests.exceptions.SSLError as e:
            logger.error("SSL issue from %s", url)
            logger.error(e)
            continue
        data = r.json()["Repertoire"]
        for study_dict in data:
            id = study_dict["study"]["study_id"]
            if id not in response_data:
                if "ireceptor" in url:
                    study_dict["publisher"] = "iReceptor"
                elif "vdjserver" in url:
                    study_dict["publisher"] = "VDJServer"
                elif "scireptor" in url:
                    study_dict["publisher"] = "sciReptor"
                response_data[id] = study_dict

    download_url = "https://vdj-staging.tacc.utexas.edu/api/v2/adc/cache/study"
    try:
        r = requests.get(download_url)
        data = r.json()["result"]
        for result in data:
            if result["study_id"] in response_data:
                response_data[result["study_id"]]["study"]["download_info"] = result
    except requests.exceptions.SSLError as e:
        logger.error("SSL issue from %s when getting download_info", url)
        logger.error(e)

    logger.info("Total number of studies: %s", len(response_data))

    with open("vdj_studies.json", "w") as f:
        json.dump(response_data, f, indent=4)

    return response_data.values()


def parse():
    studies = retrieve_study_metadata()
    count = 0
    logger.info("Parsing study metadata")
    for study in studies:
        count += 1
        logger.info(f"Parsing study {count} of {len(studies)}")
        if count % 10 == 0:
            logger.info("Parsed %s studies", count)

        output = {
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "VDJServer",
                "url": "https://vdj-staging.tacc.utexas.edu/community/",
                "versionDate": datetime.date.today().isoformat(),
                "dataset": "https://vdj-staging.tacc.utexas.edu/community/",
            },
            "@type": "Dataset",
            "url": "https://vdj-staging.tacc.utexas.edu/community/",
        }
        if study_info := study.get("study"):
            if study_id := study_info.get("study_id"):
                output["_id"] = f"vdj_{study_id}".replace(":", "_").replace(" ", "_").replace("/", "_")

            if study_title := study_info.get("study_title"):
                output["name"] = study_title

            if study_type := study_info.get("study_type"):
                measurement_technique = {}
                if id := study_type.get("id"):
                    measurement_technique["identifier"] = id
                if label := study_type.get("label"):
                    measurement_technique["name"] = label
                if len(measurement_technique) > 0:
                    measurement_technique["inDefinedTermSet"] = "NCI Thesaurus"
                    output["measurementTechnique"] = measurement_technique

            if study_description := study_info.get("study_description"):
                output["description"] = study_description

            if inclusion_exclusion_criteria := study_info.get("inclusion_exclusion_criteria"):
                if output.get("description"):
                    output["description"] += f"\n {inclusion_exclusion_criteria}"
                else:
                    output["description"] = inclusion_exclusion_criteria

            authors = []
            author_dict = {}
            if lab_name := study_info.get("lab_name"):
                author_dict["name"] = lab_name

            if lab_address := study_info.get("lab_address"):
                author_dict["affiliation"] = {"name": lab_address}

            if len(author_dict) > 0:
                authors.append(author_dict)

            if submitted_by := study_info.get("submitted_by"):
                info = submitted_by.split(", ")
                if len(info) == 2:
                    authors.append({"name": info[0], "affiliation": {"name": info[1]}})
                else:
                    authors.append({"name": submitted_by})

            if collected_by := study_info.get("collected_by"):
                info = collected_by.split(", ")
                for author in info:
                    # if author is an email address, skip
                    if "@" not in author:
                        authors.append({"name": author})

            if len(authors) > 0:
                output["author"] = authors

            if grants := study_info.get("grants"):
                if output.get("description"):
                    output["description"] += f"\n {grants}"
                else:
                    output["description"] = grants

            # if grants := study_info.get('grants'):
            #     output['funding'] = {'description': grants}

            if pub_ids := study_info.get("pub_ids"):
                if "pmid" in pub_ids:
                    output["pmids"] = pub_ids.split(":")[1].strip()
                elif "doi" in pub_ids.lower():
                    doi = re.search(r"10\.\d{4,9}\/[-._;()/:A-Z0-9]+", pub_ids, re.IGNORECASE)
                    if doi:
                        output["citation"] = {"doi": doi.group(0)}
            if adc_publish_date := study_info.get("adc_publish_date"):
                adc_publish_date = adc_publish_date.split(".")[0].strip()
                output["datePublished"] = datetime.datetime.strptime(adc_publish_date, "%Y-%m-%dT%H:%M:%S").strftime(
                    "%Y-%m-%d"
                )

            if adc_update_date := study_info.get("adc_update_date"):
                adc_update_date = adc_update_date.split(".")[0].strip()
                output["dateModified"] = datetime.datetime.strptime(adc_update_date, "%Y-%m-%dT%H:%M:%S").strftime(
                    "%Y-%m-%d"
                )

        if publisher := study.get("publisher"):
            output["sdPublisher"] = {"name": publisher}

        if download_info := study.get("download_info"):
            distribution_dict = {}

            if archive_file := download_info.get("archive_file"):
                distribution_dict["contentUrl"] = archive_file

            if download_url := download_info.get("download_url"):
                distribution_dict["contentUrl"] = download_url

            if file_size := download_info.get("file_size"):
                distribution_dict["contentSize"] = file_size

            if len(distribution_dict) > 0:
                output["distribution"] = distribution_dict

        yield output

    logger.info(f"Parsed {count} studies")
