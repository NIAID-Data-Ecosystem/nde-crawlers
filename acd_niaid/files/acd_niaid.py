import json
import logging
from datetime import datetime

import requests
import validators

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def parse():
    # schema
    # https://docs.google.com/spreadsheets/d/1XlpvoeSWSiqfw1pPAHHFTjr9Wuowj3TA-eCsyGEj0ng/edit#gid=0

    url = "https://accessclinicaldata.niaid.nih.gov/guppy/graphql"
    query = """query ($filter: JSON) {
        clinical_trials (filter: $filter, first: 10000, accessibility: accessible) {
            title,cmc_unique_id,brief_summary,data_availability_date,most_recent_update,data_available,creator,nct_number,condition,clinical_trial_website,publications,data_available_for_request,description
        }
    }"""

    r = requests.post(url, json={"query": query})

    # parse json string and convert to dictionary
    json_obj = json.loads(r.text)

    # grab list of trials from dictionary
    trials = json_obj["data"]["clinical_trials"]

    count = 0

    for trial in trials:
        if trial["nct_number"] == "NCT04280705":
            trial["isPartOf"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial",
                    "identifier": "ACTT",
                    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
                }
            ]

            trial["isRelatedTo"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial 2 (ACTT-2) - Dataset update released October 2021",
                    "identifier": "accessclinicaldata_NCT04401579",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
                {
                    "name": "Adaptive COVID-19 Treatment Trial 3 (ACTT-3) - New dataset released October 2021",
                    "identifier": "accessclinicaldata_NCT04492475",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
                {
                    "name": "Adaptive COVID-19 Treatment Trial 4 (ACTT-4) - New dataset released May 2022",
                    "identifier": "accessclinicaldata_NCT04640168",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
            ]
        if trial["nct_number"] == "NCT04401579":
            trial["isPartOf"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial",
                    "identifier": "ACTT",
                    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
                }
            ]

            trial["isRelatedTo"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial (ACTT-1) - Dataset update released August 2021",
                    "identifier": "accessclinicaldata_NCT04280705",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
                {
                    "name": "Adaptive COVID-19 Treatment Trial 3 (ACTT-3) - New dataset released October 2021",
                    "identifier": "accessclinicaldata_NCT04492475",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
                {
                    "name": "Adaptive COVID-19 Treatment Trial 4 (ACTT-4) - New dataset released May 2022",
                    "identifier": "accessclinicaldata_NCT04640168",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
            ]
        if trial["nct_number"] == "NCT04492475":
            trial["isPartOf"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial",
                    "identifier": "ACTT",
                    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
                }
            ]

            trial["isRelatedTo"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial (ACTT-1) - Dataset update released August 2021",
                    "identifier": "accessclinicaldata_NCT04280705",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
                {
                    "name": "Adaptive COVID-19 Treatment Trial 2 (ACTT-2) - Dataset update released October 2021",
                    "identifier": "accessclinicaldata_NCT04401579",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
                {
                    "name": "Adaptive COVID-19 Treatment Trial 4 (ACTT-4) - New dataset released May 2022",
                    "identifier": "accessclinicaldata_NCT04640168",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
            ]
        if trial["nct_number"] == "NCT04640168":
            trial["isPartOf"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial",
                    "identifier": "ACTT",
                    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
                }
            ]

            trial["isRelatedTo"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial (ACTT-1) - Dataset update released August 2021",
                    "identifier": "accessclinicaldata_NCT04280705",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
                {
                    "name": "Adaptive COVID-19 Treatment Trial 2 (ACTT-2) - Dataset update released October 2021",
                    "identifier": "accessclinicaldata_NCT04401579",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
                {
                    "name": "Adaptive COVID-19 Treatment Trial 3 (ACTT-3) - New dataset released October 2021",
                    "identifier": "accessclinicaldata_NCT04492475",
                    "hasPart": {"identifier": "ACTT"},
                    "@type": "Dataset",
                    "includedInDataCatalog": {"name": "AccessClinicalData@NIAID"},
                    "relationship": "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial",
                },
            ]

    for trial in trials:
        trial["name"] = trial.pop("title")
        trial["identifier"] = trial.pop("cmc_unique_id")
        trial["description"] = trial.pop("description")
        trial["abstract"] = trial.pop("brief_summary")

        trial_query = """
        query ($filter: JSON) {
            oafile (filter: $filter, first: 10000, accessibility: accessible) {
                file_name, file_size, data_format, data_type, cmc_unique_id, doc_url
            }
        }
        """

        trial_variables = {"filter": {"in": {"cmc_unique_id": [trial["identifier"]]}}}

        response = requests.post(
            url, json={"query": trial_query, "variables": trial_variables}, headers={"Content-Type": "application/json"}
        )

        test_json = json.loads(response.text)
        has_part_list = []
        for file in test_json["data"]["oafile"]:
            has_part_list.append(
                {
                    "@type": "CreativeWork",
                    "name": file["file_name"],
                    "url": file["doc_url"],
                    "encodingFormat": file["data_format"],
                }
            )
        if len(has_part_list):
            trial["hasPart"] = has_part_list
        trial["usageInfo"] = {
            "url": "https://accessclinicaldata.niaid.nih.gov/dashboard/Public/files/NIAIDDUA2021Accessclinicaldata@NIAID.pdf"
        }

        # convert date to iso format
        trial["datePublished"] = trial.pop("data_availability_date")
        if trial["datePublished"] == "Coming Soon":
            trial["datePublished"] = None
        if trial["datePublished"] is not None:
            try:
                iso_date = datetime.strptime(trial["datePublished"], "%B %Y")
            except ValueError:
                iso_date = datetime.strptime(trial["datePublished"], "%B %d, %Y")
            trial["datePublished"] = iso_date.strftime("%Y-%m-%d")

        # convert date to iso format
        trial["dateModified"] = trial.pop("most_recent_update")
        if trial["dateModified"] is not None:
            iso_date = datetime.strptime(trial["dateModified"], "%B %Y")
            trial["dateModified"] = iso_date.strftime("%Y-%m-%d")

        trial["additionalType"] = trial.pop("data_available")
        trial["funding"] = [{"funder": {"name": trial.pop("creator")}}]
        trial["nctid"] = trial.pop("nct_number")
        trial["healthCondition"] = {"name": trial.pop("condition")}
        trial["mainEntityOfPage"] = trial.pop("clinical_trial_website")

        # check if url is valid then take pubmed id
        citation_URL = trial.pop("publications")
        if citation_URL is not None and validators.url(citation_URL):
            if "pubmed" in citation_URL:
                trial["pmids"] = citation_URL.split("/")[-2]
            # To convert doi id to pubmed id, use requests library to search doi id on pubmed search engine, catch redirect and take the pubmed id from url
            elif "doi" in citation_URL:
                doi_id = citation_URL.split("/")[-1]
                r = requests.get("https://pubmed.ncbi.nlm.nih.gov/?term=" + doi_id)
                if "pubmed" in r.url:
                    trial["pmids"] = r.url.split("/")[-2]
                else:
                    trial["citation"] = None
            else:
                trial["citation"] = [{"url": citation_URL}]
        else:
            trial["citation"] = None

        if trial["data_available_for_request"] == "true":
            trial["data_available_for_request"] = "Restricted"
        if trial["data_available_for_request"] == "false":
            trial["data_available_for_request"] = "Closed"

        trial["conditionsOfAccess"] = trial.pop("data_available_for_request")

        # de-duplication of identifier
        if trial["identifier"] != trial["nctid"]:
            trial["identifier"] = [trial["identifier"], trial["nctid"]]
        else:
            trial["identifier"] = [trial["nctid"]]

        # unique _id appending identifier
        trial["_id"] = "accessclinicaldata_" + trial["identifier"][0]
        trial["includedInDataCatalog"] = {"name": "AccessClinicalData@NIAID"}
        trial["@type"] = "Dataset"
        dataset_url = "https://accessclinicaldata.niaid.nih.gov/study-viewer/clinical_trials/" + trial["identifier"][0]
        trial["url"] = dataset_url
        trial["includedInDataCatalog"]["dataset"] = dataset_url

        # getting rid of None values
        result = {k: v for k, v in trial.items() if v is not None}

        # list properties that weren't in trial
        missing_properties = {k: v for k, v in trial.items() if v is None}

        yield result

        count += 1

        if len(missing_properties.keys()) > 0:
            logger.warning("Missing type transformation: {}".format(str(missing_properties.keys())))
        logger.info("Parsed %s records", count)

    logger.info("Finished Parsing. Total Records: %s", count)
    assert count < 10000, "Records has reached 10000, check API if records exceed 10000."
