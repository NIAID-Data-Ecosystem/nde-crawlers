import logging
from datetime import datetime

import requests
import validators

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def parse():
    url = "https://accessclinicaldata.niaid.nih.gov/api/studies"
    response = requests.get(url)
    studies = response.json()
    count = 0

    for study in studies:
        if study["nct_number"] == "NCT04280705":
            study["isPartOf"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial",
                    "identifier": "ACTT",
                    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
                }
            ]
            study["isRelatedTo"] = [
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
        if study["nct_number"] == "NCT04401579":
            study["isPartOf"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial",
                    "identifier": "ACTT",
                    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
                }
            ]
            study["isRelatedTo"] = [
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
        if study["nct_number"] == "NCT04492475":
            study["isPartOf"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial",
                    "identifier": "ACTT",
                    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
                }
            ]
            study["isRelatedTo"] = [
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
        if study["nct_number"] == "NCT04640168":
            study["isPartOf"] = [
                {
                    "name": "Adaptive COVID-19 Treatment Trial",
                    "identifier": "ACTT",
                    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
                }
            ]
            study["isRelatedTo"] = [
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

        study["name"] = study.pop("title")
        study["identifier"] = study.pop("cmc_unique_id")
        study["description"] = study.pop("description")
        study["abstract"] = study.pop("brief_study_description")

        has_part_list = []
        for doc in study.get("study_documents", []):
            if doc["s3_location"] is None or doc["s3_location"].endswith(".zip"):
                continue
            has_part_list.append(
                {
                    "@type": "CreativeWork",
                    "name": doc["file_name"],
                    "url": "https://accessclinicaldata.niaid.nih.gov/api/files/" + doc["s3_location"],
                    "encodingFormat": doc["data_format"],
                }
            )
        if has_part_list:
            study["hasPart"] = has_part_list

        # URL no longer available
        # study["usageInfo"] = {
        #     "url": "https://accessclinicaldata.niaid.nih.gov/dashboard/Public/files/NIAIDDUA2021Accessclinicaldata@NIAID.pdf"
        # }

        # Convert publication dates to ISO format
        study["datePublished"] = study.pop("data_availability_date")
        if study["datePublished"] == "Coming Soon":
            study["datePublished"] = None
        if study["datePublished"] is not None:
            try:
                iso_date = datetime.strptime(study["datePublished"], "%B %Y")
            except ValueError:
                iso_date = datetime.strptime(study["datePublished"], "%B %d, %Y")
            study["datePublished"] = iso_date.strftime("%Y-%m-%d")

        study["dateModified"] = study.pop("most_recent_update")
        if study["dateModified"] is not None:
            iso_date = datetime.strptime(study["dateModified"], "%B %Y")
            study["dateModified"] = iso_date.strftime("%Y-%m-%d")

        study["additionalType"] = study.pop("data_available")
        study["funding"] = [{"funder": {"name": study.pop("creator")}}]
        study["nctid"] = study.pop("nct_number")
        study["healthCondition"] = {"name": study.pop("condition")}
        study["mainEntityOfPage"] = study.pop("clinical_trial_website")

        # Process the publications field for valid citation URLs
        citation_URL = study.pop("publications")
        if citation_URL is not None and validators.url(citation_URL):
            if "pubmed" in citation_URL:
                study["pmids"] = citation_URL.split("/")[-2]
            elif "doi" in citation_URL:
                doi_id = citation_URL.split("/")[-1]
                r = requests.get("https://pubmed.ncbi.nlm.nih.gov/?term=" + doi_id)
                if "pubmed" in r.url:
                    study["pmids"] = r.url.split("/")[-2]
                else:
                    study["citation"] = None
            else:
                study["citation"] = [{"url": citation_URL}]
        else:
            study["citation"] = None

        study["conditionsOfAccess"] = "Restricted" if study.pop("data_available_for_request") else "Closed"

        # De-duplicate identifiers if needed
        if study["identifier"] != study["nctid"]:
            study["identifier"] = [study["identifier"], study["nctid"]]
        else:
            study["identifier"] = [study["nctid"]]

        # Append a unique _id and set additional catalog fields
        study["_id"] = "accessclinicaldata_" + study["identifier"][0].lower()
        study["includedInDataCatalog"] = {"name": "AccessClinicalData@NIAID"}
        study["@type"] = "Dataset"
        dataset_url = "https://accessclinicaldata.niaid.nih.gov/study-viewer/clinical_trials/" + study["identifier"][0]
        study["url"] = dataset_url
        study["includedInDataCatalog"]["dataset"] = dataset_url

        # Remove any None values before yielding the record
        result = {k: v for k, v in study.items() if v is not None}
        missing_properties = {k: v for k, v in study.items() if v is None}

        yield result

        count += 1

        if len(missing_properties.keys()) > 0:
            logger.warning("Missing type transformation: {}".format(str(missing_properties.keys())))
        logger.info("Parsed %s records", count)

    logger.info("Finished Parsing. Total Records: %s", count)
    assert count < 10000, "Records have reached 10000, check API if records exceed 10000."
