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
        result = {}

        nct_number = study.get("nct_number")

        # Conditionally add relationship fields based on nct_number.
        if nct_number == "NCT04280705":
            result["isPartOf"] = [{
                "name": "Adaptive COVID-19 Treatment Trial",
                "identifier": "ACTT",
                "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
            }]
            result["isRelatedTo"] = [
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
        elif nct_number == "NCT04401579":
            result["isPartOf"] = [{
                "name": "Adaptive COVID-19 Treatment Trial",
                "identifier": "ACTT",
                "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
            }]
            result["isRelatedTo"] = [
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
        elif nct_number == "NCT04492475":
            result["isPartOf"] = [{
                "name": "Adaptive COVID-19 Treatment Trial",
                "identifier": "ACTT",
                "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
            }]
            result["isRelatedTo"] = [
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
        elif nct_number == "NCT04640168":
            result["isPartOf"] = [{
                "name": "Adaptive COVID-19 Treatment Trial",
                "identifier": "ACTT",
                "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
            }]
            result["isRelatedTo"] = [
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

        result["name"] = study.get("title")
        unique_id = study.get("cmc_unique_id")
        result["description"] = study.get("description")
        result["abstract"] = study.get("brief_study_description")
        result["usageInfo"] = {"url":"https://accessclinicaldata.niaid.nih.gov/api/files/NIAIDDUAAccessclinicaldata@NIAID.pdf"}

        has_part_list = []
        for doc in study.get("study_documents", []):
            s3_location = doc.get("s3_location")
            if s3_location is None or s3_location.endswith(".zip"):
                continue

            creative_work = {
                "@type": "CreativeWork",
                "name": doc.get("file_name"),
                "url": "https://accessclinicaldata.niaid.nih.gov/api/files/" + s3_location
            }
            data_format = doc.get("data_format")
            if data_format is not None:
                creative_work["encodingFormat"] = data_format

            has_part_list.append(creative_work)

        if has_part_list:
            result["hasPart"] = has_part_list

        date_published = study.get("data_availability_date")
        if date_published == "Coming Soon":
            result["datePublished"] = None
        elif date_published:
            if "T" in date_published:
                result["datePublished"] = datetime.fromisoformat(date_published).strftime("%Y-%m-%d")
            else:
                try:
                    iso_date = datetime.strptime(date_published, "%B %Y")
                except ValueError:
                    iso_date = datetime.strptime(date_published, "%B %d, %Y")
                result["datePublished"] = iso_date.strftime("%Y-%m-%d")
        else:
            result["datePublished"] = None

        date_modified = study.get("most_recent_update")
        if date_modified:
            if "T" in daPte_modified:
                result["dateModified"] = datetime.fromisoformat(date_modified).strftime("%Y-%m-%d")
            else:
                iso_date = datetime.strptime(date_modified, "%B %Y")
                result["dateModified"] = iso_date.strftime("%Y-%m-%d")
        else:
            result["dateModified"] = None

        result["additionalType"] = study.get("data_available")
        result["funding"] = [{"funder": {"name": study.get("creator")}}]
        if nct_number and nct_number != "N/A":
            result["nctid"] = nct_number
        result["healthCondition"] = {"name": study.get("condition")}
        result["mainEntityOfPage"] = study.get("clinical_trial_website")

        citation_url = study.get("publications")
        if citation_url and validators.url(citation_url):
            if "pubmed" in citation_url:
                result["pmids"] = citation_url.split("/")[-2]
            elif "doi" in citation_url:
                doi_id = citation_url.split("/")[-1]
                r = requests.get("https://pubmed.ncbi.nlm.nih.gov/?term=" + doi_id)
                if "pubmed" in r.url:
                    result["pmids"] = r.url.split("/")[-2]
                else:
                    result["citation"] = None
            else:
                result["citation"] = [{"url": citation_url}]
        else:
            result["citation"] = None

        result["conditionsOfAccess"] = "Restricted" if study.get("data_available_for_request") else "Closed"

        identifiers = [x for x in (unique_id, nct_number) if x is not None]
        seen = set()
        identifiers = [x for x in identifiers if x not in seen and not seen.add(x)]
        result["identifier"] = identifiers

        primary_id = result["identifier"][0]
        result["_id"] = "accessclinicaldata_" + primary_id.lower()
        dataset_url = "https://accessclinicaldata.niaid.nih.gov/study-viewer/clinical_trials/" + primary_id
        result["url"] = dataset_url
        result["includedInDataCatalog"] = {"name": "AccessClinicalData@NIAID", "dataset": dataset_url}
        result["@type"] = "Dataset"

        # Remove any keys with None values.
        clean_result = {k: v for k, v in result.items() if v is not None}

        yield clean_result

        count += 1
        logger.info("Parsed %s records", count)

    logger.info("Finished Parsing. Total Records: %s", count)
    assert count < 10000, "Records have reached 10000, check API if records exceed 10000."
