import logging
import re
from datetime import datetime

import dateutil.parser
import requests
import validators

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


# Strip a trailing all-uppercase abbreviation in parens (e.g. "Tuberculosis (TB)"
# -> "Tuberculosis"). Limited to 2-7 char uppercase tokens so we don't accidentally
# strip mixed-case parentheticals like "(SARS-CoV-2)" or full alt-name annotations.
_TRAILING_UPPER_ABBR = re.compile(r"\s*\(([A-Z][A-Z0-9]{1,6})\)\s*$")

# The pmid inside a PubMed article URL.
_PUBMED_PMID = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")
# A full DOI, registrant prefix included. The citations stage looks the DOI up
# verbatim, so the "10.xxxx/" prefix has to survive.
_DOI = re.compile(r"10\.\d{4,9}/[^\s\"'<>]+")

# Month-precision dates ("December 2025") carry no day. Default to the 1st rather
# than to today, which is what dateutil would otherwise fill in.
_DATE_DEFAULT = datetime(2000, 1, 1)

_ACTT_STUDY = {
    "name": "Adaptive COVID-19 Treatment Trial",
    "identifier": "ACTT",
    "url": "https://www.nih.gov/news-events/news-releases/fourth-iteration-covid-19-treatment-trial-underway",
}
# The four ACTT iterations are separate studies in the API that cross-reference
# each other. Keyed by nct number; the value is the dataset name a sibling uses
# when pointing at it.
_ACTT_ITERATIONS = {
    "NCT04280705": "Adaptive COVID-19 Treatment Trial (ACTT-1) - Dataset update released August 2021",
    "NCT04401579": "Adaptive COVID-19 Treatment Trial 2 (ACTT-2) - Dataset update released October 2021",
    "NCT04492475": "Adaptive COVID-19 Treatment Trial 3 (ACTT-3) - New dataset released October 2021",
    "NCT04640168": "Adaptive COVID-19 Treatment Trial 4 (ACTT-4) - New dataset released May 2022",
}
_ACTT_RELATIONSHIP = "Different iteration of the same study, the Adaptive COVID-19 Treatment Trial"


def _clean_health_condition_name(value):
    if not isinstance(value, str):
        return value
    m = _TRAILING_UPPER_ABBR.search(value)
    if m:
        return value[: m.start()].strip()
    return value


def _to_iso_date(value):
    """Normalize an API date to an ISO `YYYY-MM-DD` string, or None.

    The API mixes full ISO timestamps with month-precision names ("December
    2025") and the placeholder "Coming Soon". An unparseable value is logged and
    dropped, so one odd date cannot take down the whole crawl.
    """
    if not value or value == "Coming Soon":
        return None
    try:
        return dateutil.parser.parse(str(value), ignoretz=True, default=_DATE_DEFAULT).date().isoformat()
    except (dateutil.parser.ParserError, OverflowError, TypeError, ValueError):
        logger.warning("Could not parse date: %r", value)
        return None


def _actt_links(nct_number):
    """The curated isPartOf / isRelatedTo links between the four ACTT iterations."""
    if nct_number not in _ACTT_ITERATIONS:
        return {}
    return {
        "isPartOf": [{"@type": "CreativeWork", **_ACTT_STUDY}],
        "isRelatedTo": [
            {
                "@type": "Dataset",
                "name": name,
                "identifier": "accessclinicaldata_" + nct,
                "hasPart": {"@type": "CreativeWork", "identifier": "ACTT"},
                "includedInDataCatalog": {"@type": "DataCatalog", "name": "accessclinicaldata@NIAID"},
                "relationship": _ACTT_RELATIONSHIP,
            }
            for nct, name in _ACTT_ITERATIONS.items()
            if nct != nct_number
        ],
    }


def parse():
    url = "https://accessclinicaldata.niaid.nih.gov/api/studies"
    response = requests.get(url)
    studies = response.json()
    count = 0

    for study in studies:
        result = {}

        nct_number = study.get("nct_number")
        result.update(_actt_links(nct_number))

        result["name"] = study.get("title")
        unique_id = study.get("cmc_unique_id")
        result["description"] = study.get("description")
        result["abstract"] = study.get("brief_study_description")
        result["usageInfo"] = {
            "@type": "CreativeWork",
            "url": "https://accessclinicaldata.niaid.nih.gov/api/files/NIAIDDUAAccessclinicaldata@NIAID.pdf",
        }

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
            if data_format := doc.get("data_format"):
                creative_work["encodingFormat"] = data_format

            has_part_list.append(creative_work)

        if has_part_list:
            result["hasPart"] = has_part_list

        result["datePublished"] = _to_iso_date(study.get("data_availability_date"))
        result["dateModified"] = _to_iso_date(study.get("most_recent_update"))

        result["additionalType"] = study.get("data_available")
        if creator := study.get("creator"):
            result["funding"] = [{"@type": "MonetaryGrant", "funder": {"@type": "Organization", "name": creator}}]
        if nct_number and nct_number != "N/A":
            result["nctid"] = nct_number
        if condition := study.get("condition"):
            result["healthCondition"] = {"@type": "DefinedTerm", "name": _clean_health_condition_name(condition)}
        result["mainEntityOfPage"] = study.get("clinical_trial_website")

        citation_url = study.get("publications")
        if citation_url and validators.url(citation_url):
            if pmid_match := _PUBMED_PMID.search(citation_url):
                result["pmids"] = pmid_match.group(1)
            elif "doi" in citation_url and (doi_match := _DOI.search(citation_url)):
                result["citation"] = [{"@type": "ScholarlyArticle", "doi": doi_match.group()}]
            else:
                result["citation"] = [{"@type": "ScholarlyArticle", "url": citation_url}]

        result["conditionsOfAccess"] = "Restricted" if study.get("data_available_for_request") else "Closed"

        # nct_number arrives as "", None or "N/A" when the study has no trial
        # registration; none of those belong in identifier.
        identifiers = []
        for value in (unique_id, nct_number):
            if value and value.strip().upper() != "N/A" and value not in identifiers:
                identifiers.append(value)
        if not identifiers:
            logger.warning("Skipping study with no usable identifier: %r", study.get("title"))
            continue
        result["identifier"] = identifiers

        primary_id = identifiers[0]
        result["_id"] = "accessclinicaldata_" + primary_id.lower()
        dataset_url = "https://accessclinicaldata.niaid.nih.gov/study-viewer/clinical_trials/" + primary_id
        result["url"] = dataset_url
        result["includedInDataCatalog"] = {
            "@type": "DataCatalog",
            "name": "accessclinicaldata@NIAID",
            "archivedAt": dataset_url,
        }
        result["@type"] = "Dataset"

        # Drop keys the API left empty. "" is as useless to Elasticsearch as None,
        # and both show up here (clinical_trial_website, most notably).
        clean_result = {k: v for k, v in result.items() if v not in (None, "", [], {})}

        yield clean_result

        count += 1
        logger.info("Parsed %s records", count)

    logger.info("Finished Parsing. Total Records: %s", count)
    assert count < 10000, "Records have reached 10000, check API if records exceed 10000."
