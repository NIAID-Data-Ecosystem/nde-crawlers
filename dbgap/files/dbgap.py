import datetime
import json
import logging
import os

import regex as re
import requests
from download_files import download
from lxml import html

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

UNIT_MAP = {
    "year": "year",
    "years": "year",
    "yr": "year",
    "yrs": "year",
    "y": "year",
    "month": "month",
    "months": "month",
    "m": "month",
    "d": "day",
    "day": "day",
    "days": "day",
}

UNIT_RE = r"(?:years?|yrs?|y|months?|m|days?|d)"

RACE_MAP = {
    "White": ["white", "caucasian"],
    "Black or African American": ["black", "african american", "african-american"],
    "Asian": ["asian", "chinese", "japanese", "korean", "filipino", "vietnamese", "indian"],
    "Hispanic or Latino": ["hispanic", "latino", "latina"],
    "Native American": ["native american", "american indian"],
    "Pacific Islander": ["pacific islander"],
    "Middle Eastern": ["middle eastern", "arab"],
}

_NEGATION_RE = re.compile(r"\b(?:not|no|non|without|excluding|except|never)\b", re.I)

def parse_race_string(text: str) -> list:
    """
    Return a list of canonical race labels found in `text`.
    - Matches synonyms with word boundaries (so 'asian' won't match 'vegasian').
    - Skips matches if the immediately previous word is a negation word (uses _NEGATION_RE).
    - Handles combined forms like "Black/African American" (both will be detected).
    """
    if not text:
        return []

    s = text.lower()
    found = set()

    for canon, terms in RACE_MAP.items():
        candidates = list(terms) + [canon.lower()]
        for term in candidates:
            term = term.strip()
            if not term:
                continue

            for m in re.finditer(rf"\b{re.escape(term)}\b", s):
                # inspect only the single word immediately before the match
                prev = re.search(r"([a-z]+)\W*$", s[:m.start()])
                if prev and _NEGATION_RE.search(prev.group(1)):
                    continue
                found.add(canon)

    return sorted(found)

def _normalize_unit(unit_text: str | None) -> str | None:
    if not unit_text:
        return None
    return UNIT_MAP.get(unit_text.lower())


def _extract_unit(text: str) -> str | None:
    m = re.search(rf"\b({UNIT_RE})\b", text.lower())
    return _normalize_unit(m.group(1)) if m else None


def parse_age_string(age_text: str) -> dict:
    """
    Parse age text into {"minVal": int?, "maxVal": int?, "unitText": "year|month|day"?}
    """
    s = age_text.strip().lower()
    result: dict = {}

    # 1) Range: "0-20 years", "2 - 36 months", "45 to 80 year"
    m = re.search(
        rf"\b(?:age(?:d|s)?(?:\s+of)?)?\s*(\d+)\s*(?:[-–—]|to)\s*(\d+)\s*({UNIT_RE})?\b",
        s,
    )
    if m:
        result["minVal"] = int(m.group(1))
        result["maxVal"] = int(m.group(2))
        unit = _normalize_unit(m.group(3)) or _extract_unit(s)
        if unit:
            result["unitText"] = unit
        return result

    # 2) Max only: "Less than 18 years of age", "under 18 years",
    # "18 years of age and younger", "18 or younger", "at most 18 years"
    m = re.search(
        rf"\b(?:"
        rf"(?:less than|under|below|younger than|at most|up to)\s*(\d+)\s*({UNIT_RE})?"
        rf"|"
        rf"(\d+)\s*({UNIT_RE})?(?:\s+of\s+age)?\s*(?:or|and)?\s*(?:younger|below|less)"
        rf")\b",
        s,
    )
    if m:
        num = m.group(1) or m.group(3)
        unit = m.group(2) or m.group(4)
        result["maxVal"] = int(num)
        unit = _normalize_unit(unit) or _extract_unit(s)
        if unit:
            result["unitText"] = unit
        return result

    # 3) Min only: "at least 18 years", "over 18", "18 years of age or older", "18+ years"
    m = re.search(
        rf"\b(?:"
        rf"(?:at least|over|above|greater than|more than|older than)\s*(\d+)\s*({UNIT_RE})?"
        rf"|"
        rf"(\d+)\s*(?:\+\s*({UNIT_RE})?"
        rf"|"
        rf"({UNIT_RE})?(?:\s+of\s+age)?\s*(?:or|and)?\s*(?:older|above|more|greater|over))"
        rf")\b",
        s,
    )
    if m:
        num = m.group(1) or m.group(3)
        unit = m.group(2) or m.group(4) or m.group(5)
        result["minVal"] = int(num)
        unit = _normalize_unit(unit) or _extract_unit(s)
        if unit:
            result["unitText"] = unit
        return result

    return result


def remove_html_tags(text):
    # Parse the HTML string
    parsed_html = html.fromstring(text)
    # Extract the text content, effectively removing HTML tags
    clean_text = parsed_html.text_content()
    return clean_text


def add_info(sample: dict, info_url: str):
    data = requests.get(info_url)
    if data:
        data = data.json()
    else:
        return
    html_string = data.get("data", {}).get("inclusion_criteria", "")
    if html_string:
        clean_text = remove_html_tags(html_string)
        parsed_age = parse_age_string(clean_text)
        if parsed_age:
            parsed_age["@type"] = "QuantitativeValue"
            sample["developmentalStage"] = parsed_age

        parsed_race = list(parse_race_string(clean_text))
        if parsed_race:
            for race in parsed_race:
                sample.setdefault("associatedPhenotype", []).append({"@type": "DefinedTerm", "name": race})

def add_sample_quantity(sample: dict, subject_num_url: str):
    data = requests.get(subject_num_url)
    if data:
        data = data.json()
    else:
        return
    if count := data.get("data", {}).get("num_subjects"):
        sample["sampleQuantity"] = {
            "@type": "QuantitativeValue",
            "value": count,
            "unitText": "enrolled subjects",
            "unitCode": "http://purl.obolibrary.org/obo/NCIT_C207572"
        }


def parse():
    download()
    # Define the path to the json directory
    json_dir = "./json"
    # define the path to the xml directory
    xml_dir = "./xml"

    # List all files in the json directory
    json_files = [f for f in os.listdir(json_dir) if os.path.isfile(os.path.join(json_dir, f))]

    # phs002203 has a issue where all the keys are in the description field and
    # not in the configuration field we need to move them back to the configuration field
    # phs002978 has the same issue but the keys are in the studyinex field
    bad_keys = [
        "attributions",
        "consentgroups",
        "description",
        "diseases",
        "displaypublicsummary",
        "publication",
        "studyhistory",
        "studyinex",
        "studynameentrez",
        "studynamereportpage",
        "studytypes",
        "studyurls",
    ]

    def move_keys(description, description_keys, configuration):
        keys_to_move = []

        if isinstance(description, dict):
            # Collect keys to move
            for key, value in description.items():
                if key in description_keys:
                    keys_to_move.append(key)
                elif isinstance(value, dict):
                    move_keys(value, description_keys, configuration)
                elif isinstance(value, list):
                    for item in value:
                        if isinstance(item, dict):
                            move_keys(item, description_keys, configuration)

            # Move keys
            for key in keys_to_move:
                configuration[key] = description.pop(key)

    # Iterate over each file and process it
    for json_file in json_files:
        output = {
            "@type": "Dataset",
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "The Database of Genotypes and Phenotypes",
                "url": "https://www.ncbi.nlm.nih.gov/gap/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "conditionsOfAccess": "Restricted",
        }

        file_path = os.path.join(json_dir, json_file)
        with open(file_path, "r", encoding="utf-8") as file:
            data = json.load(file)
            data = data["gapexchange"]["studies"]["study"]
            assert data.get("@accession"), "Accession number is missing cannot format _id"
            if accession := data.get("@accession"):
                output["identifier"] = accession
                url = f"https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={accession}"
                output["url"] = url
                output["includedInDataCatalog"]["archivedAt"] = url
                output["_id"] = accession.split(".")[0]
                subject_num_url = f"https://dbgap.ncbi.nlm.nih.gov/beta/api/study/{accession}/subject-count/?format=json"
                info_url = f"https://dbgap.ncbi.nlm.nih.gov/beta/api/study/{accession}/info/?format=json"
                sample = {}
                add_sample_quantity(sample, subject_num_url)
                add_info(sample, info_url)
                if sample:
                    output["sample"] = sample

            if identifier := data.get("@parentstudy"):
                if identifier != accession:
                    output["isPartOf"] = {"identifier": identifier}

            if date_created := data.get("@createdate"):
                output["dateCreated"] = date_created

            if date_modified := data.get("@moddate"):
                output["dateModified"] = date_modified

            funding_list = [
                "broad institute of mit and harvard",
                "candidate gene association resource (care) - broad institute",
                "european community's seventh framework programme",
                "national institutes of health",
                "program gifts",
                "the wellcome trust",
                "medical research council (mrc) centre",
                "german federal ministry of education research",
            ]
            organization_list = [
                "consortium",
                "consortiums",
                "data repository",
                "dna repository",
                "dna sequencing",
                "genetics laboratory",
                "genome core",
                "genotyping institute",
                "immunophenotyping laboratory",
                "institute",
                "institutes",
                "institution",
                "institutions",
                "participating institutions",
                "sequence centers",
                "sequencing centers",
                "sequencing facility",
                "sequencing institute",
                "sequencing institute/company",
                "tissue core",
                "tissue donation",
                "tissue procurement",
                "tissue source",
            ]
            authors = []
            fundings = []

            # fix the bad keys
            configuration = data.get("configuration", {})
            description = configuration.get("description", {})
            move_keys(description, bad_keys, configuration)
            description = configuration.get("studyinex", {})
            move_keys(description, bad_keys, configuration)
            description = configuration.get("studyhistory", {})
            move_keys(description, bad_keys, configuration)
            configuration = {k: v for k, v in configuration.items() if v}

            authorized_access = data.get("authorizedaccess", {})
            authorized_access = {k: v for k, v in authorized_access.items() if v}

            if attributions := configuration.get("attributions", {}):
                if headers := attributions.get("header"):
                    if not isinstance(headers, list):
                        headers = [headers]
                    for header in headers:
                        title = header.get("@title")
                        name = header.get("attname")
                        institution = header.get("institution")
                        if "investigator" in title.casefold() or "fundus" in title.casefold():
                            author = {
                                "@type": "Person",
                                "name": name,
                                "affiliation": {"name": institution},
                            }
                            authors.append(author)
                        elif title.casefold() in funding_list:
                            funding = {"identifier": name, "funder": {"name": institution}}
                            fundings.append(funding)
                        elif any(
                            term in title.casefold()
                            for term in [
                                "funding",
                                "funder",
                                "funders",
                                "fund",
                                "funds",
                                "financial",
                                "grant",
                                "founding",
                            ]
                        ):
                            funding = {"identifier": name, "funder": {"name": institution}}
                            fundings.append(funding)
                        else:
                            author = {
                                "@type": ("Organization" if title.casefold() in organization_list else "Person"),
                                "name": name,
                                "affiliation": {"name": institution},
                            }
                            authors.append(author)

            if authors:
                output["author"] = authors
            if fundings:
                output["funding"] = fundings

            usage_infos = []
            # get the consent group from the authorized access
            if consent_group := configuration.get("consentgroups", {}).get("consentgroup"):
                if not isinstance(consent_group, list):
                    consent_group = [consent_group]
                for group in consent_group:
                    description = None
                    if psets := authorized_access.get("consentgroups", {}).get("participantset"):
                        if not isinstance(psets, list):
                            psets = [psets]
                        for pset in psets:
                            if pset.get("@groupnum-ref") == group.get("@groupnum"):
                                description = pset.get("uselimitation")

                    if name := group.get("@longname"):
                        usage_info = {
                            "@type": "CreativeWork",
                            "name": name,
                        }
                        if description:
                            usage_info["description"] = description
                        usage_infos.append(usage_info)

            if usage_infos:
                output["usageInfo"] = usage_infos

            xml_file_path = os.path.join(xml_dir, output["_id"] + ".xml")
            # Parse the XML file
            with open(xml_file_path, "rb") as xml_file:
                root = html.fromstring(xml_file.read())

            # Find the description elements
            for description in root.xpath("//description"):
                output["description"] = output.get("description", "") + description.text_content()

            for description in root.xpath("//studyinex"):
                output["description"] = output.get("description", "") + description.text_content()

            for description in root.xpath("//studyhistory"):
                output["description"] = output.get("description", "") + description.text_content()

            health_conditions = []
            if diseases := configuration.get("diseases", {}).get("disease"):
                if not isinstance(diseases, list):
                    diseases = [diseases]
                for disease in diseases:
                    health_condition = {
                        "inDefinedTermSet": disease.get("@vocab_source"),
                        "name": disease.get("@vocab_term"),
                    }
                    health_conditions.append(health_condition)
            if health_conditions:
                output["healthCondition"] = health_conditions

            pmids = ""
            if publications := configuration.get("publications", {}).get("publication"):
                if not isinstance(publications, list):
                    publications = [publications]
                for publication in publications:
                    if pmid := publication.get("pubmed", {}).get("@pmid"):
                        pmids += pmid if not pmids else ", " + pmid

            if pmids:
                output["pmids"] = pmids

            if name := configuration.get("studynameentrez"):
                output["name"] = name
            elif name := configuration.get("studynamereportpage"):
                output["name"] = name

            measurement_techniques = []

            if study_types := configuration.get("studytypes", {}).get("studytype"):
                if not isinstance(study_types, list):
                    study_types = [study_types]
                for study_type in study_types:
                    measurement_technique = {
                        "name": study_type,
                    }
                    measurement_techniques.append(measurement_technique)
            if measurement_techniques:
                output["measurementTechnique"] = measurement_techniques

            is_related_to = []
            if study_urls := configuration.get("studyurls", {}).get("studyurl"):
                if not isinstance(study_urls, list):
                    study_urls = [study_urls]
                for study_url in study_urls:
                    name = study_url.get("@name")
                    url = study_url.get("@url")
                    is_related_to.append({"name": name, "url": url})

            if is_related_to:
                output["isRelatedTo"] = is_related_to

            # all conditions of access are restricted by default
            # if consent_groups := authorized_access.get("consentgroups", {}).get("participantset"):
            #     if not isinstance(consent_groups, list):
            #         consent_groups = [consent_groups]
            #     for consent_group in consent_groups:
            #         if irbrequired := consent_group.get("irbrequired"):
            #             if irbrequired.casefold() == "yes":
            #                 output["conditionsOfAccess"] = "Restricted"
            #             else:
            #                 embargo_length = int(authorized_access.get("policy", {}).get("embargolength"))
            #                 display_public_summary = configuration.get("displaypublicsummary").casefold()
            #                 au_display_public_summary = (
            #                     authorized_access.get("policy", {}).get("displaypublicsummary").casefold()
            #                 )

            #                 if embargo_length > 0:
            #                     output["conditionsOfAccess"] = "Embargoed"
            #                 elif embargo_length == 0 and display_public_summary == "no":
            #                     output["conditionsOfAccess"] = "Closed"
            #                 elif (
            #                     embargo_length == 0
            #                     and display_public_summary == "yes"
            #                     and au_display_public_summary == "yes"
            #                 ):
            #                     output["conditionsOfAccess"] = "Open"
            #                 elif display_public_summary != au_display_public_summary:
            #                     output["conditionsOfAccess"] = "Closed"
            #             break

            policy = authorized_access.get("policy", {})
            acknowledgementtext = policy.get("acknowledgementtext", {})
            documentset = policy.get("documentset", {})
            if documentset:
                datausecertificate = documentset.get("datausecertificate", {})
                if license := datausecertificate.get("@filepath"):
                    output["license"] = license
            elif license := acknowledgementtext.get("para"):
                output["license"] = license
        yield output
