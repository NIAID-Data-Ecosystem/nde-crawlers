import datetime
import json
import logging
import os

from download_files import download
from lxml import html

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


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
                "@type": "Dataset",
                "name": "The Database of Genotypes and Phenotypes",
                "url": "https://www.ncbi.nlm.nih.gov/gap/",
                "versionDate": datetime.date.today().isoformat(),
            },
        }

        file_path = os.path.join(json_dir, json_file)
        with open(file_path, "r", encoding="utf-8") as file:
            data = json.load(file)
            data = data["gapexchange"]["studies"]["study"]
            assert data.get("@accession"), "Accession number is missing cannot format _id"
            if accession := data.get("@accession"):
                output["identifier"] = accession
                output["url"] = f"https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={accession}"
                output["_id"] = accession.split(".")[0]

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
                                "@type": ("Organization" if title.casefold() in organization_list else "person"),
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

            if consent_groups := authorized_access.get("consentgroups", {}).get("participantset"):
                if not isinstance(consent_groups, list):
                    consent_groups = [consent_groups]
                for consent_group in consent_groups:
                    if irbrequired := consent_group.get("irbrequired"):
                        if irbrequired.casefold() == "yes":
                            output["conditionsOfAccess"] = "Restricted"
                        else:
                            embargo_length = int(authorized_access.get("policy", {}).get("embargolength"))
                            display_public_summary = configuration.get("displaypublicsummary").casefold()
                            au_display_public_summary = (
                                authorized_access.get("policy", {}).get("displaypublicsummary").casefold()
                            )

                            if embargo_length > 0:
                                output["conditionsOfAccess"] = "Embargoed"
                            elif embargo_length == 0 and display_public_summary == "no":
                                output["conditionsOfAccess"] = "Closed"
                            elif (
                                embargo_length == 0
                                and display_public_summary == "yes"
                                and au_display_public_summary == "yes"
                            ):
                                output["conditionsOfAccess"] = "Open"
                            elif display_public_summary != au_display_public_summary:
                                output["conditionsOfAccess"] = "Closed"
                        break

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
