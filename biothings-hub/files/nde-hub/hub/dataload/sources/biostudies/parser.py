import datetime
import json
import os
import re

import dateutil

try:
    from config import logger
except ImportError:
    import logging

    logger = logging.getLogger(__name__)

known_attributes = [
    "reviewtype",
    "template",
    "rootpath",
    "rembi_pagetab conversion script version",
]
known_section_attributes = [
    "organ",
    "peak identification",
    "software",
    "cell type",
    "dataset size (in gb)",
    "method",
]
known_subsection_types = [
    "specimen",
    "study component",
    "samples",
    "minseqe score",
    "image analysis",
    "miame score",
    "annotations",
    "sample",
]
missing_subattributes = {}
missing_attributes = {}
count_links = 0


def parse_file(doc):
    global count_links
    global missing_attributes
    global missing_subattributes
    global known_attributes
    global known_section_attributes
    global known_subsection_types

    accno = doc.pop("accno")

    output = {
        "_id": f"biostudies_{accno}",
        "url": f"https://www.ebi.ac.uk/studies/{accno}",
        "includedInDataCatalog": {
            "@type": "Dataset",
            "name": "BioStudies",
            "url": "https://www.ebi.ac.uk/biostudies/",
            "versionDate": datetime.date.today().isoformat(),
        },
    }

    pattern = re.compile(r"annotations \d+")
    subsection_pattern = re.compile(r"figure \d+")

    # handle attributes
    if attributes := doc.pop("attributes", None):
        for attribute in attributes:
            key = attribute.pop("name").casefold()
            value = attribute.pop("value")

            if key == "attachto":
                output["url"] = f"https://www.ebi.ac.uk/{value}/studies/{accno}"
            elif key == "releasedate":
                output["releaseDate"] = dateutil.parser.parse(value, ignoretz=True).date().isoformat()
            elif key == "doi":
                output["doi"] = value
            elif key == "title":
                output["name"] = value
            elif key in known_attributes:
                pass
            else:
                if key not in missing_attributes:
                    missing_attributes[key] = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accno}"
            if attribute:
                logger.info(
                    f"This attribute has something other than name/value pair. name: {key}, url: https://www.ebi.ac.uk/biostudies/api/v1/studies/{accno}"
                )

    if doc.get("section"):
        doc["section"].pop("type", None)
        doc["section"].pop("accno", None)
        if section_attributes := doc["section"].pop("attributes", None):
            for attribute in section_attributes:
                key = attribute.pop("name").casefold()
                value = attribute.pop("value")

                if key == "title":
                    if title := output.get("name"):
                        assert value == title
                    else:
                        output["name"] = value
                elif (
                    key == "abstract" or key == "description" or key == "acknowledgements" or key == "funding statement"
                ):
                    output["description"] = output.get("description") + value if output.get("description") else value
                elif key == "keywords" or key == "keyword":
                    if output.get("keywords"):
                        output["keywords"].append(value)
                    else:
                        output["keywords"] = [value]
                elif key == "license":
                    output["license"] = value
                    if valquals := attribute.get("valqual"):
                        for valqual in valquals:
                            if valqual.get("name").casefold() == "url":
                                output["license"] = valqual.get("value")
                            else:
                                logger.info(
                                    f"valqual has something thats not just url https://www.ebi.ac.uk/biostudies/api/v1/studies/{accno}"
                                )
                    attribute.pop("valqual", None)
                elif key == "organism":
                    species = {"name": value}
                    if output.get("species"):
                        output["species"].append(species)
                    else:
                        output["species"] = [species]
                elif key == "method" or key == "study type" or "experimental design":
                    attribute.pop("valqual", None)
                    mt = {"name": value}
                    if output.get("measurementTechnique"):
                        output["measurementTechnique"].append(mt)
                    else:
                        output["measurementTechnique"] = [mt]
                elif key == "experimental factor":
                    vm = {"name": value}
                    if output.get("variableMeasured"):
                        output["variableMeasured"].append(vm)
                    else:
                        output["variableMeasured"] = [vm]
                elif key in known_section_attributes or pattern.match(key):
                    pass
                else:
                    if key not in missing_attributes:
                        missing_attributes[key] = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accno}"
                if attribute:
                    logger.info(
                        f"This attribute has something other than name/value pair. name: {key}, url: https://www.ebi.ac.uk/biostudies/api/v1/studies/{accno}"
                    )
        if files := doc["section"].pop("files", None):
            if isinstance(files[0], list):
                files = files[0]
            for file in files:
                file.get("path")
                if path := file.get("path"):
                    dist_url = {"contentUrl": f"https://www.ebi.ac.uk/biostudies/files/{accno}/{path}"}
                    if output.get("distribution"):
                        output["distribution"].append(dist_url)
                    else:
                        output["distribution"] = [dist_url]

        if links := doc["section"].pop("links", None):
            count_links += 1
            if isinstance(links[0], list):
                links = links[0]
            for link in links:
                irt = {}
                if url := link.get("url"):
                    irt["url"] = url
                if link.get("attributes"):
                    for attribute in link.get("attributes"):
                        if attribute.get("name").casefold() == "description":
                            irt["name"] = attribute.get("value")
                        elif attribute.get("name").casefold() == "type":
                            pass
                        else:
                            logger.info(
                                f"Links has more than just description as an attribute https://www.ebi.ac.uk/biostudies/api/v1/studies/{accno}"
                            )
                if output.get("isRelatedTo"):
                    output["isRelatedTo"].append(irt)
                else:
                    output["isRelatedTo"] = [irt]

        if subsections := doc["section"].pop("subsections", None):
            for subsection in subsections:
                if isinstance(subsection, dict):
                    sub_accno = subsection.pop("accno", None)
                    type = subsection.pop("type").casefold()
                    attributes = subsection.pop("attributes", None)
                    if type == "author":
                        author = {}
                        if attributes:
                            for attribute in attributes:
                                key = attribute.pop("name").casefold()
                                value = attribute.pop("value")
                                if key == "name":
                                    author["name"] = value
                                elif key == "orcid":
                                    author["identifier"] = value
                                elif key == "affiliation":
                                    author["affiliation"] = {"name": value}
                            if output.get("author"):
                                output["author"].append(author)
                            else:
                                output["author"] = [author]
                    elif type == "organization" or type == "organisation":
                        if sub_accno:
                            for author in output.get("author"):
                                if author.get("affiliation") and author["affiliation"].get("name") == sub_accno:
                                    for attribute in attributes:
                                        key = attribute.get("name").casefold()
                                        value = attribute.get("value")
                                        if key == "name":
                                            author["affiliation"]["name"] = value
                                        elif key == "rorid":
                                            author["affiliation"]["identifier"] = value
                    elif type == "funding":
                        funding = {}
                        if attributes:
                            for attribute in attributes:
                                key = attribute.pop("name").casefold()
                                value = attribute.pop("value")
                                if key == "agency":
                                    funding["funder"] = {"name": value}
                                elif key == "grant_id":
                                    funding["identifier"] = value
                            if output.get("funding"):
                                output["funding"].append(funding)
                            else:
                                output["funding"] = [funding]
                    elif type == "biosample":
                        if attributes:
                            for attribute in attributes:
                                if attribute.get("name") and attribute.get("value"):
                                    key = attribute.pop("name").casefold()
                                    value = attribute.pop("value")
                                if key == "organism":
                                    species = {"name": value}
                                    if output.get("species"):
                                        output["species"].append(species)
                                    else:
                                        output["species"] = [species]
                                elif key == "experimental variable" or key == "extrinsic variable":
                                    vm = {"name": value}
                                    if output.get("variableMeasured"):
                                        output["variableMeasured"].append(vm)
                                    else:
                                        output["variableMeasured"] = [vm]
                    elif type == "image acquisition" or type == "assays and data":
                        if attributes:
                            for attribute in attributes:
                                if attribute.get("name") and attribute.get("value"):
                                    key = attribute.pop("name").casefold()
                                    value = attribute.pop("value")
                                if key == "imaging method" or key == "technology" or key == "assay by molecule":
                                    mt = {"name": value}
                                    if output.get("measurementTechnique"):
                                        output["measurementTechnique"].append(mt)
                                    else:
                                        output["measurementTechnique"] = [mt]
                    elif type == "publication":
                        if subsection.get("links"):
                            count_links += 1
                        if attributes:
                            citation = {}
                            for attribute in attributes:
                                key = attribute.pop("name").casefold()
                                value = attribute.pop("value")
                                if key == "doi":
                                    citation["doi"] = value
                                elif key == "pubmed id":
                                    citation["pmid"] = value
                            if citation:
                                if output.get("citation"):
                                    output["citation"].append(citation)
                                else:
                                    output["citation"] = [citation]
                    elif type in known_subsection_types or subsection_pattern.match(type):
                        pass
                    else:
                        if type not in missing_subattributes:
                            missing_subattributes[type] = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accno}"
    return output


def parse_file_dir(input_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith(".json"):
            with open(os.path.join(input_dir, filename), "r") as infile:
                doc = json.load(infile)
                yield parse_file(doc)
    logger.info(f"Missing attributes: {missing_attributes}")
    logger.info(f"Missing subattributes: {missing_subattributes}")
    logger.info(f"Count links: {count_links}")
