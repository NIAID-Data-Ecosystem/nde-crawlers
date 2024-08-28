import datetime
import logging
import re

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


COOKIE = 'user="2|1:0|10:1722831328|4:user|32:ImR5bGFud2VsemVsQGdtYWlsLmNvbSI=|65fe32998c777a26ed86ecbbba7253015c5c4d9e565fb0980b248e84f28c3371"'


def retrieve_study_metadata():
    logger.info("Retrieving study metadata from Qiita")

    url = "https://qiita.ucsd.edu/study/list_studies/?&user=dylanwelzel@gmail.com&visibility=public&sEcho=344&query=&_=1665437236850"
    headers = {"Cookie": COOKIE}
    r = requests.get(url, headers=headers)
    return r.json()["aaData"]


def retrieve_study_samples(study_id):
    url = f"https://qiita.ucsd.edu/study/description/sample_template/columns/?study_id={study_id}"
    headers = {"Cookie": COOKIE}
    r = requests.get(url, headers=headers)

    result = []

    for column_name in r.json()["values"]:
        url = f"https://qiita.ucsd.edu/study/description/sample_template/columns/?study_id={study_id}&column={column_name}"
        headers = {"Cookie": COOKIE}
        r = requests.get(url, headers=headers)

        sample_dict = {"name": column_name, "values": []}

        for value in r.json()["values"]:
            count = r.json()["values"].count(value)
            value_dict = {"@type": "DefinedTerm", "name": value, "count": count}
            if value_dict not in sample_dict["values"]:
                sample_dict["values"].append(value_dict)

        result.append(sample_dict)

    return result


def parse():
    studies = retrieve_study_metadata()
    count = 0
    logger.info("Parsing study metadata")
    for study in studies:
        count += 1
        logger.info(f"Parsing study {count} of {len(studies)}")
        if count % 50 == 0:
            logger.info("Parsed %s studies", count)

        output = {
            "includedInDataCatalog": {
                "@type": "Dataset",
                "name": "Qiita",
                "url": "https://qiita.ucsd.edu/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Dataset",
        }

        if study_abstract := study.get("study_abstract"):
            output["description"] = study_abstract

        if study_id := study.get("study_id"):
            output["url"] = f"https://qiita.ucsd.edu/study/description/{study_id}"
            output["_id"] = f"qiita_{study_id}"

            # TODO use helper to import sample metadata to proper mapping
            # study_samples = retrieve_study_samples(study_id)
            # output['material'] = study_samples

        if study_title := study.get("study_title"):
            output["name"] = study_title

        if study_tags := study.get("study_tags"):
            output["keywords"] = study_tags

        # if owner := study.get('owner'):
        #     print(owner)

        if pi := study.get("pi"):
            name = re.findall(r">(.*?)<", pi)
            output["author"] = {"name": name[0]}

        if pubs := study.get("pubs"):
            publications = re.findall(r">(.*?)<", pubs)
            doi_list = []
            pmid_list = []
            for publication in publications:
                if "10." in publication:
                    if "http" in publication:
                        doi_list.append(publication.replace("https://doi.org/", ""))
                    else:
                        doi_list.append(publication)
                if len(publication) == 8:
                    pmid_list.append(publication)
            if len(pmid_list):
                output["pmids"] = ",".join(pmid_list)
            if len(doi_list):
                output["doi"] = doi_list

        if ebi_study_accession := study.get("ebi_study_accession"):
            output["mainEntityOfPage"] = f"https://www.ebi.ac.uk/ena/browser/view/{ebi_study_accession}"

        yield output
