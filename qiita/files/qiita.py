import datetime
import logging
import re

import requests

from api_secret import QIITA_EMAIL, QIITA_PASSWORD

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

BASE_URL = "https://qiita.ucsd.edu"


def get_authenticated_session():
    """Login to Qiita and return an authenticated requests session."""
    session = requests.Session()

    response = session.post(
        f"{BASE_URL}/auth/login/",
        data={
            "username": QIITA_EMAIL,
            "password": QIITA_PASSWORD,
        },
        headers={
            "Referer": f"{BASE_URL}/auth/login/",
            "Origin": BASE_URL,
        },
    )

    if "user" not in session.cookies:
        raise RuntimeError("Failed to authenticate with Qiita — check credentials in api_secret.py")

    logger.info("Successfully authenticated with Qiita")
    return session


session = get_authenticated_session()


def retrieve_study_metadata():
    logger.info("Retrieving study metadata from Qiita")

    url = f"{BASE_URL}/study/list_studies/?&user={QIITA_EMAIL}&visibility=public&sEcho=344&query=&_=1665437236850"
    r = session.get(url)
    return r.json()["aaData"]


def retrieve_study_samples(study_id):
    url = f"{BASE_URL}/study/description/sample_template/columns/?study_id={study_id}"
    r = session.get(url)

    result = []

    for column_name in r.json()["values"]:
        url = f"{BASE_URL}/study/description/sample_template/columns/?study_id={study_id}&column={column_name}"
        r = session.get(url)

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
                "@type": "DataCatalog",
                "name": "Qiita",
                "url": f"{BASE_URL}/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Dataset",
        }

        if study_abstract := study.get("study_abstract"):
            output["description"] = study_abstract

        if study_id := study.get("study_id"):
            output["url"] = f"{BASE_URL}/study/description/{study_id}"
            output["includedInDataCatalog"]["archivedAt"] = f"{BASE_URL}/study/description/{study_id}"
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
