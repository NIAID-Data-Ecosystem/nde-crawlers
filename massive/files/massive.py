import logging
import urllib.parse
from datetime import datetime

import requests
from bs4 import BeautifulSoup

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def fetch_dataset_ids():
    base_url = "https://massive.ucsd.edu/ProteoSAFe/QueryDatasets"
    query_params = {
        "pageSize": 1000,
        "offset": 0,
        "query": urllib.parse.quote_plus('#{"query":{},"table_sort_history":"createdMillis_dsc"}'),
    }

    dataset_ids = []
    page = 1

    while True:
        response = requests.get(base_url, params=query_params)
        if response.status_code == 200:
            logger.info(f"Processed page {page}")
            data = response.json()
            row_data = data.get("row_data", [])
            if not row_data:
                break  # No more data to process
            for item in row_data:
                dataset_ids.append(item["dataset"])
            query_params["offset"] += query_params["pageSize"]  # Increment offset for next page
            page += 1
        else:
            raise Exception(f"Failed to fetch dataset IDs, status code: {response.status_code}")

    return dataset_ids


def fetch_dataset(dataset_id):
    url = f"https://massive.ucsd.edu/ProteoSAFe/QueryDatasets?pageSize=30&offset=0&query=%7B%22title_input%22:%22{dataset_id}%22%7D"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Failed to fetch dataset {dataset_id}, status code: {response.status_code}")


def fetch_html_content(task):
    url = f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task={task}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Failed to fetch HTML content for task {task}, status code: {response.status_code}")


def parse_html_for_doi_and_license(html_content):
    soup = BeautifulSoup(html_content, "html.parser")
    p_tag = soup.find("p")

    doi = None
    license_url = None

    if p_tag:
        doi_link = p_tag.find("a", href=lambda href: href and "doi.org" in href)
        if doi_link:
            doi = doi_link.text.strip().replace("doi:", "")

        license_link = p_tag.find("a", href=lambda href: href and "creativecommons.org" in href)
        if license_link:
            license_url = license_link["href"].strip()

    logger.info(f"DOI: {doi}, License URL: {license_url}")

    return doi, license_url


def parse_dataset(json_data):
    output = {
        "includedInDataCatalog": {
            "name": "MassIVE",
            "url": "https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp",
            "@type": "Dataset",
            "versionDate": datetime.today().strftime("%Y-%m-%d"),
        },
        "@type": "Dataset",
    }

    for item in json_data.get("row_data", []):
        if identifier := item.get("dataset"):
            output["identifier"] = identifier
            output["_id"] = "massive_" + identifier.lower()

        if task := item.get("task"):
            output["url"] = "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=" + task

            # Fetch and parse HTML content for DOI and License
            html_content = fetch_html_content(task)
            doi, license_url = parse_html_for_doi_and_license(html_content)
            if doi:
                output["doi"] = doi
            if license_url:
                output["license"] = license_url

        if repo_path := item.get("repo_path"):
            path = repo_path.split("-")[-1]
            if path != "/data/massive":
                output["distribution"] = {"contentUrl": f"ftp://massive.ucsd.edu/{path}/{identifier}"}
            else:
                output["distribution"] = {"contentUrl": f"ftp://massive.ucsd.edu/v08/{identifier}"}

        if title := item.get("title"):
            output["name"] = title

        if sdPublisher := item.get("site"):
            output["sdPublisher"] = {"name": sdPublisher}

        if description := item.get("description"):
            if description != "null":
                output["description"] = description

        if keywords := item.get("keywords"):
            if "###" in keywords:
                output["keywords"] = keywords.split("###")
            else:
                output["keywords"] = keywords.split(", ")

        if date_created := item.get("create_time"):
            parsed_date = datetime.strptime(date_created, "%Y-%m-%d %H:%M:%S.%f").strftime("%Y-%m-%d")
            output["dateCreated"] = parsed_date

        if measurement_techniques := item.get("instrument_resolved"):
            mt_list = measurement_techniques.split("###")
            mt_names = [{"name": sp for sp in mt_list}]
            output["measurementTechnique"] = mt_names

        if species := item.get("species_resolved"):
            species_list = species.split("###")
            species_names = [{"name": sp.split(" (NCBITaxon:")[0]} for sp in species_list]
            output["species"] = species_names

        if author_name := item.get("pis"):
            pis_list = []
            for pi in author_name:
                pi_data = {}
                if name := pi.get("name"):
                    pi_data["name"] = name
                if email := pi.get("email"):
                    pi_data["email"] = email
                if institution := pi.get("institution"):
                    pi_data["affiliation"] = {"name": institution}
                pis_list.append(pi_data)
            output["author"] = pis_list

        if publications := item.get("publications"):
            pub_list = []
            for pub in publications:
                pub_data = {}
                if citation_id := pub.get("id"):
                    pub_data["identifier"] = citation_id
                if authors := pub.get("authors"):
                    pub_data["author"] = {"name": authors}
                if title := pub.get("title"):
                    pub_data["name"] = title
                if doi := pub.get("citation"):
                    pub_data["doi"] = doi
                if description := pub.get("abstract"):
                    pub_data["description"] = description
                pub_list.append(pub_data)
            output["citation"] = pub_list

        if privacy := item.get("privacy"):
            if privacy == "Public":
                output["conditionsOfAccess"] = "Open"
            else:
                output["conditionsOfAccess"] = "Closed"

        if "measurementTechnique" not in output:
            output["measurementTechnique"] = [
                {
                    "name": "Mass Spectrometry",
                    "url": "https://ontobee.org/ontology/MMO?iri=http://purl.obolibrary.org/obo/MMO_0000534",
                    "identifier": "MMO_0000534",
                }
            ]
        else:
            output["measurementTechnique"].append(
                {
                    "name": "Mass Spectrometry",
                    "url": "https://ontobee.org/ontology/MMO?iri=http://purl.obolibrary.org/obo/MMO_0000534",
                    "identifier": "MMO_0000534",
                }
            )

    return output


def parse():
    dataset_ids = fetch_dataset_ids()
    logger.info(f"Found {len(dataset_ids)} datasets")
    count = 0
    for dataset_id in dataset_ids:
        count += 1
        if count % 100 == 0:
            logger.info(f"Processed {count} datasets")
        try:
            dataset_json = fetch_dataset(dataset_id)
            parsed_dataset = parse_dataset(dataset_json)
            if "_id" not in parsed_dataset:
                logger.info(f"Dataset {dataset_id} has no identifier. Skipping...")
                continue
            yield parsed_dataset
        except Exception as e:
            logger.info(f"Error processing dataset {dataset_id}: {e}")
