import json
import logging
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import dateutil.parser as parser
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from sickle import Sickle

# used to test single record
# record = sickle.GetRecord(
#     identifier='doi:10.17632/s6y4yzssy6.1', metadataPrefix='datacite')
# properties = ['creator', 'title', 'publisher', 'description', 'subject',
#               'contributor', 'type', 'identifier', 'rights', 'relation', 'date']
# publishers = ['Mendeley Data', 'Elsevier BV',
#               'San Raffaele Open Research Data Repository', 'Teesside University', 'MURI/AUSMURI Project', 'University of La Laguna', 'Simon Bolivar University', 'Bicocca Open Archive Research Data', 'University of Canberra', 'University of Nottingham Ningbo China', 'Meiji University']

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

url_count = 0

# Function that returns a response using requests, used in ThreadPoolExecutor later on

session = requests.Session()
retries = Retry(total=5, backoff_factor=0.1, status_forcelist=[500, 502, 503, 504])
adapter = HTTPAdapter(max_retries=retries)
session.mount("http://", adapter)
session.mount("https://", adapter)

url_count = 0


def get_url(url):
    global url_count, session
    url_count += 1
    if url_count % 1000 == 0:
        logger.info(f"Retrieved {url_count} Metadata Sources")
    try:
        response = session.get(url, timeout=60)
        response.raise_for_status()
        return response
    except requests.exceptions.RequestException as e:
        logger.error(f"Error retrieving URL {url}: {e}")
        return None


def parse():
    sickle = Sickle("https://data.mendeley.com/oai", max_retries=3)
    records = sickle.ListRecords(metadataPrefix="oai_dc", ignore_deleted=True)

    logger.info("Retrieving Dataset ids")

    urls = []
    count = 0
    metadata_count = 0
    # Our first step is to grab all the individual dataset ids using the OAI-PMH and save them to a list
    while True:
        try:
            count += 1

            if count % 1000 == 0:
                logger.info(f"Retrieved {count} ids")

            record = records.next()
            metadata = record.metadata

            if relation := metadata.get("relation"):
                (
                    urls.append(
                        "https://data.mendeley.com/api/datasets-v2/datasets/"
                        + relation[0].split("/")[-1]
                        + "?fields=repository.*"
                    )
                )

        except StopIteration:
            logger.info("Finished Retrieving ids. Total ids: %s", count)
            # if StopIteration is raised, break from loop
            break

    logger.info("Retrieving Metadata Sources")

    # After we have the ids we use the function declared above and map the list of ids we've obtained from the OAI-PMH and ping their api for metadata.
    with ThreadPoolExecutor(max_workers=3) as pool:
        response_list = list(pool.map(get_url, urls))
        # TODO RESPONSE_LIST IS WHAT WE WANT TO CACHE

    logger.info("Finished Retrieving Metadata Sources. Total Metadata Sources: %s", url_count)

    logger.info("Parsing records")

    # Finally we handle the transformations after retrieving all the metadata.
    for response in response_list:
        if response and response.status_code == 200:
            metadata_count += 1
            if metadata_count % 1000 == 0:
                logger.info(f"Parsed {metadata_count} records")
            metadata = json.loads(response.text)
            output = {
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "Mendeley",
                    "versionDate": datetime.today().isoformat(),
                    "url": "https://data.mendeley.com/",
                },
                "@type": "Dataset",
            }
            if mendeley_id := metadata.get("id"):
                output["identifier"] = mendeley_id
                output["_id"] = "Mendeley_" + mendeley_id
            if doi := metadata.get("doi"):
                output["doi"] = doi["id"]
            if name := metadata.get("name"):
                output["name"] = name
            if description := metadata.get("description"):
                output["description"] = description
            if contributors := metadata.get("contributors"):
                authors = []
                for contributor in contributors:
                    author_obj = {}
                    if first_name := contributor.get("first_name"):
                        author_obj["givenName"] = first_name
                    if last_name := contributor.get("last_name"):
                        author_obj["familyName"] = last_name
                    if first_name and last_name:
                        author_obj["name"] = first_name + " " + last_name
                    authors.append(author_obj)
                output["author"] = authors
            if files := metadata.get("files"):
                distribution_obj = {}
                for file in files:
                    if filename := file.get("filename"):
                        distribution_obj["name"] = filename
                    if content_details := file.get("content_details"):
                        if content_type := content_details.get("content_type"):
                            distribution_obj["encodingFormat"] = content_type
                        if size := content_details.get("size"):
                            distribution_obj["contentSize"] = size
                        if created_date := content_details.get("created_date"):
                            distribution_obj["dateModified"] = parser.parse(created_date).strftime("%Y-%m-%d")
                        if download_url := content_details.get("download_url"):
                            distribution_obj["contentUrl"] = download_url
                    if metrics := content_details.get("metrics"):
                        if downloads := metrics.get("downloads"):
                            distribution_obj["aggregateRating"] = {
                                "@type": "AggregateRating",
                                "reviewCount": downloads,
                                "reviewAspect": "downloads",
                            }
                        if previews := metrics.get("previews"):
                            distribution_obj["aggregateRating"] = {
                                "@type": "AggregateRating",
                                "reviewCount": previews,
                                "reviewAspect": "downloads",
                            }
                output["distribution"] = distribution_obj

            citation_list = []
            if articles := metadata.get("articles"):
                for article in articles:
                    citation_obj = {}
                    if title := article.get("title"):
                        citation_obj["name"] = title
                    if doi := article.get("doi"):
                        citation_obj["doi"] = doi
                    if issueNumber := article.get("issueNumber"):
                        citation_obj["issn"] = issueNumber
                    if name := article.get("name"):
                        citation_obj["journalName"] = name
                    if url := article.get("url"):
                        citation_obj["url"] = url
                    if citation_obj:
                        citation_list.append(citation_obj)

            if categories := metadata.get("categories"):
                category_list = []
                for category in categories:
                    if label := category.get("label"):
                        category_list.append(label)
                output["keywords"] = category_list

            if publish_date := metadata.get("publish_date"):
                output["datePublished"] = parser.parse(publish_date).strftime("%Y-%m-%d")

            if related_links := metadata.get("related_links"):
                for link in related_links:
                    citation_obj = {}
                    if citation_type := link.get("type"):
                        citation_obj["@type"] = citation_type
                    if url := link.get("url"):
                        citation_obj["url"] = url
                    if href := link.get("href"):
                        citation_obj["url"] = href
                    if citation_obj:
                        citation_list.append(citation_obj)
            if citation_list:
                output["citation"] = citation_list

            if modified_on := metadata.get("modified_on"):
                output["dateModified"] = parser.parse(modified_on).strftime("%Y-%m-%d")

            if links := metadata.get("links"):
                output["url"] = links["view"]
                output["includedInDataCatalog"]["archivedAt"] = links["view"]
            if repository := metadata.get("repository"):
                output["sdPublisher"] = {"name": repository["name"]}

            if license_info := metadata.get("data_licence"):
                output["license"] = license_info.get("url")

            if funders := metadata.get("funders"):
                funding = []
                for funder in funders:
                    funding_entry = {}

                    funder_info = {}
                    if name := funder.get("name"):
                        funder_info["name"] = name
                    if identity := funder.get("identity"):
                        funder_info["identifier"] = identity

                    if funder_info:
                        funding_entry["funder"] = funder_info

                    if grant_id := funder.get("grant_id"):
                        funding_entry["identifier"] = grant_id

                    if funding_entry:
                        funding.append(funding_entry)

                if funding:
                    output["funding"] = funding

            yield output
