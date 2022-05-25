import json
import time
import logging
from datetime import datetime
import requests
from sickle import Sickle
from concurrent.futures import ThreadPoolExecutor

# used to test single record
# record = sickle.GetRecord(
#     identifier='doi:10.17632/s6y4yzssy6.1', metadataPrefix='datacite')
# properties = ['creator', 'title', 'publisher', 'description', 'subject',
#               'contributor', 'type', 'identifier', 'rights', 'relation', 'date']
# publishers = ['Mendeley Data', 'Elsevier BV',
#               'San Raffaele Open Research Data Repository', 'Teesside University', 'MURI/AUSMURI Project', 'University of La Laguna', 'Simon Bolivar University', 'Bicocca Open Archive Research Data', 'University of Canberra', 'University of Nottingham Ningbo China', 'Meiji University']

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')

url_count = 0


def get_url(url):
    global url_count
    url_count += 1
    if url_count % 1000 == 0:
        logger.info(f"Retrieved {url_count} Metadata Sources")
    return requests.get(url)


def parse():
    start = time.time()

    sickle = Sickle('https://data.mendeley.com/oai', max_retries=3)
    records = sickle.ListRecords(
        metadataPrefix='oai_dc', ignore_deleted=True)

    logger.info(f"Retrieving Dataset ids")

    urls = []
    count = 0
    metadata_count = 0
    while True:
        try:
            count += 1

            if count % 1000 == 0:
                logger.info(f"Retrieved {count} ids")

            record = records.next()
            metadata = record.metadata

            if relation := metadata.get('relation'):
                urls.append('https://data.mendeley.com/api/datasets-v2/datasets/' +
                            relation[0].split('/')[-1])

        except StopIteration:
            logger.info("Finished Retrieving ids. Total ids: %s", count)
            # if StopIteration is raised, break from loop
            break

    logger.info(f"Retrieving Metadata Sources")

    with ThreadPoolExecutor(max_workers=10) as pool:
        response_list = list(pool.map(get_url, urls))

    logger.info(
        "Finished Retrieving Metadata Sources. Total Metadata Sources: %s", url_count)

    logger.info(f"Parsing records")
    for response in response_list:
        if response.status_code == 200:
            metadata_count += 1
            if metadata_count % 1000 == 0:
                logger.info(f"Parsed {metadata_count} records")
            metadata = json.loads(response.text)
            output = {"includedInDataCatalog":
                      {"name": "Mendeley",
                       'versionDate': datetime.today().isoformat(),
                       'url': "https://data.mendeley.com/"},
                      "@type": "Dataset"
                      }
            if id := metadata.get('id'):
                output['identifier'] = id
                output['_id'] = 'Mendeley_' + id
            if doi := metadata.get('doi'):
                output['doi'] = doi['id']
            if name := metadata.get('name'):
                output['name'] = name
            if description := metadata.get('description'):
                output['description'] = description
            if version := metadata.get('version'):
                output['version'] = version
            if contributors := metadata.get('contributors'):
                authors = []
                for contributor in contributors:
                    author_obj = {}
                    if first_name := contributor.get('first_name'):
                        author_obj['givenName'] = first_name
                    if last_name := contributor.get('last_name'):
                        author_obj['familyName'] = last_name
                    if first_name and last_name:
                        author_obj['name'] = first_name + " " + last_name
                    authors.append(author_obj)
                output['author'] = authors
            if files := metadata.get('files'):
                distribution_obj = {}
                for file in files:
                    if filename := file.get('filename'):
                        distribution_obj['name'] = filename
                    if content_details := file.get('content_details'):
                        if content_type := content_details.get('content_type'):
                            distribution_obj['encodingFormat'] = content_type
                        if size := content_details.get('size'):
                            distribution_obj['contentSize'] = size
                        if created_date := content_details.get('created_date'):
                            distribution_obj['dateModified'] = created_date
                        if download_url := content_details.get('download_url'):
                            distribution_obj['contentUrl'] = download_url
                    if metrics := content_details.get('metrics'):
                        if downloads := metrics.get('downloads'):
                            distribution_obj['aggregateRating'] = {
                                "@type": "AggregateRating",
                                "reviewCount": downloads,
                                "reviewAspect": "downloads"
                            }
                        if previews := metrics.get('previews'):
                            distribution_obj['aggregateRating'] = {
                                "@type": "AggregateRating",
                                "reviewCount": previews,
                                "reviewAspect": "downloads"
                            }
                output['distribution'] = distribution_obj

            if articles := metadata.get('articles'):
                citation_obj = {}
                for article in articles:
                    if title := article.get('title'):
                        citation_obj['name'] = title
                    if doi := article.get('doi'):
                        citation_obj['doi'] = doi
                    if issueNumber := article.get('issueNumber'):
                        citation_obj['issn'] = issueNumber
                    if name := article.get('name'):
                        citation_obj['journalName'] = name
                    if url := article.get('url'):
                        citation_obj['url'] = url
                output['citation'] = citation_obj

            if categories := metadata.get('categories'):
                category_list = []
                for category in categories:
                    if label := category.get('label'):
                        category_list.append(label)
                output['keywords'] = category_list

            if publish_date := metadata.get('publish_date'):
                output['datePublished'] = publish_date

            if related_links := metadata.get('related_links'):
                citation_obj = {}
                for link in related_links:
                    if type := link.get('type'):
                        citation_obj['@type'] = type
                    if url := link.get('url'):
                        citation_obj['url'] = url
                output['citation'] = citation_obj

            if modified_on := metadata.get('modified_on'):
                output['dateModified'] = modified_on

            if links := metadata.get('links'):
                output['url'] = links['view']
            if repository := metadata.get('repository'):
                output['sdPublisher'] = {"name": repository['name']}

            yield output
