import logging
import xml.etree.ElementTree as ET

import requests
from parse import parse
from tenacity import RetryCallState, retry, retry_if_exception, stop_after_attempt, wait_fixed

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")


def log_retry(retry_state: RetryCallState):
    """
    Log the URL being retried and the attempt number.
    """
    url = retry_state.args[0]  # The first argument to the fetch_url function
    logger.debug(f"Retrying URL: {url} (Attempt {retry_state.attempt_number})")


def should_retry(exception):
    """
    Custom retry condition: Retry only if the exception is not a 500 error.
    """
    if isinstance(exception, requests.exceptions.HTTPError):
        # Check if the status code is 500
        if exception.response.status_code == 500:
            logger.error(f"Not retrying for 500 error: {exception.response.url}")
            return False
    return True

@retry(
    stop=stop_after_attempt(3),
    wait=wait_fixed(2),
    reraise=True,
    before_sleep=log_retry,
    retry=retry_if_exception(should_retry),
)
def fetch_url(url):
    """
    Fetch the content of a URL with retry logic.
    Retries up to 3 times with a 2-second wait between attempts,
    but does not retry for 500 errors.
    """
    response = requests.get(url, timeout=30)
    response.raise_for_status()  # Raise an error if the request fails
    return response

def get_dataset_names():
    """
    Fetch the dataset names from the XML file and return them as a list.
    """
    # URL of the XML file
    url = "https://www.ebi.ac.uk/ebisearch/ws/rest/omics"

    # Fetch the XML content with retry logic
    response = fetch_url(url)

    # Parse the XML content
    root = ET.fromstring(response.content)

    # Extract all 'id' attributes within 'domain' tags
    domain_ids = [domain.get("id") for domain in root.findall(".//domain") if domain.get("id")]

    # Add the prefixes to exclude_ids
    exclude_ids = [
        "massive",
        "geo",
        "lincs",
        "biostudies-literature",
        "project",
        "biostudies-arrayexpress",
        "dbgap",
        "paxdb",
        "iprox",
    ]

    # Filter out excluded IDs
    dataset_names = [domain_id for domain_id in domain_ids if domain_id not in exclude_ids]

    return dataset_names


def get_dataset_ids(dataset_name):
    start = 0
    hits = 1
    while start < hits:
        size = 1000
        response = fetch_url(
            f"https://www.ebi.ac.uk/ebisearch/ws/rest/{dataset_name}?query=*:*&size={size}&start={start}"
        )

        # Parse the XML content
        root = ET.fromstring(response.content)

        if start == 0:
            # Extract the total number of hits from the response
            hits = int(root.find(".//hitCount").text)
            logger.info(f"Total hits for dataset {dataset_name}: {hits}")

        # Extract all 'id' attributes within 'entry' tags
        entry_ids = [entry.get("id") for entry in root.findall(".//entry") if entry.get("id")]
        if not entry_ids:
            break
        else:
            for entry_id in entry_ids:
                if not entry_id.casefold().startswith("s-epmc"):
                    yield entry_id
            start += len(entry_ids)
            logger.info(f"Processed {start} entries for dataset: {dataset_name}")


def process_dataset_records(dataset_name, queue):
    dataset_ids = get_dataset_ids(dataset_name)
    base_url = f"https://www.omicsdi.org/ws/dataset/{dataset_name}/"

    logger.info(f"Starting processing for dataset: {dataset_name}")
    for count, dataset_id in enumerate(dataset_ids, start=1):
        # Construct the URL for each entry
        url = f"{base_url}{dataset_id}"
        try:
            response = fetch_url(url).json()
            data = parse(response, dataset_name, dataset_id, url)
            if data:
                queue.put(data)
        except Exception as e:
            logger.error(f"Error fetching URL for entry {url}: {e}")
            continue
    logger.info(f"Finished processing for dataset: {dataset_name}")
