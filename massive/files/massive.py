import json
import logging
import math
import os
import time
from datetime import datetime

from bs4 import BeautifulSoup
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

BASE_URL = "https://massive.ucsd.edu/ProteoSAFe"
QUERY_DATASETS_URL = f"{BASE_URL}/QueryDatasets"
MASSIVE_CATALOG = {
    "name": "MassIVE",
    "url": "https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp",
    "@type": "DataCatalog",
}

PAGE_SIZE = int(os.getenv("MASSIVE_PAGE_SIZE", "1000"))
REQUEST_TIMEOUT = int(os.getenv("MASSIVE_REQUEST_TIMEOUT", "120"))
MIN_EXPECTED_RECORDS = int(os.getenv("MASSIVE_MIN_EXPECTED_RECORDS", "1000"))
MIN_COMPLETION_RATIO = float(os.getenv("MASSIVE_MIN_COMPLETION_RATIO", "0.95"))
FETCH_HTML_DETAILS = os.getenv("MASSIVE_FETCH_HTML_DETAILS", "true").casefold() in {"1", "true", "yes", "y"}
HTML_DETAILS_CACHE_PATH = os.getenv("MASSIVE_HTML_DETAILS_CACHE", "/cache/massive_html_details.json")
HTML_DETAIL_SLEEP = float(os.getenv("MASSIVE_HTML_DETAIL_SLEEP", "0.2"))
HTML_CACHE_SAVE_INTERVAL = int(os.getenv("MASSIVE_HTML_CACHE_SAVE_INTERVAL", "100"))


def make_session():
    retry = Retry(
        total=6,
        connect=6,
        read=6,
        status=6,
        backoff_factor=2,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET"]),
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=4, pool_maxsize=4)
    session = requests.Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(
        {
            "Accept": "application/json, text/plain, */*",
            "User-Agent": "NDE-MassIVE-crawler/1.0",
        }
    )
    return session


def _request_json(session, url, params):
    response = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    return response.json()


def _get_query_value():
    # The leading # is part of the MassIVE UI query format. requests will URL-encode it once.
    return '#{"query":{},"table_sort_history":"createdMillis_dsc"}'


def fetch_dataset_pages(session=None, page_size=PAGE_SIZE):
    close_session = session is None
    session = session or make_session()
    query_params = {
        "pageSize": page_size,
        "offset": 0,
        "query": _get_query_value(),
    }
    page = 1
    total_rows = None

    try:
        while True:
            data = _request_json(session, QUERY_DATASETS_URL, query_params)
            row_data = data.get("row_data", [])
            if total_rows is None and data.get("total_rows") is not None:
                try:
                    total_rows = int(data["total_rows"])
                except (TypeError, ValueError):
                    logger.warning("Unexpected MassIVE total_rows value: %s", data["total_rows"])

            logger.info(
                "Processed MassIVE page %s: offset=%s rows=%s total_rows=%s",
                page,
                query_params["offset"],
                len(row_data),
                total_rows,
            )
            if not row_data:
                break

            yield row_data, total_rows

            query_params["offset"] += len(row_data)
            page += 1
            if total_rows is not None and query_params["offset"] >= total_rows:
                break
    finally:
        if close_session:
            session.close()


def fetch_dataset_ids():
    dataset_ids = []
    with make_session() as session:
        for row_data, _total_rows in fetch_dataset_pages(session=session):
            for item in row_data:
                if dataset_id := item.get("dataset"):
                    dataset_ids.append(dataset_id)

    return dataset_ids


def fetch_dataset(dataset_id, session=None):
    close_session = session is None
    session = session or make_session()
    params = {
        "pageSize": 30,
        "offset": 0,
        "query": json.dumps({"title_input": dataset_id}, separators=(",", ":")),
    }
    try:
        return _request_json(session, QUERY_DATASETS_URL, params=params)
    finally:
        if close_session:
            session.close()


def fetch_html_content(task, session=None):
    close_session = session is None
    url = f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task={task}"
    session = session or make_session()
    try:
        response = session.get(url, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        return response.text
    finally:
        if close_session:
            session.close()


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

    return doi, license_url


def load_html_details_cache(cache_path=HTML_DETAILS_CACHE_PATH):
    if not cache_path or not os.path.exists(cache_path):
        return {}

    try:
        with open(cache_path, "r") as cache_file:
            cache = json.load(cache_file)
    except Exception as e:
        logger.warning("Could not load MassIVE HTML details cache %s: %s", cache_path, e)
        return {}

    if not isinstance(cache, dict):
        logger.warning("Ignoring unexpected MassIVE HTML details cache format in %s", cache_path)
        return {}
    return cache


def save_html_details_cache(cache, cache_path=HTML_DETAILS_CACHE_PATH):
    if not cache_path:
        return

    try:
        cache_dir = os.path.dirname(cache_path)
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
        tmp_path = f"{cache_path}.{os.getpid()}.tmp"
        with open(tmp_path, "w") as cache_file:
            json.dump(cache, cache_file, sort_keys=True)
        os.replace(tmp_path, cache_path)
    except Exception as e:
        logger.warning("Could not save MassIVE HTML details cache %s: %s", cache_path, e)


def get_html_details(task, session=None, cache=None):
    if cache is not None and task in cache:
        detail = cache[task] or {}
        return detail.get("doi"), detail.get("license")

    html_content = fetch_html_content(task, session=session)
    doi, license_url = parse_html_for_doi_and_license(html_content)
    if cache is not None:
        cache[task] = {"doi": doi, "license": license_url}
    if HTML_DETAIL_SLEEP > 0:
        time.sleep(HTML_DETAIL_SLEEP)
    return doi, license_url


def parse_dataset_row(item, session=None, fetch_html_details=FETCH_HTML_DETAILS, html_details_cache=None):
    identifier = item.get("dataset")
    if not identifier:
        return None

    catalog = {**MASSIVE_CATALOG, "versionDate": datetime.today().strftime("%Y-%m-%d")}
    output = {
        "includedInDataCatalog": catalog,
        "@type": "Dataset",
        "identifier": identifier,
        "_id": identifier.lower(),
    }

    if task := item.get("task"):
        output["url"] = f"{BASE_URL}/dataset.jsp?task={task}"
        output["includedInDataCatalog"]["archivedAt"] = output["url"]

        if fetch_html_details:
            doi, license_url = get_html_details(task, session=session, cache=html_details_cache)
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

    if doi := item.get("doi"):
        output["doi"] = doi

    if license_url := item.get("license"):
        output["license"] = license_url

    if sdPublisher := item.get("site"):
        output["sdPublisher"] = {"name": sdPublisher}

    if description := item.get("description"):
        if description != "null":
            output["description"] = description

    if keywords := item.get("keywords"):
        if keywords != "null":
            if "###" in keywords:
                output["keywords"] = [keyword for keyword in keywords.split("###") if keyword]
            else:
                output["keywords"] = [keyword for keyword in keywords.split(", ") if keyword]

    if date_created := item.get("create_time"):
        try:
            parsed_date = datetime.strptime(date_created, "%Y-%m-%d %H:%M:%S.%f").strftime("%Y-%m-%d")
            output["dateCreated"] = parsed_date
        except ValueError:
            logger.warning("Could not parse create_time for %s: %s", identifier, date_created)

    if measurement_techniques := item.get("instrument_resolved"):
        mt_list = [mt for mt in measurement_techniques.split("###") if mt and mt != "null"]
        output["measurementTechnique"] = [{"name": mt} for mt in mt_list]

    if species := item.get("species_resolved"):
        species_list = [sp for sp in species.split("###") if sp and sp != "null"]
        output["species"] = [{"name": sp.split(" (NCBITaxon:")[0]} for sp in species_list]

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
            if pi_data:
                pis_list.append(pi_data)
        if pis_list:
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
            if pub_data:
                pub_list.append(pub_data)
        if pub_list:
            output["citation"] = pub_list

    if privacy := item.get("privacy"):
        if privacy == "Public":
            output["conditionsOfAccess"] = "Open"
        else:
            output["conditionsOfAccess"] = "Closed"

    default_measurement_technique = {
        "name": "Mass Spectrometry",
        "url": "https://ontobee.org/ontology/MMO?iri=http://purl.obolibrary.org/obo/MMO_0000534",
        "identifier": "MMO_0000534",
    }
    output.setdefault("measurementTechnique", [])
    output["measurementTechnique"].append(default_measurement_technique)

    return output


def parse_dataset(json_data, session=None, fetch_html_details=FETCH_HTML_DETAILS, html_details_cache=None):
    row_data = json_data.get("row_data", [])
    if not row_data:
        return {}
    return parse_dataset_row(
        row_data[0],
        session=session,
        fetch_html_details=fetch_html_details,
        html_details_cache=html_details_cache,
    )


def parse(
    fetch_html_details=FETCH_HTML_DETAILS,
    min_expected_records=MIN_EXPECTED_RECORDS,
    min_completion_ratio=MIN_COMPLETION_RATIO,
):
    count = 0
    skipped = 0
    total_rows = None
    seen_ids = set()
    html_details_cache = load_html_details_cache() if fetch_html_details else None
    saved_cache_size = len(html_details_cache) if html_details_cache is not None else 0

    try:
        with make_session() as session:
            for row_data, page_total_rows in fetch_dataset_pages(session=session):
                if page_total_rows is not None:
                    total_rows = page_total_rows
                for item in row_data:
                    parsed_dataset = parse_dataset_row(
                        item,
                        session=session,
                        fetch_html_details=fetch_html_details,
                        html_details_cache=html_details_cache,
                    )
                    if not parsed_dataset:
                        skipped += 1
                        continue

                    doc_id = parsed_dataset["_id"]
                    if doc_id in seen_ids:
                        logger.warning("Skipping duplicate MassIVE id %s", doc_id)
                        skipped += 1
                        continue

                    seen_ids.add(doc_id)
                    count += 1
                    if count % 1000 == 0:
                        logger.info("Parsed %s MassIVE datasets", count)

                    if (
                        html_details_cache is not None
                        and len(html_details_cache) - saved_cache_size >= HTML_CACHE_SAVE_INTERVAL
                    ):
                        save_html_details_cache(html_details_cache)
                        saved_cache_size = len(html_details_cache)

                    yield parsed_dataset
    finally:
        if html_details_cache is not None and len(html_details_cache) != saved_cache_size:
            save_html_details_cache(html_details_cache)

    expected_floor = min_expected_records
    if total_rows is not None:
        expected_floor = max(expected_floor, math.ceil(total_rows * min_completion_ratio))

    if count < expected_floor:
        raise RuntimeError(
            f"MassIVE crawl produced {count} records, below expected floor {expected_floor} "
            f"(total_rows={total_rows}, skipped={skipped}). Refusing to publish a partial dump."
        )

    logger.info("Finished MassIVE crawl: parsed=%s skipped=%s total_rows=%s", count, skipped, total_rows)
