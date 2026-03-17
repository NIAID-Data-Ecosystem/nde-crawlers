import datetime
import logging

import dateutil.parser
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def insert_value(d, key, value):
    if key in d:
        if isinstance(d[key], list) and value not in d[key]:
            d[key].append(value)
        if not isinstance(d[key], list) and d[key] != value:
            d[key] = [d[key], value]
    else:
        d[key] = value


def make_session_with_retries(total=3, backoff_factor=1, status_forcelist=(429,500,502,503,504), allow_post=False, pool_maxsize=10):
    allowed = frozenset(["HEAD", "GET", "OPTIONS"])
    if allow_post:
        allowed = allowed.union(["POST"])
    retry = Retry(
        total=total,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
        allowed_methods=allowed,
        raise_on_status=False
    )
    adapter = HTTPAdapter(pool_connections=10, pool_maxsize=pool_maxsize, max_retries=retry)
    session = requests.Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        "Accept": "application/json, text/plain, */*",
        "Connection": "keep-alive",
    })
    return session


def get_ids(session):
    url = "https://www.biosino.org/node/api/app/browse/search"
    page_num = 1
    payload = {
      "queryWord": "",
      "advQueryWord": "",
      "sortKey": "modifiedDate",
      "sortType": "desc",
      "pageNum": page_num,
      "pageSize": 100,
      "leftStatQueries": [
        {"type": "Level", "name": "project", "treeRefId": "0_lsTree__Level"}
      ]
    }

    response = session.post(url, json=payload)
    response.raise_for_status()  # will raise an HTTPError if the HTTP request returned an unsuccessful status code
    response = response.json()
    pages = response["data"]["pageInfo"]["totalPages"]
    logger.info(f"Total pages to crawl: {pages + 1}")
    for page in range(1, pages + 1):
        logger.info(f"Crawling page {page} of {pages}")
        payload["pageNum"] = page
        response = session.post(url, json=payload)
        response.raise_for_status()
        response = response.json()
        for project in response["data"]["pageInfo"]["content"]:
            yield project["id"]


def parse():
    with make_session_with_retries() as session:
        for project_id in get_ids(session):
            general_info = session.get(f"https://www.biosino.org/node/api/app/project/getGeneralInfo/{project_id}").json()
            general_info = general_info["data"]
            author_info = session.get(f"https://www.biosino.org/node/api/app/project/getAuthorInfo/{project_id}").json()
            author_info = author_info["data"]
            _id = general_info["projectNo"]
            url = f"https://www.biosino.org/node/project/detail/{_id}"
            output = {
                "@context": "http://schema.org/",
                "@type": "Dataset",
                "_id": "NODE_" + _id.casefold(),
                "identifier": _id,
                "url": url,
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "National Omics Data Encyclopedia",
                    "url": "https://www.biosino.org/node/home",
                    "versionDate": datetime.date.today().isoformat(),
                    "archivedAt": url,
                },
            }

            if name := general_info.get("name"):
                insert_value(output, "name", name)

            if description := general_info.get("description"):
                insert_value(output, "description", description)

            if publishes := general_info.get("publishes"):
                if not isinstance(publishes, list):
                    publishes = [publishes]
                for publish in publishes:
                    if pmid := publish.get("pmid"):
                        if "pmids" not in output:
                            output["pmids"] = pmid
                        else:
                            output["pmids"] += f", {pmid}"
                    # citation = {}
                    # if doi := publish.get("doi"):
                    #     insert_value(citation, "doi", doi)
                    # if pmid := publish.get("pmid"):
                    #     insert_value(citation, "pmid", pmid)
                    # if citation:
                    #     insert_value(output, "citation", citation)

            author = {}
            if given_name := author_info.get("firstName"):
                insert_value(author, "givenName", given_name)
            if family_name := author_info.get("lastName"):
                insert_value(author, "familyName", family_name)
            if name := author_info.get("orgName"):
                insert_value(author, "affiliation", {"name": name})
            if author:
                insert_value(output, "author", author)

            id_stat_info = session.get(f"https://www.biosino.org/node/api/app/project/getStatDetail/{project_id}").json()
            for exp_stat in id_stat_info["data"]["expStats"]:
                data_type = exp_stat["name"]
                exp_info = session.get(f"https://www.biosino.org/node/api/app/project/getExpAndSampleList?projectNo={project_id}&type=experiment&dataType={data_type}&total=0&pageNum=1&pageSize=100&sortKey=expNo&sortType=asc").json()
                pages = exp_info["data"]["expTableData"]["totalPages"]
                for page in range(1, pages + 1):
                    exp_info = session.get(f"https://www.biosino.org/node/api/app/project/getExpAndSampleList?projectNo={project_id}&type=experiment&dataType={data_type}&total=0&pageNum={page}&pageSize=100&sortKey=expNo&sortType=asc").json()
                    for exp in exp_info["data"]["expTableData"]["content"]:
                        if date_created := exp.get("createDate"):
                            try:
                                date_created = dateutil.parser.parse(date_created).date().isoformat()
                                insert_value(output, "dateCreated", date_created)
                            except (ValueError, TypeError):
                                logger.info(f"Failed to parse date: {date_created} for project {project_id}")

                        if date_modified := exp.get("updateDate"):
                            try:
                                date_modified = dateutil.parser.parse(date_modified).date().isoformat()
                                insert_value(output, "dateModified", date_modified)
                            except (ValueError, TypeError):
                                logger.info(f"Failed to parse date: {date_modified} for project {project_id}")

                        if submitter := exp.get("submitter"):
                            author = {}
                            if given_name := submitter.get("firstName"):
                                insert_value(author, "givenName", given_name)
                            if family_name := submitter.get("lastName"):
                                insert_value(author, "familyName", family_name)
                            if name := submitter.get("orgName"):
                                insert_value(author, "affiliation", {"name": name})
                            if author:
                                insert_value(output, "author", author)

                        if attributes := exp.get("attributes"):
                            if name := attributes.get("library_selection"):
                                insert_value(output, "measurementTechnique", {"name": name})
                            if name := attributes.get("library_strategy"):
                                insert_value(output, "measurementTechnique", {"name": name})
                            if name := attributes.get("platform"):
                                insert_value(output, "measurementTechnique", {"name": name})

            if date_created := output.get("dateCreated"):
                if not isinstance(date_created, list):
                    date_created = [date_created]
                date_created.sort()
                date_created = datetime.datetime.fromisoformat(date_created[0]).date().isoformat()
                output["dateCreated"] = date_created

            if date_modified := output.get("dateModified"):
                if not isinstance(date_modified, list):
                    date_modified = [date_modified]
                date_modified.sort()
                date_modified = datetime.datetime.fromisoformat(date_modified[-1]).date().isoformat()
                output["dateModified"] = date_modified

            yield output
