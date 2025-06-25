import logging
from datetime import date

import requests
from tenacity import retry, stop_after_attempt, wait_fixed

# Create a logger called "nde_crawlers"
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

PDB_API = "https://data.rcsb.org/rest/v1/core/entry"


# Define a function to print a message before retrying
def print_retry_message(retry_state):
    # Access the arguments passed to the retryable function
    args = retry_state.args
    if len(args) >= 2:  # Ensure both `url` and `params` are available
        url = args[0]
        payload = args[1]
        print(f"Retrying URL: {url}")
        print(f"Retrying Payload: {payload}")
    else:
        print("Retrying...")  # Fallback if arguments are not available


# Define a retryable function for making POST requests
@retry(stop=stop_after_attempt(5), wait=wait_fixed(2), after=print_retry_message)
def retryable_post(url, payload):
    response = requests.post(url, json=payload)
    if response.status_code == 999:
        logger.info("received end of PDB ids")
        return None
    return response.json()


@retry(stop=stop_after_attempt(5), wait=wait_fixed(2))
def retryable_get(url):
    return requests.get(url)


def paginate_through_PDB_ids(index=0):
    url = "https://www.rcsb.org/search/data"
    payload = {
        "report": "search_summary",
        "request": {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "group",
                        "logical_operator": "and",
                        "nodes": [
                            {
                                "type": "group",
                                "nodes": [
                                    {
                                        "type": "terminal",
                                        "service": "text",
                                        "parameters": {
                                            "attribute": "rcsb_entry_container_identifiers.entry_id",
                                            "operator": "exists",
                                            "negation": False,
                                        },
                                    },
                                    {
                                        "type": "terminal",
                                        "service": "text",
                                        "parameters": {
                                            "attribute": "rcsb_assembly_container_identifiers.assembly_id",
                                            "operator": "exists",
                                            "negation": False,
                                        },
                                    },
                                ],
                                "logical_operator": "and",
                            }
                        ],
                        "label": "text",
                    }
                ],
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": index, "rows": 100},
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "score", "direction": "desc"}],
                "scoring_strategy": "combined",
            },
            "request_info": {"query_id": "4da2b4cda592c491ae23704f087ac9ad"},
        },
        "getDrilldown": False,
        "attributes": None,
    }

    response = retryable_post(url, payload)
    if response is None:
        return None, None

    ids = [
        (i["identifier"], i["data"].get("organisms")) for i in response["result_set"]
    ]
    return ids, response["result_set_count"]


def get_PDB_ids():
    index = 0
    length = None
    while True:
        ids, length = paginate_through_PDB_ids(index)
        if not ids:
            logger.warning(f"breaking index {index} length {length}, no new ids found")
            break
        for _id, organisms in ids:
            yield _id, organisms
        index += 100
        logger.info(f"index: {index}")
        if length is not None and index >= length:
            break
        break


def parse():
    ids = get_PDB_ids()
    for i, (id, organisms) in enumerate(ids, start=1):
        doc = getPDBmetadata(id, organisms)
        if doc:
            yield doc

    logger.info(f"finished parsing {i} PDB metadata")


def getPDBmetadata(id, organisms):
    logger.info(f"Parsing url: {PDB_API}/{id} for id: {id}")
    try:
        resp = retryable_get(f"{PDB_API}/{id}")
    except Exception as e:
        logger.error(f"Error fetching metadata for ID {id}: {e}")
        return None
    if resp.status_code == 200:
        raw_data = resp.json()

        today = date.today().strftime("%Y-%m-%d")
        url = f"https://www.rcsb.org/structure/{raw_data['rcsb_id']}"

        md = {
            "@type": "Dataset",
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "archivedAt": url,
                "name": "Protein Data Bank",
                "url": "https://www.rcsb.org/",
                "versionDate": date.today().isoformat(),
            },
        }
        if organisms:
            species = []
            for organism in organisms:
                if organism:
                    species.append({"name": organism})
            if species:
                md["species"] = species

        md["name"] = raw_data["struct"]["title"]
        md["description"] = raw_data["struct"]["title"]
        md["_id"] = f"pdb_{raw_data['rcsb_id']}"
        md["identifier"] = raw_data["rcsb_id"]
        md["doi"] = f"10.2210/{md['_id']}/pdb"
        md["author"] = [
            {"@type": "Person", "name": author["name"]}
            for author in raw_data["audit_author"]
        ]
        if citations := raw_data.get("citation"):
            md["citedBy"] = [getCitation(citation) for citation in citations]

        if raw_data.get("exptl"):
            md["measurementTechnique"] = [
                {"name": technique["method"].lower()} for technique in raw_data["exptl"]
            ]
        if "pdbx_audit_support" in raw_data.keys():
            funding_list = []
            for funder in raw_data["pdbx_audit_support"]:
                funding = getFunding(funder)
                if funding is not None:
                    funding_list.append(funding)
            if funding_list:
                md["funding"] = funding_list

        md["datePublished"] = raw_data["rcsb_accession_info"]["deposit_date"][0:10]
        md["dateModified"] = raw_data["rcsb_accession_info"]["revision_date"][0:10]
        if raw_data.get("struct_keywords"):
            md["keywords"] = getKeywords(raw_data)
        md["url"] = url
        md["curatedBy"] = {
            "@type": "Organization",
            "url": md["url"],
            "name": "The Protein Data Bank",
            "versionDate": today,
        }
        if "rcsb_external_references" in raw_data.keys():
            md["sameAs"] = [
                link["link"] for link in raw_data["rcsb_external_references"]
            ]
        return md
    else:
        # logger.info(f"ID {id} returned an error from the API")
        logger.warning(f"ID {id} returned an error from the API")


def getCitation(citation):
    cite = {"@type": "Publication"}
    if citation.get("journal_abbrev"):
        cite["journalNameAbbrev"] = citation["journal_abbrev"]
    if title := citation.get("title"):
        cite["name"] = title

    if authors := citation.get("rcsb_authors"):
        cite["author"] = [{"@type": "Person", "name": author} for author in authors]
    if ("page_first" in citation.keys()) & ("page_last" in citation.keys()):
        cite["pagination"] = f"{citation['page_first']} - {citation['page_last']}"
    if "journal_volume" in citation.keys():
        cite["volumeNumber"] = citation["journal_volume"]
    if "year" in citation.keys():
        cite["datePublished"] = citation["year"]
    if "pdbx_database_id_doi" in citation.keys():
        cite["doi"] = citation["pdbx_database_id_doi"]
    if "pdbx_database_id_pub_med" in citation.keys():
        cite["pmid"] = citation["pdbx_database_id_pub_med"]
    return cite


def getFunding(funding):
    if name := funding.get("funding_organization"):
        funder = {"@type": "Organization", "name": name}
        obj = {"@type": "MonetaryGrant", "funder": funder}
        if "grant_number" in funding.keys():
            obj["identifier"] = funding["grant_number"]
        return obj
    else:
        return None


def getKeywords(result):
    keys = []
    if "pdbx_keywords" in result["struct_keywords"].keys():
        keys.extend(result["struct_keywords"]["pdbx_keywords"].split(","))
    if "text" in result["struct_keywords"].keys():
        keys.extend(result["struct_keywords"]["text"].split(","))

    keys = [key.strip() for key in keys]
    return list(set(keys))


if __name__ == "__main__":
    import json

    with open("d.json", "w") as d:
        docs = parse()
        for count, m in enumerate(docs, start=1):  # Start counting at 1
            json.dump(m, d)
            d.write("\n")
