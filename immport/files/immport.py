import logging
import time

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


def get_ids():
    base_url = "https://www.immport.org/shared/data/query/api/search/study?term="
    page_size = 500
    from_record = 0
    ids = []
    while True:
        immport_search_payload = {"pageSize": page_size, "fromRecord": from_record, "sortFieldDirection": "desc"}
        try:
            j = requests.get(
                base_url,
                timeout=5,
                params=immport_search_payload,
            ).json()
            ids += [h["_id"] for h in j["hits"]["hits"]]
        except requests.exceptions.Timeout:
            logger.info("Timeout searching ImmPort IDs")

        from_record += page_size
        time.sleep(1)
        logger.info(f"Retrieved ids: {len(ids)}")
        if len(j["hits"]["hits"]) < page_size:
            break
    return ids


def map_schema_json(schema_json):
    """Map the ImmPort JSON-LD export to the NDE schema.

    The JSON-LD export is already largely schema.org compatible.
    This function handles the remaining transformations:
    - Remove JSON-LD context fields
    - creator → author (with affiliation wrapping)
    - citation → citedBy
    - species/measurementTechnique: wrap strings in {"name": ...}
    """
    # Remove JSON-LD context fields not needed downstream
    schema_json.pop("@context", None)
    schema_json.pop("@id", None)

    # creator → author (wrap affiliation string in {"name": ...})
    if author := schema_json.pop("creator", None):
        for data in author:
            if isinstance(data.get("affiliation"), str):
                data["affiliation"] = {"name": data["affiliation"]}
        schema_json["author"] = author

    # citation → citedBy
    if cited_by := schema_json.pop("citation", None):
        schema_json["citedBy"] = cited_by

    # species: list of strings → list of {"name": ...}
    if species := schema_json.pop("species", None):
        if isinstance(species, list):
            species_list = [{"name": s} for s in species]
            if species_list:
                schema_json["species"] = species_list
        else:
            schema_json["species"] = {"name": species}

    # measurementTechnique: list of strings → list of {"name": ...}
    if measurement_techniques := schema_json.pop("measurementTechnique", None):
        mt_list = [{"name": mt} for mt in measurement_techniques]
        if mt_list:
            schema_json["measurementTechnique"] = mt_list

    return schema_json


def map_protocol_json(protocol_json, immport_id):
    protocol_dict = {}
    if protocols := protocol_json.get("protocols"):
        is_based_on_list = []
        for protocol in protocols:
            is_based_on_dict = {}
            if description := protocol.get("description"):
                is_based_on_dict["description"] = description
            if file_name := protocol.get("fileName"):
                is_based_on_dict["name"] = file_name
            if name := protocol.get("name"):
                is_based_on_dict["name"] = name
            if original_file_name := protocol.get("originalFileName"):
                is_based_on_dict["name"] = original_file_name
            if protocol_accession := protocol.get("protocolAccession"):
                is_based_on_dict["identifier"] = protocol_accession
            if protocol_type := protocol.get("type"):
                is_based_on_dict["additionalType"] = {"name": protocol_type}
            if bool(is_based_on_dict):
                is_based_on_dict["url"] = f"https://www.immport.org/browser/?path={immport_id}"
                is_based_on_list.append(is_based_on_dict)
        if len(is_based_on_list):
            protocol_dict["isBasedOn"] = is_based_on_list
    return protocol_dict


def parse():
    ids = get_ids()

    count = 0

    for immport_id in ids:
        count += 1
        logger.info(f"Processing ImmPort ID: {immport_id} ({count}/{len(ids)})")

        schema_url = f"https://s3.immport.org/release/metadata/Study_{immport_id}"
        try:
            schema_json = requests.get(schema_url, timeout=5).json()
            schema_json = map_schema_json(schema_json)
        except requests.exceptions.Timeout:
            logger.info(f"Timeout fetching ImmPort JSON-LD: {immport_id}")
            continue

        protocol_url = f"https://www.immport.org/shared/data/query/ui/study/design/{immport_id}"
        try:
            protocol_json = requests.get(protocol_url, timeout=5).json()
            protocol_json = map_protocol_json(protocol_json, immport_id)
        except requests.exceptions.Timeout:
            logger.info(f"Timeout fetching ImmPort protocol: {immport_id}")
            protocol_json = {}

        if protocol_json:
            schema_json.update(protocol_json)

        yield schema_json
