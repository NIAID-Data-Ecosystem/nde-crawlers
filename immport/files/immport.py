import datetime
import logging
import time

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


def _set_single_element_or_all_in_list(obj, key, value):
    """mutates obj, list of dict or dict"""
    if isinstance(obj, list):
        for element in obj:
            element[key] = value
    else:
        obj[key] = value
    return obj


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
            # catch the rather harmless exception, and do nothing
            logger.info("Timeout searching ImmPort IDs")

        from_record += page_size
        time.sleep(1)
        logger.info(f"Retrieved ids: {len(ids)}")
        if len(j["hits"]["hits"]) < page_size:
            break
    return ids


def map_schema_json(schema_json):
    if _id := schema_json.get("url", None):
        schema_json["_id"] = _id.rsplit("/", 1)[-1]
    if author := schema_json.pop("creator", None):
        for data in author:
            for key, value in data.items():
                data[key] = {"name": value} if key == "affiliation" else value
        schema_json["author"] = author
    if cited_by := schema_json.pop("citations", None):
        schema_json["citedBy"] = cited_by
    if identifier := schema_json.pop("identifiers", None):
        schema_json["identifier"] = identifier
    if species := schema_json.pop("species", None):
        if isinstance(species, list):
            species_list = []
            for single_species in species:
                species_list.append({"name": single_species})
            if len(species_list):
                schema_json["species"] = species_list
        else:
            schema_json["species"] = {"name": species}
    if measurement_techniques := schema_json.pop("measurementTechnique", None):
        mt_list = []
        for measurement_technique in measurement_techniques:
            mt_list.append({"name": measurement_technique})
        if mt_list:
            schema_json["measurementTechnique"] = mt_list
    if health_conditions := schema_json.pop("keywords", None):
        no_dups = set()
        hc = []
        for health_condition in health_conditions:
            lower = health_condition.lower()
            if lower not in no_dups:
                no_dups.add(lower)
                hc.append({"name": lower})
        schema_json["healthCondition"] = hc

    if date := schema_json.pop("curationDate", None):
        date = datetime.datetime.strptime(date, "%m/%d/%Y")  # mm/dd/YYYY is my guess
        if distribution := schema_json.pop("distribution", None):
            schema_json["distribution"] = _set_single_element_or_all_in_list(distribution, "dateModified", date)
            # if curated_by := schema_json.pop("curatedBy", None):
            # remove curatedBy
            # schema_json['curatedBy'] = _set_single_element_or_all_in_list(
            #     curated_by, 'curationDate', date
            # )
            # pass
        if schema_json.get("includedInDataCatalog"):
            schema_json["includedInDataCatalog"]["versionDate"] = date
            schema_json["includedInDataCatalog"]["name"] = "ImmPort"
        if "date" not in schema_json:
            schema_json["date"] = date

        # new mapping for api
        if description := schema_json.get("description"):
            schema_json["abstract"] = description

    # temp fix for outbreak.info
    if cb_outbreak := schema_json.get("includedInDataCatalog"):
        schema_json["curatedBy"] = cb_outbreak

    if schema_json["_id"] in ["IMMPORT_SDY1760", "IMMPORT_SDY2112"]:
        schema_json["conditionsOfAccess"] = "Restricted"
        schema_json["isAccessibleForFree"] = True
    else:
        schema_json["conditionsOfAccess"] = "Closed"
        schema_json["isAccessibleForFree"] = True

    # set includedInDataCatalog.dataset
    if schema_json.get("includedInDataCatalog"):
        schema_json["includedInDataCatalog"]["dataset"] = schema_json["url"]
        schema_json["includedInDataCatalog"]["@type"] = "DataCatalog"
    return schema_json


def map_details_json(details_json, existing_health_conditions):
    details_dict = {}
    if detailed_description := details_json.get("detailedDescription"):
        details_dict["description"] = detailed_description
    if condition_studied := details_json.get("conditionStudied"):
        health_conditions = []
        for condition in condition_studied.split(","):
            health_conditions.append({"name": condition.strip().lower()})
        if len(health_conditions):
            if existing_health_conditions:
                health_conditions.extend(existing_health_conditions)
            health_conditions = list({v["name"]: v for v in health_conditions}.values())
            details_dict["healthCondition"] = health_conditions
    keywords = []
    if program := details_json.get("program"):
        keywords.append(program)
    if endpoints := details_json.get("endpoints"):
        keywords.append(endpoints)
    if keywords:
        details_dict["keywords"] = keywords
    return details_dict


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

        schema_url = f"https://s3.immport.org/release/metadata/Study_{immport_id}.compact"
        try:
            schema_json = requests.get(schema_url, timeout=5).json()
            schema_json = map_schema_json(schema_json)
            existing_health_conditions = schema_json.get("healthCondition", None)

        except requests.exceptions.Timeout:
            # catch the rather harmless exception, and do nothing
            logger.info(f"Timeout searching ImmPort ID schema_url: {immport_id}")

        details_url = f"https://www.immport.org/shared/data/query/ui/study/summary/{immport_id}"
        try:
            details_json = requests.get(details_url, timeout=5).json()
            details_json = map_details_json(details_json, existing_health_conditions)
        except requests.exceptions.Timeout:
            # catch the rather harmless exception, and do nothing
            logger.info(f"Timeout searching ImmPort ID details_url: {immport_id}")

        protocol_url = f"https://www.immport.org/shared/data/query/ui/study/design/{immport_id}"
        try:
            protocol_json = requests.get(protocol_url, timeout=5).json()
            protocol_json = map_protocol_json(protocol_json, immport_id)
        except requests.exceptions.Timeout:
            # catch the rather harmless exception, and do nothing
            logger.info(f"Timeout searching ImmPort ID protocol_url: {immport_id}")

        if details_json:
            schema_json.update(details_json)
        if protocol_json:
            schema_json.update(protocol_json)

        yield schema_json
