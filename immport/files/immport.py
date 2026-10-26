import logging
import time

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")

CATALOG_NAME = "Immunology Database and Analysis Portal (ImmPort)"

REQUEST_TIMEOUT = 30
REQUEST_ATTEMPTS = 5
TRANSIENT_ERRORS = (
    requests.exceptions.Timeout,
    requests.exceptions.ConnectionError,
    requests.exceptions.ChunkedEncodingError,
)


def get_json(url, **kwargs):
    """GET `url` as JSON, retrying TRANSIENT_ERRORS. Re-raises the last one after REQUEST_ATTEMPTS tries."""
    for attempt in range(1, REQUEST_ATTEMPTS + 1):
        try:
            return requests.get(url, timeout=REQUEST_TIMEOUT, **kwargs).json()
        except TRANSIENT_ERRORS as e:
            logger.info(f"{type(e).__name__} fetching {url} (attempt {attempt}/{REQUEST_ATTEMPTS})")
            if attempt == REQUEST_ATTEMPTS:
                raise
            time.sleep(5 * attempt)


def get_studies():
    """Return {study accession: study search hit} for every ImmPort study.

    The search hit reports fields the JSON-LD export leaves out, such as
    actual_enrollment. A search page that keeps failing fails the crawl
    instead of leaving its studies out.
    """
    base_url = "https://www.immport.org/shared/data/query/api/search/study?term="
    page_size = 500
    from_record = 0
    studies = {}
    while True:
        immport_search_payload = {"pageSize": page_size, "fromRecord": from_record, "sortFieldDirection": "desc"}
        j = get_json(base_url, params=immport_search_payload)
        studies.update((h["_id"], h.get("_source", {})) for h in j["hits"]["hits"])

        from_record += page_size
        time.sleep(1)
        logger.info(f"Retrieved ids: {len(studies)}")
        if len(j["hits"]["hits"]) < page_size:
            break
    return studies


def _type_terms(schema_json, field):
    """Type an ontology term list ImmPort exports untyped, and its `curatedBy`."""
    terms = schema_json.get(field)
    if not terms:
        return
    for term in terms if isinstance(terms, list) else [terms]:
        if not isinstance(term, dict):
            continue
        term.setdefault("@type", "DefinedTerm")
        # "ImmPort Curation" is a curation effort by the repository, not a tool.
        if isinstance(term.get("curatedBy"), dict):
            term["curatedBy"].setdefault("@type", "Organization")


def map_schema_json(schema_json):
    """Map the ImmPort JSON-LD export to the NDE schema.

    The JSON-LD export is already largely schema.org compatible.
    This function handles the remaining transformations:
    - Remove JSON-LD context fields
    - includedInDataCatalog: use the curated ImmPort catalog name
    - creator → author (with affiliation wrapping)
    - citation → citedBy
    - species/measurementTechnique: wrap strings in {"name": ...}
    """
    # Remove JSON-LD context fields not needed downstream
    schema_json.pop("@context", None)
    schema_json.pop("@id", None)

    # The export labels ImmPort's own catalog with its tagline ("ImmPort -
    # Bioinformatics For the Future of Immunology"); NDE searches by the curated name.
    catalogs = schema_json.get("includedInDataCatalog")
    for catalog in catalogs if isinstance(catalogs, list) else [catalogs]:
        if isinstance(catalog, dict) and "immport" in str(catalog.get("name", "")).casefold():
            catalog["name"] = CATALOG_NAME

    # creator → author (wrap affiliation string in {"name": ...})
    if author := schema_json.pop("creator", None):
        for data in author:
            if isinstance(data.get("affiliation"), str):
                data["affiliation"] = {"@type": "Organization", "name": data["affiliation"]}
            # The upstream JSON-LD does not always type its creators.
            data.setdefault("@type", "Person")
        schema_json["author"] = author

    # citation → citedBy
    if cited_by := schema_json.pop("citation", None):
        schema_json["citedBy"] = cited_by

    # species: list of strings → list of {"name": ...}
    if species := schema_json.pop("species", None):
        if isinstance(species, list):
            species_list = [{"@type": "DefinedTerm", "name": s} for s in species]
            if species_list:
                schema_json["species"] = species_list
        else:
            schema_json["species"] = {"@type": "DefinedTerm", "name": species}

    # measurementTechnique: list of strings → list of {"name": ...}
    if measurement_techniques := schema_json.pop("measurementTechnique", None):
        mt_list = [{"@type": "DefinedTerm", "name": mt} for mt in measurement_techniques]
        if mt_list:
            schema_json["measurementTechnique"] = mt_list

    # healthCondition comes through from the JSON-LD export untouched, and
    # ImmPort types neither the term nor the curation record on it.
    _type_terms(schema_json, "healthCondition")

    # identifier: list of strings/dicts → list of strings (extract "value" if dict)
    if identifiers := schema_json.pop("identifier", None):
        id_list = []
        for identifier in identifiers:
            if isinstance(identifier, dict) and "value" in identifier:
                id_list.append(identifier["value"])
            elif isinstance(identifier, str):
                id_list.append(identifier)
        if id_list:
            schema_json["identifier"] = id_list

    schema_json["conditionsOfAccess"] = "Closed"

    return schema_json


def map_study_search(study):
    """Map the study search fields the JSON-LD export leaves out."""
    search_dict = {}
    # 0 marks unreported enrollment, as on the CITC consortium study (SDY1432).
    if actual_enrollment := study.get("actual_enrollment"):
        search_dict["sample"] = {
            "@type": "Sample",
            "sampleQuantity": {
                "@type": "QuantitativeValue",
                "value": actual_enrollment,
                "unitText": "enrolled subjects",
                "unitCode": "http://purl.obolibrary.org/obo/NCIT_C207572",
            },
        }
    return search_dict


def map_protocol_json(protocol_json, immport_id):
    protocol_dict = {}
    if protocols := protocol_json.get("protocols"):
        is_based_on_list = []
        for protocol in protocols:
            is_based_on_dict = {"@type": "CreativeWork"}
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
    studies = get_studies()

    count = 0

    for immport_id, study in studies.items():
        count += 1
        logger.info(f"Processing ImmPort ID: {immport_id} ({count}/{len(studies)})")

        schema_url = f"https://s3.immport.org/release/metadata/Study_{immport_id}"
        try:
            schema_json = get_json(schema_url)
            schema_json = map_schema_json(schema_json)
        except TRANSIENT_ERRORS as e:
            logger.info(f"Skipping {immport_id}, could not fetch its JSON-LD: {e}")
            continue

        schema_json.update(map_study_search(study))

        protocol_url = f"https://www.immport.org/shared/data/query/ui/study/design/{immport_id}"
        try:
            protocol_json = get_json(protocol_url)
            protocol_json = map_protocol_json(protocol_json, immport_id)
        except TRANSIENT_ERRORS as e:
            logger.info(f"Could not fetch ImmPort protocols for {immport_id}: {e}")
            protocol_json = {}

        if protocol_json:
            schema_json.update(protocol_json)

        yield schema_json
