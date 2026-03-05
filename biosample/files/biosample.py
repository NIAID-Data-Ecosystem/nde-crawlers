import datetime
import json
import logging
import re
import socket
from pathlib import Path

import dateutil.parser
import xmltodict
from Bio import Entrez
from tenacity import retry, stop_after_attempt, wait_fixed

DEFAULT_TIMEOUT = 30  # seconds
socket.setdefaulttimeout(DEFAULT_TIMEOUT)
GEO_API_KEY = "3048f6bdb7c91cc8ad7af802559ec470e609"
GEO_EMAIL = "cwu@scripps.edu"

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


Entrez.email = GEO_EMAIL
Entrez.api_key = GEO_API_KEY


def insert_value(d, key, value):
    if key in d:
        if isinstance(d[key], list) and value not in d[key]:
            d[key].append(value)
        if not isinstance(d[key], list) and d[key] != value:
            d[key] = [d[key], value]
    else:
        d[key] = value


def normalize_key(s: str) -> str:
    s = s.strip().lower()
    # if the key contains no alphabetical letter, treat it as missing
    if not re.search(r"[a-z]", s):
        return None
    s = re.sub(r"\s+", "_", s)
    return s

def parse_lat_lon(s):
    """Assumes format '54 N 20 W' -> (54.0, -20.0)."""
    parts = s.strip().split()
    if len(parts) < 4:
        return None, None
    try:
        lat_deg = float(parts[0])
        lat_dir = parts[1].upper()
        lon_deg = float(parts[2])
        lon_dir = parts[3].upper()
    except ValueError:
        return None, None
    lat = lat_deg if lat_dir == 'N' else -lat_deg
    lon = lon_deg if lon_dir == 'E' else -lon_deg
    return lat, lon


def insert_value(d, key, value):
    if key in d:
        if isinstance(d[key], list) and value not in d[key]:
            d[key].append(value)
        if not isinstance(d[key], list) and d[key] != value:
            d[key] = [d[key], value]
    else:
        d[key] = value

@retry(stop=stop_after_attempt(5), wait=wait_fixed(5), reraise=True)
def query_acc(term, retstart, retmax):
    handle = Entrez.esearch(db="biosample", term=term, usehistory="y", timeout=DEFAULT_TIMEOUT)
    record = Entrez.read(handle)
    handle.close()
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    handle = Entrez.esummary(
        db="biosample",
        query_key=query_key,
        webenv=webenv,
        retmode="json",
        retmax=retmax,
        retstart=retstart,
        timeout=DEFAULT_TIMEOUT,
    )

    records = json.load(handle)
    handle.close()
    return records["result"]


def fetch_all_samples():
    """
    Fetch all Accession numbers from NCBI Biosample using Biopython Entrez ESearch with pagination.

    """

    # First, get total count
    handle = Entrez.esearch(db="biosample", term="all[filter]", usehistory="y", timeout=DEFAULT_TIMEOUT)
    record = Entrez.read(handle)
    total = int(record["Count"])
    logger.info(f"Total records to download: {total}")
    handle.close()

    retmax = 500  # max allowed by NCBI for json output
    for retstart in range(0, total, retmax):
        if retstart > 1000:
            break  # for testing, remove this in production

        accs = query_acc("all[filter]", retstart, retmax)
        yield accs
        if (retstart + retmax) % 10000 == 0:
            logger.info(f"Fetched {retstart + retmax} of {total} records")


def parse_xml(sample_dict, output, sample_mapping, nde_mapping):
    if cod := sample_dict["BioSample"].get("@access"):
        if cod == "public":
            output["conditionsOfAccess"] = "Open"
        else:
            output["conditionsOfAccess"] = "Restricted"

    attributes = {}
    if a := sample_dict["BioSample"].get("Attributes"):
        if attr_list := a.get("Attribute"):
            if not isinstance(attr_list, list):
                attr_list = [attr_list]
            for attr in attr_list:
                if attr.get("@harmonized_name"):
                    key = attr["@harmonized_name"]
                elif attr.get("@attribute_name"):
                    key = attr["@attribute_name"]

                if value := attr.get("#text"):
                    invalid_values = [
                        "-",
                        "?",
                        "blank",
                        "ignore",
                        "missing",
                        "na",
                        "n/a",
                        "nan",
                        "none",
                        "null",
                        "not applicable",
                        "not appicable",
                        "not collected",
                        "not determined",
                        "not provided",
                        "not recorded",
                        "not recorded",
                        "restricted access",
                        "unspecified",
                        "unk",
                        "unknown",
                    ]
                    if (
                        not attributes.get(key)
                        and value.casefold() not in invalid_values
                        and not value.casefold().startswith("missing:")
                    ):
                        key = normalize_key(key)
                        if key:
                            attributes[key] = value

    # subproperty, field_value is the normalized key, value pair from XML attributes
    # sample_mapping[subproperty] is the mapping info for the output
    # value is either the subproperty itself or the field_value depending on the mapping
    geo = {}
    for subproperty, field_value in attributes.items():
        if subproperty in sample_mapping:
            if len(sample_mapping[subproperty].keys()) > 1 and "locationOfOrigin.geo.latitude" in sample_mapping[subproperty]:
                try:
                    lat, lon = parse_lat_lon(field_value)
                    if lat is not None and lon is not None:
                        geo["latitude"] = lat
                        geo["longitude"] = lon
                        insert_value(output, {"locationOfOrigin": {"geo": geo}})
                        geo = {}  # reset geo after inserting
                    continue
                except Exception as e:
                    logger.warning(f"Failed to parse latitude and longitude from value '{field_value}' for subproperty '{subproperty}': {e}")
                    continue

            for k, v in sample_mapping[subproperty].items():
                value = subproperty if v == "property" else field_value
                if k == "additionalProperty":
                    insert_value(output, k, {"@type": "PropertyValue", "propertyID": subproperty, "value": field_value})
                elif "." in k:
                    key_split = k.split(".")
                    if len(key_split) == 2:
                        parent, child = key_split
                        insert_value(output, parent, {child: value})
                    if len(key_split) == 3:
                        parent, child1, child2 = key_split
                        if child2 == "latitude" or child2 == "longitude":
                            try:
                                value = float(value)
                            except ValueError:
                                logger.warning(f"Failed to convert latitude/longitude '{value}' to float for key '{k}'")
                                continue
                        if child2 in ["latitude", "longitude", "altitude", "depth"]:
                            insert_value(geo, child2, value)
                        else:
                            insert_value(output, parent, {child1: {child2: value}})
                elif k == "sampleAvailability":
                    try:
                        value = int(value)
                        value = False if value < 1 else True
                    except ValueError:
                        value = False if value.casefold() in ["false", "no"] else True
                    insert_value(output, k, value)
                else:
                    if k in nde_mapping and nde_mapping[k][0] == "object":
                        d = {nde_mapping[k][1]: value}
                        if k == "variableMeasured" or k == "anatomicalStructure":
                            d["@type"] = "DefinedTerm"
                        insert_value(output, k, d)
                    elif k in nde_mapping and nde_mapping[k][0] == "value":
                        if k == "date":
                            try:
                                if not output.get("date"):
                                    value = dateutil.parser.parse(value, ignoretz=True).date().isoformat()
                                    insert_value(output, k, value)
                            except Exception as e:
                                logger.warning(f"Error parsing date '{value}': {e}")
                        else:
                            insert_value(output, k, value)
                    else:
                        logger.warning(f"Unmapped nde_mapping property: {k}")
        else:
            insert_value(output, "additionalProperty", {"@type": "PropertyValue", "propertyID": subproperty, "value": field_value})

    if geo:
        insert_value(output, {"locationOfOrigin": {"geo": geo}})

def parse():

    with open(Path(__file__).resolve().parent / "nde_mapping.json", "r") as f:
        nde_mapping = json.load(f)
    with open(Path(__file__).resolve().parent / "mapping_dict.json", "r") as f:
        sample_mapping = json.load(f)

    for sample_list in fetch_all_samples():
        for key, sample in sample_list.items():
            if key == "uids":
                continue
            if sample.get("accession"):
                _id = sample.get("accession")
                url = f"https://www.ncbi.nlm.nih.gov/biosample/{key}"
            else:
                continue

            # parse at the end
            xml_string = sample.pop("sampledata", None)

            output = {
                "@context": "http://schema.org/",
                "@type": "Sample",
                "_id": _id.casefold(),
                "identifier": _id,
                "url": url,
                "distribution": [{"@type": "DataDownload", "contentUrl": url}],
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "BioSample",
                    "url": "https://www.ncbi.nlm.nih.gov/biosample/",
                    "versionDate": datetime.date.today().isoformat(),
                    "archivedAt": url,
                },
            }

            if name := sample.get("title"):
                output["name"] = name

            if date := sample.get("date"):
                output["date"] = dateutil.parser.parse(date, ignoretz=True).date().isoformat()

            if date_published := sample.get("publicationdate"):
                output["datePublished"] = dateutil.parser.parse(date_published, ignoretz=True).date().isoformat()

            if date_modified := sample.get("modificationdate"):
                output["dateModified"] = dateutil.parser.parse(date_modified, ignoretz=True).date().isoformat()

            if name := sample.get("organization"):
                output["author"] = {"affiliation": {"@type": "Organization", "name": name}}

            if sample.get("taxonomy") or sample.get("organism"):
                species = {}
                if sample.get("taxonomy"):
                    species["identifier"] = sample.get("taxonomy")
                if sample.get("organism"):
                    species["name"] = sample.get("organism")

            if ids := sample.get("identifiers"):
                alternate_identifiers = []
                parts = ids.split("; ")
                for part in parts:
                    if ":" in part:
                        key, value = [x.strip() for x in part.split(":", 1)]
                        if key != "BioSample":
                            alternate_identifiers.append(value)
                if alternate_identifiers:
                    output["alternateIdentifier"] = alternate_identifiers

            sample_dict = xmltodict.parse(xml_string) if xml_string else {}
            if sample_dict:
                parse_xml(sample_dict, output, sample_mapping, nde_mapping)
            yield output
