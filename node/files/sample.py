import datetime
import logging
import re

import dateutil.parser

from .node import insert_value

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

def parse_ontology_term(text):
    """
    Parse a string like "water body [ENVO_00000063]" into:
      {"name": "water body", "identifier": "ENVO_00000063"}
    If there are no brackets, returns:
      {"name": "water body"}
    """
    if text is None:
        return None
    m = re.match(r"^(.+?)\s*\[([^\]]+)\]\s*$", text.strip())
    if m:
        return {"name": m.group(1).strip(), "identifier": m.group(2).strip()}
    return {"name": text.strip()}

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

def parse_latitude(text):
    """
    Parse a latitude string like "26.60695 N" or "26.60695 S" into a float.
    N is positive, S is negative. Returns None if not parseable.
    """
    if text is None:
        return None
    m = re.match(r"^\s*(\d+(?:\.\d+)?)\s*([NS])\s*$", str(text).strip(), re.IGNORECASE)
    if m:
        value = float(m.group(1))
        if m.group(2).upper() == "S":
            value = -value
        return value
    return None


def parse_longitude(text):
    """
    Parse a longitude string like "118.32417 W" or "118.32417 E" into a float.
    E is positive, W is negative. Returns None if not parseable.
    """
    if text is None:
        return None
    m = re.match(r"^\s*(\d+(?:\.\d+)?)\s*([EW])\s*$", str(text).strip(), re.IGNORECASE)
    if m:
        value = float(m.group(1))
        if m.group(2).upper() == "W":
            value = -value
        return value
    return None


def parse_temperature(text):
    """
    Parse a temperature string like "-80℃", "10C", or "10F" into:
      {"value": -80, "unitText": "Celsius"}
      {"value": 10, "unitText": "Fahrenheit"}
    Returns None if not parseable.
    """
    if text is None:
        return None
    m = re.match(r"^\s*(-?\d+(?:\.\d+)?)\s*(?:°|℃|℉)?\s*([CF])\s*$", str(text).strip(), re.IGNORECASE)
    if not m:
        # handle ℃ / ℉ symbols without a trailing letter
        m = re.match(r"^\s*(-?\d+(?:\.\d+)?)\s*(℃|℉)\s*$", str(text).strip())
    if m:
        raw_unit = m.group(2)
        if raw_unit in ("℃",) or raw_unit.upper() == "C":
            unit = "Celsius"
        else:
            unit = "Fahrenheit"
        return {"value": float(m.group(1)), "unitText": unit}
    return None


def parse_weight(text):
    """
    Parse a weight string like "68.8kg" into:
      {"name": "weight", "value": 68.8, "unitText": "kg"}
    Returns None if not parseable.
    """
    if text is None:
        return None
    m = re.match(r"^\s*(\d+(?:\.\d+)?)\s*([a-zA-Z]+)\s*$", str(text).strip())
    if m:
        return {"name": "weight", "value": float(m.group(1)), "unitText": m.group(2)}
    else:
        return {"name": text.strip()}


def parse_height(text):
    """
    Parse a height string like "170cm" into:
      {"name": "height", "value": 170.0, "unitText": "cm"}
    Returns None if not parseable.
    """
    if text is None:
        return None
    m = re.match(r"^\s*(\d+(?:\.\d+)?)\s*([a-zA-Z]+)\s*$", str(text).strip())
    if m:
        return {"name": "height", "value": float(m.group(1)), "unitText": m.group(2)}
    return {"name": text.strip()}


def parse_value_unit(text):
    """
    Parse a generic value+unit string like "5g" into:
      {"value": 5.0, "unitText": "g"}
    Returns None if not parseable.
    """
    if text is None:
        return None
    m = re.match(r"^\s*(-?\d+(?:\.\d+)?)\s*([a-zA-Z]+)\s*$", str(text).strip())
    if m:
        return {"value": float(m.group(1)), "unitText": m.group(2)}
    return {"name": text.strip()}


def find_age(text):
    """
    Find an age integer in a string.
    Matches patterns like:
      - "age: 60"
      - "age 60"
      - "age-60"
      - "male;age: 20"
      - "60 years old"
    Returns int age or None.
    """
    if text is None:
        return None
    s = str(text).lower()
    # common "age" followed by optional "(years)" then separator/space then number
    m = re.search(
        r"\bage\b\s*(?:\(\s*years?\s*\)\s*)?(?:\s*[:=\-]\s*|\s+)(\d{1,3})\b",
        s,
    )
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    # fallback: "<number> years old" or "<number> yrs"
    m2 = re.search(r"\b(\d{1,3})\s*(?:years?|yrs?|y)\b", s)
    if m2:
        try:
            return int(m2.group(1))
        except ValueError:
            return None
    return None

def get_samples(session, table_info, project_id):
    if table_info.get("data") and table_info["data"].get("sapTables"):
        for sample_table in table_info["data"]["sapTables"]:
            data_type = sample_table["type"]
            sap_info = session.get(f"https://www.biosino.org/node/api/app/project/getExpAndSampleList?projectNo={project_id}&type=sample&dataType={data_type}&total=0&pageNum=1&pageSize=100&sortKey=sapNo&sortType=asc").json()
            pages = sap_info["data"]["sapTableData"]["totalPages"]
            for page in range(1, pages + 1):
                print(f"Crawling sample info page {page} of {pages} for project {project_id}")
                sap_info = session.get(f"https://www.biosino.org/node/api/app/project/getExpAndSampleList?projectNo={project_id}&type=sample&dataType={data_type}&total=0&pageNum={page}&pageSize=100&sortKey=sapNo&sortType=asc").json()
                for sap in sap_info["data"]["sapTableData"]["content"]:
                    yield sap

def parse_sample(sample):
    if _id := sample.get("sapNo"):
        url = f"https://www.biosino.org/node/sample/detail/{_id}"
    else:
        logger.warning(f"Sample missing sapNo: {sample}")
        return None

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

    if name := sample.get("name"):
        insert_value(output, "name", name)

    if identifier := sample.get("taxId"):
        insert_value(output.setdefault("species", {}), "identifier", identifier)
        insert_value(output.setdefault("infectiousAgent", {}), "identifier", identifier)

    if name := sample.get("organism"):
        insert_value(output.setdefault("species", {}), "name", name)
        insert_value(output.setdefault("infectiousAgent", {}), "name", name)

    if date_created := sample.get("createDate"):
        try:
            date_created = dateutil.parser.parse(date_created).date().isoformat()
            insert_value(output, "dateCreated", date_created)
        except (ValueError, TypeError):
            logger.info(f"Failed to parse dateCreated: {date_created} for sample {_id}")

    if date_modified := sample.get("updateDate"):
        try:
            date_modified = dateutil.parser.parse(date_modified).date().isoformat()
            insert_value(output, "dateModified", date_modified)
        except (ValueError, TypeError):
            logger.info(f"Failed to parse dateModified: {date_modified} for sample {_id}")

    if author := sample.get("submitter"):
        if given_name := author.get("firstName"):
            insert_value(output.setdefault("author", {}), "givenName", given_name)
        if family_name := author.get("lastName"):
            insert_value(output.setdefault("author", {}), "familyName", family_name)
        if name := author.get("orgName"):
            insert_value(output.setdefault("author", {}).setdefault("affiliation", {}), "name", name)

    if name := sample.get("subjectType"):
        insert_value(output, "sampleType", {"name": name})

    if description := sample.get("description"):
        insert_value(output, "description", description)

    if name := sample.get("tissue"):
        insert_value(output, "anatomicalStructure", {"name": name})

    if identifier := sample.get("usedIds"):
        if not isinstance(identifier, list):
            identifier = [identifier]
        for id in identifier:
            insert_value(output, "isRelatedTo", {"identifier": id})

    if date_published := sample.get("publicDate"):
        try:
            date_published = dateutil.parser.parse(date_published).date().isoformat()
            insert_value(output, "datePublished", date_published)
        except (ValueError, TypeError):
            logger.info(f"Failed to parse datePublished: {date_published} for sample {_id}")

    if sample_process := sample.get("protocol"):
        insert_value(output, "sampleProcess", sample_process)


    if attributes := sample.get("attributes"):
        for key, value in attributes.items():
            if not value:
                continue

            if key == "biomaterial_provider":
                insert_value(output.setdefault("author", {}), "name", value)
            elif key == "host":
                insert_value(output, "species", {"name": value})
            elif key == "environmental_package":
                insert_value(output, "environmentalSystem", {"name": value})
            elif key == "host_sex":
                insert_value(output, "sex", value)
            elif key == "host_age":
                if age := find_age(attributes.get("host_age")):
                    insert_value(output, "developmentalStage", {"value": age, "unitText": "year"})
                else:
                    insert_value(output, "developmentalStage", {"name": value})
            elif key == "env_biome":
                insert_value(output, "environmentalSystem", parse_ontology_term(value))
            elif key == "env_feature":
                insert_value(output, "environmentalSystem", parse_ontology_term(value))
            elif key == "env_material":
                insert_value(output, "environmentalSystem", parse_ontology_term(value))
            elif key == "geographic_location":
                insert_value(output.setdefault("locationOfOrigin", {}), "name", value)
            elif key == "latitude":
                if latitude := parse_latitude(value):
                    insert_value(output.setdefault("locationOfOrigin", {}).setdefault("geo", {}), "latitude", latitude)
            elif key == "longitude":
                if longitude := parse_longitude(value):
                    insert_value(output.setdefault("locationOfOrigin", {}).setdefault("geo", {}), "longitude", longitude)
            elif key == "sample_collection_date":
                try:
                    value = dateutil.parser.parse(value).date().isoformat()
                    insert_value(output, "dateCollected", value)
                except (ValueError, TypeError):
                    logger.info(f"Failed to parse dateCollected: {value} for sample {_id}")
            elif key == "oxygenation_status_of_sample":
                insert_value(output, "sampleState", value)
            elif key == "developmental_stage":
                insert_value(output, "developmentalStage", {"name": value})
            elif key == "sex":
                insert_value(output, "sex", value)
            elif key == "age":
                if age := find_age(value):
                    insert_value(output, "developmentalStage", {"value": age, "unitText": "year"})
                else:
                    insert_value(output, "developmentalStage", {"name": value})
            elif key == "disease_name":
                insert_value(output, "healthCondition", {"name": value})
            elif key == "cell_line_name":
                insert_value(output, "cellType", {"name": value})
            elif key == "cell_type":
                insert_value(output, "cellType", {"name": value})
            elif key == "breed":
                insert_value(output, "associatedGenotype", value)
            elif key == "molecular_types_extracted_from_samples":
                insert_value(output, "sampleType", {"name": value})
            elif key == "strain_name":
                insert_value(output, "associatedGenotype", value)
            elif key == "sample_storage_temperature":
                if temp := parse_temperature(value):
                    insert_value(output, "sampleStorageTemperature", temp)
                else:
                    insert_value(output, "sampleStorageTemperature", {"name": value})
            elif key == "elevation":
                insert_value(output.setdefault("locationOfOrigin", {}).setdefault("geo", {}), "elevation", value)
            elif key == "weight":
                insert_value(output, "associatedPhenotype", parse_weight(value))
            elif key == "host_weight":
                insert_value(output, "associatedPhenotype", parse_weight(value))
            elif key == "host_height":
                insert_value(output, "associatedPhenotype", parse_height(value))
            elif key == "genotype":
                insert_value(output, "associatedGenotype", value)
            elif key == "bioproject_accession":
                insert_value(output, "isBasisFor", {"identifier": value})
            elif key == "sample_collection_device":
                insert_value(output, "instrument", {"name": value})
            elif key == "isolation_source":
                insert_value(output, "keywords", value)
            elif key == "disease_pathological_type":
                insert_value(output, "healthCondition", {"name": value})
            elif key == "mutation_type":
                insert_value(output, "associatedGenotype", value)
            elif key == "cultivar":
                insert_value(output, "associatedGenotype", value)
            elif key == "lat_lon":
                geo = {}
                try:
                    lat, lon = parse_lat_lon(value)
                    if lat is not None and lon is not None:
                        geo["latitude"] = lat
                        geo["longitude"] = lon
                        insert_value(output, "locationOfOrigin", {"geo": geo})
                    continue
                except Exception as e:
                    logger.warning(f"Failed to parse latitude and longitude from value '{value}': {e}")
                    continue
            elif key == "latitude_start" or key == "longitude_start":
                if key == "latitude_start":
                    try:
                        geo = {}
                        latitude = float(value)
                        geo["latitude"] = latitude
                        if attributes.get("longitude_start"):
                            longitude = float(attributes.get("longitude_start"))
                            geo["longitude"] = longitude
                        insert_value(output, "locationOfOrigin", {"geo": geo})
                    except Exception as e:
                        logger.warning(f"Failed to parse latitude_start or longitude_start from value '{value}': {e}")
                        continue
            elif key == "latitude_end" or key == "longitude_end":
                if key == "latitude_end":
                    try:
                        geo = {}
                        latitude = float(value)
                        geo["latitude"] = latitude
                        if attributes.get("longitude_end"):
                            longitude = float(attributes.get("longitude_end"))
                            geo["longitude"] = longitude
                        insert_value(output, "locationOfOrigin", {"geo": geo})
                    except Exception as e:
                        logger.warning(f"Failed to parse latitude_end or longitude_end from value '{value}': {e}")
                        continue
            elif key == "altitude":
                insert_value(output.setdefault("locationOfOrigin", {}).setdefault("geo", {}), "altitude", value)
            elif key == "sample_storage_location":
                insert_value(output, "itemLocation", {"name": value})
            elif key == "sample_volume_or_weight_for_dna_extraction":
                insert_value(output, "sampleQuantity", parse_value_unit(value))
            elif key == "natural_host":
                insert_value(output, "species", parse_ontology_term(value))
            elif key == "sample_collection_method":
                insert_value(output, "collectionMethod", value)
            elif key == "height":
                insert_value(output, "associatedPhenotype", parse_height(value))
            elif key == "phenotype":
                insert_value(output, "associatedPhenotype", {"name": value})
            elif key == "age_at_surgery":
                if age := find_age(value):
                    insert_value(output, "developmentalStage", {"value": age, "unitText": "year"})
                else:
                    insert_value(output, "developmentalStage", {"name": value})
            elif key == "pathogen_type":
                insert_value(output, "infectiousAgent", {"name": value})
            elif key == "pathogen_subtype":
                insert_value(output, "infectiousAgent", {"name": value})
            else:
                insert_value(output, "additionalProperty", {"@type": "PropertyValue", "propertyID": key, "value": value})


    if cinf := sample.get("calc_info"):
        for key, value in cinf.items():
            if not value:
                continue
            if key == "update_time":
                try:
                    value = dateutil.parser.parse(value).date().isoformat()
                    insert_value(output, "dateModified", value)
                except (ValueError, TypeError):
                    logger.info(f"Failed to parse dateModified from calc_info: {value} for sample {_id}")
            elif key == "exp_type":
                insert_value(output, "keywords", value)
            else:
                insert_value(output, "additionalProperty", {"@type": "PropertyValue", "propertyID": key, "value": value})

    if cattributes := sample.get("customAttr"):
        for key, value in cattributes.items():
            if not value:
                continue
            if key == "phenotype":
                insert_value(output, "associatedPhenotype", {"name": value})
            elif key == "host_species":
                insert_value(output, "species", {"name": value})
            elif key == "breed":
                insert_value(output, "associatedGenotype", value)
            elif key == "natural_host":
                insert_value(output, "species", {"name": value})
            elif key == "sourceProject":
                insert_value(output, "keywords", value)
            else:
                insert_value(output, "additionalProperty", {"@type": "PropertyValue", "propertyID": key, "value": value})

    return output



def add_aggregate_element(aggregate_element, sample):
    keys = [
        "additionalType",
        "associatedGenotype",
        "associatedPhenotype",
        "anatomicalStructure",
        "anatomicalSystem",
        "cellType",
        "developmentalStage",
        "experimentalPurpose",
        "sampleAvailability",
        "sampleQuantity",
        "sampleState",
        "sampleStorageTemperature",
        "sampleType",
        "sex",
        "url",
    ]

    for key in keys:
        if key in sample:
            insert_value(aggregate_element, key, sample[key])







