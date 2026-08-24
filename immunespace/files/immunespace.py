import csv
import datetime
import io
import logging
import re

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")

SOURCE_ORGANIZATION_URL = (
    "https://raw.githubusercontent.com/NIAID-Data-Ecosystem/nde-metadata-corrections/"
    "main/collections_corrections_production/HIPC_correction.json"
)
STUDIES_URL = "https://immunespace.org/api_kb/get_study_id_dropdown"
SIGNATURES_URL = "https://immunespace.org/api_kb/get_full_signatures_data/?"
SIGNATURE_RESULTS_URL = "https://immunespace.org/query/results/?ordering_tab=signatures_tab"
STUDY_URL_TEMPLATE = "https://immunespace.org/query/study/{study_id}"
REQUEST_TIMEOUT = 120

# Columns a signature record cannot be built without. Losing one of these is a
# real schema break and should stop the crawl.
SIGNATURE_FIELDS = {
    "Signature ID",
    "Response Component Type",
    "Description",
    "Up",
    "Down",
    "Changed",
    "Arm ID",
    "Study ID",
    "Disease",
}

# Columns that only enrich a record. Everything downstream already guards for
# their absence, so ImmuneSpace dropping one should be logged, not fatal --
# "Disease Stage" disappeared in August 2026 and took the whole crawl with it.
OPTIONAL_SIGNATURE_FIELDS = {"Disease Stage", "Material"}

# "Material" arrived in the same August 2026 change that removed "Disease Stage",
# and mixes three different kinds of value in one column: what the sample was
# taken from, what the subjects were given, and what infected them. Each belongs
# somewhere different, so they are sorted rather than dumped into one field.
# Anything unrecognized falls through to keywords, where it is at least
# searchable, and is logged so a new value gets noticed.
MATERIAL_CELL_TYPES = {
    "macrophage",
    "peripheral blood mononuclear cell",
    "t cell",
}

MATERIAL_ANATOMICAL_STRUCTURES = {
    "bone marrow",
    "colon",
    "ileum",
    "inguinal lymph node",
    "jejunum",
    "lung",
    "lymph node",
    "mesenteric lymph node",
    "pulmonary lymph node",
    "spleen",
    "thymus",
    "tonsil",
}

MATERIAL_SAMPLE_TYPES = {
    "blood",
    "blood plasma",
    "blood serum",
}

MATERIAL_INFECTIOUS_AGENTS = {
    "chikungunya virus",
    "dengue virus",
    "mycobacterium tuberculosis variant bovis bcg",
    "sars-cov-2",
}

# Carries no information about the specimen, so it is dropped rather than indexed.
MATERIAL_PLACEHOLDERS = {
    "pool of specimens",
    "specimen type: other",
}

CATALOG = {
    "@type": "DataCatalog",
    "name": "ImmuneSpace",
    "identifier": "ImmuneSpace",
    "url": "https://immunespace.org/",
    "archivedAt": SIGNATURE_RESULTS_URL,
}

RESPONSE_COMPONENTS = {
    "gene": {
        "@type": "DefinedTerm",
        "additional_type": "Gene",
        "measured_property": {
            "@type": "Property",
            "identifier": "SIO_001078",
            "inDefinedTermSet": "SIO",
            "name": "Gene Expression",
            "url": "http://semanticscience.org/resource/SIO_001078",
        },
        "observation_about_type": "Gene",
        "observation_type": {
            "@type": "DefinedTerm",
            "alternateName": "Differential Gene Expression",
            "identifier": "SIO_001078",
            "inDefinedTermSet": "SIO",
            "name": "differential gene expression ratio",
            "url": "http://semanticscience.org/resource/SIO_001078",
        },
    },
    "protein": {
        "@type": "DefinedTerm",
        "additional_type": "Protein",
        "measured_property": {"@type": "Property", "name": "Protein Abundance"},
        "observation_about_type": "Protein",
        "observation_type": {"@type": "DefinedTerm", "name": "Differential Protein Abundance"},
    },
    "metabolite": {
        "@type": "DefinedTerm",
        "additional_type": "Metabolite",
        "measured_property": {"@type": "Property", "name": "Metabolite Abundance"},
        "observation_about_type": "ChemicalSubstance",
        "observation_type": {"@type": "DefinedTerm", "name": "Differential Metabolite Abundance"},
    },
    "cell": {
        "@type": "DefinedTerm",
        "additional_type": "Cell",
        "measured_property": {"@type": "Property", "name": "Cell Abundance"},
        "observation_about_type": "DefinedTerm",
        "observation_type": {"@type": "DefinedTerm", "name": "Differential Cell Abundance"},
    },
    "weight": {
        "@type": "DefinedTerm",
        "additional_type": "Weight",
        "measured_property": {"@type": "Property", "name": "Body Weight"},
        "observation_about_type": "DefinedTerm",
        "observation_type": {"@type": "DefinedTerm", "name": "Differential Body Weight"},
    },
}

DIRECTIONS = {
    "Up": {
        "@type": "SemanticTriple",
        "gene_label": "Upregulated",
        "other_label": "Increased",
        "qualifier": "increase",
        "gene_identifier": "https://w3id.org/biolink/vocab/DirectionQualifierEnum#upregulated",
    },
    "Down": {
        "@type": "SemanticTriple",
        "gene_label": "Downregulated",
        "other_label": "Decreased",
        "qualifier": "decrease",
        "gene_identifier": "https://w3id.org/biolink/vocab/DirectionQualifierEnum#downregulated",
    },
    "Changed": {
        "@type": "SemanticTriple",
        "gene_label": "Differentially expressed",
        "other_label": "Changed",
        "qualifier": "changed",
        "gene_identifier": "https://semanticscience.org/resource/SIO_001078",
    },
}

NCIT_RELATIONSHIP_PROPERTY = "http://purl.obolibrary.org/obo/NCIT_C25648"
CURIE_RE = re.compile(r"^(?P<prefix>[A-Za-z][A-Za-z0-9_.-]*):(?P<identifier>[A-Za-z0-9_.-]+)$")
STUDY_ID_RE = re.compile(r"(SDY\d+)", re.IGNORECASE)


def _clean(value):
    if value is None:
        return None
    value = re.sub(r"\s+", " ", str(value).replace("\xa0", " ")).strip()
    return value or None


def _slug(value, max_len=160):
    value = _clean(value) or "unknown"
    value = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").casefold()
    return (value or "unknown")[:max_len]


def _append_unique(values, value):
    value = _clean(value)
    if value and value.casefold() not in {item.casefold() for item in values}:
        values.append(value)


def _single_or_list(values):
    if not values:
        return None
    return values[0] if len(values) == 1 else values


def _classify_material(value):
    """Sort a `Material` value into the group it belongs to, or None to drop it.

    A vaccine is recognized by name rather than by list, since ImmuneSpace keeps
    adding them; the rest are small, stable vocabularies.
    """
    key = (_clean(value) or "").casefold()
    if not key or key in MATERIAL_PLACEHOLDERS:
        return None
    if key in MATERIAL_CELL_TYPES:
        return "cell_types"
    if key in MATERIAL_ANATOMICAL_STRUCTURES:
        return "anatomical_structures"
    if key in MATERIAL_SAMPLE_TYPES:
        return "sample_types"
    if key in MATERIAL_INFECTIOUS_AGENTS:
        return "infectious_agents"
    # A vaccine has no dedicated NDE field, so it stays a keyword rather than
    # being forced into one that does not mean it.
    if "vaccine" in key or "vax" in key:
        return "material_keywords"
    logger.info("Unrecognized ImmuneSpace Material %r, keeping it as a keyword", value)
    return "material_keywords"


def _normalize_response_type(value):
    value = (_clean(value) or "").casefold()
    if value == "genes":
        return "gene"
    if value not in RESPONSE_COMPONENTS:
        raise ValueError(f"Unsupported ImmuneSpace Response Component Type: {value!r}")
    return value


def _study_id(signature_id, explicit_study_id=None):
    if study_id := _clean(explicit_study_id):
        return study_id.upper()
    match = STUDY_ID_RE.search(signature_id or "")
    if not match:
        raise ValueError(f"Could not determine Study ID for ImmuneSpace signature {signature_id!r}")
    return match.group(1).upper()


def _split_components(value):
    return [_clean(item) for item in str(value or "").split(";") if _clean(item)]


def _catalog(crawl_date):
    return {**CATALOG, "versionDate": crawl_date}


def _curie_fields(value):
    match = CURIE_RE.fullmatch(value)
    if not match:
        return {}
    prefix = match.group("prefix").upper()
    identifier = match.group("identifier")
    return {
        "identifier": f"{prefix}:{identifier}",
        "inDefinedTermSet": prefix,
        "url": f"http://purl.obolibrary.org/obo/{prefix}_{identifier}",
    }


def fetch_source_organization(requester=requests):
    """Fetch the sourceOrganization used by the existing ImmuneSpace records."""
    logger.info("Fetching sourceOrganization from: %s", SOURCE_ORGANIZATION_URL)
    response = requester.get(SOURCE_ORGANIZATION_URL, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    return response.json().get("sourceOrganization")


def _fetch_signature_rows(requester=requests):
    logger.info("Making request: %s", SIGNATURES_URL)
    response = requester.get(SIGNATURES_URL, headers={"Accept": "text/csv"}, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    text = response.content.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    fieldnames = set(reader.fieldnames or [])
    missing_fields = SIGNATURE_FIELDS - fieldnames
    if missing_fields:
        raise ValueError(f"ImmuneSpace signatures CSV is missing fields: {sorted(missing_fields)}")
    if missing_optional := OPTIONAL_SIGNATURE_FIELDS - fieldnames:
        logger.warning("ImmuneSpace signatures CSV no longer provides: %s", sorted(missing_optional))
    if new_fields := fieldnames - SIGNATURE_FIELDS - OPTIONAL_SIGNATURE_FIELDS:
        logger.info("ImmuneSpace signatures CSV has new unmapped fields: %s", sorted(new_fields))
    return reader


def _group_signature_rows(rows):
    groups = {}
    for row in rows:
        signature_id = _clean(row.get("Signature ID"))
        description = _clean(row.get("Description"))
        if not signature_id or not description:
            raise ValueError("ImmuneSpace signature rows require Signature ID and Description")

        response_type = _normalize_response_type(row.get("Response Component Type"))
        group = groups.setdefault(
            signature_id,
            {
                "anatomical_structures": [],
                "arms": [],
                "cell_types": [],
                "components": {direction: [] for direction in DIRECTIONS},
                "description": description,
                "disease_stages": [],
                "diseases": [],
                "infectious_agents": [],
                "material_keywords": [],
                "response_type": response_type,
                "sample_types": [],
                "studies": [],
            },
        )
        if group["description"] != description or group["response_type"] != response_type:
            raise ValueError(f"Conflicting metadata for ImmuneSpace signature {signature_id}")

        _append_unique(group["arms"], row.get("Arm ID"))
        _append_unique(group["studies"], _study_id(signature_id, row.get("Study ID")))
        _append_unique(group["diseases"], row.get("Disease"))
        _append_unique(group["disease_stages"], row.get("Disease Stage"))
        if bucket := _classify_material(row.get("Material")):
            _append_unique(group[bucket], row.get("Material"))
        for direction in DIRECTIONS:
            for component in _split_components(row.get(direction)):
                _append_unique(group["components"][direction], component)

    return groups


def _subject_dataset(group, crawl_date):
    primary_study = group["studies"][0]
    identifiers = [*group["studies"], *group["arms"]]
    return {
        "@type": "Dataset",
        "identifier": _single_or_list(identifiers),
        "includedInDataCatalog": _catalog(crawl_date),
        "name": primary_study,
        "_id": primary_study.casefold(),
        "url": STUDY_URL_TEMPLATE.format(study_id=primary_study),
    }


def _observation_about(component, response_config):
    return {
        "@type": response_config["observation_about_type"],
        "additionalType": response_config["additional_type"],
        "name": component,
        **_curie_fields(component),
    }


def _variable_measured(component, response_config):
    return {
        "@type": "StatisticalVariable",
        "additionalType": response_config["additional_type"],
        "constraintProperty": ["schema:healthCondition", "nde:sample"],
        "name": component,
        "populationType": "schema:Patient",
        "statType": "Direction",
        **_curie_fields(component),
    }


def _semantic_mapping(observation_about, subject_of, direction_config, direction_label, description, response_type):
    subject_identifier = observation_about.get("identifier")
    triple_subject = {
        "@type": "PropertyValue",
        "name": RESPONSE_COMPONENTS[response_type]["additional_type"],
        "propertyID": f"schema.org/{observation_about['@type']}",
        "value": observation_about["name"],
    }
    if subject_identifier:
        triple_subject["identifier"] = subject_identifier
    if observation_about.get("url"):
        triple_subject["url"] = observation_about["url"]

    triple_predicate = {
        "@type": "PropertyValue",
        "name": "Relationship",
        "propertyID": NCIT_RELATIONSHIP_PROPERTY,
        "value": direction_label,
    }
    if response_type == "gene":
        triple_predicate["identifier"] = direction_config["gene_identifier"]

    primary_study = subject_of["identifier"]
    if isinstance(primary_study, list):
        primary_study = primary_study[0]
    return {
        "@type": "SemanticTriple",
        "tripleSubject": triple_subject,
        "triplePredicate": triple_predicate,
        "tripleObject": {
            "@type": "PropertyValue",
            "identifier": primary_study,
            "name": "Dataset",
            "propertyID": "schema.org/Dataset",
            "url": subject_of["url"],
            "value": primary_study,
        },
        "tripleSubjectQualifier": [
            {
                "@type": "PropertyValue",
                "name": "aspect",
                "value": RESPONSE_COMPONENTS[response_type]["measured_property"]["name"].casefold(),
            },
            {"@type": "PropertyValue", "name": "direction", "value": direction_config["qualifier"]},
        ],
        "triplePredicateQualifier": [{"@type": "PropertyValue", "name": "comparison", "value": description}],
    }


def _build_signature_doc(signature_id, group, direction, component, source_organization, crawl_date):
    response_type = group["response_type"]
    response_config = RESPONSE_COMPONENTS[response_type]
    direction_config = DIRECTIONS[direction]
    direction_label = direction_config["gene_label"] if response_type == "gene" else direction_config["other_label"]
    subject_of = _subject_dataset(group, crawl_date)
    observation_about = _observation_about(component, response_config)
    primary_study = group["studies"][0]
    identifier = f"{signature_id}:{direction.casefold()}:{component}"
    doc = {
        "@context": "http://schema.org/",
        "@type": "Inference",
        "_id": (
            f"immunespace_signature_{_slug(signature_id.removeprefix('sig:'))}_"
            f"{direction.casefold()}_{_slug(component)}"
        ),
        "identifier": identifier,
        "name": f"{component} is {direction_label.casefold()} in {group['description']} ({primary_study})",
        "description": group["description"],
        "url": SIGNATURE_RESULTS_URL,
        "date": crawl_date,
        "includedInDataCatalog": _catalog(crawl_date),
        "subjectOf": subject_of,
        "observationType": dict(response_config["observation_type"]),
        "observationAbout": observation_about,
        "measuredProperty": dict(response_config["measured_property"]),
        "measurementQualifier": group["description"],
        "variableMeasured": _variable_measured(component, response_config),
        "semanticMapping": _semantic_mapping(
            observation_about,
            subject_of,
            direction_config,
            direction_label,
            group["description"],
            response_type,
        ),
        "keywords": [
            "Immune signature",
            response_config["additional_type"],
            direction_label,
            *group["diseases"],
            *group["disease_stages"],
            *group["material_keywords"],
        ],
        "species": {"@type": "DefinedTerm", "name": "Homo sapiens"},
    }
    if source_organization:
        doc["sourceOrganization"] = source_organization
    if group["diseases"]:
        conditions = [{"@type": "DefinedTerm", "name": disease} for disease in group["diseases"]]
        doc["healthCondition"] = _single_or_list(conditions)
    if group["infectious_agents"]:
        agents = [{"@type": "DefinedTerm", "name": agent} for agent in group["infectious_agents"]]
        doc["infectiousAgent"] = _single_or_list(agents)

    # The sample block is assembled from whatever the signature described about
    # the specimen: the disease stage it was taken at, and the material it was.
    sample = {"@type": "Sample"}
    if group["disease_stages"]:
        properties = [
            {
                "@type": "PropertyValue",
                "name": "disease stage",
                "propertyID": "immunespace:DiseaseStage",
                "value": stage,
            }
            for stage in group["disease_stages"]
        ]
        sample["additionalProperty"] = _single_or_list(properties)
    for bucket, field in (
        ("sample_types", "sampleType"),
        ("cell_types", "cellType"),
        ("anatomical_structures", "anatomicalStructure"),
    ):
        if group[bucket]:
            sample[field] = _single_or_list([{"@type": "DefinedTerm", "name": name} for name in group[bucket]])
    if len(sample) > 1:
        doc["sample"] = sample
    return doc


def parse_datasets(source_organization=None, requester=requests):
    """Yield the existing ImmuneSpace Dataset records."""
    if source_organization is None:
        source_organization = fetch_source_organization(requester=requester)

    logger.info("Making request: %s", STUDIES_URL)
    response = requester.get(STUDIES_URL, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    records = response.json()
    logger.info("Parsing ImmuneSpace Dataset records...")
    count = 0
    for count, record in enumerate(records, start=1):
        if count % 100 == 0:
            logger.info("Parsed %d ImmuneSpace Dataset records", count)

        output = {
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "ImmuneSpace",
                "url": "https://immunespace.org",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": f"https://immunespace.org/query/study/{record['value']}",
            },
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": record["value"],
            "url": f"https://immunespace.org/query/study/{record['value']}",
            "sourceOrganization": source_organization,
            "species": {"@type": "DefinedTerm", "name": "Homo sapiens"},
        }
        yield output

    logger.info("Finished parsing %d ImmuneSpace Dataset records", count)


def parse_signatures(source_organization=None, requester=requests, rows=None, crawl_date=None):
    """Yield feature-level ImmuneSpace signature records as NDE Inferences."""
    if source_organization is None:
        source_organization = fetch_source_organization(requester=requester)
    if rows is None:
        rows = _fetch_signature_rows(requester=requester)
    crawl_date = crawl_date or datetime.date.today().isoformat()

    groups = _group_signature_rows(rows)
    count = 0
    seen_ids = set()
    for signature_id, group in groups.items():
        for direction in DIRECTIONS:
            for component in group["components"][direction]:
                doc = _build_signature_doc(signature_id, group, direction, component, source_organization, crawl_date)
                if doc["_id"] in seen_ids:
                    raise ValueError(f"Duplicate ImmuneSpace Inference _id: {doc['_id']}")
                seen_ids.add(doc["_id"])
                count += 1
                yield doc

    logger.info("Finished parsing %d ImmuneSpace signature Inferences", count)


def parse(requester=requests):
    """Yield both Dataset and Inference records for the shared ImmuneSpace download."""
    source_organization = fetch_source_organization(requester=requester)
    yield from parse_datasets(source_organization=source_organization, requester=requester)
    yield from parse_signatures(source_organization=source_organization, requester=requester)
