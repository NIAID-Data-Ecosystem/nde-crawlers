import csv
import datetime
import io
import json
import logging
import os
import re
import time
from typing import Iterable
from urllib.parse import quote, urlencode, urljoin

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

API_ROOT = "https://www.ebi.ac.uk/gxa"
EXPERIMENTS_URL = f"{API_ROOT}/json/experiments"
EXPERIMENT_JSON_URL_TEMPLATE = f"{API_ROOT}/json/experiments/{{accession}}"
EXPERIMENT_RESOURCES_URL_TEMPLATE = f"{API_ROOT}/json/experiments/{{accession}}/resources/DATA"
EXPERIMENT_RESULTS_URL_TEMPLATE = f"{API_ROOT}/experiments/{{accession}}/Results"
EXPERIMENT_DOWNLOAD_URL_TEMPLATE = f"{API_ROOT}/experiments-content/{{accession}}/download/{{raw_type}}"

TIMEOUT = int(os.getenv("GXA_TIMEOUT", "120"))
SLEEP = float(os.getenv("GXA_SLEEP", "0.25"))
P_VALUE_CUTOFF = os.getenv("GXA_P_VALUE_CUTOFF", "0.05")
MAX_EXPERIMENTS = int(os.getenv("GXA_MAX_EXPERIMENTS", "0")) or None
MAX_RECORDS = int(os.getenv("GXA_MAX_RECORDS", "0")) or None
SPECIES_FILTER = os.getenv("GXA_SPECIES", "").strip().casefold()
FETCH_DESIGN_TERMS = os.getenv("GXA_FETCH_DESIGN_TERMS", "true").casefold() not in {"0", "false", "no"}

RESULT_COLUMN_RE = re.compile(r"^(?P<contrast>.+?)\s*\.(?P<metric>foldChange|pValue)$")
DESIGN_VALUE_COLUMN_RE = re.compile(r"^(?P<kind>Sample Characteristic|Factor Value)\[(?P<property>.+)]$")
DESIGN_TERM_COLUMN_RE = re.compile(
    r"^(?P<kind>Sample Characteristic Ontology Term|Factor Value Ontology Term)\[(?P<property>.+)]$"
)

INVALID_VALUES = {
    "",
    "-",
    "--",
    "na",
    "n/a",
    "none",
    "null",
    "not applicable",
    "not available",
    "unknown",
    "unspecified",
}

HEALTH_CONDITION_SKIP_VALUES = {
    "control",
    "healthy",
    "healthy control",
    "normal",
    "normal control",
    "not applicable",
    "not available",
    "untreated",
    "wild type",
}

INFECTIOUS_AGENT_SKIP_VALUES = {
    "control",
    "mock",
    "none",
    "normal",
    "not applicable",
    "not available",
    "untreated",
}

PROPERTY_NAME_ALIASES = {
    "development stage": "developmental stage",
    "developmentalstage": "developmental stage",
    "organismpart": "organism part",
    "sample with known storage state": "specimen with known storage state",
}

SAMPLE_PROPERTIES_WITH_SPECIFIC_FIELDS = {
    "cell type",
    "developmental stage",
    "disease",
    "infect",
    "infectious agent",
    "organism",
    "organism part",
    "pathogen",
    "sex",
    "specimen with known storage state",
}

CATALOG = {
    "@type": "DataCatalog",
    "name": "Gene Expression Atlas",
    "alternateName": "GXA",
    "identifier": "GXA",
    "url": API_ROOT + "/",
}

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "identifier": "topic_0203",
        "inDefinedTermSet": "EDAM",
        "name": "Gene expression",
        "url": "http://edamontology.org/topic_0203",
    }
]

OBSERVATION_TYPE = {
    "@type": "DefinedTerm",
    "alternateName": "Differential Gene Expression",
    "identifier": "SIO_001078",
    "inDefinedTermSet": "SIO",
    "name": "differential gene expression ratio",
    "url": "http://semanticscience.org/resource/SIO_001078",
}

MEASURED_PROPERTY = {
    "@type": "Property",
    "identifier": "SIO_001078",
    "inDefinedTermSet": "SIO",
    "name": "Gene Expression",
    "url": "http://semanticscience.org/resource/SIO_001078",
}

NCIT_RELATIONSHIP_PROPERTY = "http://purl.obolibrary.org/obo/NCIT_C25648"
FOLD_CHANGE_UNIT_CODE = "https://www.wikidata.org/wiki/Q17014303"

DIRECTION_TERMS = {
    "up": {
        "label": "Upregulated",
        "qualifier": "increase",
        "identifier": "https://w3id.org/biolink/vocab/DirectionQualifierEnum#upregulated",
    },
    "down": {
        "label": "Downregulated",
        "qualifier": "decrease",
        "identifier": "https://w3id.org/biolink/vocab/DirectionQualifierEnum#downregulated",
    },
}

COVID_19_INFECTIOUS_AGENT = {
    "@type": "DefinedTerm",
    "alternateName": ["SARS-CoV-2", "SARS-ncov-2"],
    "identifier": "2697049",
    "inDefinedTermSet": "NCBITaxon",
    "name": "severe acute respiratory syndrome coronavirus 2",
    "url": "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2697049",
}


def _split_env_values(name):
    value = os.getenv(name, "")
    if not value.strip():
        return []
    return [part.strip() for part in re.split(r"[,\s]+", value) if part.strip()]


ACCESSIONS = _split_env_values("GXA_ACCESSIONS")


def clean_value(value):
    if value is None:
        return None
    value = str(value).replace("\xa0", " ")
    value = re.sub(r"\s+", " ", value).strip()
    if value.casefold() in INVALID_VALUES:
        return None
    return value


def _norm(value):
    value = clean_value(value)
    return value.casefold() if value else None


def _property_key(value):
    key = _norm(value)
    return PROPERTY_NAME_ALIASES.get(key, key)


def _property_lookup_names(value):
    name = clean_value(value)
    names = []
    if name:
        names.append(name)
    key = _property_key(name)
    if key and key not in {item.casefold() for item in names}:
        names.append(key)
    return names


def _single_or_list(values):
    values = _unique(values)
    if not values:
        return None
    return values[0] if len(values) == 1 else values


def _unique(values):
    seen = set()
    unique_values = []
    for value in values:
        if value in (None, "", [], {}):
            continue
        if isinstance(value, dict):
            key = tuple(sorted((k, str(v)) for k, v in value.items()))
        else:
            key = str(value).casefold()
        if key not in seen:
            unique_values.append(value)
            seen.add(key)
    return unique_values


def _to_float(value):
    value = clean_value(value)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _parse_gxa_date(value):
    value = clean_value(value)
    if not value:
        return None
    for fmt in ("%d-%m-%Y", "%Y-%m-%d"):
        try:
            return datetime.datetime.strptime(value, fmt).date().isoformat()
        except ValueError:
            pass
    return None


def _slug(value, max_len=120):
    value = clean_value(value) or "unknown"
    value = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").casefold()
    return (value or "unknown")[:max_len]


def _absolute_url(path):
    if not path:
        return None
    return urljoin(API_ROOT + "/", path)


def _request(session, url, response_type="json", required=True, retries=3):
    headers = {
        "Accept": "application/json" if response_type == "json" else "text/plain,application/json",
        "User-Agent": "nde-gxa-crawler/0.1",
    }
    for attempt in range(1, retries + 1):
        try:
            response = session.get(url, headers=headers, timeout=TIMEOUT)
            if response.status_code in {429, 500, 502, 503, 504} and attempt < retries:
                logger.warning("GXA API returned %s for %s; retrying", response.status_code, url)
                time.sleep(attempt)
                continue
            if response.status_code in {400, 404} and not required:
                logger.warning("Skipping unavailable GXA URL %s: HTTP %s", url, response.status_code)
                return None
            response.raise_for_status()
            return response.json() if response_type == "json" else response.text
        except (requests.RequestException, ValueError) as exc:
            if attempt == retries:
                if required:
                    raise
                logger.warning("Skipping failed GXA URL %s: %s", url, exc)
                return None
            logger.warning("GXA request failed for %s: %s; retrying", url, exc)
            time.sleep(attempt)
    return None


def _iter_experiment_summaries(session):
    if ACCESSIONS:
        for accession in ACCESSIONS:
            yield {"experimentAccession": accession}
        return

    data = _request(session, EXPERIMENTS_URL)
    experiments = data.get("experiments", []) if isinstance(data, dict) else []
    yielded = 0
    for experiment in experiments:
        raw_type = clean_value(experiment.get("rawExperimentType")) or ""
        experiment_type = clean_value(experiment.get("experimentType")) or ""
        if experiment_type != "Differential" and "DIFFERENTIAL" not in raw_type:
            continue
        if SPECIES_FILTER and SPECIES_FILTER not in (experiment.get("species") or "").casefold():
            continue

        yield experiment
        yielded += 1
        if MAX_EXPERIMENTS and yielded >= MAX_EXPERIMENTS:
            return


def _fetch_experiment(session, accession):
    url = EXPERIMENT_JSON_URL_TEMPLATE.format(accession=accession)
    return _request(session, url, required=False)


def _raw_experiment_type(summary, experiment_payload):
    raw_type = clean_value(summary.get("rawExperimentType"))
    if raw_type:
        return raw_type
    raw_type = clean_value((experiment_payload.get("experiment") or {}).get("type"))
    return raw_type.upper() if raw_type else None


def _download_results(session, accession, raw_type):
    params = {
        "geneQuery": "[]",
        "unit": "FOLD_CHANGE",
        "cutoff": P_VALUE_CUTOFF,
        "heatmapMatrixSize": "50",
        "selectedColumnIds": "",
        "type": raw_type,
    }
    url = EXPERIMENT_DOWNLOAD_URL_TEMPLATE.format(accession=accession, raw_type=raw_type)
    return _request(session, f"{url}?{urlencode(params)}", response_type="text", required=False)


def _fetch_design_terms(session, accession):
    if not FETCH_DESIGN_TERMS:
        return {}

    resources_url = EXPERIMENT_RESOURCES_URL_TEMPLATE.format(accession=accession)
    resources = _request(session, resources_url, required=False)
    if not isinstance(resources, list):
        return {}

    design_path = None
    for resource in resources:
        description = clean_value(resource.get("description")) or ""
        if resource.get("type") == "icon-experiment-design" or "experiment design" in description.casefold():
            design_path = resource.get("url")
            break
    if not design_path:
        return {}

    design_text = _request(session, _absolute_url(design_path), response_type="text", required=False)
    return _parse_design_terms(design_text) if design_text else {}


def _parse_design_terms(design_text):
    reader = csv.reader(io.StringIO(design_text), delimiter="\t")
    header = None
    for row in reader:
        if row and not row[0].startswith("#"):
            header = row
            break
    if not header:
        return {}

    term_columns = {}
    value_columns = []
    for index, column in enumerate(header):
        if match := DESIGN_TERM_COLUMN_RE.match(column):
            kind = match.group("kind").replace(" Ontology Term", "")
            term_columns[(kind, match.group("property").casefold())] = index
        elif match := DESIGN_VALUE_COLUMN_RE.match(column):
            value_columns.append((index, match.group("kind"), match.group("property")))

    terms = {}
    for row in reader:
        for value_index, kind, property_name in value_columns:
            if value_index >= len(row):
                continue
            value = clean_value(row[value_index])
            if not value:
                continue
            term_index = term_columns.get((kind, property_name.casefold()))
            if term_index is None or term_index >= len(row):
                continue
            url = clean_value(row[term_index])
            if not url:
                continue
            terms.setdefault((property_name.casefold(), value.casefold()), url)
    return terms


def _iter_tsv_results(tsv_text):
    reader = csv.reader(io.StringIO(tsv_text), delimiter="\t")
    header = None
    for row in reader:
        if not row:
            continue
        if row[0].startswith("#"):
            continue
        header = row
        break
    if not header:
        return

    contrast_columns = {}
    for index, column in enumerate(header[2:], start=2):
        match = RESULT_COLUMN_RE.match(column.strip())
        if not match:
            continue
        contrast = match.group("contrast").strip()
        metric = match.group("metric")
        contrast_columns.setdefault(contrast, {})[metric] = index

    for row in reader:
        if len(row) < 2:
            continue
        gene_id = clean_value(row[0])
        gene_name = clean_value(row[1]) or gene_id
        if not gene_id:
            continue

        for contrast_name, columns in contrast_columns.items():
            fold_change = _to_float(_row_value(row, columns.get("foldChange")))
            if fold_change is None:
                continue
            p_value = _to_float(_row_value(row, columns.get("pValue")))
            yield {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "contrast_name": contrast_name,
                "fold_change": fold_change,
                "p_value": p_value,
            }


def _row_value(row, index):
    if index is None or index >= len(row):
        return None
    return row[index]


def _ontology_parts(url):
    url = clean_value(url)
    if not url:
        return None, None
    if "ncbi.nlm.nih.gov/Taxonomy" in url and "id=" in url:
        return url.rsplit("id=", 1)[-1], "NCBITaxon"
    identifier = url.rstrip("/").rsplit("/", 1)[-1].rsplit("#", 1)[-1]
    if not identifier:
        return None, None
    if "_" in identifier:
        prefix = identifier.split("_", 1)[0]
        return identifier, prefix
    return identifier, None


def _defined_term(name, url=None, term_type="DefinedTerm"):
    name = clean_value(name)
    if not name:
        return None
    term = {"@type": term_type, "name": name}
    if url := clean_value(url):
        term["url"] = url
        identifier, term_set = _ontology_parts(url)
        if identifier:
            term["identifier"] = identifier
        if term_set:
            term["inDefinedTermSet"] = term_set
    return term


def _term_from_design(property_name, value, design_terms, term_type="DefinedTerm"):
    value = clean_value(value)
    if not value:
        return None
    url = None
    for name in _property_lookup_names(property_name):
        url = design_terms.get((name.casefold(), value.casefold()))
        if url:
            break
    return _defined_term(value, url=url, term_type=term_type)


def _contrast_properties(contrast):
    summary = contrast.get("contrastSummary") or {}
    return [prop for prop in summary.get("properties") or [] if isinstance(prop, dict)]


def _contrast_property_values(properties, property_name):
    values = []
    target_key = _property_key(property_name)
    for prop in properties:
        if _property_key(prop.get("propertyName")) != target_key:
            continue
        values.extend([clean_value(prop.get("testValue")), clean_value(prop.get("referenceValue"))])
    return _unique(values)


def _factor_condition_text(properties, side):
    parts = []
    for prop in properties:
        prop_name = clean_value(prop.get("propertyName"))
        if not prop_name:
            continue
        is_factor = prop.get("contrastPropertyType") == "FACTOR"
        prop_key = _property_key(prop_name)
        if not is_factor and prop_key not in {"disease", "organism part", "cell type", "genotype"}:
            continue
        value = clean_value(prop.get(f"{side}Value"))
        if value:
            parts.append(f"{prop_name}: {value}")
    return "; ".join(parts) if parts else None


def _primary_condition(properties, side):
    parts = []
    for prop in properties:
        if prop.get("contrastPropertyType") != "FACTOR":
            continue
        value = clean_value(prop.get(f"{side}Value"))
        if value:
            parts.append(value)
    return ", ".join(parts) if parts else None


def _population_type(species):
    return "schema:Patient" if clean_value(species) == "Homo sapiens" else "schema:BioSample"


def _statistical_variable(gene_name, gene_id, condition, species):
    variable = {
        "@type": "StatisticalVariable",
        "additionalType": "Gene",
        "constraintProperty": ["schema:healthCondition", "nde:sample"],
        "identifier": gene_id,
        "name": gene_name or gene_id,
        "populationType": _population_type(species),
        "statType": "Fold Change",
    }
    if condition:
        variable["description"] = f"{gene_name or gene_id} expression in {condition}"
    return variable


def _gene_term(gene_id, gene_name):
    term = {
        "@type": "Gene",
        "identifier": gene_id,
        "name": gene_name or gene_id,
        "url": f"{API_ROOT}/genes/{quote(gene_id)}",
    }
    if gene_name and gene_name != gene_id:
        term["alternateName"] = gene_id
    return term


def _result_url(accession, gene_id, gene_name):
    value = gene_name or gene_id
    category = "symbol" if gene_name and gene_name != gene_id else "gene_id"
    gene_query = quote(json.dumps([{"value": value, "category": category}], separators=(",", ":")))
    return f"{EXPERIMENT_RESULTS_URL_TEMPLATE.format(accession=accession)}?geneQuery={gene_query}"


def _dataset_subject(accession, description):
    url = EXPERIMENT_RESULTS_URL_TEMPLATE.format(accession=accession)
    dataset = {
        "@type": "Dataset",
        "identifier": accession,
        "includedInDataCatalog": dict(CATALOG),
        "name": description or accession,
        "url": url,
    }
    return dataset


def _direction(fold_change):
    if fold_change > 0:
        return DIRECTION_TERMS["up"]
    if fold_change < 0:
        return DIRECTION_TERMS["down"]
    return {
        "label": "Differentially expressed",
        "qualifier": "changed",
        "identifier": "https://semanticscience.org/resource/SIO_001078",
    }


def _semantic_mapping(gene, direction, subject_of, measurement_qualifier):
    mapping = {
        "@type": "SemanticTriple",
        "tripleSubject": {
            "@type": "PropertyValue",
            "identifier": gene.get("identifier"),
            "name": "Gene",
            "propertyID": "schema.org/Gene",
            "url": gene.get("url"),
            "value": gene.get("name"),
        },
        "triplePredicate": {
            "@type": "PropertyValue",
            "identifier": direction["identifier"],
            "name": "Relationship",
            "propertyID": NCIT_RELATIONSHIP_PROPERTY,
            "value": direction["label"],
        },
        "tripleObject": {
            "@type": "PropertyValue",
            "identifier": subject_of.get("identifier"),
            "name": "Dataset",
            "propertyID": "schema.org/Dataset",
            "url": subject_of.get("url"),
            "value": subject_of.get("identifier"),
        },
        "tripleSubjectQualifier": [
            {"@type": "PropertyValue", "name": "aspect", "value": "abundance"},
            {"@type": "PropertyValue", "name": "direction", "value": direction["qualifier"]},
        ],
    }
    if measurement_qualifier:
        mapping["triplePredicateQualifier"] = [
            {"@type": "PropertyValue", "name": "comparison", "value": measurement_qualifier}
        ]
    return mapping


def _health_conditions(properties, design_terms):
    health_conditions = []
    for value in _contrast_property_values(properties, "disease"):
        if value.casefold() in HEALTH_CONDITION_SKIP_VALUES:
            continue
        term = _term_from_design("disease", value, design_terms)
        if term:
            health_conditions.append(term)
    return _single_or_list(health_conditions)


def _species_terms(properties, design_terms, experiment_species):
    species_values = _contrast_property_values(properties, "organism") or [experiment_species]
    species = []
    for value in species_values:
        term = _term_from_design("organism", value, design_terms, term_type="Taxon")
        if term:
            species.append(term)
    return _single_or_list(species)


def _infectious_agents(properties, design_terms, health_conditions):
    agents = []
    for property_name in ("infect", "infectious agent", "pathogen"):
        for value in _contrast_property_values(properties, property_name):
            if value.casefold() in INFECTIOUS_AGENT_SKIP_VALUES:
                continue
            term = _term_from_design(property_name, value, design_terms, term_type="Taxon")
            if term:
                agents.append(term)

    conditions = health_conditions if isinstance(health_conditions, list) else [health_conditions]
    for condition in conditions:
        name = condition.get("name") if isinstance(condition, dict) else condition
        if clean_value(name) and "covid" in name.casefold():
            agents.append(dict(COVID_19_INFECTIOUS_AGENT))
    return _single_or_list(agents)


def _sample_additional_properties(properties, design_terms):
    additional_properties = []
    for prop in properties:
        property_name = clean_value(prop.get("propertyName"))
        if not property_name or prop.get("contrastPropertyType") != "SAMPLE":
            continue
        property_key = _property_key(property_name)
        if property_key in SAMPLE_PROPERTIES_WITH_SPECIFIC_FIELDS:
            continue

        values = []
        test_value = clean_value(prop.get("testValue"))
        reference_value = clean_value(prop.get("referenceValue"))
        if test_value:
            values.append(("test condition", test_value))
        if reference_value and reference_value != test_value:
            values.append(("reference condition", reference_value))

        for label, value in values:
            item = {
                "@type": "PropertyValue",
                "description": "GXA SAMPLE contrast property",
                "name": property_name if len(values) == 1 else f"{property_name} ({label})",
                "propertyID": f"gxa:{_slug(property_key or property_name)}",
                "value": value,
            }
            if term := _term_from_design(property_name, value, design_terms):
                if identifier := term.get("identifier"):
                    item["identifier"] = identifier
                if url := term.get("url"):
                    item["url"] = url
            additional_properties.append(item)

    return _single_or_list(additional_properties)


def _sample(properties, design_terms):
    sample = {"@type": "Sample"}
    organism_parts = _contrast_property_values(properties, "organism part")
    if organism_parts:
        terms = [_term_from_design("organism part", value, design_terms) for value in organism_parts]
        value = _single_or_list(terms)
        if value:
            sample["anatomicalStructure"] = value

    cell_types = _contrast_property_values(properties, "cell type")
    if cell_types:
        terms = [_term_from_design("cell type", value, design_terms) for value in cell_types]
        value = _single_or_list(terms)
        if value:
            sample["cellType"] = value

    developmental_stages = _contrast_property_values(properties, "developmental stage")
    if developmental_stages:
        terms = [_term_from_design("developmental stage", value, design_terms) for value in developmental_stages]
        value = _single_or_list(terms)
        if value:
            sample["developmentalStage"] = value

    sex_values = [value for value in _contrast_property_values(properties, "sex") if "," not in value]
    if sex_values:
        sample["sex"] = _single_or_list(sex_values)

    specimen_states = _contrast_property_values(properties, "specimen with known storage state")
    if specimen_states:
        sample["sampleState"] = _single_or_list(specimen_states)

    if additional_properties := _sample_additional_properties(properties, design_terms):
        sample["additionalProperty"] = additional_properties

    return sample if sample.keys() - {"@type"} else None


def _measurement_technique(summary, raw_type):
    techniques = []
    for technique in summary.get("technologyType") or []:
        if term := _defined_term(technique):
            techniques.append(term)
    if not techniques and raw_type:
        techniques.append(_defined_term(raw_type.replace("_", " ").title()))
    return _single_or_list(techniques)


def _analytical_method(raw_type):
    raw_type = raw_type or ""
    if "RNASEQ" in raw_type:
        return {
            "@type": "ComputationalTool",
            "name": "iRAP Pipeline",
            "url": "https://github.com/nunofonseca/irap",
        }
    if "MICROARRAY" in raw_type:
        return {
            "@type": "ComputationalTool",
            "name": "Expression Atlas microarray differential expression analysis",
        }
    if "PROTEOMICS" in raw_type:
        return {
            "@type": "ComputationalTool",
            "name": "Expression Atlas proteomics differential expression analysis",
        }
    return None


def _build_doc(result, summary, experiment_payload, contrast, raw_type, design_terms):
    accession = summary["experimentAccession"]
    experiment = experiment_payload.get("experiment") or {}
    experiment_species = clean_value(experiment.get("species")) or clean_value(summary.get("species"))
    description = clean_value(experiment.get("description")) or clean_value(summary.get("experimentDescription"))
    contrast_id = clean_value(contrast.get("id")) or _slug(result["contrast_name"])
    properties = _contrast_properties(contrast)
    direction = _direction(result["fold_change"])
    gene = _gene_term(result["gene_id"], result["gene_name"])
    subject_of = _dataset_subject(accession, description)
    test_condition = _factor_condition_text(properties, "test")
    reference_condition = _factor_condition_text(properties, "reference")
    measurement_qualifier = (
        clean_value((contrast.get("contrastSummary") or {}).get("contrastDescription")) or result["contrast_name"]
    )

    doc = {
        "@context": "http://schema.org/",
        "@type": "Inference",
        "_id": f"gxa_{_slug(accession)}_{_slug(contrast_id)}_{_slug(result['gene_id'])}",
        "identifier": f"{accession}:{contrast_id}:{result['gene_id']}",
        "name": f"{gene['name']} is {direction['label'].casefold()} in {measurement_qualifier} ({accession})",
        "url": _result_url(accession, result["gene_id"], result["gene_name"]),
        "includedInDataCatalog": {**CATALOG, "archivedAt": EXPERIMENT_RESULTS_URL_TEMPLATE.format(accession=accession)},
        "subjectOf": subject_of,
        "observationType": dict(OBSERVATION_TYPE),
        "observationAbout": gene,
        "measuredProperty": dict(MEASURED_PROPERTY),
        "measurementQualifier": measurement_qualifier,
        "measurementDenominator": _statistical_variable(
            result["gene_name"], result["gene_id"], reference_condition, experiment_species
        ),
        "variableMeasured": _statistical_variable(
            result["gene_name"], result["gene_id"], test_condition, experiment_species
        ),
        "value": result["fold_change"],
        "unitText": "Log2 fold change",
        "unitCode": FOLD_CHANGE_UNIT_CODE,
        "semanticMapping": _semantic_mapping(gene, direction, subject_of, measurement_qualifier),
        "topicCategory": [dict(term) for term in TOPIC_CATEGORY],
        "keywords": _unique([raw_type, "differential gene expression", direction["label"]]),
    }

    if result["p_value"] is not None:
        doc["marginOfError"] = {
            "@type": "QuantitativeValue",
            "name": "adjusted p-value",
            "value": result["p_value"],
        }
        doc["description"] = (
            f"{gene['name']} has log2 fold change {result['fold_change']} "
            f"with adjusted p-value {result['p_value']} for {measurement_qualifier}."
        )

    if observation_date := _parse_gxa_date(summary.get("lastUpdate")):
        doc["date"] = observation_date
        doc["observationDate"] = observation_date
        doc["dateModified"] = observation_date
    else:
        doc["date"] = datetime.date.today().isoformat()

    if date_created := _parse_gxa_date(summary.get("loadDate")):
        doc["dateCreated"] = date_created

    if health_conditions := _health_conditions(properties, design_terms):
        doc["healthCondition"] = health_conditions

    if infectious_agents := _infectious_agents(properties, design_terms, doc.get("healthCondition")):
        doc["infectiousAgent"] = infectious_agents

    if species := _species_terms(properties, design_terms, experiment_species):
        doc["species"] = species

    if sample := _sample(properties, design_terms):
        doc["sample"] = sample

    if primary_test_condition := _primary_condition(properties, "test"):
        doc["variableMeasured"]["value"] = primary_test_condition
    if primary_reference_condition := _primary_condition(properties, "reference"):
        doc["measurementDenominator"]["value"] = primary_reference_condition

    if technique := _measurement_technique(summary, raw_type):
        doc["measurementTechnique"] = technique
        doc["measurementMethod"] = technique

    if method := _analytical_method(raw_type):
        doc["analyticalMethod"] = method

    return doc


def parse(experiment_accessions: Iterable[str] | None = None):
    session = requests.Session()
    yielded_records = 0
    summaries = (
        ({"experimentAccession": accession} for accession in experiment_accessions)
        if experiment_accessions is not None
        else _iter_experiment_summaries(session)
    )

    for summary in summaries:
        accession = clean_value(summary.get("experimentAccession"))
        if not accession:
            continue

        logger.info("Fetching GXA experiment %s", accession)
        experiment_payload = _fetch_experiment(session, accession)
        if not isinstance(experiment_payload, dict):
            continue

        raw_type = _raw_experiment_type(summary, experiment_payload)
        if not raw_type or "DIFFERENTIAL" not in raw_type:
            logger.info("Skipping non-differential GXA experiment %s with type %s", accession, raw_type)
            continue

        contrast_by_name = {
            clean_value(contrast.get("displayName")): contrast
            for contrast in experiment_payload.get("columnHeaders") or []
            if clean_value(contrast.get("displayName"))
        }
        if not contrast_by_name:
            logger.warning("No contrasts found for GXA experiment %s", accession)
            continue

        tsv_text = _download_results(session, accession, raw_type)
        if not tsv_text:
            continue

        design_terms = _fetch_design_terms(session, accession)
        experiment_records = 0
        for result in _iter_tsv_results(tsv_text):
            contrast = contrast_by_name.get(result["contrast_name"])
            if not contrast:
                logger.warning("No contrast metadata for %s in GXA experiment %s", result["contrast_name"], accession)
                continue
            yield _build_doc(result, summary, experiment_payload, contrast, raw_type, design_terms)
            yielded_records += 1
            experiment_records += 1
            if MAX_RECORDS and yielded_records >= MAX_RECORDS:
                logger.info("Reached GXA_MAX_RECORDS=%s", MAX_RECORDS)
                return

        logger.info("Yielded %s inference records for GXA experiment %s", experiment_records, accession)
        time.sleep(SLEEP)
