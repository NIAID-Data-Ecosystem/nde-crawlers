import copy
import datetime
import json
import logging
import re

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

# All AIRR endpoints to crawl
AIRR_ENDPOINTS = [
    "https://vdjserver.org/airr/v1/repertoire",
    "https://roche-airr.ireceptor.org/airr/v1/repertoire",
    "https://ipa1.ireceptor.org/airr/v1/repertoire",
    "https://ipa4.ireceptor.org/airr/v1/repertoire",
    "https://ipa5.ireceptor.org/airr/v1/repertoire",
    "https://covid19-3.ireceptor.org/airr/v1/repertoire",
    "https://covid19-1.ireceptor.org/airr/v1/repertoire",
    "https://t1d-1.ireceptor.org/airr/v1/repertoire",
    "https://ipa3.ireceptor.org/airr/v1/repertoire",
    "https://ipa6.ireceptor.org/airr/v1/repertoire",
    "https://ipa2.ireceptor.org/airr/v1/repertoire",
    "https://covid19-2.ireceptor.org/airr/v1/repertoire",
    "https://covid19-4.ireceptor.org/airr/v1/repertoire",
    "https://scireptor.dkfz.de/airr/v1/repertoire",
    "https://agschwab.uni-muenster.de/airr/v1/repertoire",
]


def _ensure_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _add_unique_dict(target_list, seen_keys, data):
    cleaned = {k: v for k, v in data.items() if v not in (None, "", [], {}, set())}
    if not cleaned:
        return
    key = tuple(sorted((k, tuple(v) if isinstance(v, list) else v) for k, v in cleaned.items()))
    if key in seen_keys:
        return
    seen_keys.add(key)
    target_list.append(cleaned)


def _add_unique_value(target_list, seen_values, value):
    if value in (None, ""):
        return
    norm = value.strip() if isinstance(value, str) else value
    if norm in (None, ""):
        return
    if norm in seen_values:
        return
    seen_values.add(norm)
    target_list.append(norm)


def _iter_string_values(value):
    if value in (None, "", [], {}, set()):
        return
    if isinstance(value, list):
        for item in value:
            yield from _iter_string_values(item)
        return
    if isinstance(value, dict):
        for key in ("label", "name", "value", "id"):
            if key in value and value[key] not in (None, ""):
                yield from _iter_string_values(value[key])
                return
        return
    if isinstance(value, str):
        parts = [part.strip() for part in re.split(r"[;,]", value) if part.strip()]
        if parts:
            for part in parts:
                yield part
            return
        yield value.strip()
        return
    yield value


def _sanitize_identifier(value, fallback="entry"):
    """Normalize identifiers so they are filesystem/URL friendly."""
    if value in (None, ""):
        return fallback
    if not isinstance(value, str):
        value = str(value)
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value)
    cleaned = cleaned.strip("_")
    return cleaned or fallback


def retrieve_all_study_ids():
    logger.info(f"Retrieving all study IDs from {len(AIRR_ENDPOINTS)} AIRR endpoints")
    all_ids = set()

    for endpoint in AIRR_ENDPOINTS:
        logger.info(f"Crawling endpoint: {endpoint}")
        endpoint_ids = retrieve_study_ids_from_endpoint(endpoint)
        all_ids.update(endpoint_ids)
        logger.info(f"Found {len(endpoint_ids)} study IDs from {endpoint}")

    logger.info(f"Retrieved {len(all_ids)} total unique study IDs")
    return all_ids


def retrieve_study_ids_from_endpoint(endpoint_url):
    """Retrieve study IDs from a single AIRR endpoint"""
    endpoint_ids = set()
    size, page = 100, 0

    while True:
        q = {"from": page * size, "size": size}
        try:
            resp = requests.post(endpoint_url, json=q, timeout=30)
            resp.raise_for_status()
        except requests.exceptions.RequestException as e:
            logger.error(f"Error paging study IDs from {endpoint_url}: {e}")
            break

        try:
            hits = resp.json().get("Repertoire", [])
        except ValueError as e:
            logger.error(f"Invalid JSON response from {endpoint_url}: {e}")
            break

        if not hits:
            break

        for rec in hits:
            sid = rec.get("study", {}).get("study_id")
            if sid:
                endpoint_ids.add(sid)

        if len(hits) < size:
            break
        page += 1

    return endpoint_ids


def fetch_study_records(endpoint_url, study_id):
    records = []
    size, page = 100, 0

    while True:
        q = {
            "filters": {"op": "=", "content": {"field": "study.study_id", "value": study_id}},
            "from": page * size,
            "size": size,
            "include_fields": "airr-core",
        }
        try:
            resp = requests.post(endpoint_url, json=q, timeout=30)
            resp.raise_for_status()
        except requests.exceptions.RequestException as e:
            logger.warning(f"Failed to fetch study {study_id} page {page} from {endpoint_url}: {e}")
            break

        try:
            hits = resp.json().get("Repertoire", [])
        except ValueError as e:
            logger.error(f"Invalid JSON for study {study_id} from {endpoint_url}: {e}")
            break

        if not hits:
            break

        records.extend(hits)
        if len(hits) < size:
            break
        page += 1

    return records


def merge_samples(repertoires):
    samples = []
    seen_keys = set()

    for rep in repertoires:
        rep_samples = rep.get("sample") or []
        if isinstance(rep_samples, dict):
            rep_samples = [rep_samples]
        for sample in rep_samples:
            if not isinstance(sample, dict):
                continue
            key = sample.get("sample_id") or sample.get("sequencing_run_id") or sample.get("vdjserver_uuid")
            if not key:
                key = json.dumps(sample, sort_keys=True, default=str)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            samples.append(copy.deepcopy(sample))

    return samples


def merge_subjects(repertoires):
    subjects = []
    seen_keys = set()

    for rep in repertoires:
        rep_subjects = rep.get("subject")
        if not rep_subjects:
            continue
        rep_subjects = _ensure_list(rep_subjects)
        for subject in rep_subjects:
            if not isinstance(subject, dict):
                continue
            key = subject.get("subject_id") or subject.get("vdjserver_uuid")
            if not key:
                key = json.dumps(subject, sort_keys=True, default=str)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            subjects.append(copy.deepcopy(subject))

    return subjects


def merge_data_processing(repertoires):
    data_processing = []
    seen_keys = set()

    for rep in repertoires:
        rep_dp = _ensure_list(rep.get("data_processing"))
        for entry in rep_dp:
            if not isinstance(entry, dict):
                continue
            key = json.dumps(entry, sort_keys=True, default=str)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            data_processing.append(copy.deepcopy(entry))

    return data_processing


def aggregate_repertoire_records(repertoires):
    if not repertoires:
        return {"core": {}, "samples": [], "subjects": []}

    core = copy.deepcopy(repertoires[0])
    samples = merge_samples(repertoires)
    subjects = merge_subjects(repertoires)
    data_processing = merge_data_processing(repertoires)

    if samples:
        core["sample"] = samples
    if subjects:
        core["subject"] = subjects[0] if len(subjects) == 1 else subjects
    if data_processing:
        core["data_processing"] = data_processing

    return {
        "core": core,
        "samples": samples,
        "subjects": subjects,
    }


def retrieve_study_metadata():
    ids = retrieve_all_study_ids()
    meta = {}

    for sid in ids:
        logger.info(f"Fetching core record for study {sid}")

        # Try each endpoint until we find the study
        study_data = None
        for endpoint in AIRR_ENDPOINTS:
            records = fetch_study_records(endpoint, sid)
            if not records:
                continue
            study_data = aggregate_repertoire_records(records)
            study_data["repertoires"] = [copy.deepcopy(rep) for rep in records]
            study_data["endpoint"] = endpoint
            study_data["disease_facet"] = get_disease_facet_data(sid, endpoint)
            study_data["species_facet"] = get_species_facet_data(sid, endpoint)
            break

        if study_data:
            meta[sid] = study_data
        else:
            logger.error(f"Could not find study {sid} in any endpoint")

    logger.info(f"Metadata collected for {len(meta)} studies")
    return meta


def get_disease_facet_data(sid, endpoint):
    return _facet(sid, "subject.diagnosis.disease_diagnosis", endpoint)


def get_species_facet_data(sid, endpoint):
    return _facet(sid, "subject.species", endpoint)


def _facet(sid, field, endpoint):
    q = {"filters": {"op": "=", "content": {"field": "study.study_id", "value": sid}}, "facets": field}
    try:
        r = requests.post(endpoint, json=q, timeout=30)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"Facet {field} failed for {sid} on {endpoint}: {e}")
        return {}


def _build_sample_collection_record(
    samples,
    dataset_record,
    dataset_identifier,
    dataset_url,
    dataset_name,
    repertoire_id=None,
    subjects=None,
):
    if not samples:
        return None

    dataset_identifier = dataset_identifier or _sanitize_identifier(dataset_record.get("_id"), fallback="dataset")
    dataset_record_id = dataset_record.get("_id") or _sanitize_identifier(dataset_identifier, fallback="dataset")

    if repertoire_id:
        rep_slug = _sanitize_identifier(repertoire_id, fallback="repertoire")
        record_id = f"{dataset_record_id}_samples_{rep_slug}"
        identifier = f"{dataset_identifier}_sample_{rep_slug}"
        display_name_source = dataset_name or dataset_identifier or "Sample"
        display_name = f"{display_name_source} ({repertoire_id})"
    else:
        record_id = f"{dataset_record_id}_samples"
        identifier = f"{dataset_identifier}_sample"
        display_name = dataset_name or dataset_identifier or "Sample"

    sample_record = {
        "_id": record_id,
        "@type": "Sample",
        "identifier": identifier,
        "name": display_name,
        "url": dataset_url,
        "isPartOf": {
            "@type": "Dataset",
            "identifier": dataset_identifier,
            "url": dataset_url,
        },
    }

    if dataset_name:
        sample_record["isPartOf"]["name"] = dataset_name

    if dataset_record.get("includedInDataCatalog"):
        sample_record["includedInDataCatalog"] = copy.deepcopy(dataset_record["includedInDataCatalog"])
    if dataset_record.get("conditionsOfAccess"):
        sample_record["conditionsOfAccess"] = dataset_record["conditionsOfAccess"]
    if dataset_record.get("isAccessibleForFree") is not None:
        sample_record["isAccessibleForFree"] = dataset_record["isAccessibleForFree"]
    if dataset_record.get("variableMeasured"):
        sample_record["variableMeasured"] = copy.deepcopy(dataset_record["variableMeasured"])
    if dataset_record.get("license"):
        sample_record["license"] = dataset_record["license"]
    if dataset_record.get("usageInfo"):
        sample_record["usageInfo"] = copy.deepcopy(dataset_record["usageInfo"])
    if dataset_record.get("measurementTechnique"):
        sample_record["measurementTechnique"] = copy.deepcopy(dataset_record["measurementTechnique"])

    sample_ids = set()
    anatomical_structures = []
    anatomical_seen = set()
    collectors = []
    collectors_seen = set()
    collection_methods = []
    collection_methods_seen = set()
    sample_processes = []
    sample_processes_seen = set()
    sample_states = []
    sample_states_seen = set()
    sample_types = []
    sample_types_seen = set()
    measurement_techniques = []
    measurement_techniques_seen = set()
    phenotypes = []
    phenotypes_seen = set()
    species_entries = []
    species_seen = set()
    cell_types = []
    cell_types_seen = set()
    temporal_coverages = []
    temporal_seen = set()
    health_conditions = []
    health_conditions_seen = set()
    distribution_entries = []
    distribution_seen = set()
    instruments = []
    instruments_seen = set()
    dates_processed = []
    dates_seen = set()
    alternate_identifiers = []
    alternate_ids_seen = set()
    sample_quantities = []
    quantity_seen = set()
    extra_variables = []
    extra_variables_seen = set()
    developmental_stages = []
    developmental_seen = set()
    locations_of_origin = []
    locations_seen = set()
    spatial_coverages = []
    spatial_seen = set()
    associated_genotypes = []
    associated_genotypes_seen = set()
    subject_sexes = []
    subject_sexes_seen = set()

    def add_anatomical(identifier=None, name=None):
        _add_unique_dict(anatomical_structures, anatomical_seen, {"identifier": identifier, "name": name})

    def add_collector(name):
        if name:
            _add_unique_dict(collectors, collectors_seen, {"name": name})

    def add_collection_method(value):
        if value:
            _add_unique_value(collection_methods, collection_methods_seen, value)

    def add_sample_process(value):
        if value:
            _add_unique_value(sample_processes, sample_processes_seen, value)

    def add_sample_state(value):
        if value:
            _add_unique_value(sample_states, sample_states_seen, value)

    def add_sample_type(value):
        if value:
            _add_unique_dict(sample_types, sample_types_seen, {"name": value})

    def add_measurement(name=None, description=None):
        entry = {"name": name, "description": description}
        _add_unique_dict(measurement_techniques, measurement_techniques_seen, entry)

    def add_phenotype(value):
        if value:
            _add_unique_dict(phenotypes, phenotypes_seen, {"name": value})

    def add_species(identifier=None, name=None):
        _add_unique_dict(species_entries, species_seen, {"identifier": identifier, "name": name})

    def add_cell_type(identifier=None, name=None):
        _add_unique_dict(cell_types, cell_types_seen, {"identifier": identifier, "name": name})

    def add_temporal(name=None, duration=None):
        _add_unique_dict(temporal_coverages, temporal_seen, {"name": name, "duration": duration})

    def add_health_condition(name, identifier=None):
        entry = {"name": name, "identifier": identifier}
        _add_unique_dict(health_conditions, health_conditions_seen, entry)

    def add_distribution(entry):
        _add_unique_dict(distribution_entries, distribution_seen, entry)

    def add_instrument(value):
        if not value:
            return
        if isinstance(value, str):
            value = value.strip()
            if not value:
                return
        entry = {"name": value}
        _add_unique_dict(instruments, instruments_seen, entry)

    def normalize_date(value):
        # Normalize various external date formats to YYYY-MM-DD for consistency.
        if isinstance(value, datetime.datetime):
            return value.date().isoformat()
        if isinstance(value, datetime.date):
            return value.isoformat()
        if not isinstance(value, str):
            return None
        text = value.strip()
        if not text:
            return None

        candidates = [text]
        if " " in text:
            candidates.append(text.split(" ", 1)[0])

        formats = [
            "%Y-%m-%d",
            "%Y/%m/%d",
            "%m/%d/%Y",
            "%m/%d/%y",
            "%m-%d-%Y",
            "%m-%d-%y",
            "%Y-%m-%dT%H:%M:%S",
            "%Y-%m-%dT%H:%M:%S.%f",
            "%Y-%m-%d %H:%M:%S",
            "%m/%d/%Y %H:%M",
            "%m/%d/%y %H:%M",
        ]


        for candidate in candidates:
            trimmed = candidate.strip()
            if not trimmed:
                continue
            try:
                dt = datetime.datetime.fromisoformat(trimmed.replace("Z", ""))
                return dt.date().isoformat()
            except ValueError:
                pass
            for fmt in formats:
                try:
                    dt = datetime.datetime.strptime(trimmed, fmt)
                    return dt.date().isoformat()
                except ValueError:
                    continue
        return None

    def add_date(value):
        normalized = normalize_date(value)
        if normalized:
            _add_unique_value(dates_processed, dates_seen, normalized)
        elif value:
            _add_unique_value(dates_processed, dates_seen, value)

    def add_alt_identifier(value):
        if value:
            _add_unique_value(alternate_identifiers, alternate_ids_seen, value)

    def add_quantity(value, unit=None, label=None):
        if value in (None, ""):
            return
        # Ensure quantity values stay aligned with ES integer mapping.
        def _coerce_int(raw):
            if isinstance(raw, bool):
                return None
            if isinstance(raw, int):
                return raw
            if isinstance(raw, float):
                return int(raw) if raw.is_integer() else None
            text = str(raw).strip()
            if not text:
                return None
            if text.isdigit() or (text.startswith("-") and text[1:].isdigit()):
                return int(text)
            try:
                float_value = float(text)
            except (TypeError, ValueError):
                return None
            return int(float_value) if float_value.is_integer() else None

        coerced_value = _coerce_int(value)
        if coerced_value is None:
            logger.warning(
                "Skipping non-integer sampleQuantity value %r for label %r",
                value,
                label or "sampleQuantity",
            )
            return

        entry = {
            "value": coerced_value,
            "unitText": unit,
            "name": label,
        }
        _add_unique_dict(sample_quantities, quantity_seen, entry)

    def add_variable_measured(name=None, description=None):
        entry = {"name": name, "description": description}
        _add_unique_dict(extra_variables, extra_variables_seen, entry)

    def add_variable_measured_value(name, raw_value):
        if raw_value in (None, "", [], {}, set()):
            return
        for value in _iter_string_values(raw_value):
            if value in (None, ""):
                continue
            add_variable_measured(name=name, description=str(value))

    def add_developmental_stage(name=None, min_value=None, max_value=None, unit_code=None, unit_text=None):
        entry = {
            "name": name,
            "minValue": min_value,
            "maxValue": max_value,
            "unitCode": unit_code,
            "unitText": unit_text,
        }
        _add_unique_dict(developmental_stages, developmental_seen, entry)

    def add_location_of_origin(value):
        if value in (None, ""):
            return

        if isinstance(value, dict):
            entry = {k: v for k, v in value.items() if v not in (None, "", [], {}, set())}
            if not entry:
                return
            if not entry.get("name"):
                for candidate_key in ("label", "value", "identifier"):
                    candidate = entry.get(candidate_key)
                    if isinstance(candidate, str) and candidate.strip():
                        entry["name"] = candidate.strip()
                        break
            name_value = entry.get("name")
            if isinstance(name_value, str):
                entry["name"] = name_value.strip()
            if not entry.get("name"):
                return
        else:
            text = value.strip() if isinstance(value, str) else str(value).strip()
            if not text:
                return
            entry = {"name": text}

        key = json.dumps(entry, sort_keys=True, default=str)
        if key in locations_seen:
            return
        locations_seen.add(key)
        locations_of_origin.append(entry)

    def add_spatial_coverage(name=None):
        if name in (None, ""):
            return
        _add_unique_dict(spatial_coverages, spatial_seen, {"name": name})

    def add_associated_genotype(value):
        if value in (None, "", [], {}, set()):
            return
        for entry in _iter_string_values(value):
            if entry in (None, ""):
                continue
            _add_unique_value(associated_genotypes, associated_genotypes_seen, entry.strip() if isinstance(entry, str) else entry)

    def add_subject_sex(value):
        if value in (None, ""):
            return
        _add_unique_value(subject_sexes, subject_sexes_seen, value.strip() if isinstance(value, str) else value)

    if repertoire_id:
        add_alt_identifier(str(repertoire_id))

    for subject in subjects or []:
        if not isinstance(subject, dict):
            continue

        age = subject.get("age")
        age_min = subject.get("age_min")
        age_max = subject.get("age_max")
        age_unit = subject.get("age_unit")
        unit_code = None
        unit_text = None
        if isinstance(age_unit, dict):
            unit_code = age_unit.get("id")
            unit_text = age_unit.get("label")
        elif isinstance(age_unit, str):
            unit_text = age_unit.strip()
        if any(value not in (None, "") for value in (age, age_min, age_max, unit_code, unit_text)):
            add_developmental_stage(
                name=age,
                min_value=age_min,
                max_value=age_max,
                unit_code=unit_code,
                unit_text=unit_text,
            )

        add_subject_sex(subject.get("sex"))

        add_variable_measured_value("ethnicity", subject.get("ethnicity"))
        add_variable_measured_value("race", subject.get("race"))

        strain_name = subject.get("strain_name")
        if isinstance(strain_name, str) and strain_name.strip():
            add_associated_genotype(strain_name.strip())

        for species_block in _ensure_list(subject.get("species")):
            if isinstance(species_block, dict):
                add_species(species_block.get("id"), species_block.get("label"))
            elif isinstance(species_block, str):
                add_species(name=species_block.strip())

        genotype = subject.get("genotype") or {}
        mhc_set = genotype.get("mhc_genotype_set")
        if isinstance(mhc_set, dict):
            add_associated_genotype(mhc_set.get("mhc_genotype_set_id"))
            for entry in _ensure_list(mhc_set.get("mhc_genotype_list")):
                if isinstance(entry, dict):
                    add_associated_genotype(entry.get("mhc_genotype_id"))
            for entry in _ensure_list(mhc_set.get("mhc_genotype_class_list")):
                if isinstance(entry, dict):
                    add_associated_genotype(entry.get("mhc_genotype_id"))

        for diagnosis in _ensure_list(subject.get("diagnosis")):
            if not isinstance(diagnosis, dict):
                continue
            disease_block = diagnosis.get("disease_diagnosis") or {}
            if isinstance(disease_block, dict):
                name = disease_block.get("label") or disease_block.get("name")
                identifier = disease_block.get("id")
                if name or identifier:
                    add_health_condition(name or identifier, identifier)
            else:
                add_health_condition(str(disease_block), None)

            add_variable_measured_value("disease severity", diagnosis.get("disease_stage") or diagnosis.get("diseaseStage"))
            add_variable_measured_value("intervention", diagnosis.get("intervention"))
            add_variable_measured_value("prior therapies", diagnosis.get("prior_therapies") or diagnosis.get("priorTherapies"))

            immune_block = diagnosis.get("ImmuneCODE")
            if isinstance(immune_block, dict):
                collection_region = immune_block.get("collection_region")
                for region in _iter_string_values(collection_region):
                    if region in (None, ""):
                        continue
                    region_text = region.strip() if isinstance(region, str) else region
                    add_location_of_origin(region_text)
                    add_spatial_coverage(region_text)

    for sample in samples:
        sample_id = sample.get("sample_id")
        if sample_id not in (None, ""):
            sample_ids.add(str(sample_id))

        for value in _ensure_list(sample.get("anatomic_site")):
            if isinstance(value, str):
                add_anatomical(name=value.strip())

        for tissue in _ensure_list(sample.get("tissue")):
            if isinstance(tissue, dict):
                add_anatomical(identifier=tissue.get("id"), name=(tissue.get("label") or tissue.get("value")))
            elif isinstance(tissue, str):
                add_anatomical(name=tissue.strip())

        for provider in _ensure_list(sample.get("biomaterial_provider")):
            if isinstance(provider, str):
                add_collector(provider.strip())

        for facility in _ensure_list(sample.get("sequencing_facility")):
            if isinstance(facility, str):
                add_collector(facility.strip())

        for method in _ensure_list(sample.get("cell_isolation")):
            if isinstance(method, str):
                add_collection_method(method.strip())

        for process in _ensure_list(sample.get("cell_processing_protocol")):
            if isinstance(process, str):
                add_sample_process(process.strip())

        for process in _ensure_list(sample.get("cell_storage")):
            if isinstance(process, str):
                add_sample_process(process.strip())

        for process in _ensure_list(sample.get("tissue_processing")):
            if isinstance(process, str):
                add_sample_process(process.strip())

        if sample.get("single_cell") not in (None, ""):
            add_sample_process(f"Single cell: {sample['single_cell']}")

        for linkage in _ensure_list(sample.get("physical_linkage")):
            if isinstance(linkage, str):
                add_sample_process(linkage.strip())

        for state in _ensure_list(sample.get("cell_quality")):
            if isinstance(state, str):
                add_sample_state(state.strip())

        for state in _ensure_list(sample.get("template_quality")):
            if isinstance(state, str):
                add_sample_state(state.strip())

        for state in _ensure_list(sample.get("complete_sequences")):
            if isinstance(state, str):
                add_sample_state(state.strip())

        for stype in _ensure_list(sample.get("sample_type")):
            if isinstance(stype, str):
                add_sample_type(stype.strip())

        for template in _ensure_list(sample.get("template_class")):
            if isinstance(template, str):
                add_sample_type(template.strip())

        for phenotype in _ensure_list(sample.get("cell_phenotype")):
            if isinstance(phenotype, str):
                add_phenotype(phenotype.strip())

        for species in _ensure_list(sample.get("cell_species")):
            if isinstance(species, dict):
                add_species(species.get("id"), species.get("label"))
            elif isinstance(species, str):
                add_species(name=species.strip())

        for subset in _ensure_list(sample.get("cell_subset")):
            if isinstance(subset, dict):
                add_cell_type(subset.get("id"), subset.get("label"))
            elif isinstance(subset, str):
                add_cell_type(name=subset.strip())

        reference = sample.get("collection_time_point_reference")
        relative = sample.get("collection_time_point_relative")
        unit = sample.get("collection_time_point_relative_unit")
        duration = None
        if relative not in (None, ""):
            rel_str = str(relative).strip()
            unit_label = None
            if isinstance(unit, dict):
                unit_label = unit.get("label") or unit.get("id")
            elif isinstance(unit, str):
                unit_label = unit.strip()
            duration = f"{rel_str} {unit_label}".strip() if unit_label else rel_str
        if reference or duration:
            ref_value = reference.strip() if isinstance(reference, str) else reference
            add_temporal(ref_value, duration)

        for disease in _ensure_list(sample.get("disease_state_sample")):
            if isinstance(disease, str):
                add_health_condition(disease.strip())

        sequencing_files = sample.get("sequencing_files")
        for file_entry in _ensure_list(sequencing_files):
            if not isinstance(file_entry, dict):
                continue
            main_entry = {}
            content_url = file_entry.get("filename") or file_entry.get("file_url") or file_entry.get("content_url")
            if content_url:
                main_entry["contentUrl"] = content_url
            if file_entry.get("file_type"):
                main_entry["encodingFormat"] = file_entry["file_type"]
            if file_entry.get("sequencing_data_id"):
                data_id = str(file_entry["sequencing_data_id"])
                main_entry["identifier"] = data_id
                add_alt_identifier(data_id)
            desc_parts = []
            if file_entry.get("read_direction"):
                desc_parts.append(f"direction: {file_entry['read_direction']}")
            if file_entry.get("read_length"):
                desc_parts.append(f"read_length: {file_entry['read_length']}")
            if file_entry.get("paired_read_length"):
                desc_parts.append(f"paired_read_length: {file_entry['paired_read_length']}")
            if file_entry.get("paired_read_direction"):
                desc_parts.append(f"paired_direction: {file_entry['paired_read_direction']}")
            if desc_parts:
                main_entry["description"] = "; ".join(desc_parts)
            add_distribution(main_entry)

            paired_filename = file_entry.get("paired_filename")
            if paired_filename:
                paired_entry = {"contentUrl": paired_filename}
                if file_entry.get("file_type"):
                    paired_entry["encodingFormat"] = file_entry["file_type"]
                if file_entry.get("sequencing_data_id"):
                    paired_entry["identifier"] = f"{file_entry['sequencing_data_id']}-paired"
                paired_desc = []
                if file_entry.get("paired_read_direction"):
                    paired_desc.append(f"direction: {file_entry['paired_read_direction']}")
                if file_entry.get("paired_read_length"):
                    paired_desc.append(f"read_length: {file_entry['paired_read_length']}")
                if paired_desc:
                    paired_entry["description"] = "; ".join(paired_desc)
                add_distribution(paired_entry)

        for kit in _ensure_list(sample.get("library_generation_kit_version")):
            if isinstance(kit, str):
                add_measurement(name=kit.strip())

        for method in _ensure_list(sample.get("library_generation_method")):
            if isinstance(method, str):
                add_measurement(name=method.strip())

        for protocol in _ensure_list(sample.get("library_generation_protocol")):
            if isinstance(protocol, str):
                add_measurement(description=protocol.strip())

        for pcr in _ensure_list(sample.get("pcr_target")):
            if isinstance(pcr, dict):
                locus = pcr.get("pcr_target_locus")
                desc_parts = []
                fwd = pcr.get("forward_pcr_primer_target_location")
                if fwd:
                    desc_parts.append(f"forward primer: {fwd}")
                rev = pcr.get("reverse_pcr_primer_target_location")
                if rev:
                    desc_parts.append(f"reverse primer: {rev}")
                description = "; ".join(desc_parts) if desc_parts else None
                add_variable_measured(name=locus, description=description)

        cell_number = sample.get("cell_number")
        if cell_number not in (None, ""):
            add_quantity(cell_number, label="cell_number")

        cells_per_reaction = sample.get("cells_per_reaction")
        if cells_per_reaction not in (None, ""):
            add_quantity(cells_per_reaction, label="cells_per_reaction")

        template_amount = sample.get("template_amount")
        template_unit = sample.get("template_amount_unit")
        if template_amount not in (None, ""):
            unit_label = None
            if isinstance(template_unit, dict):
                unit_label = template_unit.get("label") or template_unit.get("id")
            elif isinstance(template_unit, str):
                unit_label = template_unit.strip()
            add_quantity(template_amount, unit=unit_label, label="template_amount")

        for instrument in _ensure_list(sample.get("sequencing_platform")):
            if isinstance(instrument, str):
                add_instrument(instrument.strip())

        for instrument in _ensure_list(sample.get("sequencing_kit")):
            if isinstance(instrument, str):
                add_instrument(instrument.strip())

        run_date = sample.get("sequencing_run_date")
        if run_date:
            add_date(run_date)

        for alt in ("sequencing_run_id", "vdjserver_uuid", "sample_processing_id"):
            alt_value = sample.get(alt)
            if alt_value not in (None, ""):
                add_alt_identifier(str(alt_value))

        vdjserver_block = sample.get("vdjserver")
        if isinstance(vdjserver_block, dict):
            alt_value = vdjserver_block.get("vdjserver_uuid")
            if alt_value not in (None, ""):
                add_alt_identifier(str(alt_value))

    sample_list = sorted(sample_ids)
    sample_record["sampleList"] = sample_list
    sample_record["collectionSize"] = {
        "value": len(sample_list) if sample_list else len(samples)
    }

    if developmental_stages:
        sample_record["developmentalStage"] = (
            developmental_stages[0] if len(developmental_stages) == 1 else developmental_stages
        )
    if locations_of_origin:
        sample_record["locationOfOrigin"] = (
            locations_of_origin[0] if len(locations_of_origin) == 1 else locations_of_origin
        )
    if spatial_coverages:
        sample_record["spatialCoverage"] = spatial_coverages
    if subject_sexes:
        sample_record["sex"] = subject_sexes[0] if len(subject_sexes) == 1 else subject_sexes
    if associated_genotypes:
        sample_record["associatedGenotype"] = associated_genotypes

    if anatomical_structures:
        sample_record["anatomicalStructure"] = anatomical_structures
    if collectors:
        sample_record["collector"] = collectors
    if collection_methods:
        sample_record["collectionMethod"] = collection_methods
    if sample_processes:
        sample_record["sampleProcess"] = sample_processes
    if sample_states:
        sample_record["sampleState"] = sample_states
    if sample_types:
        sample_record["sampleType"] = sample_types
    if measurement_techniques:
        existing_measurements = sample_record.get("measurementTechnique", [])
        if isinstance(existing_measurements, dict):
            existing_measurements = [existing_measurements]
        elif not isinstance(existing_measurements, list):
            existing_measurements = []
        existing_keys = {
            tuple(sorted((k, v) for k, v in ((key, val) for key, val in item.items() if val not in (None, ""))))
            for item in existing_measurements
            if isinstance(item, dict)
        }
        for entry in measurement_techniques:
            key = tuple(sorted((k, v) for k, v in entry.items()))
            if key not in existing_keys:
                existing_measurements.append(entry)
                existing_keys.add(key)
        if existing_measurements:
            sample_record["measurementTechnique"] = existing_measurements
    if phenotypes:
        sample_record["associatedPhenotype"] = phenotypes
    if species_entries:
        sample_record["species"] = species_entries
    if cell_types:
        sample_record["cellType"] = cell_types
    if temporal_coverages:
        sample_record["temporalCoverage"] = temporal_coverages
    if health_conditions:
        sample_record.setdefault("healthCondition", [])
        existing = {
            tuple(sorted((k, v) for k, v in entry.items()))
            for entry in sample_record["healthCondition"]
            if isinstance(entry, dict)
        }
        for entry in health_conditions:
            key = tuple(sorted((k, v) for k, v in entry.items()))
            if key not in existing:
                sample_record["healthCondition"].append(entry)
                existing.add(key)
    if distribution_entries:
        sample_record.setdefault("distribution", [])
        sample_record["distribution"].extend(distribution_entries)
    if instruments:
        sample_record["instrument"] = instruments
    if dates_processed:
        sample_record["dateProcessed"] = dates_processed
    if alternate_identifiers:
        sample_record["alternateIdentifier"] = alternate_identifiers
    if sample_quantities:
        sample_record["sampleQuantity"] = sample_quantities
    if extra_variables:
        existing_variables = sample_record.get("variableMeasured", [])
        if not isinstance(existing_variables, list):
            existing_variables = []
        existing_keys = {
            tuple(sorted((k, v) for k, v in entry.items()))
            for entry in existing_variables
            if isinstance(entry, dict)
        }
        for entry in extra_variables:
            key = tuple(sorted((k, v) for k, v in entry.items()))
            if key not in existing_keys:
                existing_variables.append(entry)
                existing_keys.add(key)
        sample_record["variableMeasured"] = existing_variables

    return sample_record


def build_sample_collection_records(sid, metadata, dataset_record):
    dataset_identifier = dataset_record.get("identifier") or sid
    dataset_url = dataset_record.get("url")
    dataset_name = dataset_record.get("name") or sid

    all_subjects = metadata.get("subjects")
    if all_subjects is None:
        all_subjects = metadata.get("core", {}).get("subject") or []
    if isinstance(all_subjects, dict):
        all_subjects = [all_subjects]
    all_subjects = [s for s in all_subjects if isinstance(s, dict)]

    sample_groups = []
    repertoires = metadata.get("repertoires") or []
    if repertoires:
        for rep in repertoires:
            rep_samples = rep.get("sample") or []
            if isinstance(rep_samples, dict):
                rep_samples = [rep_samples]
            rep_samples = [s for s in rep_samples if isinstance(s, dict)]
            if not rep_samples:
                continue
            rep_subjects = [s for s in _ensure_list(rep.get("subject")) if isinstance(s, dict)]
            if not rep_subjects:
                rep_subjects = all_subjects
            sample_groups.append((rep_samples, rep.get("repertoire_id"), rep_subjects))
    else:
        samples = metadata.get("samples")
        if samples is None:
            samples = metadata.get("core", {}).get("sample") or []
        if isinstance(samples, dict):
            samples = [samples]
        samples = [s for s in samples if isinstance(s, dict)]
        if samples:
            sample_groups.append((samples, None, all_subjects))

    records = []
    for samples, repertoire_id, subjects in sample_groups:
        record = _build_sample_collection_record(
            samples,
            dataset_record,
            dataset_identifier,
            dataset_url,
            dataset_name,
            repertoire_id=repertoire_id,
            subjects=subjects,
        )
        if record:
            records.append(record)

    return records


def aggregate_sample_properties_to_dataset(dataset_record, sample_collections):
    """
    Aggregate sample-level properties to the dataset record for search purposes.
    Properties like species, healthCondition, cellType, etc. are collected from all
    sample collections and added to the dataset if not already present.
    """
    if not sample_collections:
        return

    # Properties to aggregate from samples to dataset
    properties_to_aggregate = {
        "cellType": set(),
        "anatomicalStructure": set(),
        "sex": set(),
        "developmentalStage": set(),
        "associatedGenotype": set(),
        "spatialCoverage": set(),
        "temporalCoverage": set(),
    }

    # Collect unique values from all sample collections
    for sample_coll in sample_collections:
        for prop_name in properties_to_aggregate.keys():
            prop_value = sample_coll.get(prop_name)
            if prop_value:
                # Handle both single values and arrays
                values_list = prop_value if isinstance(prop_value, list) else [prop_value]
                for val in values_list:
                    if isinstance(val, dict):
                        # For dict values, create a hashable key
                        key = json.dumps(val, sort_keys=True, default=str)
                        properties_to_aggregate[prop_name].add(key)
                    elif val not in (None, ""):
                        properties_to_aggregate[prop_name].add(val if isinstance(val, str) else json.dumps(val, sort_keys=True, default=str))

    # Add aggregated properties to dataset if they have values and aren't already more complete
    for prop_name, values_set in properties_to_aggregate.items():
        if not values_set:
            continue

        # Convert back from JSON strings to dicts where needed
        aggregated_values = []
        for val_str in values_set:
            try:
                # Try to parse as JSON for dict values
                parsed = json.loads(val_str) if isinstance(val_str, str) and (val_str.startswith("{") or val_str.startswith("[")) else val_str
                aggregated_values.append(parsed)
            except (json.JSONDecodeError, ValueError):
                aggregated_values.append(val_str)

        # Check if dataset already has this property with sufficient data
        existing_value = dataset_record.get(prop_name)
        if existing_value:
            # If dataset already has values, merge with sample values
            existing_list = existing_value if isinstance(existing_value, list) else [existing_value]
            existing_keys = {json.dumps(v, sort_keys=True, default=str) for v in existing_list}

            for new_val in aggregated_values:
                new_key = json.dumps(new_val, sort_keys=True, default=str)
                if new_key not in existing_keys:
                    existing_list.append(new_val)
                    existing_keys.add(new_key)

            if len(existing_list) > 1:
                dataset_record[prop_name] = existing_list
            elif len(existing_list) == 1:
                dataset_record[prop_name] = existing_list[0]
        else:
            # Add new aggregated values to dataset
            if len(aggregated_values) == 1:
                dataset_record[prop_name] = aggregated_values[0]
            elif len(aggregated_values) > 1:
                dataset_record[prop_name] = aggregated_values


def parse():
    all_meta = retrieve_study_metadata()
    total = len(all_meta)

    for idx, (sid, md) in enumerate(all_meta.items(), start=1):
        logger.info(f"Parsing {idx}/{total}: {sid}")
        core = md["core"]
        out = {}

        samples = md.get("samples")
        if samples is None:
            samples = core.get("sample") or []
        if isinstance(samples, dict):
            samples = [samples]
        samples = [s for s in samples if isinstance(s, dict)]
        logger.info("Study %s parsed %d samples", sid, len(samples))

        subjects = md.get("subjects")
        if subjects is None:
            subject_block = core.get("subject")
            if isinstance(subject_block, dict):
                subjects = [subject_block]
            elif isinstance(subject_block, list):
                subjects = [s for s in subject_block if isinstance(s, dict)]
            else:
                subjects = []
        else:
            subjects = [s for s in subjects if isinstance(s, dict)]

        stud_url = f"https://vdjserver.org/community?study_id={sid}"
        out["_id"] = f"vdj_{sid}".replace(":", "_").replace(" ", "_").replace("/", "_")
        out["url"] = stud_url
        out["identifier"] = sid
        out["@type"] = "Dataset"
        out["includedInDataCatalog"] = {
            "@type": "DataCatalog",
            "name": "VDJServer",
            "url": stud_url,
            "versionDate": datetime.date.today().isoformat(),
            "archivedAt": stud_url,
        }
        if core.get("study", {}).get("study_title"):
            out["name"] = core["study"]["study_title"]

        stype = core.get("study", {}).get("study_type", {})
        if stype.get("id") or stype.get("label"):
            out["measurementTechnique"] = {}
            if stype.get("id"):
                out["measurementTechnique"]["identifier"] = stype["id"]
            if stype.get("label"):
                out["measurementTechnique"]["name"] = stype["label"]

        def flatten_text(x):
            if isinstance(x, list):
                return "".join(x)
            return x or ""

        raw_desc = core.get("study", {}).get("study_description", "")
        raw_inc_exc = core.get("study", {}).get("inclusion_exclusion_criteria", "")

        desc = flatten_text(raw_desc)
        inc_exc = flatten_text(raw_inc_exc)

        out["description"] = " ".join(filter(None, [desc, inc_exc])).strip()
        authors = []
        lab = {}
        if core.get("study", {}).get("lab_name"):
            lab["name"] = core["study"]["lab_name"]
        addr = core.get("study", {}).get("lab_address")
        if addr:
            lab.setdefault("affiliation", {})["name"] = addr
        if lab:
            authors.append(lab)
        sub = core.get("study", {}).get("submitted_by")
        if sub:
            split = [s.strip() for s in re.split(r"[;,]", sub) if s.strip()]
            if split:
                sb = {"name": split[0]}
                if len(split) > 1:
                    sb.setdefault("affiliation", {})["name"] = split[1]
                authors.append(sb)
        coll = core.get("study", {}).get("collected_by", [])
        if isinstance(coll, str):
            coll = [coll]
        if coll:
            for name in coll:
                authors.append({"name": name})
        if authors:
            out["author"] = authors
        grants = core.get("study", {}).get("grants") or []
        if isinstance(grants, str):
            grants = [grants]
        if grants:
            extra = " ".join(grants)
            out["description"] += " " + extra if out.get("description") else extra
        if out.get("description") == "":
            del out["description"]
        pids = core.get("study", {}).get("pub_ids", [])
        if isinstance(pids, str):
            pids = [pids]
        pmids, dois = [], []
        if pids:
            for pid in pids:
                if pid.lower().startswith("pmid:") or re.match(r"^\d+$", pid):
                    pmids.append(pid.replace("PMID:", ""))
                elif "doi" in pid.lower():
                    doi = pid.lower().replace("DOI:", "").strip()
                    dois.append(doi)
        if pmids:
            out["pmids"] = ", ".join(pmids)
        if dois:
            out.setdefault("citation", {})["doi"] = dois if len(dois) > 1 else dois[0]
        if core.get("study", {}).get("adc_publish_date"):
            out["datePublished"] = core["study"]["adc_publish_date"].split("T")[0]
        if core.get("study", {}).get("adc_update_date"):
            out["dateModified"] = core["study"]["adc_update_date"].split("T")[0]
        if core.get("study", {}).get("publisher"):
            out.setdefault("sdPublisher", {})["name"] = core["study"]["publisher"]
        dists = []
        di = core.get("download_info", {})
        if di.get("archive_file"):
            dists.append({"contentUrl": di["archive_file"]})
        if di.get("download_url"):
            dists.append({"contentUrl": di["download_url"]})
        if di.get("file_size"):

            for d in dists:
                d["contentSize"] = di["file_size"]
        if dists:
            out["distribution"] = dists

        data_processings = core.get("data_processing", [])
        if isinstance(data_processings, dict):
            data_processings = [data_processings]

        for dp in data_processings:
            for fn in dp.get("data_processing_files", []):
                out.setdefault("distribution", []).append({"contentUrl": fn})
        grants = core.get("study", {}).get("grants") or []

        if isinstance(grants, str):
            grants = [grants]

        # for grant in grants:
        # out.setdefault("funding", []).append({"description": grant})
        kw_n, kw_id = [], []
        for sample in samples:
            for fld in ("anatomic_site", "disease_state_sample", "template_class"):
                v = sample.get(fld)
                if v:
                    kw_n.append(v)
            for nested in ("cell_subset", "tissue"):
                nested_obj = sample.get(nested) or {}
                label = nested_obj.get("label")
                if label:
                    kw_n.append(label)
                nid = nested_obj.get("id")
                if nid:
                    kw_id.append(nid)

        diagnosis_entries = []
        for subject in subjects:
            for diag in _ensure_list(subject.get("diagnosis")):
                if isinstance(diag, dict):
                    diagnosis_entries.append(diag)
        for diag in diagnosis_entries:
            stage = diag.get("disease_stage") or diag.get("diseaseStage")
            if stage:
                kw_n.append(stage)
                break

        all_keywords = kw_n + kw_id
        unique_keywords = list(dict.fromkeys(all_keywords))
        if unique_keywords:
            out["keywords"] = unique_keywords

        species_seen = set()
        health_conditions_seen = set()

        species_facet = md.get("species_facet") or {}
        disease_facet = md.get("disease_facet") or {}

        for b in species_facet.get("Facet", []):
            logger.debug("Processing species facet: %s", b)
            species_data = b.get("subject.species", {})
            species_id = species_data.get("id")
            species_name = species_data.get("label")

            unique_key = species_name or species_id
            if unique_key and unique_key not in species_seen:
                species_seen.add(unique_key)
                sp = {}
                if species_id:
                    sp["identifier"] = species_id
                if species_name:
                    sp["name"] = species_name
                out.setdefault("species", []).append(sp)

        for sample in samples:
            cs = sample.get("cell_species", {})
            if not isinstance(cs, dict):
                continue
            label = cs.get("label")
            if label and label not in species_seen:
                species_seen.add(label)
                out.setdefault("species", []).append({"name": label})

        for b in disease_facet.get("Facet", []):
            logger.debug("Processing disease facet: %s", b)
            disease_data = b.get("subject.diagnosis.disease_diagnosis") or {}
            if not isinstance(disease_data, dict):
                logger.debug("Skipping facet with unexpected disease data: %s", disease_data)
                continue
            disease_name = disease_data.get("label") or disease_data.get("id")

            if disease_name and disease_name not in health_conditions_seen:
                health_conditions_seen.add(disease_name)
                disease_entry = {"name": disease_name}
                if disease_data.get("id"):
                    disease_entry["identifier"] = disease_data["id"]
                out.setdefault("healthCondition", []).append(disease_entry)

        out["conditionsOfAccess"] = "Open"
        out["isAccessibleForFree"] = True
        out["variableMeasured"] = [
            {
                "identifier": "C20971",
                "name": "V(D)J Recombination",
                "inDefinedTermSet": "NCIT",
                "url": "https://evsexplore.semantics.cancer.gov/evsexplore/concept/ncit/C20971",
            },
            {
                "identifier": "data_2977",
                "name": "Nucleic acid sequence",
                "inDefinedTermSet": "EDAM",
                "url": "https://edamontology.org/data_2977",
            },
        ]

        variable_measured = out["variableMeasured"]
        variable_seen = set()
        for entry in variable_measured:
            if isinstance(entry, dict):
                cleaned = {k: v for k, v in entry.items() if v not in (None, "", [], {}, set())}
                key = tuple(sorted((k, tuple(v) if isinstance(v, list) else v) for k, v in cleaned.items()))
                variable_seen.add(key)

        associated_genotypes = out.setdefault("associatedGenotype", [])
        associated_seen = set()
        cleaned_associated_genotypes = []
        for entry in associated_genotypes:
            candidate = None
            if isinstance(entry, str):
                candidate = entry.strip()
            elif isinstance(entry, dict):
                for key in ("name", "identifier", "value"):
                    val = entry.get(key)
                    if isinstance(val, str) and val.strip():
                        candidate = val.strip()
                        break
                if candidate is None:
                    for val in entry.values():
                        if isinstance(val, str) and val.strip():
                            candidate = val.strip()
                            break
            elif entry is not None:
                candidate = str(entry).strip()
            if not candidate:
                continue
            if candidate in associated_seen:
                continue
            associated_seen.add(candidate)
            cleaned_associated_genotypes.append(candidate)
        associated_genotypes.clear()
        associated_genotypes.extend(cleaned_associated_genotypes)

        sexes = set()

        def add_variable_measurement_entry(name, raw_value):
            if raw_value in (None, "", [], {}, set()):
                return
            for value in _iter_string_values(raw_value):
                if value in (None, ""):
                    continue
                entry = {"name": name, "description": str(value)}
                _add_unique_dict(variable_measured, variable_seen, entry)

        def add_associated_genotype_entry(raw_value):
            if raw_value in (None, "", [], {}, set()):
                return
            for value in _iter_string_values(raw_value):
                if value in (None, ""):
                    continue
                text = str(value)
                _add_unique_value(associated_genotypes, associated_seen, text)

        for subject in subjects:
            sex_value = subject.get("sex")
            if isinstance(sex_value, str) and sex_value.strip():
                sexes.add(sex_value.strip())

            add_variable_measurement_entry("ethnicity", subject.get("ethnicity"))
            add_variable_measurement_entry("race", subject.get("race"))

            strain_name = subject.get("strain_name")
            if isinstance(strain_name, str) and strain_name.strip():
                add_associated_genotype_entry(strain_name.strip())

            for species_block in _ensure_list(subject.get("species")):
                if not isinstance(species_block, dict):
                    continue
                species_id = species_block.get("id")
                species_name = species_block.get("label")
                unique_key = species_name or species_id
                if unique_key and unique_key not in species_seen:
                    species_seen.add(unique_key)
                    species_entry = {}
                    if species_id:
                        species_entry["identifier"] = species_id
                    if species_name:
                        species_entry["name"] = species_name
                    out.setdefault("species", []).append(species_entry)

            genotype = subject.get("genotype") or {}
            mhc_set = genotype.get("mhc_genotype_set")
            if isinstance(mhc_set, dict):
                add_associated_genotype_entry(mhc_set.get("mhc_genotype_set_id"))
                for entry in _ensure_list(mhc_set.get("mhc_genotype_list")):
                    if isinstance(entry, dict):
                        add_associated_genotype_entry(entry.get("mhc_genotype_id"))
                for entry in _ensure_list(mhc_set.get("mhc_genotype_class_list")):
                    if isinstance(entry, dict):
                        add_associated_genotype_entry(entry.get("mhc_genotype_id"))

            for diagnosis in _ensure_list(subject.get("diagnosis")):
                if not isinstance(diagnosis, dict):
                    continue
                add_variable_measurement_entry("disease severity", diagnosis.get("disease_stage") or diagnosis.get("diseaseStage"))
                add_variable_measurement_entry("intervention", diagnosis.get("intervention"))
                add_variable_measurement_entry("prior therapies", diagnosis.get("prior_therapies") or diagnosis.get("priorTherapies"))

        if sexes:
            sex_values = sorted(sexes)
            out["sex"] = sex_values[0] if len(sex_values) == 1 else sex_values
        elif "sex" in out:
            out.pop("sex")

        if not associated_genotypes and "associatedGenotype" in out:
            out.pop("associatedGenotype")

        out["license"] = "https://creativecommons.org/public-domain/"
        out["usageInfo"] = {
            "@type": "CreativeWork",
            "name": "AIRR Data Commons usage info",
            "description": "The software is under the GNU Affero General Public License, the data is Public Domain",
        }

        sample_collections = build_sample_collection_records(sid, md, out)

        # Aggregate sample properties to dataset level for search purposes
        aggregate_sample_properties_to_dataset(out, sample_collections)

        yield out

        for sample_collection in sample_collections:
            yield sample_collection

    logger.info(f"Finished parsing {total} studies")
