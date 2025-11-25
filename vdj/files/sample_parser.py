import copy
import datetime
import json
import logging
import re
from decimal import Decimal, InvalidOperation

from parser_utils import _add_unique_dict, _add_unique_value, _ensure_list, _iter_string_values, _sanitize_identifier

logger = logging.getLogger("nde-logger")


def _build_sample_record(
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

    sample_objects = []
    sample_seen_keys = set()
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
    quantity_warning_tracker = set()
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

    numeric_pattern = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

    def _normalize_unit(value):
        if isinstance(value, dict):
            for key in ("label", "name", "id", "identifier"):
                candidate = value.get(key)
                if isinstance(candidate, str):
                    candidate = candidate.strip()
                if candidate:
                    return candidate
            return None
        if isinstance(value, str):
            text = value.strip()
            return text or None
        return value

    def _coerce_numeric(raw):
        if raw in (None, ""):
            return None
        if isinstance(raw, bool):
            return None
        if isinstance(raw, (int, float)):
            return int(raw) if isinstance(raw, int) or (isinstance(raw, float) and raw.is_integer()) else float(raw)
        if isinstance(raw, (list, tuple)) and raw:
            return _coerce_numeric(raw[0])
        text = str(raw).strip()
        if not text:
            return None
        cleaned = text.replace(",", "")
        match = numeric_pattern.search(cleaned)
        candidate = match.group(0) if match else cleaned
        try:
            decimal_value = Decimal(candidate)
        except InvalidOperation:
            return None
        if decimal_value == decimal_value.to_integral_value():
            return int(decimal_value)
        return float(decimal_value)

    def add_quantity(value, unit=None, label=None):
        if value in (None, "", [], {}, set()):
            return

        derived_unit = _normalize_unit(unit)
        derived_label = label.strip() if isinstance(label, str) else label
        raw_value = value

        if isinstance(value, dict):
            raw_value = value.get("value")
            if raw_value in (None, ""):
                raw_value = value.get("amount")
            if raw_value in (None, ""):
                raw_value = value.get("minValue")
            if raw_value in (None, ""):
                raw_value = value.get("maxValue")
            derived_unit = derived_unit or _normalize_unit(value.get("unitText") or value.get("unit") or value.get("unitLabel"))
            candidate_label = value.get("name") or value.get("label")
            if isinstance(candidate_label, str) and candidate_label.strip():
                derived_label = derived_label or candidate_label.strip()
        elif isinstance(value, (list, tuple)) and value:
            raw_value = value[0]

        numeric_value = _coerce_numeric(raw_value)
        if numeric_value is None:
            warning_key = (str(value), derived_label or label or "sampleQuantity")
            if warning_key not in quantity_warning_tracker:
                quantity_warning_tracker.add(warning_key)
                logger.warning(
                    "Skipping non-numeric sampleQuantity value %r for label %r",
                    value,
                    derived_label or label or "sampleQuantity",
                )
            return

        entry = {"value": numeric_value}
        if derived_unit:
            entry["unitText"] = derived_unit
        if derived_label:
            entry["name"] = derived_label
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

    def _sample_fingerprint(entry):
        for candidate in ("sample_id", "sample_processing_id", "vdjserver_uuid", "identifier", "name"):
            value = entry.get(candidate)
            if value not in (None, ""):
                return f"id::{value}"
        try:
            return f"hash::{json.dumps(entry, sort_keys=True, default=str)}"
        except TypeError:
            return f"hash::{str(entry)}"

    def add_sample_object(entry):
        if not isinstance(entry, dict):
            return
        sample_copy = copy.deepcopy(entry)
        fingerprint = _sample_fingerprint(sample_copy)
        if fingerprint in sample_seen_keys:
            return
        sample_seen_keys.add(fingerprint)
        sample_objects.append(sample_copy)

    for sample in samples:
        add_sample_object(sample)

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

    def _sample_sort_key(entry):
        if isinstance(entry, dict):
            for candidate in ("sample_id", "sample_processing_id", "vdjserver_uuid", "identifier", "name"):
                value = entry.get(candidate)
                if value not in (None, ""):
                    return str(value)
            try:
                return json.dumps(entry, sort_keys=True, default=str)
            except TypeError:
                return str(entry)
        return str(entry)

    sample_objects.sort(key=_sample_sort_key)
    if sample_objects:
        # Preserve the normalized raw samples temporarily so dataset.sample.sampleList can
        # enumerate every underlying sample rather than just the aggregate Sample record.
        sample_record["sample"] = sample_objects
    sample_record["collectionSize"] = {
        "value": len(sample_objects) if sample_objects else len(samples)
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


def build_sample_records(sid, metadata, dataset_record):
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
        record = _build_sample_record(
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
