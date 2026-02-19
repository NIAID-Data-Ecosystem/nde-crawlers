import copy
import datetime
import json
import logging
import re

from parser_utils import _add_unique_dict, _add_unique_value, _ensure_list, _iter_string_values

logger = logging.getLogger("nde-logger")

SAMPLE_ONLY_DATASET_FIELDS = {
    "cellType",
    "anatomicalStructure",
    "anatomicalSystem",
    "sex",
    "developmentalStage",
    "associatedGenotype",
    "associatedPhenotype",
    "sampleAvailability",
    "sampleQuantity",
    "sampleType",
}

SAMPLE_ALLOWED_FIELDS = {
    "@type",
    "aggregateElement",
    "itemListElement",
    "numberOfItems",
}


def _remove_sample_only_dataset_fields(dataset_record):
    """Strip sample-only properties from dataset root to enforce schema policy."""
    stray_fields = [field for field in SAMPLE_ONLY_DATASET_FIELDS if field in dataset_record]
    if not stray_fields:
        return

    dataset_id = dataset_record.get("_id") or dataset_record.get("identifier")
    logger.warning(
        "Dataset %s contained sample-only fields on the root object; removing: %s",
        dataset_id,
        ", ".join(sorted(stray_fields)),
    )
    for field in stray_fields:
        dataset_record.pop(field, None)


def build_dataset_record(sid, metadata, core, samples, subjects):
    out = {}
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
        out["measurementTechnique"] = {"@type": "DefinedTerm"}
        if stype.get("id"):
            out["measurementTechnique"]["identifier"] = stype["id"]
            if isinstance(stype["id"], str) and ":" in stype["id"]:
                out["measurementTechnique"]["inDefinedTermSet"] = stype["id"].split(":", 1)[0]
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
    kw_n = []

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

    unique_keywords = list(dict.fromkeys(kw_n))
    if unique_keywords:
        out["keywords"] = unique_keywords

    species_seen = set()
    health_conditions_seen = set()

    species_facet = metadata.get("species_facet") or {}
    disease_facet = metadata.get("disease_facet") or {}

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

    def add_variable_measurement_entry(name, raw_value):
        if raw_value in (None, "", [], {}, set()):
            return
        for value in _iter_string_values(raw_value):
            if value in (None, ""):
                continue
            entry = {"name": name, "description": str(value)}
            _add_unique_dict(variable_measured, variable_seen, entry)

    for subject in subjects:
        add_variable_measurement_entry("ethnicity", subject.get("ethnicity"))
        add_variable_measurement_entry("race", subject.get("race"))

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
            for entry in _ensure_list(mhc_set.get("mhc_genotype_list")):
                if isinstance(entry, dict):
                    add_variable_measurement_entry("mhc_genotype", entry.get("mhc_genotype_id"))
            for entry in _ensure_list(mhc_set.get("mhc_genotype_class_list")):
                if isinstance(entry, dict):
                    add_variable_measurement_entry("mhc_genotype_class", entry.get("mhc_genotype_id"))

        for diagnosis in _ensure_list(subject.get("diagnosis")):
            if not isinstance(diagnosis, dict):
                continue
            add_variable_measurement_entry("disease severity", diagnosis.get("disease_stage") or diagnosis.get("diseaseStage"))
            add_variable_measurement_entry("intervention", diagnosis.get("intervention"))
            add_variable_measurement_entry("prior therapies", diagnosis.get("prior_therapies") or diagnosis.get("priorTherapies"))

    out["license"] = "https://creativecommons.org/public-domain/"
    out["usageInfo"] = {
        "@type": "CreativeWork",
        "name": "AIRR Data Commons usage info",
        "description": "The software is under the GNU Affero General Public License, the data is Public Domain",
    }

    _remove_sample_only_dataset_fields(out)

    return out


def build_dataset_sample_objects(sample_collections):
    if not sample_collections:
        return {}

    valid_collections = [sc for sc in sample_collections if isinstance(sc, dict)]
    if not valid_collections:
        return {}

    excluded_keys = {
        "_id",
        "@type",
        "isBasisFor",
        "sampleList",
        "itemListElement",
        "sample",
        "collectionSize",
        "numberOfItems",
        "identifier",
        "url",
        "includedInDataCatalog",
    }
    multi_value_fields = {
        "associatedPhenotype",
        "associatedGenotype",
        "cellType",
        "anatomicalSystem",
        "anatomicalStructure",
        "sex",
        "developmentalStage",
        "sampleQuantity",
        "sampleType",
        "includedInDataCatalog",
    }

    list_trackers = {}
    sample_collection = {"@type": "SampleCollection"}
    aggregate_data = {}
    collection_size_total = 0
    collection_size_found = False
    EMPTY_VALUES = (None, "", [], {}, set())
    sample_list_context = {"parent_url": None}

    def _is_empty(value):
        return value in EMPTY_VALUES

    def _update_sample_list_context(collection):
        if not isinstance(collection, dict):
            sample_list_context["parent_url"] = None
            return
        sample_list_context["parent_url"] = collection.get("url")

    def _sample_entry_fingerprint(entry):
        if not isinstance(entry, dict):
            return str(entry)
        for candidate in (
            "_id",
            "sample_id",
            "sample_processing_id",
            "vdjserver_uuid",
            "identifier",
            "name",
            "url",
        ):
            value = entry.get(candidate)
            if not _is_empty(value):
                return f"id::{value}"
        try:
            return f"hash::{json.dumps(entry, sort_keys=True, default=str)}"
        except TypeError:
            return f"hash::{str(entry)}"

    def _is_number(value):
        return isinstance(value, (int, float)) and not isinstance(value, bool)

    def _sanitize_sample_list_entry(entry):
        def _stringify(value):
            if isinstance(value, str):
                text = value.strip()
                return text or None
            if value in (None, ""):
                return None
            return str(value)

        if isinstance(entry, dict):
            sanitized = {}

            def _extract_identifier(*fields):
                for candidate in fields:
                    candidate_value = entry.get(candidate)
                    text = _stringify(candidate_value)
                    if text:
                        return text
                return None

            primary_id = _stringify(entry.get("_id")) or _extract_identifier(
                "sample_id",
                "sample_processing_id",
                "vdjserver_uuid",
                "identifier",
                "name",
            )
            if not primary_id:
                fingerprint = _sample_entry_fingerprint(entry)
                primary_id = _stringify(fingerprint)

            identifier_value = _stringify(entry.get("identifier")) or _extract_identifier(
                "name",
                "sample_id",
                "sample_processing_id",
                "vdjserver_uuid",
            )
            if not identifier_value:
                identifier_value = primary_id

            if primary_id:
                sanitized["_id"] = primary_id
            if identifier_value:
                sanitized["identifier"] = identifier_value

            url_value = entry.get("url") or sample_list_context.get("parent_url")
            url_text = _stringify(url_value)
            if url_text:
                sanitized["url"] = url_text
            return sanitized if sanitized else None

        text_value = _stringify(entry)
        if not text_value:
            return None
        sanitized = {"_id": text_value, "identifier": text_value}
        url_text = _stringify(sample_list_context.get("parent_url"))
        if url_text:
            sanitized["url"] = url_text
        return sanitized

    def _normalize_unique_string_list(values):
        normalized = []
        seen = set()

        def _consume(val):
            if isinstance(val, (list, tuple, set)):
                for nested in val:
                    _consume(nested)
                return
            candidate = None
            if isinstance(val, dict):
                candidate = val.get("name") or val.get("value") or val.get("identifier")
            elif val is not None:
                candidate = val
            if candidate is None:
                return
            text = str(candidate).strip().lower()
            if not text or text in seen:
                return
            seen.add(text)
            normalized.append(text)

        for value in _ensure_list(values):
            _consume(value)

        return normalized

    def _collect_numeric_bounds(entry):
        values = []
        for key in ("value", "minValue", "maxValue"):
            val = entry.get(key)
            if _is_number(val):
                values.append(val)
        if not values:
            return None, None
        return min(values), max(values)

    def _aggregate_sample_quantities(entries):
        aggregated = {}
        for raw_entry in _ensure_list(entries):
            if not isinstance(raw_entry, dict):
                continue
            name = raw_entry.get("name")
            if not name:
                continue
            min_val, max_val = _collect_numeric_bounds(raw_entry)
            if min_val is None and max_val is None:
                continue
            bucket = aggregated.setdefault(
                name,
                {
                    "name": name,
                    "minValue": None,
                    "maxValue": None,
                    "unitText": raw_entry.get("unitText"),
                    "unitCode": raw_entry.get("unitCode"),
                },
            )
            if bucket["minValue"] is None or (min_val is not None and min_val < bucket["minValue"]):
                bucket["minValue"] = min_val
            if bucket["maxValue"] is None or (max_val is not None and max_val > bucket["maxValue"]):
                bucket["maxValue"] = max_val
            if not bucket.get("unitText") and raw_entry.get("unitText"):
                bucket["unitText"] = raw_entry["unitText"]
            if not bucket.get("unitCode") and raw_entry.get("unitCode"):
                bucket["unitCode"] = raw_entry["unitCode"]
        results = []
        for bucket in aggregated.values():
            if bucket["minValue"] is None and bucket["maxValue"] is None:
                continue
            entry = {"name": bucket["name"]}
            if bucket["minValue"] is not None:
                entry["minValue"] = bucket["minValue"]
            if bucket["maxValue"] is not None:
                entry["maxValue"] = bucket["maxValue"]
            if bucket.get("unitText"):
                entry["unitText"] = bucket["unitText"]
            if bucket.get("unitCode"):
                entry["unitCode"] = bucket["unitCode"]
            results.append(entry)
        return results

    def _aggregate_developmental_stages(entries):
        aggregated = {}
        for raw_entry in _ensure_list(entries):
            if not isinstance(raw_entry, dict):
                continue
            key = (raw_entry.get("unitCode"), raw_entry.get("unitText"))
            min_val, max_val = _collect_numeric_bounds(raw_entry)
            if min_val is None and max_val is None:
                continue
            bucket = aggregated.setdefault(
                key,
                {
                    "unitCode": raw_entry.get("unitCode"),
                    "unitText": raw_entry.get("unitText"),
                    "minValue": None,
                    "maxValue": None,
                },
            )
            if bucket["minValue"] is None or (min_val is not None and min_val < bucket["minValue"]):
                bucket["minValue"] = min_val
            if bucket["maxValue"] is None or (max_val is not None and max_val > bucket["maxValue"]):
                bucket["maxValue"] = max_val
            if not bucket.get("unitText") and raw_entry.get("unitText"):
                bucket["unitText"] = raw_entry["unitText"]
            if not bucket.get("unitCode") and raw_entry.get("unitCode"):
                bucket["unitCode"] = raw_entry["unitCode"]
        results = []
        for bucket in aggregated.values():
            if bucket["minValue"] is None and bucket["maxValue"] is None:
                continue
            entry = {}
            if bucket.get("unitCode"):
                entry["unitCode"] = bucket["unitCode"]
            if bucket.get("unitText"):
                entry["unitText"] = bucket["unitText"]
            if bucket["minValue"] is not None:
                entry["minValue"] = bucket["minValue"]
            if bucket["maxValue"] is not None:
                entry["maxValue"] = bucket["maxValue"]
            results.append(entry)
        return results

    def _set_number_of_items(value):
        sample_collection["numberOfItems"] = {"value": value, "unitText": "sample"}

    def _merge_list_field(field, values, target):
        if _is_empty(values):
            return
        entries = _ensure_list(values)
        trackers = list_trackers.setdefault((id(target), field), {"dict": set(), "value": set(), "sample_keys": set()})
        current = target.get(field)
        if not isinstance(current, list):
            normalized = []
            if isinstance(current, dict):
                if field == "itemListElement":
                    sanitized_target = _sanitize_sample_list_entry(current)
                    if sanitized_target:
                        fingerprint = _sample_entry_fingerprint(sanitized_target)
                        if fingerprint not in trackers["sample_keys"]:
                            trackers["sample_keys"].add(fingerprint)
                            normalized.append(sanitized_target)
                else:
                    _add_unique_dict(normalized, trackers["dict"], current)
            elif not _is_empty(current):
                _add_unique_value(normalized, trackers["value"], current)
            target[field] = normalized
            current = normalized
        for entry in entries:
            if _is_empty(entry):
                continue
            if isinstance(entry, list):
                _merge_list_field(field, entry, target)
                continue
            if isinstance(entry, dict):
                if field == "itemListElement":
                    sanitized_entry = _sanitize_sample_list_entry(entry)
                    if not sanitized_entry:
                        continue
                    fingerprint = _sample_entry_fingerprint(sanitized_entry)
                    if fingerprint in trackers["sample_keys"]:
                        continue
                    trackers["sample_keys"].add(fingerprint)
                    current.append(sanitized_entry)
                else:
                    _add_unique_dict(current, trackers["dict"], entry)
            else:
                _add_unique_value(current, trackers["value"], entry)

    def _merge_scalar_field(field, value, target):
        if _is_empty(value):
            return
        if field in multi_value_fields:
            _merge_list_field(field, value, target)
            return
        if field == "sampleAvailability":
            target[field] = bool(target.get(field)) or bool(value)
            return

        existing = target.get(field)
        if existing is None:
            target[field] = copy.deepcopy(value)
            return
        if isinstance(existing, list):
            _merge_list_field(field, value, target)
            return
        if isinstance(existing, dict) and isinstance(value, dict):
            for key, val in value.items():
                if key not in existing and not _is_empty(val):
                    existing[key] = copy.deepcopy(val)
            return
        if isinstance(existing, bool) and isinstance(value, bool):
            target[field] = existing or value
            return
        if existing == value:
            return
        _merge_list_field(field, [existing, value], target)

    def _coerce_collection_size(size_value):
        if isinstance(size_value, dict):
            for candidate in ("value", "minValue", "maxValue"):
                val = size_value.get(candidate)
                if isinstance(val, (int, float)) and not isinstance(val, bool):
                    return val
            return None
        if isinstance(size_value, (int, float)) and not isinstance(size_value, bool):
            return size_value
        return None

    for sample_coll in valid_collections:
        _update_sample_list_context(sample_coll)
        explicit_samples = sample_coll.get("sample")
        has_explicit_samples = bool(_ensure_list(explicit_samples))
        _merge_list_field("itemListElement", sample_coll.get("sampleList"), sample_collection)
        _merge_list_field("itemListElement", explicit_samples, sample_collection)

        # Always expose the generated Sample record itself so dataset.sample.itemListElement references the
        # actual Sample documents (e.g., vdj_* identifiers) instead of only the raw source labels.
        if not has_explicit_samples:
            derived_sample_entry = {
                "_id": sample_coll.get("_id"),
                "identifier": sample_coll.get("identifier"),
                "url": sample_coll.get("url"),
                "name": sample_coll.get("name"),
            }
            derived_sample_entry = {k: v for k, v in derived_sample_entry.items() if v not in EMPTY_VALUES}
            if derived_sample_entry:
                _merge_list_field("itemListElement", [derived_sample_entry], sample_collection)

        size_value = _coerce_collection_size(sample_coll.get("numberOfItems") or sample_coll.get("collectionSize"))
        if size_value is not None:
            collection_size_total += size_value
            collection_size_found = True

        for field, value in sample_coll.items():
            if field in excluded_keys:
                continue
            _merge_scalar_field(field, value, aggregate_data)

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

    if sample_collection.get("itemListElement"):
        item_list_entries = sample_collection["itemListElement"]
        if all(isinstance(entry, dict) for entry in item_list_entries):
            item_list_entries.sort(key=_sample_sort_key)
        else:
            sample_collection["itemListElement"] = sorted(item_list_entries)

    if sample_collection.get("itemListElement"):
        _set_number_of_items(len(sample_collection["itemListElement"]))
    elif collection_size_found:
        _set_number_of_items(collection_size_total)
    else:
        _set_number_of_items(len(valid_collections))

    if aggregate_data.get("sex"):
        normalized_sex = _normalize_unique_string_list(aggregate_data["sex"])
        if normalized_sex:
            aggregate_data["sex"] = normalized_sex
        else:
            aggregate_data.pop("sex", None)

    if aggregate_data.get("sampleQuantity"):
        aggregated_quantities = _aggregate_sample_quantities(aggregate_data["sampleQuantity"])
        if aggregated_quantities:
            aggregate_data["sampleQuantity"] = aggregated_quantities
        else:
            aggregate_data.pop("sampleQuantity", None)

    if aggregate_data.get("developmentalStage"):
        aggregated_stages = _aggregate_developmental_stages(aggregate_data["developmentalStage"])
        if aggregated_stages:
            aggregate_data["developmentalStage"] = aggregated_stages
        else:
            aggregate_data.pop("developmentalStage", None)

    if aggregate_data:
        aggregate_element = {"@type": "Sample"}
        aggregate_element.update(aggregate_data)
        sample_collection["aggregateElement"] = aggregate_element

    for key in list(sample_collection.keys()):
        if key not in SAMPLE_ALLOWED_FIELDS:
            sample_collection.pop(key, None)

    if not sample_collection.get("aggregateElement") and not sample_collection.get("itemListElement"):
        return {}

    return sample_collection
