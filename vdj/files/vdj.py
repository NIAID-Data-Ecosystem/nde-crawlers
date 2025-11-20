import copy
import json
import logging

import requests
from dataset_parser import build_dataset_record, build_dataset_sample_objects
from parser_utils import _ensure_list
from sample_parser import build_sample_records

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

def parse():
    all_meta = retrieve_study_metadata()
    total = len(all_meta)

    for idx, (sid, md) in enumerate(all_meta.items(), start=1):
        logger.info(f"Parsing {idx}/{total}: {sid}")
        core = md["core"]

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

        dataset_record = build_dataset_record(sid, md, core, samples, subjects)
        sample_records = build_sample_records(sid, md, dataset_record)

        sample_objects = build_dataset_sample_objects(sample_records)
        if sample_objects:
            dataset_record["sample"] = sample_objects

        yield dataset_record

        for sample_collection in sample_records:
            yield sample_collection

    logger.info(f"Finished parsing {total} studies")
