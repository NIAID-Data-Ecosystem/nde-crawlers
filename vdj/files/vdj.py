import datetime
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


def retrieve_study_metadata():
    ids = retrieve_all_study_ids()
    meta = {}

    for sid in ids:
        logger.info(f"Fetching core record for study {sid}")

        # Try each endpoint until we find the study
        study_data = None
        for endpoint in AIRR_ENDPOINTS:
            core_q = {
                "filters": {"op": "=", "content": {"field": "study.study_id", "value": sid}},
                "size": 1,
                "include_fields": "airr-core",
            }
            try:
                r = requests.post(endpoint, json=core_q, timeout=30)
                r.raise_for_status()
                reps = r.json().get("Repertoire", [])
                if reps:
                    study_data = {
                        "core": reps[0],
                        "endpoint": endpoint,
                        "disease_facet": get_disease_facet_data(sid, endpoint),
                        "species_facet": get_species_facet_data(sid, endpoint),
                    }
                    break
            except requests.exceptions.RequestException as e:
                logger.warning(f"Failed to fetch {sid} from {endpoint}: {e}")
                continue

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
        samples = core.get("sample") or []

        if isinstance(samples, dict):
            samples = [samples]

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

        diagnosis = core.get("subject", {}).get("diagnosis", [])
        if isinstance(diagnosis, dict):
            diagnosis = [diagnosis]
        for diag in diagnosis:
            if diag.get("disease_stage"):
                kw_n.append(diag["disease_stage"])
                break

        all_keywords = kw_n + kw_id
        unique_keywords = list(dict.fromkeys(all_keywords))
        if unique_keywords:
            out["keywords"] = unique_keywords

        species_seen = set()
        health_conditions_seen = set()

        for b in md["species_facet"].get("Facet", []):
            print("Processing species facet:", b)
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

        cs = sample.get("cell_species", {})
        if cs.get("label") and cs["label"] not in species_seen:
            species_seen.add(cs["label"])
            out.setdefault("species", []).append({"name": cs["label"]})

        for b in md["disease_facet"].get("Facet", []):
            print("Processing disease facet:", b)
            disease_data = b.get("subject.diagnosis.disease_diagnosis", {})
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

        out["license"] = "https://creativecommons.org/public-domain/"
        out["usageInfo"] = {
            "@type": "CreativeWork",
            "name": "AIRR Data Commons usage info",
            "description": "The software is under the GNU Affero General Public License, the data is Public Domain",
        }

        yield out

    logger.info(f"Finished parsing {total} studies")
