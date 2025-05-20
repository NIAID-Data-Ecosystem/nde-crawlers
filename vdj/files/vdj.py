import datetime
import logging
import re

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def retrieve_all_study_ids():
    logger.info("Retrieving all study IDs from VDJ Server API")
    study_url = "https://vdjserver.org/airr/v1/repertoire"
    all_ids = set()
    size, page = 100, 0

    while True:
        q = {"from": page * size, "size": size}
        try:
            resp = requests.post(study_url, json=q)
            resp.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logger.error(f"Error paging study IDs: {e}")
            break

        hits = resp.json().get("Repertoire", [])
        if not hits:
            break
        for rec in hits:
            sid = rec.get("study", {}).get("study_id")
            if sid:
                all_ids.add(sid)
        if len(hits) < size:
            break
        page += 1

    logger.info(f"Retrieved {len(all_ids)} study IDs")
    return all_ids


def retrieve_study_metadata():
    base = "https://vdjserver.org/airr/v1/repertoire"
    ids = retrieve_all_study_ids()
    meta = {}

    for sid in ids:
        logger.info(f"Fetching core record for study {sid}")
        core_q = {
            "filters": {"op": "=", "content": {"field": "study.study_id", "value": sid}},
            "size": 1,
            "include_fields": "airr-core",
        }
        try:
            r = requests.post(base, json=core_q)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logger.error(f"Core fetch failed for {sid}: {e}")
            continue

        reps = r.json().get("Repertoire", [])
        if not reps:
            continue

        meta[sid] = {
            "core": reps[0],
            "disease_facet": get_disease_facet_data(sid),
            "species_facet": get_species_facet_data(sid),
        }

    logger.info(f"Metadata collected for {len(meta)} studies")
    return meta


def get_disease_facet_data(sid):
    return _facet(sid, "subject.diagnosis.disease_diagnosis")


def get_species_facet_data(sid):
    return _facet(sid, "subject.species")


def _facet(sid, field):
    q = {"filters": {"op": "=", "content": {"field": "study.study_id", "value": sid}}, "facets": field}
    try:
        r = requests.post("https://vdjserver.org/airr/v1/repertoire", json=q)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.HTTPError as e:
        logger.error(f"Facet {field} failed for {sid}: {e}")
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
            "dataset": stud_url,
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

        for b in md["species_facet"].get("facets", {}).get("subject.species", []):
            sp = {}
            species_id = b.get("id")
            species_name = b.get("label")
            if species_id:
                sp["identifier"] = species_id
            if species_name:
                sp["name"] = species_name
            if species_id or species_name:
                out.setdefault("species", []).append(sp)

        cs = sample.get("cell_species", {})
        if cs.get("label"):
            out.setdefault("species", []).append({"name": cs["label"]})
        for b in md["disease_facet"].get("facets", {}).get("subject.diagnosis.disease_diagnosis", []):
            dnm = b.get("key") or b.get("value")
            if dnm:
                out.setdefault("healthCondition", []).append({"name": dnm})

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
            "name": "VDJ Server usage info",
            "description": "The software is under the GNU Affero General Public License, the data is Public Domain",
        }

        yield out

    logger.info(f"Finished parsing {total} studies")
