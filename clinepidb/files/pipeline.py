import datetime
import logging
import re

import requests
import xmltodict

logger = logging.getLogger("clinepidb-logger")


_NCIT_LABEL_CACHE = {}


def _looks_like_url(value):
    return isinstance(value, str) and value.startswith(("http://", "https://"))


def _ncit_code_from_obo_purl(url):
    if not _looks_like_url(url):
        return None
    token = url.strip().rstrip("/").split("/")[-1]
    match = re.search(r"NCIT_C(\d+)$", token, flags=re.IGNORECASE)
    if not match:
        return None
    return "C" + match.group(1)


def _resolve_ncit_label_from_obo_purl(url):
    if not _looks_like_url(url):
        return None
    if url in _NCIT_LABEL_CACHE:
        return _NCIT_LABEL_CACHE[url]

    code = _ncit_code_from_obo_purl(url)
    if not code:
        _NCIT_LABEL_CACHE[url] = None
        return None

    try:
        resp = requests.get(
            f"https://api-evsrest.nci.nih.gov/api/v1/concept/ncit/{code}",
            timeout=10,
        )
        if resp.status_code != 200:
            _NCIT_LABEL_CACHE[url] = None
            return None
        payload = resp.json() if resp.content else {}
        label = payload.get("name")
        if isinstance(label, str) and label.strip():
            _NCIT_LABEL_CACHE[url] = label.strip()
            return _NCIT_LABEL_CACHE[url]
    except Exception:
        pass

    _NCIT_LABEL_CACHE[url] = None
    return None


def _split_csv(value):
    if value in (None, ""):
        return []
    if isinstance(value, (list, tuple, set)):
        parts = []
        for entry in value:
            parts.extend(_split_csv(entry))
        return parts
    if not isinstance(value, str):
        value = str(value)
    return [part.strip() for part in value.split(",") if part and part.strip()]


def _coerce_int(value):
    if value in (None, ""):
        return None
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if not isinstance(value, str):
        value = str(value)
    text = value.strip()
    if not text or text.upper() in {"N/A", "NA", "NONE", "NULL"}:
        return None
    text = text.replace(",", "")
    try:
        return int(float(text))
    except Exception:
        return None


def _get_study_characteristic_value(record, key):
    try:
        table = record.get("tables", {}).get("StudyCharacteristicTable")
        if isinstance(table, list) and table:
            return table[0].get(key)
    except Exception:
        return None
    return None


def _add_sample_quantity(sample, raw_value, unit_text, unit_code=None):
    value = _coerce_int(raw_value)
    if value is None:
        return

    resolved_unit_text = unit_text
    # If unitText is actually an NCIT URL, resolve it to the preferred label.
    if _looks_like_url(unit_text):
        resolved = _resolve_ncit_label_from_obo_purl(unit_text)
        if resolved:
            resolved_unit_text = resolved
        if unit_code is None:
            unit_code = unit_text
    # Fallback: if unitCode is an NCIT URL and unitText is missing/URL-like, resolve from unitCode.
    elif _looks_like_url(unit_code) and (unit_text in (None, "") or _looks_like_url(unit_text)):
        resolved = _resolve_ncit_label_from_obo_purl(unit_code)
        if resolved:
            resolved_unit_text = resolved

    if unit_text == "Clinical Study Participants":
        entry = {"@type": "QuantitativeValue", "value": value, "unitText": "Study Subjects", "unitCode": "http://purl.obolibrary.org/obo/NCIT_C41189"}
    else:
        entry = {"@type": "QuantitativeValue", "value": value, "unitText": resolved_unit_text}
        if unit_code:
            entry["unitCode"] = unit_code

    existing = sample.get("sampleQuantity")
    if existing is None:
        sample["sampleQuantity"] = [entry]
        return
    if isinstance(existing, dict):
        sample["sampleQuantity"] = [existing]
        existing = sample["sampleQuantity"]
    if isinstance(existing, list) and entry not in existing:
        existing.append(entry)


def get_variableMeasured(entity):
    variableMeasured = []
    # Skip entities with displayName "Sample"
    if entity.get("displayName") == "Sample":
        return variableMeasured

    # Collect variable display names
    variables = entity.get("variables", [])
    for variable in variables:
        displayName = variable.get("displayName")
        if displayName:
            variableMeasured.append({"name": displayName})

    # Recursively traverse child entities
    children = entity.get("children", [])
    for child in children:
        variableMeasured.extend(get_variableMeasured(child))

    return variableMeasured


def record_generator():
    logger.info("record generator function started")
    xml_url = "https://clinepidb.org/ce/sitemap-SitemapDatasets.xml"
    sitemap_request = requests.get(xml_url)
    dict_data = xmltodict.parse(sitemap_request.content)
    data_urls = [hit["loc"] for hit in dict_data["urlset"]["url"]]
    logger.info("%s datasets found" % len(data_urls))
    cond_acc_dict = {
        "public": "Open",
        "protected": "Restricted",
        "controlled": "Restricted",
        "prerelease": "Restricted",
    }

    for url in data_urls:
        dataset_id = url.split("/")[-1]

        cookies = {
            "JSESSIONID": "34D7C5C6EDE1FBA34B0956E13896A19D",
        }

        headers = {
            "Accept": "*/*",
            "Accept-Language": "en-US,en;q=0.9",
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "Origin": "https://clinepidb.org",
            "Pragma": "no-cache",
            "Referer": "https://clinepidb.org/ce/app/workspace/analyses/%s/new/details" % dataset_id,
            "content-type": "application/json",
        }

        data = (
            '{"attributes":["dataset_id","eda_study_id","study_access", "disease","sex","Sample_Type","Study_Type","Study_Design", \
                "HH_count","Variable_count","approved", "Part_count","percent_approved","total","total_number","Obser_count","Samp_count","Country","Years","AdditionalData","WHO", \
                "ProjectName","contact","institution","description","eupath_release","email","summary"], \
                "primaryKey":[{"name":"dataset_id","value":%s}], \
                "tables":["StudyCharacteristicTable","Publications","DownloadVersion","Contacts","AccessRequestStats","AccessRequest","HyperLinks"]}'
            % dataset_id
        )

        response = requests.post(
            "https://clinepidb.org/ce/service/record-types/dataset/records", headers=headers, cookies=cookies, data=data
        )
        record = response.json()

        # set initial available values
        record["identifier"] = record["id"][0].pop("value")
        record["name"] = record.pop("displayName")

        url = "https://clinepidb.org/ce/app/record/dataset/" + record["identifier"]

        # add keys
        record.update(
            {
                "_id": "clinepidb_" + record["identifier"],
                "@type": "Dataset",
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "ClinEpiDB",
                    "url": "https://clinepidb.org/ce/app/",
                    "versionDate": datetime.date.today().isoformat(),
                    "archivedAt": url,
                },
                "url": url,
            }
        )

        try:
            # attributes
            measurement_technique = record["attributes"].pop("Study_Design")
            if measurement_technique:
                record["measurementTechnique"] = {"name": measurement_technique}
            record["description"] = record["attributes"].pop("description")
            record["description"] += record["attributes"].pop("summary")
            keywords = []
            if record["attributes"].get("WHO"):
                keywords += record["attributes"].pop("WHO").split(",")
            # clinepidb removed Participant_Type from attributes
            # record["keywords"] += record["attributes"].pop("Participant_Type")
            # record["keywords"] += record["attributes"].pop("Study_Type")
            if record["attributes"].get("Study_Type"):
                keywords += record["attributes"].pop("Study_Type").split(",")
            if keywords:
                record["keywords"] = keywords
            record["spatialCoverage"] = {"name": record["attributes"].pop("Country")}

            ##version_date =  record['attributes'].pop('eupath_release').split(",")[0]
            # iso_version_date = datetime.datetime.strptime(version_date, "%Y-%m-%d").date().isoformat()
            # record['version'] = {str( record['attributes'].pop('eupath_release'))}

            cond_acc = cond_acc_dict[record["attributes"]["study_access"].lower()]
            record["conditionsOfAccess"] = cond_acc

            temp_cov = record["attributes"].pop("Years")
            record["temporalCoverage"] = {
                "@type": "TemporalInterval",
                "startDate": (temp_cov.split(" ")[0].replace(",", "")),
                "endDate": (temp_cov.split(" ")[-1]),
            }
            if not temp_cov:
                record["temporalCoverage"] = {
                    "@type": "TemporalInterval",
                    "startDate": record["tables"]["StudyCharacteristicTable"][0]["Years"][0].replace(",", ""),
                    "endDate": record["tables"]["StudyCharacteristicTable"][0]["Years"][-1],
                }

            # tables.StudyCharacteriticTable -- get these variables if the variables above are missing
            health_cond = record["attributes"].pop("disease") or record["tables"]["StudyCharacteristicTable"][0].pop(
                "disease"
            )
            if health_cond:
                record["healthCondition"] = {"name": health_cond}

            # # Sample metadata (minimal Sample object inside Dataset)
            # base_sample_type = {
            #     "@type": "DefinedTerm",
            #     "name": "Study Subject",
            #     "url": "http://purl.obolibrary.org/obo/NCIT_C41189",
            #     "inDefinedTermSet": "NCIT",
            #     "termCode": "NCIT_C41189",
            # }
            sample = {"@type": "Sample", "sampleType": []}

            # Sample_Type is often a comma-separated list (e.g. "Blood sample, Urine sample")
            sample_type_raw = _get_study_characteristic_value(record, "Sample_Type") or record.get("attributes", {}).get(
                "Sample_Type"
            )
            sample_type_parts = _split_csv(sample_type_raw)
            if sample_type_parts:
                normalized = []
                seen = set()
                for part in sample_type_parts:
                    key = part.casefold()
                    if key in seen:
                        continue
                    seen.add(key)
                    normalized.append({"@type": "DefinedTerm", "name": part})
                # Skip sentinel values like "No samples" while still keeping the base Study Subject term
                normalized = [
                    hit
                    for hit in normalized
                    if hit.get("name", "").casefold() not in {"no samples", "no sample"}
                ]
                if normalized:
                    sample["sampleType"].extend(normalized)

            sex_value = record.get("attributes", {}).get("sex") or _get_study_characteristic_value(record, "sex")
            sex_parts = _split_csv(sex_value)
            if sex_parts:
                sample["sex"] = sex_parts[0] if len(sex_parts) == 1 else sorted(set(sex_parts))

            participant_type = _get_study_characteristic_value(record, "Participant_Type")
            participant_parts = _split_csv(participant_type)
            if participant_parts:
                sample["developmentalStage"] = [
                    {"@type": "DefinedTerm", "name": part} for part in sorted(set(participant_parts))
                ]

            # If participant metadata exists, tag sampleType with Clinical Study Participants.
            table_part_count = _get_study_characteristic_value(record, "Part_count")
            if participant_parts or table_part_count not in (None, ""):
                clinical_participants_term = {
                    "@type": "DefinedTerm",
                    "name": "Study Subject",
                    "url": "http://purl.obolibrary.org/obo/NCIT_C41189",
                    "inDefinedTermSet": "NCIT",
                    "termCode": "NCIT_C41189",
                }
                existing_types = sample.get("sampleType")
                if isinstance(existing_types, list):
                    existing_keys = {
                        (t.get("name", "").casefold(), t.get("termCode"), t.get("url"))
                        for t in existing_types
                        if isinstance(t, dict)
                    }
                    key = (
                        clinical_participants_term.get("name", "").casefold(),
                        clinical_participants_term.get("termCode"),
                        clinical_participants_term.get("url"),
                    )
                    if key not in existing_keys:
                        existing_types.append(clinical_participants_term)

            hh_count = record.get("attributes", {}).get("HH_count") or _get_study_characteristic_value(record, "HH_count")
            obser_count = record.get("attributes", {}).get("Obser_count") or _get_study_characteristic_value(record, "Obser_count")
            samp_count = record.get("attributes", {}).get("Samp_count") or _get_study_characteristic_value(record, "Samp_count")
            part_count = record.get("attributes", {}).get("Part_count") or table_part_count

            _add_sample_quantity(
                sample,
                hh_count,
                unit_text="Households",
                unit_code="http://purl.obolibrary.org/obo/NCIT_C41194",
            )
            _add_sample_quantity(
                sample,
                obser_count,
                unit_text="Observations",
                unit_code="http://purl.obolibrary.org/obo/NCIT_C93430",
            )
            _add_sample_quantity(sample, samp_count, unit_text="Samples")
            _add_sample_quantity(
                sample,
                part_count,
                unit_text="Clinical Study Participants",
                unit_code="http://purl.obolibrary.org/obo/NCIT_C70668",
            )
            # if sampleType is empty pop it
            if not sample["sampleType"]:
                sample.pop("sampleType")
            record["sample"] = sample

            # tables.Publications -- helper function to create citations
            pmid_list = [hit.pop("pmid") for hit in record["tables"]["Publications"]]

            if pmid_list:
                if len(pmid_list) == 1:
                    record["pmids"] = pmid_list[0]
                else:
                    record["pmids"] = ",".join(pmid_list)

            # tables.DowloadVersion
            date_dict = {
                "JAN": "01",
                "FEB": "02",
                "MAR": "03",
                "APR": "04",
                "MAY": "05",
                "JUN": "06",
                "JUL": "07",
                "AUG": "08",
                "SEP": "09",
                "OCT": "10",
                "NOV": "11",
                "DEC": "12",
            }
            distr_dic = {}
            recent_date = None
            for hit in record["tables"]["DownloadVersion"]:
                date = hit["release_date"]
                month = date.split("-")[1]
                if month in date_dict.keys():
                    new_date = date.split("-")[0] + "-" + date_dict[month] + "-" + date.split("-")[2]
                    distr_dic[new_date] = hit
            recent_date = sorted(
                [datetime.datetime.strptime(date, "%Y-%m-%d").date().isoformat() for date in distr_dic]
            )[-1]
            record["distribution"] = {"dateModified": recent_date, "name": distr_dic[recent_date].pop("dataset_name")}

            # tables.Contacts -- multiple
            record["author"] = [
                {"name": hit["contact_name"], "affiliation": {"name": hit["affiliation"]}}
                for hit in record["tables"]["Contacts"]
            ]

            # tables.HyperLinks
            hl_list = []
            for dict_ in record["tables"]["HyperLinks"]:
                # TODO DO WE NEED TO CAPTURE THE URL?
                # url_ = dict_.pop("url")
                dict_.pop("url")
                hl_url = dict_["hyper_link"].pop("url")
                if "NCT" in hl_url:
                    hl_list.append({"identifier": hl_url.split("/")[-1], "url": hl_url})
                else:
                    hl_list.append({"url": hl_url})

            record["isBasedOn"] = hl_list

            # Get variableMeasured
            eda_study_id = record["attributes"]["eda_study_id"]
            study_url = "https://clinepidb.org/eda/studies/%s" % eda_study_id
            cookies = {
                "Authorization": "eyJhbGciOiJFUzUxMiJ9.eyJzdWIiOiIxMTY1MjQyMjIzIiwiaXNfZ3Vlc3QiOnRydWUsImlzcyI6Imh0dHBzOi8vZXVwYXRoZGIub3JnL29hdXRoIiwiYXVkIjoiYXBpQ29tcG9uZW50U2l0ZSIsImF6cCI6ImFwaUNvbXBvbmVudFNpdGUiLCJhdXRoX3RpbWUiOjE3MzI1NTUyOTksImlhdCI6MTczMjU1NTI5OSwiZXhwIjoxODI3MTYzMjk5fQ.AVyvEqT-DWNiMlsjyRZ1KsG3KXELjc6y3PGMpRVJdre0ODsI9tVRn6UEzT2MKhJjKzo1FHYADVDe-SJUc-T6vmy6APn-1JsAitDe0MvbgenWwBoGYgXNx7T8pQ9uosJmyGBE9-yDq4Jmrk87WEsYUylQ2nHQD6fM022jRXVMG_GS7Hk1"
            }
            study_response = requests.get(study_url, headers=headers, cookies=cookies)
            study_data = study_response.json()
            variableMeasured = get_variableMeasured(study_data["study"]["rootEntity"])
            record["variableMeasured"] = variableMeasured

        except Exception as e:
            logger.error("identifier:" + record["identifier"] + ": " + str(e))

        # Remove excess data
        record.pop("attributes")
        record.pop("tables")
        record.pop("id")
        record.pop("recordClassName")
        record.pop("tableErrors")

        yield record
