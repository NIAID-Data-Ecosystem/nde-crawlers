import datetime
import logging

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


_BROWSE_CACHE = {}


def _ensure_list(value):
    if value in (None, ""):
        return []
    if isinstance(value, list):
        return value
    return [value]


def _first_scalar(value):
    for entry in _ensure_list(value):
        if entry in (None, ""):
            continue
        if isinstance(entry, str):
            entry = entry.strip()
            if not entry:
                continue
        return entry
    return None


def _is_unknown_text(value):
    if value in (None, ""):
        return True
    if not isinstance(value, str):
        return False
    return value.strip().casefold() == "unknown"


def _first_text(value):
    """Return the first non-empty, non-'Unknown' scalar text value."""
    text = _first_scalar(value)
    if isinstance(text, str):
        text = text.strip()
        if not text or _is_unknown_text(text):
            return None
        return text
    if _is_unknown_text(text):
        return None
    return text


def _fetch_browse_json(entity_type, uuid, timeout=20):
    cache_key = (entity_type, uuid)
    if cache_key in _BROWSE_CACHE:
        return _BROWSE_CACHE[cache_key]

    url = f"https://portal.hubmapconsortium.org/browse/{entity_type}/{uuid}.json"
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
        if resp.status_code != 200:
            _BROWSE_CACHE[cache_key] = None
            return None
        payload = resp.json()
        _BROWSE_CACHE[cache_key] = payload
        return payload
    except Exception:
        _BROWSE_CACHE[cache_key] = None
        return None


def get_ids():
    logger.info("Getting datasets")
    url = "https://search.api.hubmapconsortium.org/v3/portal/search"
    payload = {
        "query": {
            "bool": {
                "must": [
                    {
                        "bool": {
                            "must_not": [
                                {"exists": {"field": "next_revision_uuid"}},
                                {"exists": {"field": "sub_status"}},
                            ]
                        }
                    },
                    {
                        "bool": {
                            "must_not": [
                                {"exists": {"field": "next_revision_uuid"}},
                                {"exists": {"field": "sub_status"}},
                            ]
                        }
                    },
                ]
            }
        },
        "post_filter": {"term": {"entity_type.keyword": "Dataset"}},
        "_source": False,
        "size": 10000,
    }
    r = requests.post(url, json=payload).json()
    ids = [x["_id"] for x in r["hits"]["hits"]]
    logger.info("Total datasets: %s", len(ids))
    count = 0
    datasets = []
    for id in ids:
        count += 1
        if count % 50 == 0:
            logger.info("Retrieved %s out of %s datasets", count, len(ids))
        dataset = requests.get(
            f"https://portal.hubmapconsortium.org/browse/dataset/{id}.json",
            headers={"Accept": "application/json"},
            timeout=30,
        )
        datasets.append(dataset.json())
    logger.info("Retrieved %s out of %s datasets", count, len(ids))
    return datasets


def get_sample_ids():
    logger.info("Getting samples")
    url = "https://search.api.hubmapconsortium.org/v3/portal/search"
    payload = {
        "query": {
            "bool": {
                "must": [
                    {
                        "bool": {
                            "must_not": [
                                {"exists": {"field": "next_revision_uuid"}},
                                {"exists": {"field": "sub_status"}},
                            ]
                        }
                    },
                    {
                        "bool": {
                            "must_not": [
                                {"exists": {"field": "next_revision_uuid"}},
                                {"exists": {"field": "sub_status"}},
                            ]
                        }
                    },
                ]
            }
        },
        "post_filter": {"term": {"entity_type.keyword": "Sample"}},
        "_source": False,
        "size": 10000,
    }
    r = requests.post(url, json=payload, timeout=30).json()
    ids = [x["_id"] for x in r.get("hits", {}).get("hits", [])]
    logger.info("Total samples: %s", len(ids))
    return ids


def _build_donor_description(donor_mapped):
    if not isinstance(donor_mapped, dict) or not donor_mapped:
        return None

    parts = []
    age_value = _first_scalar(donor_mapped.get("age_value"))
    age_unit = _first_scalar(donor_mapped.get("age_unit"))
    if age_value is not None and age_unit:
        parts.append(f"{age_value} {age_unit} old")

    race_value = _first_text(donor_mapped.get("race"))
    if race_value is not None:
        parts.append(str(race_value))

    sex_value = _first_text(donor_mapped.get("sex"))
    if sex_value is not None:
        parts.append(str(sex_value))

    cause_of_death = _first_text(donor_mapped.get("cause_of_death"))
    if cause_of_death is not None:
        parts.append(str(cause_of_death))

    mechanism_of_injury = _first_text(donor_mapped.get("mechanism_of_injury"))
    if mechanism_of_injury is not None:
        parts.append(str(mechanism_of_injury))

    text = ", ".join([p for p in parts if p])
    return text or None


def parse_samples():
    sample_ids = get_sample_ids()
    count = 0
    for sid in sample_ids:
        count += 1
        if count % 200 == 0:
            logger.info("Parsed %s out of %s samples", count, len(sample_ids))

        metadata = _fetch_browse_json("sample", sid, timeout=30)
        if not isinstance(metadata, dict):
            continue

        output = {
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "HuBMAP",
                "url": "https://hubmapconsortium.org/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Sample",
        }

        if uuid := metadata.get("uuid"):
            output["_id"] = "HUBMAP_SAMPLE_" + uuid
            output["identifier"] = metadata.get("hubmap_id") or uuid
            output["url"] = f"https://portal.hubmapconsortium.org/browse/sample/{uuid}"
            output["includedInDataCatalog"]["archivedAt"] = output["url"]

        if name := metadata.get("hubmap_id"):
            output["name"] = name

        if created_timestamp := metadata.get("created_timestamp"):
            output["dateCreated"] = datetime.datetime.utcfromtimestamp(created_timestamp / 1000).strftime("%Y-%m-%d")

        if data_access_level := metadata.get("data_access_level"):
            output["isAccessibleForFree"] = data_access_level == "public"

        sample_category = (
            metadata.get("mapped_sample_category")
            or metadata.get("sample_category")
            or metadata.get("display_subtype")
        )
        if sample_category:
            output["sampleType"] = [{"@type": "DefinedTerm", "name": str(sample_category)}]

        organs = []
        for organ in _ensure_list(metadata.get("origin_samples_unique_mapped_organs")):
            organ_value = _first_scalar(organ)
            if organ_value:
                organs.append(organ_value)
        if not organs:
            for origin in _ensure_list(metadata.get("origin_samples")):
                if isinstance(origin, dict):
                    organ_value = origin.get("mapped_organ") or origin.get("organ")
                    organ_value = _first_scalar(organ_value)
                    if organ_value:
                        organs.append(organ_value)

        if organs:
            seen = set()
            anatomical = []
            for organ in organs:
                key = str(organ).casefold()
                if key in seen:
                    continue
                seen.add(key)
                anatomical.append({"@type": "DefinedTerm", "name": str(organ)})
            if anatomical:
                output["anatomicalStructure"] = anatomical

        donor_desc = None
        donor = metadata.get("donor") if isinstance(metadata, dict) else None
        donor_mapped = donor.get("mapped_metadata", {}) if isinstance(donor, dict) else {}
        if isinstance(donor_mapped, dict) and donor_mapped:
            sex_value = _first_text(donor_mapped.get("sex"))
            if sex_value:
                output["sex"] = sex_value

            age_value = _first_scalar(donor_mapped.get("age_value"))
            age_unit = _first_scalar(donor_mapped.get("age_unit"))
            if age_value is not None and age_unit:
                output["developmentalStage"] = {
                    "@type": "QuantitativeValue",
                    "value": age_value,
                    "unitText": age_unit,
                }

            phenotypes = []
            race_value = _first_text(donor_mapped.get("race"))
            if race_value:
                phenotypes.append({"@type": "DefinedTerm", "name": str(race_value)})
            blood_type = _first_text(donor_mapped.get("abo_blood_group_system"))
            if blood_type:
                phenotypes.append({"@type": "DefinedTerm", "name": str(blood_type)})
            if phenotypes:
                output["associatedPhenotype"] = phenotypes

            health_conditions = []
            for condition in _ensure_list(donor_mapped.get("medical_history")):
                condition_value = _first_text(condition)
                if condition_value:
                    health_conditions.append({"name": str(condition_value)})

            mechanism_of_injury = _first_text(donor_mapped.get("mechanism_of_injury"))
            if mechanism_of_injury:
                health_conditions.append(
                    {
                        "@type": "DefinedTerm",
                        "name": str(mechanism_of_injury),
                        "url": "http://purl.obolibrary.org/obo/MONDO_0021178",
                        "inDefinedTermSet": "MONDO",
                    }
                )

            if health_conditions:
                output["healthCondition"] = health_conditions

            variable_measured = []
            if donor_mapped.get("body_mass_index_value") not in (None, ""):
                variable_measured.append(
                    {
                        "@type": "DefinedTerm",
                        "name": "Body mass index",
                        "url": "http://purl.obolibrary.org/obo/CMO_0000105",
                        "inDefinedTermSet": "CMO",
                    }
                )
            if donor_mapped.get("height_value") not in (None, ""):
                variable_measured.append(
                    {
                        "@type": "DefinedTerm",
                        "name": "Height",
                        "url": "http://purl.obolibrary.org/obo/CMO_0000106",
                        "inDefinedTermSet": "CMO",
                    }
                )
            if donor_mapped.get("weight_value") not in (None, ""):
                variable_measured.append(
                    {
                        "@type": "DefinedTerm",
                        "name": "Weight",
                        "url": "http://purl.obolibrary.org/obo/CMO_0000012",
                        "inDefinedTermSet": "CMO",
                    }
                )
            if _first_text(donor_mapped.get("cause_of_death")):
                variable_measured.append(
                    {
                        "@type": "DefinedTerm",
                        "name": "Cause of death",
                        "url": "http://purl.obolibrary.org/obo/NCIT_C81239",
                        "inDefinedTermSet": "NCIT",
                    }
                )
            if variable_measured:
                output["variableMeasured"] = variable_measured

            donor_desc = _build_donor_description(donor_mapped)

        sample_metadata = metadata.get("metadata")
        sample_desc = _first_text(metadata.get("description"))
        if not sample_desc and isinstance(sample_metadata, dict):
            sample_desc = (
                _first_text(sample_metadata.get("description"))
                or _first_text(sample_metadata.get("sample_description"))
                or _first_text(sample_metadata.get("sampleDescription"))
            )

        # TSV rule: description = "{sample.description} from {donor.description}" only when donor description is known.
        combined_desc = None
        if sample_desc and donor_desc:
            combined_desc = f"{sample_desc} from {donor_desc}"
        elif sample_desc:
            combined_desc = sample_desc
        elif donor_desc:
            combined_desc = donor_desc

        notes = None
        if isinstance(sample_metadata, dict):
            notes = sample_metadata.get("notes")
        notes = _first_text(notes)
        if notes:
            combined_desc = f"{combined_desc}. Notes: {notes}" if combined_desc else f"Notes: {notes}"

        if combined_desc:
            output["description"] = combined_desc

        # HuBMAP lineage (all via isBasisFor):
        # - Pre-processed sample isBasisFor post-processed sample
        # - Post-processed sample isBasisFor Dataset

        basis_for = []
        for desc in _ensure_list(metadata.get("immediate_descendants")):
            if not isinstance(desc, dict) or not desc.get("uuid"):
                continue
            if desc.get("entity_type") == "Sample":
                sid2 = desc["uuid"]
                entry = {
                    "@type": "Sample",
                    "identifier": sid2,
                    "url": f"https://portal.hubmapconsortium.org/browse/sample/{sid2}",
                }
                if desc_name := desc.get("hubmap_id"):
                    entry["name"] = desc_name
                basis_for.append(entry)
            elif desc.get("entity_type") == "Dataset":
                did = desc["uuid"]
                entry = {
                    "@type": "Dataset",
                    "identifier": did,
                    "url": f"https://portal.hubmapconsortium.org/browse/dataset/{did}",
                }
                if desc_name := desc.get("hubmap_id"):
                    entry["name"] = desc_name
                basis_for.append(entry)

        if basis_for:
            output["isBasisFor"] = basis_for

        yield output

    logger.info("Finished parsing %s samples", count)


def parse_datasets():
    datasets = get_ids()
    count = 0
    for metadata in datasets:
        count += 1
        if count % 50 == 0:
            logger.info("Parsed %s out of %s datasets", count, len(datasets))

        homo_sapiens = {
            "identifier": "9606",
            "inDefinedTermSet": "UniProt",
            "url": "https://www.uniprot.org/taxonomy/9606",
            "originalName": "homo sapiens",
            "isCurated": True,
            "curatedBy": {
                "name": "PubTator",
                "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
                "dateModified": "2023-10-05",
            },
            "name": "Homo sapiens",
            "commonName": "Human",
            "displayName": "Human | Homo sapiens",
            "alternateName": [
                "Human",
                "Homo sapiens Linnaeus, 1758",
                "human",
                "Home sapiens",
                "Homo sampiens",
                "Homo sapeins",
                "Homo sapian",
                "Homo sapians",
                "Homo sapien",
                "Homo sapience",
                "Homo sapiense",
                "Homo sapients",
                "Homo sapines",
                "Homo spaiens",
                "Homo spiens",
                "Humo sapiens",
            ],
            "classification": "host",
        }

        output = {
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "HuBMAP",
                "url": "https://hubmapconsortium.org/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Dataset",
            "species": homo_sapiens,
        }

        keywords = []
        if anatomy_1 := metadata.get("anatomy_0"):
            for term in anatomy_1:
                keywords.append(term)
        if anatomy_2 := metadata.get("anatomy_1"):
            for term in anatomy_2:
                keywords.append(term)
        if anatomy_3 := metadata.get("anatomy_2"):
            for term in anatomy_3:
                keywords.append(term)
        if display_subtype := metadata.get("display_subtype"):
            keywords.append(display_subtype)
        if len(keywords):
            output["keywords"] = keywords

        author_list = []
        if contacts := metadata.get("contacts"):
            for contact in contacts:
                author = {}
                if affiliation := contact.get("affiliation"):
                    author["affiliation"] = {"name": affiliation}
                if first_name := contact.get("first_name"):
                    author["givenName"] = first_name
                if last_name := contact.get("last_name"):
                    author["familyName"] = last_name
                if middle_name := contact.get("middle_name_or_initial"):
                    author["givenName"] += f" {middle_name}"
                if name := contact.get("name"):
                    author["name"] = name
                if orcid_id := contact.get("orcid_id"):
                    author["identifier"] = orcid_id
                if bool(author):
                    author_list.append(author)

        if contributors := metadata.get("contributors"):
            for contributor in contributors:
                author = {}
                if affiliation := contributor.get("affiliation"):
                    author["affiliation"] = {"name": affiliation}
                if first_name := contributor.get("first_name"):
                    author["givenName"] = first_name
                if last_name := contributor.get("last_name"):
                    author["familyName"] = last_name
                if middle_name := contributor.get("middle_name_or_initial"):
                    author["givenName"] += f" {middle_name}"
                if name := contributor.get("name"):
                    author["name"] = name
                if orcid_id := contributor.get("orcid_id"):
                    author["identifier"] = orcid_id
                if bool(author):
                    author_list.append(author)

        if len(author_list):
            output["author"] = author_list

        if created_timestamp := metadata.get("created_timestamp"):
            output["dateCreated"] = datetime.datetime.utcfromtimestamp(created_timestamp / 1000).strftime("%Y-%m-%d")

        if data_access_level := metadata.get("data_access_level"):
            if data_access_level == "public":
                output["isAccessibleForFree"] = True
            else:
                output["isAccessibleForFree"] = False

        measurement_technique = {}
        if data_types := metadata.get("data_types"):
            measurement_technique["name"] = data_types[0]
        if dataset_info := metadata.get("dataset_info"):
            measurement_technique["description"] = dataset_info
        if "name" in measurement_technique:
            output["measurementTechnique"] = measurement_technique

        if doi_url := metadata.get("doi_url"):
            output["doi"] = doi_url.removeprefix("https://doi.org/")

        if uuid := metadata.get("uuid"):
            url = f"https://portal.hubmapconsortium.org/browse/dataset/{uuid}"
            output["url"] = url
            output["includedInDataCatalog"]["archivedAt"] = url
            output["_id"] = "HUBMAP_" + uuid

        if files := metadata.get("files"):
            distribution_list = []
            for file in files:
                distribution_dict = {}
                if file_description := file.get("description"):
                    distribution_dict["name"] = file_description
                # if edam_term := file.get('edam_term'):
                #     distribution_dict['encodingFormat'] = edam_term
                if mapped_description := file.get("mapped_description"):
                    if "name" in distribution_dict:
                        distribution_dict["description"] = mapped_description
                    else:
                        distribution_dict["name"] = mapped_description
                if rel_path := file.get("rel_path"):
                    distribution_dict["contentUrl"] = f"https://assets.hubmapconsortium.org/{uuid}/{rel_path}"
                if size := file.get("size"):
                    distribution_dict["contentSize"] = size
                if file_type := file.get("type"):
                    if file_type != "unknown":
                        distribution_dict["encodingFormat"] = file_type
                if bool(distribution_dict):
                    distribution_list.append(distribution_dict)
            if len(distribution_list):
                output["distribution"] = distribution_list

        if hubmap_id := metadata.get("hubmap_id"):
            output["name"] = hubmap_id

        if version := metadata.get("version"):
            output["version"] = version

        if last_modified_timestamp := metadata.get("last_modified_timestamp"):
            output["dateModified"] = datetime.datetime.utcfromtimestamp(last_modified_timestamp / 1000).strftime(
                "%Y-%m-%d"
            )

        if published_timestamp := metadata.get("published_timestamp"):
            output["datePublished"] = datetime.datetime.utcfromtimestamp(published_timestamp / 1000).strftime(
                "%Y-%m-%d"
            )

        if dataset_metadata := metadata.get("metadata"):
            if origin := dataset_metadata.get("dag_provenance_list"):
                output["isBasedOn"] = {"name": origin[0]["origin"]}
            if nested_metadata := dataset_metadata.get("metadata"):
                if protocols_io_doi := nested_metadata.get("protocols_io_doi"):
                    if "isBasedOn" in output:
                        output["isBasedOn"]["doi"] = protocols_io_doi
                    else:
                        output["isBasedOn"] = {"doi": protocols_io_doi}
                if section_prep_protocols_io_doi := nested_metadata.get("section_prep_protocols_io_doi"):
                    if "isBasedOn" in output:
                        output["isBasedOn"]["doi"] = section_prep_protocols_io_doi
                    else:
                        output["isBasedOn"] = {"doi": section_prep_protocols_io_doi}
                # if analyte_class := nested_metadata.get('analyte_class'):
                #     output['variableMeasured'] = analyte_class
                # if assay_category := nested_metadata.get('assay_category'):
                #     output['measurementTechnique']['name'] = assay_category
                # if antibodies_path := nested_metadata.get('antibodies_path'):
                #     print(antibodies_path)
                #     print(output['url'])
                #     if 'distribution' in output:
                #         print(output['distribution'])

        # Minimal Sample info on Dataset records (derived from donor + one related sample).
        # This supports UI layouts that conditionally render sample sections even when info is sparse.
        try:
            sample_entry = {"@type": "Sample"}

            donor = metadata.get("donor") if isinstance(metadata, dict) else None
            donor_mapped = donor.get("mapped_metadata", {}) if isinstance(donor, dict) else {}
            if isinstance(donor_mapped, dict) and donor_mapped:
                sex_value = _first_scalar(donor_mapped.get("sex"))
                if sex_value:
                    sample_entry["sex"] = sex_value

                age_value = _first_scalar(donor_mapped.get("age_value"))
                age_unit = _first_scalar(donor_mapped.get("age_unit"))
                if age_value is not None and age_unit:
                    sample_entry["developmentalStage"] = {
                        "@type": "QuantitativeValue",
                        "value": age_value,
                        "unitText": age_unit,
                    }

                phenotypes = []
                race_value = _first_scalar(donor_mapped.get("race"))
                if race_value:
                    phenotypes.append({"@type": "DefinedTerm", "name": race_value})
                blood_type = _first_scalar(donor_mapped.get("abo_blood_group_system"))
                if blood_type:
                    phenotypes.append({"@type": "DefinedTerm", "name": blood_type})
                if phenotypes:
                    sample_entry["associatedPhenotype"] = phenotypes

            sample_ids = []
            for anc in _ensure_list(metadata.get("ancestors")):
                if isinstance(anc, dict) and anc.get("entity_type") == "Sample" and anc.get("uuid"):
                    sample_ids.append(anc["uuid"])

            sample_payload = _fetch_browse_json("sample", sample_ids[0]) if sample_ids else None
            if isinstance(sample_payload, dict):
                sample_category = (
                    sample_payload.get("mapped_sample_category")
                    or sample_payload.get("sample_category")
                    or sample_payload.get("display_subtype")
                )
                if sample_category:
                    sample_entry["sampleType"] = [{"@type": "DefinedTerm", "name": str(sample_category)}]

                organs = []
                for organ in _ensure_list(sample_payload.get("origin_samples_unique_mapped_organs")):
                    organ_value = _first_scalar(organ)
                    if organ_value:
                        organs.append(organ_value)
                if not organs:
                    for origin in _ensure_list(sample_payload.get("origin_samples")):
                        if isinstance(origin, dict):
                            organ_value = origin.get("mapped_organ") or origin.get("organ")
                            organ_value = _first_scalar(organ_value)
                            if organ_value:
                                organs.append(organ_value)

                if organs:
                    # De-dupe while preserving order
                    seen = set()
                    anatomical = []
                    for organ in organs:
                        key = str(organ).casefold()
                        if key in seen:
                            continue
                        seen.add(key)
                        anatomical.append({"@type": "DefinedTerm", "name": str(organ)})
                    if anatomical:
                        sample_entry["anatomicalStructure"] = anatomical

            if len(sample_entry.keys()) > 1:
                output["sample"] = sample_entry
        except Exception:
            pass

        if title := metadata.get("title"):
            output["description"] = title

        yield output

    logger.info("Finished parsing %s datasets", count)


def parse():
    yield from parse_datasets()
    yield from parse_samples()
