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


def _search_entity_ids(entity_type, page_size=5000):
    """uuids of every current (latest revision, no sub_status) entity of `entity_type`.

    A single search returns at most 10,000 hits, so this pages with search_after.
    """
    url = "https://search.api.hubmapconsortium.org/v3/portal/search"
    payload = {
        "query": {
            "bool": {
                "must_not": [
                    {"exists": {"field": "next_revision_uuid"}},
                    {"exists": {"field": "sub_status"}},
                ],
                "filter": [{"term": {"entity_type.keyword": entity_type}}],
            }
        },
        "_source": False,
        "size": page_size,
        "sort": [{"uuid.keyword": "asc"}],
    }
    ids = []
    while True:
        resp = requests.post(url, json=payload, timeout=60)
        resp.raise_for_status()
        hits = resp.json()["hits"]["hits"]
        ids.extend(hit["_id"] for hit in hits)
        if len(hits) < page_size:
            return ids
        payload["search_after"] = hits[-1]["sort"]


def get_ids():
    logger.info("Getting datasets")
    ids = _search_entity_ids("Dataset")
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
    ids = _search_entity_ids("Sample")
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


def _sample_text(metadata):
    """A sample's own (description, notes)."""
    sample_metadata = metadata.get("metadata")
    description = _first_text(metadata.get("description"))
    notes = None
    if isinstance(sample_metadata, dict):
        description = (
            description
            or _first_text(sample_metadata.get("description"))
            or _first_text(sample_metadata.get("sample_description"))
            or _first_text(sample_metadata.get("sampleDescription"))
        )
        notes = _first_text(sample_metadata.get("notes"))
    return description, notes


def _anatomical_structures(metadata):
    """Organ DefinedTerms from an entity's origin samples, de-duped in order."""
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

    seen = set()
    anatomical = []
    for organ in organs:
        key = str(organ).casefold()
        if key in seen:
            continue
        seen.add(key)
        anatomical.append({"@type": "DefinedTerm", "name": str(organ)})
    return anatomical


def _index_source_samples(datasets):
    """Map each sample uuid to isBasisFor entries for the crawled Datasets that list it in `source_samples`."""
    datasets_by_sample = {}
    for metadata in datasets:
        if not isinstance(metadata, dict) or not (uuid := metadata.get("uuid")):
            continue
        entry = {
            "@type": "Dataset",
            "identifier": metadata.get("hubmap_id") or uuid,
            "url": f"https://portal.hubmapconsortium.org/browse/dataset/{uuid}",
        }
        if name := metadata.get("hubmap_id"):
            entry["name"] = name
        for source in _ensure_list(metadata.get("source_samples")):
            if isinstance(source, dict) and source.get("uuid"):
                entries = datasets_by_sample.setdefault(source["uuid"], [])
                if entry not in entries:
                    entries.append(entry)
    return datasets_by_sample


def _build_sample_collection(metadata):
    """Flatten a Dataset's provenance into a SampleCollection of its source samples.

    HuBMAP chains Donor -> Organ -> Block -> Section/Suspension -> Dataset -> processed Dataset.
    `source_samples` is HuBMAP's resolution of that chain to the samples that were assayed, so it
    also covers processed datasets, and it is what numberOfItems counts. The organ and donor
    upstream of those samples are summarized in aggregateElement.
    """
    items = []
    sample_types = []
    seen_samples = set()
    seen_types = set()
    for source in _ensure_list(metadata.get("source_samples")):
        if not isinstance(source, dict) or not source.get("uuid") or source["uuid"] in seen_samples:
            continue
        seen_samples.add(source["uuid"])
        items.append(
            {
                "@type": "Sample",
                "identifier": source.get("hubmap_id") or source["uuid"],
                "url": f"https://portal.hubmapconsortium.org/browse/sample/{source['uuid']}",
            }
        )
        sample_category = _first_text(source.get("mapped_sample_category") or source.get("sample_category"))
        if sample_category and str(sample_category).casefold() not in seen_types:
            seen_types.add(str(sample_category).casefold())
            sample_types.append({"@type": "DefinedTerm", "name": str(sample_category)})

    aggregate = {"@type": "Sample", "additionalType": "ExperimentalRunSample"}
    donor = metadata.get("donor")
    donor_mapped = donor.get("mapped_metadata") if isinstance(donor, dict) else None
    if isinstance(donor_mapped, dict) and donor_mapped:
        sex_value = _first_text(donor_mapped.get("sex"))
        if sex_value:
            aggregate["sex"] = sex_value

        age_value = _first_scalar(donor_mapped.get("age_value"))
        age_unit = _first_scalar(donor_mapped.get("age_unit"))
        if age_value is not None and age_unit:
            aggregate["developmentalStage"] = {
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
            aggregate["associatedPhenotype"] = phenotypes

    if sample_types:
        aggregate["sampleType"] = sample_types
    if anatomical := _anatomical_structures(metadata):
        aggregate["anatomicalStructure"] = anatomical

    collection = {"@type": "SampleCollection"}
    if items:
        collection["itemListElement"] = items
        collection["numberOfItems"] = {"@type": "QuantitativeValue", "value": len(items), "unitText": "sample"}
    if len(aggregate) > 2:
        collection["aggregateElement"] = aggregate
    return collection if len(collection) > 1 else None


def parse_samples(datasets_by_sample):
    # Only samples a crawled Dataset was assayed on become records; their upstream samples are folded in.
    sample_ids = [sid for sid in get_sample_ids() if sid in datasets_by_sample]
    logger.info("Assayed samples: %s", len(sample_ids))
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
                "name": "Human Biomolecular Atlas Program (HuBMAP)",
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

        if anatomical := _anatomical_structures(metadata):
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
                    health_conditions.append({"@type": "DefinedTerm", "name": str(condition_value)})

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

        sample_desc, notes = _sample_text(metadata)

        # TSV rule: description = "{sample.description} from {donor.description}" only when donor description is known.
        combined_desc = None
        if sample_desc and donor_desc:
            combined_desc = f"{sample_desc} from {donor_desc}"
        elif sample_desc:
            combined_desc = sample_desc
        elif donor_desc:
            combined_desc = donor_desc

        if notes:
            combined_desc = f"{combined_desc}. Notes: {notes}" if combined_desc else f"Notes: {notes}"

        # Upstream samples (organ, block, ...) aren't records of their own, so their descriptions and notes
        # are carried on each assayed sample, labelled with their sample category.
        seen_text = {text.rstrip(" .").casefold() for text in (sample_desc, notes) if text}
        for ancestor in _ensure_list(metadata.get("ancestors")):
            if not isinstance(ancestor, dict) or ancestor.get("entity_type") != "Sample" or not ancestor.get("uuid"):
                continue
            ancestor_metadata = _fetch_browse_json("sample", ancestor["uuid"], timeout=30)
            if not isinstance(ancestor_metadata, dict):
                continue
            ancestor_desc, ancestor_notes = _sample_text(ancestor_metadata)
            label = (
                ancestor_metadata.get("mapped_sample_category") or ancestor_metadata.get("sample_category") or "Sample"
            )
            for text, prefix in ((ancestor_desc, f"{label}: "), (ancestor_notes, f"{label} notes: ")):
                key = text.rstrip(" .").casefold() if text else None
                if key and key not in seen_text:
                    seen_text.add(key)
                    upstream = prefix + text.rstrip(" .")
                    combined_desc = f"{combined_desc.rstrip(' .')}. {upstream}" if combined_desc else upstream

        if combined_desc:
            output["description"] = combined_desc

        # Every record is an assayed sample. It isBasisFor each crawled Dataset that lists it in source_samples,
        # processed datasets included, which matches the Dataset's sample collection.
        output["additionalType"] = "ExperimentalRunSample"
        output["isBasisFor"] = datasets_by_sample[sid]

        yield output

    logger.info("Finished parsing %s samples", count)


def parse_datasets(datasets):
    count = 0
    for metadata in datasets:
        count += 1
        if count % 50 == 0:
            logger.info("Parsed %s out of %s datasets", count, len(datasets))

        homo_sapiens = {
            "@type": "DefinedTerm",
            "identifier": "9606",
            "inDefinedTermSet": "UniProt",
            "url": "https://www.uniprot.org/taxonomy/9606",
            "originalName": "homo sapiens",
            "isCurated": True,
            "curatedBy": {
                "@type": "SoftwareApplication",
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
                "name": "Human Biomolecular Atlas Program (HuBMAP)",
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
                author = {"@type": "Person"}
                if affiliation := contact.get("affiliation"):
                    author["affiliation"] = {"@type": "Organization", "name": affiliation}
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
                author = {"@type": "Person"}
                if affiliation := contributor.get("affiliation"):
                    author["affiliation"] = {"@type": "Organization", "name": affiliation}
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

        measurement_technique = {"@type": "DefinedTerm"}
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
                distribution_dict = {"@type": "DataDownload"}
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
            # Sample isBasisFor entries carry this identifier; the portal finds a Dataset's samples by it.
            output["identifier"] = hubmap_id
            output["name"] = hubmap_id

        # `version` is deliberately not emitted: check_schema rejects it on every
        # record, so the HuBMAP dataset version is dropped here.

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
                output["isBasedOn"] = {"@type": "CreativeWork", "name": origin[0]["origin"]}
            if nested_metadata := dataset_metadata.get("metadata"):
                if protocols_io_doi := nested_metadata.get("protocols_io_doi"):
                    if "isBasedOn" in output:
                        output["isBasedOn"]["doi"] = protocols_io_doi
                    else:
                        output["isBasedOn"] = {"@type": "ScholarlyArticle", "doi": protocols_io_doi}
                if section_prep_protocols_io_doi := nested_metadata.get("section_prep_protocols_io_doi"):
                    if "isBasedOn" in output:
                        output["isBasedOn"]["doi"] = section_prep_protocols_io_doi
                    else:
                        output["isBasedOn"] = {"@type": "ScholarlyArticle", "doi": section_prep_protocols_io_doi}
                # if analyte_class := nested_metadata.get('analyte_class'):
                #     output['variableMeasured'] = analyte_class
                # if assay_category := nested_metadata.get('assay_category'):
                #     output['measurementTechnique']['name'] = assay_category
                # if antibodies_path := nested_metadata.get('antibodies_path'):
                #     print(antibodies_path)
                #     print(output['url'])
                #     if 'distribution' in output:
                #         print(output['distribution'])

        # The Dataset lists its source samples directly (with a count) instead of leaving them to be
        # found through HuBMAP's Organ -> Block -> Section chain.
        if sample_collection := _build_sample_collection(metadata):
            output["sample"] = sample_collection

        if title := metadata.get("title"):
            output["description"] = title

        yield output

    logger.info("Finished parsing %s datasets", count)


def parse():
    datasets = get_ids()
    yield from parse_datasets(datasets)
    yield from parse_samples(_index_source_samples(datasets))
