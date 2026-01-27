import datetime
import os
import traceback

import dateutil.parser

from .sample_charateristics import parse_sex

try:
    from config import logger
except ImportError:
    import logging

    logger = logging.getLogger(__name__)


def insert_value(d, key, value):
    if key in d:
        if isinstance(d[key], list) and value not in d[key]:
            d[key].append(value)
        if not isinstance(d[key], list) and d[key] != value:
            d[key] = [d[key], value]
    else:
        d[key] = value


def get_full_name(name):
    parts = name.split(",")
    while len(parts) < 3:
        parts.append("")
    if len(parts) > 3:
        parts = [parts[0], "".join(parts[1:-1]), parts[-1]]
    first, middle, last = parts
    full_name = " ".join([part.strip() for part in [first, middle, last] if part.strip()])
    return full_name


def parse_soft_series(filepath):
    """
    Parse a GEO SOFT series file into a dictionary.
    Each key is the SOFT field (e.g., '!Series_title'), value is a list if repeated, or a string.
    """
    result = {}
    with open(filepath, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("^"):
                continue
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()
                # Store as list if key repeats
                if key in result:
                    if isinstance(result[key], list):
                        result[key].append(value)
                    else:
                        result[key] = [result[key], value]
                else:
                    result[key] = value
    return result


def parse_sample_characteristics(output, value):
    sample_mapping = {
        "tissue": ("anatomicalStructure", "field_value"),
        "tissue id": ("anatomicalStructure", "field_value"),
        "organ": ("anatomicalStructure", "field_value"),
        "tissue type": ("anatomicalStructure", "field_value"),
        "genotype": ("associatedGenotype", "field_value"),
        "genotype/variation": ("associatedGenotype", "field_value"),
        "ctnnb1 genotype": ("associatedGenotype", "field_value"),
        "mouse line": ("associatedGenotype", "field_value"),
        "phenotype": ("associatedPhenotype", "field_value"),
        "cell type": ("cellType", "field_value"),
        "cell line": ("cellType", "field_value"),
        "cell subset": ("cellType", "field_value"),
        "age": ("developmentalStage", "field_value"),
        "time of_larval_development": ("developmentalStage", "field_value"),
        "age group": ("developmentalStage", "field_value"),
        "age (months)": ("developmentalStage", "field_value"),
        "gestational age": ("developmentalStage", "field_value"),
        "pistil development stage": ("developmentalStage", "field_value"),
        "developmental stage": ("developmentalStage", "field_value"),
        "sequencing run": ("experimentalPurpose", "subproperty"),
        "disease": ("healthCondition", "field_value"),
        "diagnosis": ("healthCondition", "field_value"),
        "disease state": ("healthCondition", "field_value"),
        "condition": ("healthCondition", "field_value"),
        "injury": ("healthCondition", "field_value"),
        "treatment": ("sampleProcess", "field_value"),
        "sex": ("sex", "field_value"),
        "cultivar": ("species", "field_value"),
        "strain": ("species", "field_value"),
        "time": ("temporalCoverage", "field_value"),
        "timepoint": ("temporalCoverage", "field_value"),
        "time of_treatment": ("temporalCoverage", "field_value"),
        "time point": ("temporalCoverage", "field_value"),
        "treatment/time point": ("temporalCoverage", "field_value"),
        "time(days)": ("temporalCoverage", "field_value"),
    }
    values = value if isinstance(value, list) else [value]
    for v in values:
        # Split by ":", strip whitespace, and only split on the first ":"
        parts = [p.strip() for p in v.split(":", 1)]
        if len(parts) != 2:
            logger.warning(f"Invalid sample characteristic format: {v}")
            continue

        subproperty, field_value = parts[0], parts[1]

        sex = parse_sex(subproperty, field_value)
        if isinstance(sex, tuple):
            if sex[0] and isinstance(sex[0], list):
                [insert_value(output, "sex", s) for s in sex[0]]
            else:
                insert_value(output, "sex", sex[0])
            if sex[1]:
                insert_value(output, "developmentalStage", sex[1])
        if sex and isinstance(sex, list):
            [insert_value(output, "sex", s) for s in sex]
        elif sex:
            insert_value(output, "sex", sex)
        else:
            continue

        if subproperty.lower() in sample_mapping:
            mapping = sample_mapping[subproperty.lower()]
            if mapping[0] in [
                "species",
                "healthCondition",
                "developmentalStage",
                "cellType",
                "anatomicalStructure",
                "associatedPhenotype",
            ]:
                if not mapping[0] in output:
                    output[mapping[0]] = []
                if isinstance(field_value, list):
                    for fv in field_value:
                        d = {"name": fv}
                        if "anatomicalStructure" == mapping[0]:
                            d["@type"] = "DefinedTerm"
                        if not d in output[mapping[0]]:
                            output[mapping[0]].append(d)
                else:
                    d = {"name": field_value}
                    if "anatomicalStructure" == mapping[0]:
                        d["@type"] = "DefinedTerm"
                    if not d in output[mapping[0]]:
                        output[mapping[0]].append(d)
            elif mapping[0] == "temporalCoverage":
                if isinstance(field_value, list):
                    output[mapping[0]] = [{"duration": fv} for fv in field_value]
                else:
                    output[mapping[0]] = [{"duration": field_value}]
            else:
                if mapping[0] == "sampleProcess":
                    if "sampleProcess" in output and output["sampleProcess"]:
                        # If sampleProcess already exists, append new info
                        output["sampleProcess"] = output["sampleProcess"] + " " + field_value
                    else:
                        output["sampleProcess"] = field_value
                else:
                    if not mapping[0] in output:
                        output[mapping[0]] = []

                    if mapping[1] == "subproperty":
                        output[mapping[0]].append(subproperty)
                    else:
                        output[mapping[0]].append(field_value)
        else:
            if not "variableMeasured" in output:
                output["variableMeasured"] = []
            d = {"name": subproperty, "@type": "DefinedTerm"}
            if not d in output["variableMeasured"]:
                output["variableMeasured"].append(d)


def parse_gsm(data_folder):
    """
    Parse a GEO SOFT platform file into a dictionary.
    Each key is the SOFT field (e.g., '!Platform_title'), value is a list if repeated, or a string.
    """

    gse_dir = os.path.join(data_folder, "gsm")
    records = get_records(gse_dir)
    for item in records:
        if not item.get("!Sample_geo_accession"):
            continue
        _id = item.get("!Sample_geo_accession")
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + _id
        output = {
            "@context": "http://schema.org/",
            "@type": "Sample",
            "_id": _id.casefold(),
            "identifier": _id,
            "url": url,
            "distribution": [{"@type": "dataDownload", "contentUrl": url}],
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "NCBI GEO",
                "url": "https://www.ncbi.nlm.nih.gov/geo/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
        }

        if name := item.get("!Sample_title"):
            output["name"] = name

        if date_published := item.get("!Sample_status"):
            date_str = date_published.replace("Public on ", "")
            try:
                dt = dateutil.parser.parse(date_str, ignoretz=True).date().isoformat()
                output["datePublished"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_str}': {e}")

        if date_created := item.get("!Sample_submission_date"):
            try:
                dt = dateutil.parser.parse(date_created, ignoretz=True).date().isoformat()
                output["dateCreated"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_created}': {e}")

        if date_modified := item.get("!Sample_last_update_date"):
            try:
                dt = dateutil.parser.parse(date_modified, ignoretz=True).date().isoformat()
                output["dateModified"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_modified}': {e}")

        if sample_type := item.get("!Sample_type"):
            if isinstance(sample_type, list):
                output["sampleType"] = [{"name": s, "@type": "DefinedTerm"} for s in sample_type]
            else:
                output["sampleType"] = [{"name": sample_type, "@type": "DefinedTerm"}]

        if description := item.get("!Sample_description"):
            if isinstance(description, list):
                output["description"] = " ".join(description)
            else:
                output["description"] = description

        if sample_process := item.get("!Sample_data_processing"):
            if isinstance(sample_process, list):
                output["sampleProcess"] = " ".join(sample_process)
            else:
                output["sampleProcess"] = sample_process

        if author := item.get("!Sample_contact_name"):
            output["author"] = {"@type": "Person", "name": get_full_name(author)}
            if affiliation := item.get("!Sample_contact_institute"):
                output["author"]["affilliation"] = {"name": affiliation}

        if instrument := item.get("!Sample_instrument_model"):
            if isinstance(instrument, list):
                output["instrument"] = [{"name": i} for i in instrument]
            else:
                output["instrument"] = {"name": instrument}

        mts = []

        if mt := item.get("!Sample_library_selection"):
            if isinstance(mt, list):
                for m in mt:
                    mts.append({"name": m})
            else:
                mts.append({"name": mt})

        if mt := item.get("!Sample_library_strategy"):
            if isinstance(mt, list):
                for m in mt:
                    mts.append({"name": m})
            else:
                mts.append({"name": mt})

        if mts:
            output["measurementTechnique"] = mts

        if sample_type := item.get("!Sample_library_source"):
            if "sampleType" not in output:
                output["sampleType"] = []
            if isinstance(sample_type, list):
                output["sampleType"].extend([{"name": s, "@type": "DefinedTerm"} for s in sample_type])
            else:
                output["sampleType"].append({"name": sample_type, "@type": "DefinedTerm"})

        if same_as := item.get("!Sample_relation"):
            output["sameAs"] = same_as

        if is_basis_for := item.get("!Sample_series_id"):
            if isinstance(is_basis_for, list):
                output["isBasisFor"] = [
                    {
                        "@type": "Dataset",
                        "identifier": sid,
                        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + sid,
                    }
                    for sid in is_basis_for
                ]
            else:
                output["isBasisFor"] = {
                    "@type": "Dataset",
                    "identifier": is_basis_for,
                    "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + is_basis_for,
                }

        for key, value in item.items():
            if key.startswith("!Sample_characteristics"):
                parse_sample_characteristics(output, value)

            if "protocol" in key.lower() or key.startswith("!Sample_molecule"):
                if isinstance(value, list):
                    sample_process = " ".join(value)
                else:
                    sample_process = value
                if "sampleProcess" in output and output["sampleProcess"]:
                    # If sampleProcess already exists, append new info
                    output["sampleProcess"] = output["sampleProcess"] + " " + sample_process
                else:
                    output["sampleProcess"] = sample_process

            if key.startswith("!Sample_supplementary_file"):
                if isinstance(value, list):
                    for v in value:
                        insert_value(output, "distribution", {"@type": "DataDownload", "contentUrl": v})
                else:
                    insert_value(output, "distribution", {"@type": "DataDownload", "contentUrl": value})

        yield output


def get_records(data_folder):
    for root, dirs, _ in os.walk(data_folder):
        for dir in dirs:
            logger.info(dir)
            dirpath = os.path.join(root, dir)
            for file in os.listdir(dirpath):
                logger.info(file)
                if file.endswith(".txt"):
                    fpath = os.path.join(dirpath, file)
                    try:
                        item = parse_soft_series(fpath)
                        yield item
                    except Exception as e:
                        logger.error(f"Error parsing {fpath}: {e}")


def find_gsm_file(data_folder, acc):
    # Extract prefix (GSE/GSM) and numeric part
    prefix = acc[:3]
    num = acc[3:]
    # Pad numeric part to at least 3 digits for nnn, or 6 for full
    padded = num.zfill(6)
    subdir = prefix + padded[:3] + "nnn"
    subdir_path = os.path.join(data_folder, subdir)
    file_path = os.path.join(subdir_path, f"{acc}.txt")
    if os.path.exists(file_path):
        return file_path
    else:
        logger.warning(f"GSM file not found: {file_path}")
        return None


def parse_series_sample_characteristics(sample_elements, value):
    sample_mapping = {
        "tissue": ("anatomicalStructure", "field_value"),
        "tissue id": ("anatomicalStructure", "field_value"),
        "organ": ("anatomicalStructure", "field_value"),
        "tissue type": ("anatomicalStructure", "field_value"),
        "genotype": ("associatedGenotype", "field_value"),
        "genotype/variation": ("associatedGenotype", "field_value"),
        "ctnnb1 genotype": ("associatedGenotype", "field_value"),
        "mouse line": ("associatedGenotype", "field_value"),
        "phenotype": ("associatedPhenotype", "field_value"),
        "cell type": ("cellType", "field_value"),
        "cell line": ("cellType", "field_value"),
        "cell subset": ("cellType", "field_value"),
        "age": ("developmentalStage", "field_value"),
        "time of_larval_development": ("developmentalStage", "field_value"),
        "age group": ("developmentalStage", "field_value"),
        "age (months)": ("developmentalStage", "field_value"),
        "gestational age": ("developmentalStage", "field_value"),
        "pistil development stage": ("developmentalStage", "field_value"),
        "developmental stage": ("developmentalStage", "field_value"),
    }
    values = value if isinstance(value, list) else [value]
    for v in values:
        # Split by ":", strip whitespace, and only split on the first ":"
        parts = [p.strip() for p in v.split(":", 1)]
        if len(parts) != 2:
            logger.warning(f"Invalid sample characteristic format: {v}")
            continue

        subproperty, field_value = parts[0].casefold(), parts[1].casefold()

        if subproperty in sample_mapping:
            mapping = sample_mapping[subproperty]
            if mapping[0] in [
                "developmentalStage",
                "cellType",
                "anatomicalStructure",
                "associatedPhenotype",
            ]:

                if isinstance(field_value, list):
                    for fv in field_value:
                        d = {"name": fv}
                        if "anatomicalStructure" == mapping[0]:
                            d["@type"] = "DefinedTerm"
                        insert_value(sample_elements, mapping[0], d)
                else:
                    d = {"name": field_value}
                    if "anatomicalStructure" == mapping[0]:
                        d["@type"] = "DefinedTerm"
                    insert_value(sample_elements, mapping[0], d)
            else:
                if mapping[1] == "subproperty":
                    insert_value(sample_elements, mapping[0], subproperty)
                else:
                    if isinstance(field_value, list):
                        for fv in field_value:
                            insert_value(sample_elements, mapping[0], fv)
                    else:
                        insert_value(sample_elements, mapping[0], field_value)


def parse_series_sample(item, sample):
    """
    Given a GEO Series item and an output dictionary, parse for the sample field in the dataset
    """
    if not item.get("!Sample_geo_accession"):
        return
    sample_elements = sample["aggregateElement"]
    if sample_type := item.get("!Sample_type"):
        if isinstance(sample_type, list):
            for s in sample_type:
                insert_value(sample_elements, "sampleType", {"name": s, "@type": "DefinedTerm"})
        else:
            insert_value(sample_elements, "sampleType", {"name": sample_type, "@type": "DefinedTerm"})

    if sample_type := item.get("!Sample_library_source"):
        if isinstance(sample_type, list):
            for s in sample_type:
                insert_value(sample_elements, "sampleType", {"name": s, "@type": "DefinedTerm"})
        else:
            insert_value(sample_elements, "sampleType", {"name": sample_type, "@type": "DefinedTerm"})

    for key, value in item.items():
        if key.startswith("!Sample_characteristics"):
            parse_series_sample_characteristics(sample_elements, value)


def parse_gse(data_folder):
    """
    Parse a GEO SOFT platform file into a dictionary.
    Each key is the SOFT field (e.g., '!Platform_title'), value is a list if repeated, or a string.
    """

    gse_dir = os.path.join(data_folder, "gse")
    records = get_records(gse_dir)
    for item in records:
        if not item.get("!Series_geo_accession"):
            continue
        _id = item.get("!Series_geo_accession")
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + _id
        output = {
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": _id.casefold(),
            "identifier": _id,
            "url": url,
            "distribution": {"@type": "dataDownload", "contentUrl": url},
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "NCBI GEO",
                "url": "https://www.ncbi.nlm.nih.gov/geo/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
        }

        if gsm_ids := item.get("!Series_sample_id"):
            sample = {
                "@type": "SampleCollection",
                "itemListElement": [],
                "aggregateElement": {},
                "numberOfItems": {"value": 0, "unitText": "sample"},
            }
            gsm_dir = os.path.join(data_folder, "gsm")
            if not isinstance(gsm_ids, list):
                gsm_ids = [gsm_ids]
            for gsm_id in gsm_ids:
                sample["itemListElement"].append(
                    {
                        "@type": "Sample",
                        "identifier": gsm_id,
                        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + gsm_id,
                        "_id": gsm_id.casefold(),
                    }
                )
                sample["numberOfItems"]["value"] += 1
                gsm_file = find_gsm_file(gsm_dir, gsm_id)
                try:
                    sample_item = parse_soft_series(gsm_file)
                    parse_series_sample(sample_item, sample)
                except Exception as e:
                    logger.error(f"Error parsing GSM file {gsm_file}: {e}")
                    logger.error(traceback.format_exc())
                    continue

            if sample.get("itemListElement"):
                # logger.info(f"gse_id: {output['_id']}")
                # logger.info(f"sample: {sample}")
                output["sample"] = sample

        if name := item.get("!Series_title"):
            output["name"] = name

        if date_published := item.get("!Series_status"):
            date_str = date_published.replace("Public on ", "")
            try:
                dt = dateutil.parser.parse(date_str, ignoretz=True).date().isoformat()
                output["datePublished"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_str}': {e}")

        if date_created := item.get("!Series_submission_date"):
            try:
                dt = dateutil.parser.parse(date_created, ignoretz=True).date().isoformat()
                output["dateCreated"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_created}': {e}")

        if date_modified := item.get("!Series_last_update_date"):
            try:
                dt = dateutil.parser.parse(date_modified, ignoretz=True).date().isoformat()
                output["dateModified"] = dt
            except Exception as e:
                logger.warning(f"Error parsing date '{date_modified}': {e}")

        if pmids := item.get("!Series_pubmed_id"):
            if isinstance(pmids, list):
                output["pmids"] = ",".join(pmids)
            else:
                output["pmids"] = pmids

        if descriptions := item.get("!Series_summary"):
            if isinstance(descriptions, list):
                output["description"] = " ".join(descriptions)
            else:
                output["description"] = descriptions

        if description := item.get("!Series_overall_design"):
            if isinstance(description, list):
                output["description"] = (output.get("description", "") + " " + " ".join(description)).strip()
            else:
                output["description"] = (output.get("description", "") + " " + description).strip()

        if mt := item.get("!Series_type"):
            if isinstance(mt, list):
                output["measurementTechnique"] = [{"name": m} for m in mt]
            else:
                output["measurementTechnique"] = {"name": mt}

        if authors := item.get("!Series_contributor"):
            if isinstance(authors, list):
                output["author"] = [{"@type": "Person", "name": get_full_name(a)} for a in authors]
            else:
                output["author"] = [{"@type": "Person", "name": get_full_name(authors)}]

        if publishers := item.get("!Series_contact_institute"):
            if isinstance(publishers, list):
                output["sdPublisher"] = [{"@type": "Organization", "name": p} for p in publishers]
            else:
                output["sdPublisher"] = [{"@type": "Organization", "name": publishers}]

        if species := item.get("!Series_platform_organism"):
            if isinstance(species, list):
                output["species"] = [{"name": s} for s in species]
            else:
                output["species"] = [{"name": species}]

        yield output
