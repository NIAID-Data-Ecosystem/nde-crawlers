import datetime
import json
import logging
import os
import re
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from time import sleep
from urllib.error import ContentTooShortError

import pandas as pd
import requests
from sql_database import NDEDatabase

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

INVALID_SAMPLE_VALUES = {
    "-", "?", "blank", "ignore", "missing", "na", "n/a", "nan", "none",
    "null", "not applicable", "not collected", "not determined",
    "not provided", "not recorded", "restricted access",
    "unspecified", "unk", "unknown",
}

# Maps SRA sample attribute tags to aggregateElement field names
AGGREGATE_ATTR_MAP = {
    "sex": "sex",
    "Sex": "sex",
    "cell_type": "cellType",
    "cell type": "cellType",
    "cell_line": "cellType",
    "cell line": "cellType",
    "cell line name": "cellType",
    "dev_stage": "developmentalStage",
    "development_stage": "developmentalStage",
    "developmental_stage": "developmentalStage",
    "genotype": "associatedGenotype",
    "strain": "associatedGenotype",
    "phenotype": "associatedPhenotype",
    "tissue": "sampleType",
    "tissue_type": "sampleType",
    "source_name": "sampleType",
    "sample_type": "sampleType",
    "isolate": "sampleType",
}

# Maps SRA sample attribute tags to dataset-level field names
DATASET_ATTR_MAP = {
    "disease": "healthCondition",
    "host_disease": "healthCondition",
    "disease_state": "healthCondition",
}


class NCBI_SRA(NDEDatabase):
    # override variables
    SQL_DB = "ncbi_sra.db"
    EXPIRE = datetime.timedelta(days=365)

    # Used for testing small chunks of data
    DATA_LIMIT = 5

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ncbi_api_key = os.environ.get("NCBI_API_KEY")
        if not self.ncbi_api_key:
            try:
                from config_local import GEO_API_KEY
                self.ncbi_api_key = GEO_API_KEY
            except ImportError:
                pass

    def _eutils_params(self, params):
        """Add standard eutils parameters including API key if available."""
        params["tool"] = "nde-crawlers"
        params["email"] = "nde@scripps.edu"
        if self.ncbi_api_key:
            params["api_key"] = self.ncbi_api_key
        return params

    def _request_delay(self):
        """Sleep to respect NCBI rate limits. 3 req/s without key, 10 req/s with key."""
        if self.ncbi_api_key:
            sleep(0.34)
        else:
            sleep(1.1)

    def _download_file(self, url, dest, max_retries=10):
        """Download a file with progress logging and retry logic."""
        retry_count = 0
        while True:
            try:
                logger.info(f"Starting download: {url}")
                with requests.get(url, stream=True, timeout=300) as r:
                    r.raise_for_status()
                    total = int(r.headers.get("content-length", 0))
                    if total:
                        logger.info(f"File size: {total / (1024**3):.2f} GB")
                    else:
                        logger.info("File size: unknown")
                    downloaded = 0
                    last_log = 0
                    log_interval = 500 * 1024 * 1024  # 500 MB
                    dl_start = time.time()
                    with open(dest, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                            f.write(chunk)
                            downloaded += len(chunk)
                            if downloaded - last_log >= log_interval:
                                elapsed = time.time() - dl_start
                                speed = downloaded / elapsed / (1024**2) if elapsed > 0 else 0
                                if total:
                                    pct = downloaded / total * 100
                                    logger.info(
                                        f"Download progress: {downloaded / (1024**3):.2f} / "
                                        f"{total / (1024**3):.2f} GB ({pct:.1f}%) - {speed:.1f} MB/s"
                                    )
                                else:
                                    logger.info(
                                        f"Download progress: {downloaded / (1024**3):.2f} GB - {speed:.1f} MB/s"
                                    )
                                last_log = downloaded
                elapsed = time.time() - dl_start
                logger.info(f"Download complete: {downloaded / (1024**3):.2f} GB in {elapsed:.0f}s")
                return
            except (ContentTooShortError, requests.exceptions.RequestException) as e:
                retry_count += 1
                if retry_count > max_retries:
                    raise
                logger.error(f"Download error: {e}. Retry {retry_count}/{max_retries}")
                sleep(retry_count * 2)

    @staticmethod
    def _xml_text(element, path, attr=None):
        """Safely extract text or attribute from an XML element."""
        try:
            el = element.find(path)
            if el is None:
                return ""
            if attr:
                return el.attrib.get(attr, "")
            return el.text or ""
        except (AttributeError, KeyError):
            return ""

    def _parse_efetch_xml(self, xml_text):
        """Parse efetch XML into list of run-level metadata dicts."""
        root = ET.fromstring(xml_text)
        results = []

        for pkg in root.findall("EXPERIMENT_PACKAGE"):
            d = {}
            d["experiment_accession"] = self._xml_text(pkg, "./EXPERIMENT", "accession")
            d["experiment_title"] = self._xml_text(pkg, "./EXPERIMENT/TITLE")

            lib = pkg.find("./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR")
            if lib is not None:
                d["library_strategy"] = self._xml_text(lib, "LIBRARY_STRATEGY")

            platform_el = pkg.find("./EXPERIMENT/PLATFORM")
            if platform_el is not None and len(platform_el) > 0:
                if len(platform_el[0]) > 0:
                    model_text = platform_el[0][0].text
                    if model_text:
                        d["model"] = model_text.strip()

            d["organisation"] = self._xml_text(pkg, "./Organization/Name")
            d["organisation_contact_email"] = self._xml_text(
                pkg, "./Organization/Contact", "email"
            )

            d["sample_accession"] = self._xml_text(pkg, "./SAMPLE", "accession")
            d["sample_title"] = self._xml_text(pkg, "./SAMPLE/TITLE")
            d["taxon_name"] = self._xml_text(
                pkg, "./SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME"
            )

            # Extract BioSample accession (SAMN ID) - resolvable unlike SRS IDs
            for ext_id in pkg.findall("./SAMPLE/IDENTIFIERS/EXTERNAL_ID"):
                if ext_id.attrib.get("namespace") == "BioSample" and ext_id.text:
                    d["biosample_accession"] = ext_id.text.strip()
                    break

            # Sample attributes (captures cell line, strain, HapMap fields, etc.)
            sample_attrs = pkg.find("./SAMPLE/SAMPLE_ATTRIBUTES")
            if sample_attrs is not None:
                for attr_el in sample_attrs:
                    tag = attr_el.findtext("TAG")
                    value = attr_el.findtext("VALUE")
                    if tag and value:
                        d[tag] = value

            d["study_title"] = self._xml_text(pkg, "./STUDY/DESCRIPTOR/STUDY_TITLE")
            d["experiment_desc"] = self._xml_text(
                pkg, "./EXPERIMENT/DESIGN/DESIGN_DESCRIPTION"
            )
            d["study_abstract"] = self._xml_text(
                pkg, "./STUDY/DESCRIPTOR/STUDY_ABSTRACT"
            )

            for run_el in pkg.findall("./RUN_SET/RUN"):
                d2 = dict(d)
                d2["run"] = run_el.attrib.get("accession", "")
                d2["published"] = run_el.attrib.get("published", "")
                results.append(d2)

        return results

    def _fetch_study_metadata(self, study_acc, bio_project=None):
        """Fetch run-level metadata for a study from NCBI eutils.

        Tries the study accession first, falls back to BioProject if available.
        """
        search_terms = [f"{study_acc}[Accession]"]
        if bio_project and bio_project != "-":
            search_terms.append(f"{bio_project}[BioProject]")

        for term in search_terms:
            try:
                self._request_delay()
                res = requests.get(
                    ESEARCH_URL,
                    params=self._eutils_params(
                        {
                            "db": "sra",
                            "term": term,
                            "retmax": 0,
                            "usehistory": "y",
                        }
                    ),
                    timeout=60,
                )
                res.raise_for_status()

                search_root = ET.fromstring(res.text)
                if search_root.find("ERROR") is not None:
                    logger.warning(
                        f"esearch error for {study_acc} ({term}): "
                        f"{search_root.find('ERROR').text}"
                    )
                    continue

                count_el = search_root.find("Count")
                if count_el is None or int(count_el.text) == 0:
                    logger.debug(f"No results for {study_acc} with term: {term}")
                    continue

                count = int(count_el.text)
                webenv = search_root.find("WebEnv").text
                query_key = search_root.find("QueryKey").text

                all_data = []
                retstart = 0
                retmax = 500

                while retstart < count:
                    self._request_delay()
                    res = requests.get(
                        EFETCH_URL,
                        params=self._eutils_params(
                            {
                                "db": "sra",
                                "WebEnv": webenv,
                                "query_key": query_key,
                                "retstart": retstart,
                                "retmax": retmax,
                            }
                        ),
                        timeout=120,
                    )
                    res.raise_for_status()
                    batch = self._parse_efetch_xml(res.text)
                    all_data.extend(batch)
                    retstart += retmax

                if all_data:
                    return all_data

            except requests.exceptions.RequestException as e:
                logger.warning(f"Request error for {study_acc} ({term}): {e}")
            except ET.ParseError as e:
                logger.warning(f"XML parse error for {study_acc}: {e}")

        return None

    # API
    def query_sra(self, study_info):
        study_acc = study_info[0]
        logger.debug(f"Current Study: {study_acc}")
        if study_acc == "-":
            return (study_acc, json.dumps(None))

        bio_project = study_info[7] if len(study_info) > 7 else None

        try:
            data = self._fetch_study_metadata(study_acc, bio_project)

            if not data:
                return (study_acc, json.dumps(None))

            ftp_info = {}
            if study_info[0] != "-":
                ftp_info["ftp_acc"] = study_info[0]
            if study_info[1] != "-":
                ftp_info["ftp_type"] = study_info[1]
            if study_info[3] != "-":
                ftp_info["ftp_updated"] = study_info[3]
            if study_info[4] != "-":
                ftp_info["ftp_published"] = study_info[4]
            if study_info[5] != "-":
                ftp_info["ftp_experiment"] = study_info[5]
            if study_info[6] != "-":
                ftp_info["ftp_sample"] = study_info[6]
            if study_info[7] != "-":
                ftp_info["ftp_bioProject"] = study_info[7]
            if study_info[8] != "-":
                ftp_info["ftp_replacedBy"] = study_info[8]
            data.append(ftp_info)

            return (study_acc, json.dumps(data))
        except Exception as e:
            logger.error(f"Error for {study_acc}: {e}")
            return (study_acc, json.dumps(None))

    def load_cache(self):
        fileloc = "https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab"
        self._download_file(fileloc, "SRA_Accessions.tab")

        logger.info("Retrieving Studies from SRA_Accessions.tab")

        df = pd.read_csv(
            r"SRA_Accessions.tab",
            sep="\t",
            usecols=[
                "Accession",
                "Type",
                "Status",
                "Updated",
                "Published",
                "Experiment",
                "Sample",
                "BioProject",
                "ReplacedBy",
            ],
        )

        only_live = df[df["Status"] == "live"]
        filtered = only_live[only_live["Type"] == "STUDY"]
        accession_list = filtered[
            ["Accession", "Type", "Status", "Updated", "Published", "Experiment", "Sample", "BioProject", "ReplacedBy"]
        ].values.tolist()

        # Used for testing small chunks of data
        if self.DATA_LIMIT:
            accession_list = accession_list[: self.DATA_LIMIT]

        logger.info("Total Studies Found: {}".format(len(accession_list)))
        count = 0

        logger.info("Retrieving Individual Study Metadata from API")

        start = time.time()
        with ThreadPoolExecutor(max_workers=3) as pool:
            data = pool.map(self.query_sra, accession_list)
            for item in data:
                yield item
                count += 1
                if count % 1000 == 0:
                    end = time.time()
                    logger.info("Retrieved 1000 Studies in {} seconds".format(end - start))
                    start = time.time()

        logger.info("Removing SRA_Accessions.tab")
        # os.remove("SRA_Accessions.tab")
        logger.info("Removed SRA_Accessions.tab")

    def parse(self, studies):
        logger.info("Parsing Individual Study Metadata")

        count = 0
        too_big = 0

        for study in studies:
            count += 1
            if count % 1000 == 0:
                logger.info(f"Studies Parsed: {count}")

            study_metadata = json.loads(study[1])

            if study_metadata is None:
                continue
            if len(study_metadata) == 0:
                continue

            logger.debug(f"Study: {study[0]} - {len(study_metadata) - 1} runs")

            if len(study_metadata) > 18000:
                logger.info(f"{study[0]} has more than 18000 records: {len(study_metadata)}")
                study_metadata = study_metadata[-18000:]
                too_big += 1

            ftp_info = study_metadata.pop()
            output = {
                "@context": "https://schema.org/",
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "NCBI SRA",
                    "url": "https://www.ncbi.nlm.nih.gov/sra/",
                    "versionDate": datetime.date.today().strftime("%Y-%m-%d"),
                },
                "@type": "Dataset",
            }
            # Top Level Study Data
            if accession := ftp_info.get("ftp_acc"):
                output["_id"] = "NCBI_SRA_" + accession
                output["url"] = "https://www.ncbi.nlm.nih.gov/sra/" + accession
                output["includedInDataCatalog"]["archivedAt"] = output["url"]
            if updated := ftp_info.get("ftp_updated"):
                output["dateModified"] = datetime.datetime.strptime(updated, "%Y-%m-%dT%H:%M:%SZ").strftime("%Y-%m-%d")
            if published := ftp_info.get("ftp_published"):
                output["datePublished"] = datetime.datetime.strptime(published, "%Y-%m-%dT%H:%M:%SZ").strftime(
                    "%Y-%m-%d"
                )
            if replaced_by := ftp_info.get("ftp_replacedBy"):
                output["sameAs"] = replaced_by

            is_based_on = []
            distribution_list = []
            species_list = []
            seen_species = set()
            measurement_technique_list = []
            seen_measurements_techniques = set()
            author_list = []
            sample_collection_items = []
            seen_samples = set()
            aggregate_data = {
                "sex": set(),
                "cellType": set(),
                "developmentalStage": set(),
                "associatedGenotype": set(),
                "associatedPhenotype": set(),
                "sampleType": set(),
            }
            health_condition_set = set()

            if bio_project := ftp_info.get("ftp_bioProject"):
                bio_project_dict = {}
                bio_project_dict["identifier"] = bio_project
                bio_project_dict["additionalType"] = {
                    "name": "BioProject",
                    "url": "http://purl.obolibrary.org/obo/NCIT_C45293",
                }
                is_based_on.append(bio_project_dict)

            for run_metadata in study_metadata:
                if study_title := run_metadata.get("study_title"):
                    output["name"] = study_title
                if study_abstract := run_metadata.get("study_abstract"):
                    output["description"] = study_abstract

                author_dict = {}
                if organisation := run_metadata.get("organisation"):
                    author_dict["name"] = organisation
                if organisation_contact_email := run_metadata.get("organisation_contact_email"):
                    author_dict["email"] = organisation_contact_email
                if bool(author_dict) and author_dict not in author_list:
                    author_list.append(author_dict)

                # distribution
                distribution_dict = {}
                if gcp_url := run_metadata.get("GCP_url"):
                    distribution_dict["contentUrl"] = gcp_url
                if run_metadata.get("GCP_free_egress"):
                    output["isAccessibleForFree"] = True
                if bool(distribution_dict) and distribution_dict not in distribution_list:
                    distribution_list.append(distribution_dict)

                distribution_dict = {}
                if aws_url := run_metadata.get("AWS_url"):
                    distribution_dict["contentUrl"] = aws_url
                if run_metadata.get("AWS_free_egress"):
                    output["isAccessibleForFree"] = True
                if bool(distribution_dict) and distribution_dict not in distribution_list:
                    distribution_list.append(distribution_dict)

                # species
                if taxon_name := run_metadata.get("taxon_name"):
                    if taxon_name not in seen_species:
                        species_dict = {"name": taxon_name}
                        seen_species.add(taxon_name)
                        species_list.append(species_dict)

                # measurement techniques
                if measurement_technique := run_metadata.get("library_strategy"):
                    if measurement_technique not in seen_measurements_techniques:
                        measurement_technique_dict = {"name": measurement_technique}
                        seen_measurements_techniques.add(measurement_technique)
                        measurement_technique_list.append(measurement_technique_dict)

                # runs
                run_dict = {}
                if run_accession := run_metadata.get("run"):
                    run_dict["identifier"] = run_accession
                    run_dict["additionalType"] = {"name": "Run", "url": "http://purl.obolibrary.org/obo/NCIT_C47911"}
                    run_dict["url"] = "https://www.ncbi.nlm.nih.gov/sra/" + run_accession
                if published := run_metadata.get("published"):
                    try:
                        run_dict["datePublished"] = datetime.datetime.strptime(
                            published, "%Y-%m-%d %H:%M:%S"
                        ).strftime("%Y-%m-%d")
                    except ValueError:
                        pass
                if bool(run_dict) and run_dict not in is_based_on:
                    is_based_on.append(run_dict)

                # experiments
                experiment_dict = {}
                if experiment_accession := run_metadata.get("experiment_accession"):
                    experiment_dict["identifier"] = experiment_accession
                    experiment_dict["url"] = "https://www.ncbi.nlm.nih.gov/sra/" + experiment_accession
                if experiment_title := run_metadata.get("experiment_title"):
                    experiment_dict["name"] = experiment_title
                if experiment_desc := run_metadata.get("experiment_desc"):
                    experiment_dict["description"] = experiment_desc
                if bool(experiment_dict) and experiment_dict not in is_based_on:
                    experiment_dict["additionalType"] = {
                        "name": "Experiment",
                        "url": "http://purl.obolibrary.org/obo/NCIT_C42790",
                    }
                    is_based_on.append(experiment_dict)

                # Collect sample info for SampleCollection
                biosample_id = run_metadata.get("biosample_accession")
                sample_acc = run_metadata.get("sample_accession")
                sample_key = biosample_id or sample_acc
                if sample_key and sample_key not in seen_samples:
                    seen_samples.add(sample_key)
                    if biosample_id:
                        sample_item = {
                            "@type": "Sample",
                            "identifier": biosample_id,
                            "url": f"https://www.ncbi.nlm.nih.gov/biosample/{biosample_id}",
                            "_id": biosample_id.casefold(),
                        }
                        if sample_acc:
                            sample_item["additionalIdentifier"] = sample_acc
                    else:
                        sample_item = {
                            "@type": "Sample",
                            "identifier": sample_acc,
                            "url": f"https://www.ncbi.nlm.nih.gov/sra/{sample_acc}",
                            "_id": sample_acc.casefold(),
                        }
                    if sample_title := run_metadata.get("sample_title"):
                        sample_item["name"] = sample_title
                    sample_collection_items.append(sample_item)

                # instruments
                instrument_dict = {}
                if model := run_metadata.get("model"):
                    instrument_dict["name"] = model
                if bool(instrument_dict) and instrument_dict not in is_based_on:
                    instrument_dict["additionalType"] = {
                        "name": "Instrument",
                        "url": "http://purl.obolibrary.org/obo/NCIT_C16742",
                    }

                # cells
                cell_dict = {}
                if cell_line := run_metadata.get("cell line"):
                    cell_dict["name"] = cell_line
                if cell_line_name := run_metadata.get("cell line name"):
                    cell_dict["name"] = cell_line_name
                if cell_strain := run_metadata.get("strain"):
                    cell_dict["identifier"] = cell_strain
                if bool(cell_dict) and cell_dict not in is_based_on:
                    cell_dict["additionalType"] = {"name": "Cell", "url": "http://purl.obolibrary.org/obo/NCIT_C12508"}
                    is_based_on.append(cell_dict)

                # hapmap
                hapmap_dict = {}
                if hapmap_id := run_metadata.get("HapMap sample ID"):
                    hapmap_dict["identifier"] = hapmap_id
                if cell_line := run_metadata.get("Cell line"):
                    hapmap_dict["name"] = cell_line
                if sex := run_metadata.get("sex"):
                    hapmap_dict["gender"] = sex
                if bool(hapmap_dict) and hapmap_dict not in is_based_on:
                    hapmap_dict["additionalType"] = {
                        "name": "HapMap",
                        "url": "http://purl.obolibrary.org/obo/NCIT_C70979",
                    }
                    is_based_on.append(hapmap_dict)

                # Collect aggregate values from sample attributes
                for attr_tag, agg_field in AGGREGATE_ATTR_MAP.items():
                    if val := run_metadata.get(attr_tag):
                        val_clean = val.strip()
                        if val_clean.casefold() not in INVALID_SAMPLE_VALUES and not val_clean.casefold().startswith("missing:"):
                            aggregate_data[agg_field].add(val_clean)

                # Collect dataset-level fields from sample attributes
                for attr_tag, ds_field in DATASET_ATTR_MAP.items():
                    if val := run_metadata.get(attr_tag):
                        val_clean = val.strip()
                        if val_clean.casefold() not in INVALID_SAMPLE_VALUES and not val_clean.casefold().startswith("missing:"):
                            health_condition_set.add(val_clean)

            def sanitize_json_string(json_string):
                # Replace single quotes with double quotes
                json_string = json_string.replace("'", '"')

                # Remove any trailing commas (common mistake in JSON)
                json_string = re.sub(r",\s*}", "}", json_string)
                json_string = re.sub(r",\s*]", "]", json_string)

                return json_string

            def remove_duplicates_and_limit(items, limit=100):
                seen = set()
                unique_items = []
                for item in items:
                    sanitized_item = sanitize_json_string(json.dumps(item))
                    if sanitized_item not in seen:
                        seen.add(sanitized_item)
                        unique_items.append(item)
                    if len(unique_items) == limit:
                        break
                return unique_items

            if len(is_based_on):
                try:
                    is_based_on = remove_duplicates_and_limit(is_based_on, limit=100)
                except Exception as e:
                    logger.error(f"Error while processing isBasedOn: {e}")
                    is_based_on = []

                if len(is_based_on) > 100:
                    is_based_on = is_based_on[:100]

                output["isBasedOn"] = is_based_on
            if len(author_list):
                output["author"] = author_list
            if len(species_list):
                output["species"] = species_list
            if len(distribution_list):
                output["distribution"] = distribution_list
            if len(measurement_technique_list):
                output["measurementTechnique"] = measurement_technique_list
            if health_condition_set:
                output["healthCondition"] = [{"name": v} for v in sorted(health_condition_set)]

            if sample_collection_items:
                sample_collection = {
                    "@type": "SampleCollection",
                    "itemListElement": sample_collection_items,
                    "numberOfItems": {
                        "value": len(sample_collection_items),
                        "unitText": "sample",
                    },
                }

                # Build aggregateElement from collected sample attributes
                aggregate_element = {}
                if aggregate_data["sex"]:
                    aggregate_element["sex"] = sorted(aggregate_data["sex"])
                if aggregate_data["cellType"]:
                    aggregate_element["cellType"] = [
                        {"name": v} for v in sorted(aggregate_data["cellType"])
                    ]
                if aggregate_data["developmentalStage"]:
                    aggregate_element["developmentalStage"] = [
                        {"name": v} for v in sorted(aggregate_data["developmentalStage"])
                    ]
                if aggregate_data["associatedGenotype"]:
                    aggregate_element["associatedGenotype"] = sorted(aggregate_data["associatedGenotype"])
                if aggregate_data["associatedPhenotype"]:
                    aggregate_element["associatedPhenotype"] = [
                        {"name": v} for v in sorted(aggregate_data["associatedPhenotype"])
                    ]
                if aggregate_data["sampleType"]:
                    aggregate_element["sampleType"] = [
                        {"name": v} for v in sorted(aggregate_data["sampleType"])
                    ]

                if aggregate_element:
                    sample_collection["aggregateElement"] = aggregate_element

                output["sample"] = sample_collection

            yield output

        logger.info(f"Finished Parsing {count} Studies")
        logger.info(f"Total large documents: {too_big}")

    def update_cache(self):
        RETRIEVE_METHOD = False
        start = time.time()

        last_updated = self.retreive_last_updated()

        fileloc = "https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab"
        self._download_file(fileloc, "SRA_Accessions.tab")

        logger.info("Retrieving Studies from SRA_Accessions.tab")

        df = pd.read_csv(
            r"SRA_Accessions.tab",
            sep="\t",
            usecols=[
                "Accession",
                "Type",
                "Status",
                "Updated",
                "Published",
                "Experiment",
                "Sample",
                "BioProject",
                "ReplacedBy",
            ],
        )
        only_live = df[df["Status"] == "live"]
        filtered = only_live[only_live["Type"] == "STUDY"]

        if RETRIEVE_METHOD == True:
            last_record = self.retrieve_last_record()
            logger.info(f"Last Record: {last_record}")
            accession_list = filtered[
                [
                    "Accession",
                    "Type",
                    "Status",
                    "Updated",
                    "Published",
                    "Experiment",
                    "Sample",
                    "BioProject",
                    "ReplacedBy",
                ]
            ].values.tolist()
            # ignore studies before last_record
            for i, study in enumerate(accession_list):
                if study[0] == last_record:
                    logger.info(f"Found last record at {i} out of {len(accession_list)}")
                    accession_list = accession_list[i:]
                    break
        else:
            new_studies = filtered[filtered["Updated"] > last_updated]
            accession_list = new_studies[
                [
                    "Accession",
                    "Type",
                    "Status",
                    "Updated",
                    "Published",
                    "Experiment",
                    "Sample",
                    "BioProject",
                    "ReplacedBy",
                ]
            ].values.tolist()
            logger.info("Total Studies Found: {}".format(len(accession_list)))

        if self.DATA_LIMIT:
            accession_list = accession_list[: self.DATA_LIMIT]

        count = 0

        logger.info("Retrieving Individual Study Metadata from API")

        with ThreadPoolExecutor(max_workers=3) as pool:
            data = pool.map(self.query_sra, accession_list)
            for item in data:
                yield item
                count += 1
                if count % 1000 == 0:
                    end = time.time()
                    logger.info("Retrieved {} Studies in {} seconds".format(count, end - start))
                    start = time.time()

        logger.info("Removing SRA_Accessions.tab")
        # os.remove("SRA_Accessions.tab")
        logger.info("Removed SRA_Accessions.tab")
        self.insert_last_updated()
