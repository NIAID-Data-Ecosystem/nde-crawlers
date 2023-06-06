import datetime
import json
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import wget
from kingfisher import annotate
from sql_database import NDEDatabase

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


class NCBI_SRA(NDEDatabase):
    # override variables
    SQL_DB = "ncbi_sra2.db"
    EXPIRE = datetime.timedelta(days=90)

    # Used for testing small chunks of data
    DATA_LIMIT = None

    # API
    def query_sra(self, study_info):
        study_acc = study_info[0]
        logger.info(f"Current Study: {study_acc}")
        if study_acc == "-":
            return (study_acc, json.dumps(None))
        try:
            annotate(
                run_identifiers=None,
                run_identifiers_file=None,
                bioproject_accession=study_acc,
                output_file=f"{study_acc}.json",
                output_format="json",
                all_columns=True,
            )
            with open(f"{study_acc}.json") as f:
                data = json.load(f)
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
                os.remove(f"{study_acc}.json")
                return (study_acc, json.dumps(data))
        except KeyError as e:
            logger.error(f"KeyError for {study_acc}: {e}")
            return (study_acc, json.dumps(None))
        except Exception as e:
            if "HTTP Failure" in str(e):
                logger.error(f"HTTP Failure for {study_acc}: {e}")
                return (study_acc, json.dumps(None))
            else:
                logger.error(f"Unknown Error for {study_acc}: {e}")
                return (study_acc, json.dumps(None))

    def load_cache(self):
        logger.info("Starting FTP Download")

        fileloc = "https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab"
        retry_count = 0
        while True:
            try:
                logger.info("Starting FTP Download")
                wget.download(fileloc, out="SRA_Accessions.tab")
                break
            except ContentTooShortError as e:
                retry_count += 1
                if retry_count > 10:
                    raise e
                else:
                    logger.error("Error downloading file. %s. Retrying...", e)
                    logger.info("Try count: %s of 10", retry_count)
                    sleep(retry_count * 2)

        logger.info("FTP Download Complete")

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
                logger.info(f"Yielding {item[0]}")
                yield item
                count += 1
                if count % 1000 == 0:
                    end = time.time()
                    logger.info("Retrieved 1000 Studies in {} seconds".format(end - start))
                    start = time.time()

        logger.info("Removing SRA_Accessions.tab")
        os.remove("SRA_Accessions.tab")
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

            logger.info(f"Study: {study[0]}")

            if study_metadata is None:
                logger.info(f"No Metadata for {study[0]}")
                continue
            if len(study_metadata) == 0:
                logger.info(f"No Metadata for {study[0]}")
                continue

            logger.info(f"Total runs: {len(study_metadata) - 1}")

            if len(study_metadata) > 18000:
                logger.info(f"{study[0]} has more than 18000 records: {len(study_metadata)}")
                study_metadata = study_metadata[-18000:]
                too_big += 1

            ftp_info = study_metadata.pop()
            output = {
                "@context": "https://schema.org/",
                "includedInDataCatalog": {
                    "@type": "Dataset",
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
            author_list = []

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
                if gcp_free_egress := run_metadata.get("GCP_free_egress"):
                    output["isAccessibleForFree"] = True
                if bool(distribution_dict) and distribution_dict not in distribution_list:
                    distribution_list.append(distribution_dict)

                distribution_dict = {}
                if aws_url := run_metadata.get("AWS_url"):
                    distribution_dict["contentUrl"] = aws_url
                if aws_free_egress := run_metadata.get("AWS_free_egress"):
                    output["isAccessibleForFree"] = True
                if bool(distribution_dict) and distribution_dict not in distribution_list:
                    distribution_list.append(distribution_dict)

                # species
                species_dict = {}
                if organism_taxid := run_metadata.get("organism_taxid"):
                    species_dict["identifier"] = organism_taxid
                if organism_name := run_metadata.get("organism_name"):
                    species_dict["name"] = organism_name
                if bool(species_dict) and species_dict not in species_list:
                    species_dict["additionalType"] = {
                        "name": "Species",
                        "url": "http://purl.obolibrary.org/obo/NCIT_C45293",
                    }
                    species_list.append(species_dict)

                # runs
                run_dict = {}
                if run_accession := run_metadata.get("run"):
                    run_dict["identifier"] = run_accession
                    run_dict["additionalType"] = {"name": "Run", "url": "http://purl.obolibrary.org/obo/NCIT_C47911"}
                    run_dict["url"] = "https://www.ncbi.nlm.nih.gov/sra/" + run_accession
                if published := run_metadata.get("published"):
                    run_dict["datePublished"] = datetime.datetime.strptime(published, "%Y-%m-%d %H:%M:%S").strftime(
                        "%Y-%m-%d"
                    )
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

                sample_dict = {}
                if sample_accession := run_metadata.get("sample_accession"):
                    sample_dict["identifier"] = sample_accession
                    sample_dict["url"] = "https://www.ncbi.nlm.nih.gov/sra/" + sample_accession
                if sample_comment := run_metadata.get("sample comment"):
                    sample_dict["description"] = sample_comment
                if sample_title := run_metadata.get("sample_title"):
                    sample_dict["name"] = sample_title
                if bool(sample_dict) and sample_dict not in is_based_on:
                    sample_dict["additionalType"] = {
                        "name": "Sample",
                        "url": "http://purl.obolibrary.org/obo/NCIT_C70699",
                    }
                    is_based_on.append(sample_dict)

                # instruments
                instrument_dict = {}
                if instrument := run_metadata.get("instrument"):
                    instrument_dict["name"] = instrument
                if instrument_model := run_metadata.get("instrument_model"):
                    instrument_dict["identifier"] = instrument_model
                if instrument_model_description := run_metadata.get("instrument_model_desc"):
                    instrument_dict["description"] = instrument_model_description
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

            if len(is_based_on):
                output["isBasedOn"] = is_based_on
            if len(species_list):
                output["species"] = species_list
            if len(distribution_list):
                output["distribution"] = distribution_list

            yield output
            logger.info(f"Yielded: {study[0]}")

        logger.info(f"Finished Parsing {count} Studies")
        logger.info(f"Total large documents: {too_big}")

    def update_cache(self):
        RETRIEVE_METHOD = False
        start = time.time()

        last_updated = self.retreive_last_updated()

        fileloc = "https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab"
        retry_count = 0
        while True:
            try:
                logger.info("Starting FTP Download")
                wget.download(fileloc, out="SRA_Accessions.tab")
                break
            except ContentTooShortError as e:
                retry_count += 1
                if retry_count > 10:
                    raise e
                else:
                    logger.error("Error downloading file. %s. Retrying...", e)
                    logger.info("Try count: %s of 10", retry_count)
                    sleep(retry_count * 2)

        logger.info("FTP Download Complete")

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

        count = 0

        logger.info("Retrieving Individual Study Metadata from API")

        with ThreadPoolExecutor(max_workers=3) as pool:
            data = pool.map(self.query_sra, accession_list)
            for item in data:
                logger.info(f"Yielding {item[0]}")
                yield item
                count += 1
                if count % 1000 == 0:
                    end = time.time()
                    logger.info("Retrieved {} Studies in {} seconds".format(count, end - start))
                    start = time.time()

        logger.info("Removing SRA_Accessions.tab")
        os.remove("SRA_Accessions.tab")
        logger.info("Removed SRA_Accessions.tab")
        self.insert_last_updated()
