import datetime
import json
import logging
import re
import socket

import dateutil.parser
import xmltodict
from Bio import Entrez
from tenacity import retry, stop_after_attempt, wait_fixed

DEFAULT_TIMEOUT = 30  # seconds
socket.setdefaulttimeout(DEFAULT_TIMEOUT)
GEO_API_KEY = "3048f6bdb7c91cc8ad7af802559ec470e609"
GEO_EMAIL = "cwu@scripps.edu"

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


Entrez.email = GEO_EMAIL
Entrez.api_key = GEO_API_KEY


def insert_value(d, key, value):
    if key in d:
        if isinstance(d[key], list):
            d[key].append(value)
        else:
            d[key] = [d[key], value]
    else:
        d[key] = value


@retry(stop=stop_after_attempt(3), wait=wait_fixed(2), reraise=True)
def query_acc(term, retstart, retmax):
    handle = Entrez.esearch(db="biosample", term=term, usehistory="y", timeout=DEFAULT_TIMEOUT)
    record = Entrez.read(handle)
    handle.close()
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    handle = Entrez.esummary(
        db="biosample",
        query_key=query_key,
        webenv=webenv,
        retmode="json",
        retmax=retmax,
        retstart=retstart,
        timeout=DEFAULT_TIMEOUT,
    )

    records = json.load(handle)
    handle.close()
    return records["result"]


def fetch_all_samples():
    """
    Fetch all Accession numbers from NCBI Biosample using Biopython Entrez ESearch with pagination.

    """

    # First, get total count
    handle = Entrez.esearch(db="biosample", term="all[filter]", usehistory="y", timeout=DEFAULT_TIMEOUT)
    record = Entrez.read(handle)
    total = int(record["Count"])
    logger.info(f"Total records to download: {total}")
    handle.close()

    retmax = 500  # max allowed by NCBI for json output
    for retstart in range(0, total, retmax):
        accs = query_acc("all[filter]", retstart, retmax)
        yield accs
        if (retstart + retmax) % 10000 == 0:
            logger.info(f"Fetched {retstart + retmax} of {total} records")


def parse_xml(sample_dict, output):
    if cod := sample_dict["BioSample"].get("@access"):
        if cod == "public":
            output["conditionsOfAccess"] = "Open"
        else:
            output["conditionsOfAccess"] = "Restricted"

    attributes = {}
    if a := sample_dict["BioSample"].get("Attributes"):
        if attr_list := a.get("Attribute"):
            if not isinstance(attr_list, list):
                attr_list = [attr_list]
            for attr in attr_list:
                if attr.get("@harmonized_name"):
                    key = attr["@harmonized_name"]
                elif attr.get("@attribute_name"):
                    key = attr["@attribute_name"]

                if value := attr.get("#text"):
                    invalid_values = [
                        "na",
                        "not applicable",
                        "not collected",
                        "not provided",
                        "not available",
                        "-",
                        "none",
                        "n/a",
                        "missing",
                        "unknown",
                        "restricted access",
                        "unspecified",
                        "not determined",
                        "not recorded",
                        "unk",
                        "blank",
                        "null"
                    ]
                    if (
                        not attributes.get(key)
                        and value.casefold() not in invalid_values
                        and not value.casefold().startswith("missing:")
                    ):
                        attributes[key] = value

    if location_of_origin := attributes.get("add_host_info"):
        l = {"name": location_of_origin}
        insert_value(output, "locationOfOrigin", l)

    if collection_method := attributes.get("add_recov_method"):
        insert_value(output, "collectionMethod", collection_method)

    if env_system := attributes.get("adjacent_environment"):
        insert_value(output, "environmentalSystem", env_system)

    if dev_stage := attributes.get("age"):
        d = {"name": dev_stage}
        insert_value(output, "developmentalStage", d)

    if dev_stage := attributes.get("Age (years)"):
        d = {"value": dev_stage}
        insert_value(output, "developmentalStage", d)

    if alt_name := attributes.get("alias"):
        insert_value(output, "alternateName", alt_name)

    if alt_id := attributes.get("alternate_ID"):
        insert_value(output, "alternateIdentifier", alt_id)

    if location_of_origin := attributes.get("altitude"):
        l = {"name": location_of_origin}
        insert_value(output, "locationOfOrigin", l)

    if sample_type := attributes.get("analyte_type"):
        insert_value(output, "sampleType", {"name": sample_type})

    if env_system := attributes.get("animal_env"):
        insert_value(output, "environmentalSystem", env_system)

    if species := attributes.get("ArrayExpress-SPECIES"):
        s = {"name": species}
        insert_value(output, "species", s)

    if name := attributes.get("author"):
        author = {"name": name}
        insert_value(output, "author", author)

    if attributes.get("bacteria_available_from"):
        output["sampleAvailability"] = True

    if sample_state := attributes.get("biol_stat"):
        insert_value(output, "sampleState", sample_state)

    if location_of_origin := attributes.get("biological material altitude"):
        l = {"name": location_of_origin}
        insert_value(output, "locationOfOrigin", l)

    if location_of_origin := attributes.get("biological material coordinates uncertainty"):
        l = {"name": location_of_origin}
        insert_value(output, "locationOfOrigin", l)

    if location_of_origin := attributes.get("biological material geographic location"):
        l = {"name": location_of_origin}
        insert_value(output, "locationOfOrigin", l)

    l = {}
    if latitude := attributes.get("biological material latitude"):
        try:
            if "," in latitude:
                latitude = latitude.replace(",", ".")
            l["geo"] = {"latitude": float(latitude)}
        except Exception as e:
            logger.warning(f"Failed to parse latitude: {latitude} Error: {e}")
    if longitude := attributes.get("biological material longitude"):
        try:
            if "," in longitude:
                longitude = longitude.replace(",", ".")
            l["geo"]["longitude"] = float(longitude)
        except Exception as e:
            logger.warning(f"Failed to parse longitude: {longitude} Error: {e}")

    if l:
        insert_value(output, "locationOfOrigin", l)

    if sample_process := attributes.get("biological material preprocessing"):
        insert_value(output, "sampleProcess", sample_process)

    if contributor := attributes.get("biomaterial_provider"):
        contributor = {"name": contributor}
        insert_value(output, "contributor", contributor)

    if item_location := attributes.get("biospecimen_repository"):
        item_location = {"name": item_location}
        insert_value(output, "itemLocation", item_location)
        output["sampleAvailability"] = True

    if alt_id := attributes.get("biospecimen_repository_sample_id"):
        insert_value(output, "alternateIdentifier", alt_id)

    if alt_id := attributes.get("Catalog_number"):
        insert_value(output, "alternateIdentifier", alt_id)

    if alt_id := attributes.get("CBP Individual ID"):
        insert_value(output, "alternateIdentifier", alt_id)

    if alt_id := attributes.get("CBP Sample ID"):
        insert_value(output, "alternateIdentifier", alt_id)

    if associated_genotype := attributes.get("cell subpopulation"):
        insert_value(output, "associatedGenotype", associated_genotype)

    if cell_type := attributes.get("cell_line"):
        insert_value(output, "cellType", {"name": cell_type})

    if cell_type := attributes.get("cell_subtype"):
        insert_value(output, "cellType", {"name": cell_type})

    if cell_type := attributes.get("cell_type"):
        insert_value(output, "cellType", {"name": cell_type})

    if collector := attributes.get("collected_by"):
        collector = {"name": collector}
        insert_value(output, "collector", collector)

    if collector := attributes.get("collecting institution"):
        collector = {"name": collector}
        insert_value(output, "collector", collector)

    if date_collected := attributes.get("collection_date"):
        try:
            output["dateCollected"] = dateutil.parser.parse(date_collected, ignoretz=True).date().isoformat()
        except Exception as e:
            if "/" in date_collected:
                try:
                    date_collected = date_collected.split("/")[0]
                    output["dateCollected"] = dateutil.parser.parse(date_collected, ignoretz=True).date().isoformat()
                except Exception as e:
                    logger.warning(f"Failed to parse collection_date: {date_collected} Error: {e}")
            else:
                logger.warning(f"Failed to parse collection_date: {date_collected} Error: {e}")

    if location_of_origin := attributes.get("collection_depth"):
        l = {"name": location_of_origin}
        insert_value(output, "locationOfOrigin", l)

    if collection_method := attributes.get("collection_method"):
        insert_value(output, "collectionMethod", collection_method)

    if sample_quantity := attributes.get("collection_volume"):
        sample_quantity = {"unitText": sample_quantity}
        insert_value(output, "sampleQuantity", sample_quantity)

    if collector_name := attributes.get("collector name"):
        collector_name = collector_name.split(", ")
        for name in collector_name:
            collector = {"name": name}
            insert_value(output, "collector", collector)

    if species := attributes.get("common name"):
        s = {"name": species}
        insert_value(output, "species", s)

    if species := attributes.get("Common name"):
        s = {"name": species}
        insert_value(output, "species", s)

    if species := attributes.get("common_name"):
        s = {"name": species}
        insert_value(output, "species", s)

    if sample_process := attributes.get("concentration_method"):
        insert_value(output, "sampleProcess", sample_process)

    if date_processed := attributes.get("cult_isol_date"):
        try:
            output["dateProcessed"] = dateutil.parser.parse(date_processed, ignoretz=True).date().isoformat()
        except Exception as e:
            logger.warning(f"Failed to parse cult_isol_date: {date_processed} Error: {e}")

    if is_based_on := attributes.get("derived from assembly"):
        is_based_on = {"identifier": is_based_on}
        insert_value(output, "isBasedOn", is_based_on)

    if is_based_on := attributes.get("derived_from"):
        is_based_on = is_based_on.split(", ")
        for id in is_based_on:
            id = {"identifier": id}
            insert_value(output, "isBasedOn", id)

    if description := attributes.get("description"):
        if not output.get("description"):
            output["description"] = description
        else:
            output["description"] += description

    if description := attributes.get("description2"):
        if not output.get("description"):
            output["description"] = description
        else:
            output["description"] += description

    if description := attributes.get("design_description"):
        if not output.get("description"):
            output["description"] = description
        else:
            output["description"] += description

    if hc := attributes.get("disease"):
        hc = {"name": hc}
        insert_value(output, "healthCondition", hc)

    if attributes.get("disease_stage"):
        insert_value(output, "variableMeasured", {"name": "disease stage"})

    if contributor := attributes.get("Dissected by"):
        contributor = contributor.split(", ")
        for c in contributor:
            c = {"name": c}
            insert_value(output, "contributor", c)

    if sample_process := attributes.get("DNA_extraction_method"):
        insert_value(output, "sampleProcess", sample_process)

    if env_system := attributes.get("env_broad_scale"):
        insert_value(output, "environmentalSystem", env_system)

    if env_system := attributes.get("env_local_scale"):
        insert_value(output, "environmentalSystem", env_system)

    if env_system := attributes.get("env_medium"):
        insert_value(output, "environmentalSystem", env_system)

    if attributes.get("environmental_sample", "").casefold() == "yes" or attributes.get("environmental-sample", 0) == 1:
        insert_value(output, "sampleType", {"name": "environmental"})

    if attributes.get("ethnicity"):
        insert_value(output, "variableMeasured", {"name": "ethnicity"})

    if identifier := attributes.get("External Id"):
        insert_value(output, "identifier", identifier)

    if alternate_id := attributes.get("External Sample ID"):
        insert_value(output, "alternateIdentifier", alternate_id)

    if contributor := attributes.get("Extracted by"):
        contributor = contributor.split(", ")
        for c in contributor:
            c = {"name": c}
            insert_value(output, "contributor", c)

    if contributor := attributes.get("Extraction node"):
        contributor = contributor.split(", ")
        for c in contributor:
            c = {"name": c}
            insert_value(output, "contributor", c)

    if sample_process := attributes.get("extraction_method"):
        insert_value(output, "sampleProcess", sample_process)

    if location_of_origin := attributes.get("food_origin"):
        l = {"name": location_of_origin}
        insert_value(output, "locationOfOrigin", l)

    if sample_process := attributes.get("food_processing_method"):
        insert_value(output, "sampleProcess", sample_process)

    if associated_genotype := attributes.get("genotype"):
        insert_value(output, "associatedGenotype", associated_genotype)

    if location_of_origin := attributes.get("geo_loc_name"):
        l = {"name": location_of_origin}
        insert_value(output, "locationOfOrigin", l)

    if alternate_id := attributes.get("gisaid_accession"):
        insert_value(output, "alternateIdentifier", alternate_id)

    if env_system := attributes.get("GOLD Ecosystem Classification"):
        insert_value(output, "environmentalSystem", env_system)

    if collector := attributes.get("Greenhouse collection"):
        collector = {"name": collector}
        insert_value(output, "collector", collector)

    if alternate_id := attributes.get("GSC Sample ID"):
        insert_value(output, "alternateIdentifier", alternate_id)

    if env_system := attributes.get("habitat"):
        insert_value(output, "environmentalSystem", env_system)

    if hc := attributes.get("health_state"):
        hc = {"name": hc}
        insert_value(output, "healthCondition", hc)

    if sample_process := attributes.get("histology"):
        insert_value(output, "sampleProcess", sample_process)
        insert_value(output, "sampleType", {"name": "histology"})

    if author := attributes.get("Holding_institution"):
        author = {"name": author}
        insert_value(output, "author", author)

    if species := attributes.get("host"):
        s = {"name": species}
        insert_value(output, "species", s)

    if env_system := attributes.get("host habitat"):
        insert_value(output, "environmentalSystem", env_system)

    if dev_stage := attributes.get("host_age"):
        d = {"name": dev_stage}
        insert_value(output, "developmentalStage", d)

    if anatomical_system := attributes.get("host_body_habitat"):
        insert_value(output, "anatomicalSystem", {"name": anatomical_system})

    if species := attributes.get("host_common_name"):
        s = {"name": species}
        insert_value(output, "species", s)

    if species := attributes.get("host_description"):
        s = {"name": species}
        insert_value(output, "species", s)

    if hc := attributes.get("host_disease"):
        hc = {"name": hc}
        insert_value(output, "healthCondition", hc)

    if attributes.get("host_disease_outcome"):
        insert_value(output, "variableMeasured", {"name": "disease outcome"})

    if attributes.get("host_disease_stage"):
        insert_value(output, "variableMeasured", {"name": "disease stage"})

    if associated_genotype := attributes.get("host_genotype"):
        insert_value(output, "associatedGenotype", associated_genotype)

    if hc := attributes.get("host_health_state"):
        hc = {"name": hc}
        insert_value(output, "healthCondition", hc)

    if dev_stage := attributes.get("host_life_stage"):
        d = {"name": dev_stage}
        insert_value(output, "developmentalStage", d)

    if associated_phenotype := attributes.get("host_phenotype"):
        insert_value(output, "associatedPhenotype", {"name": associated_phenotype})

    if is_part_of := attributes.get("host_subject_id"):
        is_part_of = {"identifier": is_part_of}
        insert_value(output, "isPartOf", is_part_of)

    if species := attributes.get("host_taxid"):
        s = {"identifier": species}
        insert_value(output, "species", s)

    if sex := attributes.get("host_sex"):
        output["sex"] = sex

    if species := attributes.get("host_taxid"):
        s = {"identifier": species}
        insert_value(output, "species", s)

    if anatomical_system := attributes.get("host_tissue_sampled"):
        insert_value(output, "anatomicalSystem", {"name": anatomical_system})

    if identifier := attributes.get("identifer"):
        insert_value(output, "identifier", identifier)

    if contributor := attributes.get("identified_by"):
        contributor = contributor.split(", ")
        aff = attributes.get("identifier_affiliation")
        for c in contributor:
            c = {"name": c}
            if aff:
                c["affiliation"] = aff
            insert_value(output, "contributor", c)

    if instrument := attributes.get("instrument"):
        instrument = {"name": instrument}
        insert_value(output, "instrument", instrument)

    if instrument := attributes.get("instrument_model"):
        instrument = {"name": instrument}
        insert_value(output, "instrument", instrument)

    if experimental_purpose := attributes.get("investigation_type"):
        insert_value(output, "experimentalPurpose", experimental_purpose)

    if species := attributes.get("lab_host"):
        s = {"name": species}
        insert_value(output, "species", s)

    if latitude := attributes.get("latitude"):
        try:
            if "," in latitude:
                latitude = latitude.replace(",", ".")
            l = {"geo": {"latitude": float(latitude)}}
            if longitude := attributes.get("longitude"):
                if "," in longitude:
                    longitude = longitude.replace(",", ".")
                l["geo"]["longitude"] = float(longitude)
            insert_value(output, "locationOfOrigin", l)
        except Exception as e:
            logger.warning(f"Failed to parse latitude/longitude: {location_of_origin}, {longitude} Error: {e}")

    if experimental_purpose := attributes.get("library type"):
        insert_value(output, "experimentalPurpose", experimental_purpose)

    if mt := attributes.get("library_protocol"):
        mt = {"name": mt}
        insert_value(output, "measurementTechnique", mt)

    if mt := attributes.get("library_selection"):
        mt = {"name": mt}
        insert_value(output, "measurementTechnique", mt)

    if mt := attributes.get("library_source"):
        mt = {"name": mt}
        insert_value(output, "measurementTechnique", mt)

    if mt := attributes.get("library_strategy"):
        mt = {"name": mt}
        insert_value(output, "measurementTechnique", mt)

    if dev_stage := attributes.get("life_stage"):
        d = {"name": dev_stage}
        insert_value(output, "developmentalStage", d)

    if env_system := attributes.get("local environmental context"):
        insert_value(output, "environmentalSystem", env_system)

    if env_system := attributes.get("metagenome_source"):
        insert_value(output, "environmentalSystem", env_system)
        insert_value(output, "sampleType", {"name": "metagenome"})

    if env_system := attributes.get("metagenomic source"):
        insert_value(output, "environmentalSystem", env_system)
        insert_value(output, "sampleType", {"name": "metagenome"})

    if attributes.get("mouse") or attributes.get("mouse_strain"):
        insert_value(output, "species", {"name": "Mus musculus"})

    if species := attributes.get("organism"):
        insert_value(output, "species", {"name": species})

    if species := attributes.get("organism common name"):
        insert_value(output, "species", {"name": species})

    if location_of_origin := attributes.get("Original geographic location"):
        for loc in location_of_origin.split(", "):
            l = {"name": loc}
            insert_value(output, "locationOfOrigin", l)

    if hc := attributes.get("parasites"):
        hc = {"name": hc}
        insert_value(output, "healthCondition", hc)

    if anatomical_structure := attributes.get("plant anatomical entity"):
        insert_value(output, "anatomicalStructure", {"name": anatomical_structure})

    if dev_stage := attributes.get("plant developmental stage"):
        insert_value(output, "developmentalStage", {"name": dev_stage})

    if anatomical_structure := attributes.get("plant_struc"):
        insert_value(output, "anatomicalStructure", {"name": anatomical_structure})

    if instrument := attributes.get("platform"):
        insert_value(output, "instrument", {"name": instrument})

    if experimental_purpose := attributes.get("purpose_of_sampling"):
        insert_value(output, "experimentalPurpose", experimental_purpose)

    if experimental_purpose := attributes.get("purpose_of_sequencing"):
        insert_value(output, "experimentalPurpose", experimental_purpose)

    if experimental_purpose := attributes.get("purpose_of_ww_sampling"):
        insert_value(output, "experimentalPurpose", experimental_purpose)

    if experimental_purpose := attributes.get("purpose_of_ww_sequencing"):
        insert_value(output, "experimentalPurpose", experimental_purpose)

    if experimental_purpose := attributes.get("samp_capt_status"):
        insert_value(output, "experimentalPurpose", experimental_purpose)

    if collection_method := attributes.get("samp_collect_device"):
        insert_value(output, "collectionMethod", collection_method)

    if collection_method := attributes.get("samp_collect_method"):
        insert_value(output, "collectionMethod", collection_method)

    if env_system := attributes.get("samp_collect_point"):
        insert_value(output, "environmentalSystem", env_system)

    if sample_type := attributes.get("samp_mat_type"):
        insert_value(output, "sampleType", {"name": sample_type})

    if sample_process := attributes.get("samp_pooling"):
        insert_value(output, "sampleProcess", sample_process)

    if sample_quantity := attributes.get("samp_size"):
        sample_quantity = {"unitText": sample_quantity}
        insert_value(output, "sampleQuantity", sample_quantity)

    if item_location := attributes.get("samp_store_loc"):
        item_location = {"name": item_location}
        insert_value(output, "itemLocation", item_location)

    if sample_process := attributes.get("soil_type_meth"):
        insert_value(output, "sampleProcess", sample_process)

    if species := attributes.get("species"):
        s = {"name": species}
        insert_value(output, "species", s)

    if mt := attributes.get("study_design"):
        mt = {"name": mt}
        insert_value(output, "measurementTechnique", mt)

    if is_basis_for := attributes.get("study_name"):
        is_basis_for = {"name": is_basis_for}
        insert_value(output, "isBasisFor", is_basis_for)

    if species := attributes.get("sub_species"):
        s = {"name": species}
        insert_value(output, "species", s)

    if identifier := attributes.get("submitted_sample_id"):
        insert_value(output, "identifier", identifier)

    if is_part_of := attributes.get("submitted_subject_id"):
        is_part_of = {"identifier": is_part_of}
        insert_value(output, "isPartOf", is_part_of)

    if contributor := attributes.get("submitter_handle"):
        contributor = {"name": contributor}
        insert_value(output, "contributor", contributor)

    if experimental_purpose := attributes.get("target analysis type"):
        insert_value(output, "experimentalPurpose", experimental_purpose)

    if mt := attributes.get("target gene"):
        mt = {"name": mt}
        insert_value(output, "measurementTechnique", mt)

    if mt := attributes.get("target_gene"):
        mt = {"name": mt}
        insert_value(output, "measurementTechnique", mt)

    if species := attributes.get("Tax ID"):
        s = {"identifier": species}
        insert_value(output, "species", s)

    if species := attributes.get("taxonomic classification"):
        s = {"name": species}
        insert_value(output, "species", s)

    if temporal_coverage := attributes.get("time"):
        temporal_coverage = {"name": temporal_coverage}
        insert_value(output, "temporalCoverage", temporal_coverage)

    if temporal_coverage := attributes.get("Time (days)"):
        temporal_coverage = {"duration": temporal_coverage}
        insert_value(output, "temporalCoverage", temporal_coverage)

    if anatomical_system := attributes.get("tissue"):
        insert_value(output, "anatomicalSystem", {"name": anatomical_system})

    if sample_state := attributes.get("tissue_health"):
        insert_value(output, "sampleState", sample_state)

    if cell_type := attributes.get("transferred cells"):
        insert_value(output, "cellType", {"name": cell_type})

    if temporal_coverage := attributes.get("treatment endpoint"):
        temporal_coverage = {"name": temporal_coverage}
        insert_value(output, "temporalCoverage", temporal_coverage)

    if temporal_coverage := attributes.get("treatment starting_point"):
        temporal_coverage = {"name": temporal_coverage}
        insert_value(output, "temporalCoverage", temporal_coverage)

    if sample_process := attributes.get("virus_enrich_appr"):
        insert_value(output, "sampleProcess", sample_process)

    if sample_process := attributes.get("ww_processing_protocol"):
        insert_value(output, "sampleProcess", sample_process)

    if env_system := attributes.get("ww_sample_matrix"):
        insert_value(output, "environmentalSystem", env_system)

    if env_system := attributes.get("ww_sample_site"):
        insert_value(output, "environmentalSystem", env_system)

    geo_loc = {}
    mat_loc = {}
    for key, value in attributes.items():
        if key.casefold().startswith("filename"):
            value = {"name": value}
            if encoding_format := attributes.get("filetype"):
                value["encodingFormat"] = encoding_format
            insert_value(output, "distribution", value)

        if key.casefold().startswith("geographic location"):
            if "(region and locality)" in key:
                geo_loc["name"] = value
            try:
                if "(latitude)" in key:
                    if "," in value:
                        value = value.replace(",", ".")
                    geo_loc["geo"] = {"latitude": float(value)}
            except Exception as e:
                logger.warning(f"Failed to parse latitude: {latitude} Error: {e}")
            try:
                if "(longitude)" in key:
                    if "," in value:
                        value = value.replace(",", ".")
                    if not geo_loc.get("geo"):
                        geo_loc["geo"] = {}
                    geo_loc["geo"]["longitude"] = float(value)
            except Exception as e:
                logger.warning(f"Failed to parse longitude: {longitude} Error: {e}")

        if key.casefold().startswith("material source"):
            try:
                if "latitude" in key:
                    if "," in value:
                        value = value.replace(",", ".")
                    mat_loc["geo"] = {"latitude": float(value)}
            except Exception as e:
                logger.warning(f"Failed to parse latitude: {latitude} Error: {e}")

            try:
                if "longitude" in key:
                    if "," in value:
                        value = value.replace(",", ".")
                    if not mat_loc.get("geo"):
                        mat_loc["geo"] = {}
                    mat_loc["geo"]["longitude"] = float(value)
            except Exception as e:
                logger.warning(f"Failed to parse longitude: {longitude} Error: {e}")

            if "ID" in key:
                mat_loc["identifier"] = value
            elif "geographic location" in key:
                mat_loc["name"] = value

        if key.casefold() == "name":
            insert_value(output, "name", value)

        if key.casefold().startswith("soil_"):
            insert_value(output, "environmentalSystem", "soil")

        if key.casefold() == "subject_id":
            is_part_of = {"identifier": value}
            insert_value(output, "isPartOf", is_part_of)

        if key.casefold() == "title":
            insert_value(output, "name", value)

        if re.match(r"ww_surv_target_\d+$", key.casefold()):
            insert_value(output, "infectiousAgent", {"name": value})

        if re.match(r"ww_surv_target_\d+_protocol$", key.casefold()):
            insert_value(output, "sampleProcess", value)

    if geo_loc:
        insert_value(output, "locationOfOrigin", geo_loc)
    if mat_loc:
        insert_value(output, "locationOfOrigin", mat_loc)


def parse():
    for sample_list in fetch_all_samples():
        for key, sample in sample_list.items():
            if key == "uids":
                continue
            if sample.get("accession"):
                _id = sample.get("accession")
                url = f"https://www.ncbi.nlm.nih.gov/biosample/{key}"
            else:
                continue

            # parse at the end
            xml_string = sample.pop("sampledata", None)

            output = {
                "@context": "http://schema.org/",
                "@type": "Sample",
                "_id": _id.casefold(),
                "identifier": _id,
                "url": url,
                "distribution": [{"@type": "dataDownload", "contentUrl": url}],
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "BioSample",
                    "url": "https://www.ncbi.nlm.nih.gov/biosample/",
                    "versionDate": datetime.date.today().isoformat(),
                    "archivedAt": url,
                },
            }

            if name := sample.get("title"):
                output["name"] = name

            if date := sample.get("date"):
                output["date"] = dateutil.parser.parse(date, ignoretz=True).date().isoformat()

            if date_published := sample.get("publicationdate"):
                output["datePublished"] = dateutil.parser.parse(date_published, ignoretz=True).date().isoformat()

            if date_modified := sample.get("modificationdate"):
                output["dateModified"] = dateutil.parser.parse(date_modified, ignoretz=True).date().isoformat()

            if name := sample.get("organization"):
                output["author"] = {"affiliation": {"@type": "Organization", "name": name}}

            if sample.get("taxonomy") or sample.get("organism"):
                species = {}
                if sample.get("taxonomy"):
                    species["identifier"] = sample.get("taxonomy")
                if sample.get("organism"):
                    species["name"] = sample.get("organism")

            if ids := sample.get("identifiers"):
                alternate_identifiers = []
                parts = ids.split("; ")
                for part in parts:
                    if ":" in part:
                        key, value = [x.strip() for x in part.split(":", 1)]
                        if key != "BioSample":
                            alternate_identifiers.append(value)
                if alternate_identifiers:
                    output["alternateIdentifier"] = alternate_identifiers

            sample_dict = xmltodict.parse(xml_string) if xml_string else {}
            if sample_dict:
                parse_xml(sample_dict, output)
            yield output
