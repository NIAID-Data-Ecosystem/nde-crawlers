import datetime
import logging
import re

import ndex2.client as nc
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

technique_lookup = {
    "Pathway Figure OCR": {"name": "image processing", "url": "http://purl.obolibrary.org/obo/NCIT_C17606"},
    "https://www.biorxiv.org/content/10.1101/2020.05.29.124503v1.full": {
        "name": "image processing",
        "url": "http://purl.obolibrary.org/obo/NCIT_C17606",
    },
    "OSLOM clustering": {"name": "cluster analysis", "url": "http://purl.obolibrary.org/obo/NCIT_C63918"},
    "Database Integration": {"name": "data integration and warehousing", "url": "http://edamontology.org/topic_3366"},
    "Literature Curation": {"name": "curation", "url": "http://purl.obolibrary.org/obo/NCIT_C48292"},
    "NCMine (Tadaka &  Kinoshita, 2016)": {
        "name": "cluster analysis",
        "url": "http://purl.obolibrary.org/obo/NCIT_C63918",
    },
    "STRING": {"name": "string", "url": "http://purl.obolibrary.org/obo/MI_1014"},
    "coexpression analysis": {"name": "coexpression", "url": "http://purl.obolibrary.org/obo/MI_2231"},
    "FAIRE-seq": {"name": "FAIRE-seq", "url": "http://www.ebi.ac.uk/efo/EFO_0004428"},
    "STRING v11 (Szklarczyk et al 2019)": {"name": "string", "url": "http://purl.obolibrary.org/obo/MI_1014"},
    "NCMine": {"name": "cluster analysis", "url": "http://purl.obolibrary.org/obo/NCIT_C63918"},
    "R": {"name": "computational method", "url": "http://www.bioassayontology.org/bao#BAO_0002094"},
    "RCy3": {"name": "computational method", "url": "http://www.bioassayontology.org/bao#BAO_0002094"},
    "SPARQL": {"name": "computational method", "url": "http://www.bioassayontology.org/bao#BAO_0002094"},
    "Cytoscape": {"name": "cytoscape network analysis", "url": "http://www.bioassayontology.org/bao#BAO_0002362"},
    "computationally inferred PPI interactome using an interolog and domain-based approach": {
        "name": "computational method",
        "url": "http://www.bioassayontology.org/bao#BAO_0002094",
    },
    "CRISPR screen": {"name": "CRISPR/Cas9 method", "url": "http://www.bioassayontology.org/bao#BAO_0010249"},
    "Bayesian network (SiGN-BN)": {"name": "Bayesian approach", "url": "http://purl.obolibrary.org/obo/NCIT_C142403"},
    "Bayesian": {"name": "Bayesian approach", "url": "http://purl.obolibrary.org/obo/NCIT_C142403"},
    "Cytoscape v3.9": {"name": "cytoscape network analysis", "url": "http://www.bioassayontology.org/bao#BAO_0002362"},
    "CRISPR screening": {"name": "CRISPR/Cas9 method", "url": "http://www.bioassayontology.org/bao#BAO_0010249"},
}


VALID_NETWORK_SETS = {
    "bdba6a7a-488a-11ec-b3be-0ac135e8bacf",  # INDRA-GO
    "453c1c63-5c10-11e9-9f06-0ac135e8bacf",  # WikiPathways (Homo sapiens)
    "0e7da362-e594-11ef-8e41-005056ae3c32",  # Gene Ontology (GO-CAMs)
    "224d4de6-e23f-11ea-99da-0ac135e8bacf",  # HiDef
    "78db519f-eb1d-11eb-b666-0ac135e8bacf",  # Swaney (Science 2021)
    "0bd21f41-cd02-11ea-aaef-0ac135e8bacf",  # Zhang (bioinformatics 2018)
    "2c200d22-e2fe-11ea-99da-0ac135e8bacf",  # Patrick (cell systems 2018)
}


def get_valid_network_ids():
    """
    Iterates over each network set in VALID_NETWORK_SETS, calls the NDEx API,
    and collects network IDs from the 'networks' array in each response.

    Returns:
        A list of network IDs.
    """
    base_url = "https://www.ndexbio.org/v2/networkset/"
    valid_networks = []

    for network_set_id in VALID_NETWORK_SETS:
        if network_set_id == "78db519f-eb1d-11eb-b666-0ac135e8bacf":
            # Zhang (bioinformatics 2018) is just one network
            valid_networks.append("78db519f-eb1d-11eb-b666-0ac135e8bacf")
            continue
        url = base_url + network_set_id
        logger.info(f"Fetching network set: {network_set_id} from {url}")
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                valid_networks.extend(data.get("networks", []))
                logger.info(f"Fetched {len(data.get('networks', []))} networks from {url}")
            else:
                logger.error(f"Failed to fetch {url}: HTTP {response.status_code}")
        except Exception as e:
            logger.error(f"Exception occurred while fetching {url}: {e}")

    return valid_networks


def strip_html_tags(text: str) -> str:
    """Return text with any HTML tags removed."""
    return re.sub(r"<[^>]*>", "", text)


def convert_to_pmid(identifier: str) -> str:
    """
    Given an identifier (PMID, DOI, or URL), use the NCBI ID Converter API
    to convert it to a PMID if possible.

    If the identifier is already numeric (assumed to be a PMID), it's returned as is.
    """
    # If already a numeric PMID, return it
    if identifier.isdigit():
        return identifier

    base_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
    email = os.environ.get("NDEX_EMAIL")
    params = {"ids": identifier, "format": "json", "tool": "niaid-data-ecosystem", "email": email}
    try:
        logger.info(f"Converting identifier {identifier} to PMID")
        response = requests.get(base_url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            records = data.get("records", [])
            logger.info(f"Received {len(records)} records for identifier {identifier}")
            if records:
                logger.info(f"Records: {records}")
                # The API returns a list of records; we take the first one.
                record = records[0]
                pmid = record.get("pmid")
                if pmid:
                    logger.info(f"Converted {identifier} to PMID: {pmid}")
                    return pmid
    except Exception as e:
        print(f"Error converting identifier {identifier}: {e}")
    return None


def process_reference(ref: str) -> tuple:
    """
    Processes a reference string to extract identifiers and converts them
    to PMIDs (if possible). Also extracts any extra URL or non-PubMed identifiers.
    Returns a tuple of PMIDs (as a comma-separated string) and any extra URL found.
    """
    pmid_set = set()
    extra_url = None
    logger.info(f"Processing reference: {ref}")

    # Capture a PubMed link (supporting both NCBI and pubmed domains).
    pubmed_url_match = re.search(
        r'href=["\']https?://(?:www\.)?(?:ncbi\.nlm\.nih\.gov/pubmed|pubmed\.ncbi\.nlm\.nih\.gov)/(\d+)(?:/?)["\']',
        ref,
        re.IGNORECASE,
    )
    if pubmed_url_match:
        found_pmid = pubmed_url_match.group(1)
        logger.info(f"Found PubMed URL PMID: {found_pmid}")
        pmid_set.add(found_pmid)

    # Capture a PMC link from the href attribute.
    pmc_url_match = re.search(
        r'href=["\']https?://(?:www\.)?ncbi\.nlm\.nih\.gov/pmc/articles/PMC(\d+)(?:/?)["\']', ref, re.IGNORECASE
    )
    if pmc_url_match:
        found_pmcid = pmc_url_match.group(1)
        logger.info(f"Found PMC URL PMCID: {found_pmcid}")
        pmid_from_pmcid = convert_to_pmid(found_pmcid)
        if pmid_from_pmcid:
            pmid_set.add(pmid_from_pmcid)

    # Strip HTML tags
    clean_ref = re.sub(r"<[^>]+>", "", ref)

    # Look for a textual PMID or pubmed identifier.
    pmid_match = re.search(r"(?:PMID|pubmed)[:\s]*(\d+)", clean_ref, re.IGNORECASE)
    if pmid_match:
        found_pmid = pmid_match.group(1)
        logger.info(f"Found textual PMID: {found_pmid}")
        pmid_set.add(found_pmid)

    # Look for a textual PMCID.
    pmcid_match = re.search(r"PMCID[:\s]*PMC(\d+)", clean_ref, re.IGNORECASE)
    if pmcid_match:
        logger.info(f"Found textual PMCID: {pmcid_match.group(1)}")
        pmcid = pmcid_match.group(1)
        pmid_from_pmcid = convert_to_pmid(pmcid)
        if pmid_from_pmcid:
            pmid_set.add(pmid_from_pmcid)

    # Look for a DOI (allowing for URL-encoded characters and optional dx.).
    doi_match = re.search(
        r"(?:https?://(?:dx\.)?(?:doi\.org)/)?(10\.\d{4,9}/[-._;()/:A-Z0-9%]+)", clean_ref, re.IGNORECASE
    )
    if doi_match:
        found_doi = doi_match.group(1)
        logger.info(f"Found DOI: {found_doi}")
        pmid_from_doi = convert_to_pmid(found_doi)
        if pmid_from_doi:
            pmid_set.add(pmid_from_doi)

    # Look for any other URLs in the cleaned reference.
    url_match = re.search(r'(https?://[^\s"<]+)', clean_ref, re.IGNORECASE)
    if url_match:
        url = url_match.group(1)
        # Exclude URLs for DOI, PMC, or PubMed.
        if not any(x in url.lower() for x in ["doi", "pmc", "pubmed"]):
            logger.info(f"Found extra URL: {url}")
            extra_url = url

    return ",".join(sorted(pmid_set)), extra_url


def process_networks(networks, valid_network_ids):
    for network in networks.get("networks", []):
        properties_dict = {}
        if properties := network.get("properties"):
            properties_dict = {prop["predicateString"]: prop["value"] for prop in properties}

        # Helper to search for keys in both network and its properties
        def get_value(*keys):
            for key in keys:
                if value := network.get(key):
                    return value
                elif value := properties_dict.get(key):
                    return value
            return None

        def add_species(raw_value, seen, species_list):
            if not raw_value:
                return
            for token in raw_value.split(","):
                token = token.strip()
                if not token or token.isdigit():
                    continue
                if "human" in token.lower():
                    token = "human"
                if token not in seen:
                    species_list.append({"name": token})
                    seen.add(token)

        output = {
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "NDEx",
                "url": "https://www.ndexbio.org/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Dataset",
        }

        if external_id := get_value("externalId"):
            output["identifier"] = external_id
            output["url"] = f"https://www.ndexbio.org/viewer/networks/{external_id}"
            output["includedInDataCatalog"]["dataset"] = f"https://www.ndexbio.org/viewer/networks/{external_id}"
            output["_id"] = f"ndex_{external_id}"
        else:
            logger.warning("Network missing externalId")
            logger.warning(network)
            continue

        if get_value("isDeleted"):
            logger.warning(f"Network {external_id} is deleted")
            continue

        if name := get_value("name"):
            output["name"] = name
        if description := get_value("description"):
            output["description"] = description

        if creation_time := get_value("creationTime"):
            output["dateCreated"] = datetime.datetime.fromtimestamp(int(creation_time) / 1000).isoformat()

        if modification_time := get_value("modificationTime"):
            output["dateModified"] = datetime.datetime.fromtimestamp(int(modification_time) / 1000).isoformat()
        elif properties_modification_time := get_value("lastmodifieddate", "ndex:modificationTime"):
            output["dateModified"] = properties_modification_time

        if visibility := get_value("visibility"):
            output["conditionsOfAccess"] = "Open" if visibility == "PUBLIC" else "Restricted"
        elif properties_visibility := get_value("visibility"):
            if properties_visibility is False:
                output["conditionsOfAccess"] = "Open"

        distribution = {}
        if cx_file_size := get_value("cxFileSize"):
            distribution["contentSize"] = cx_file_size
        if cx_format := get_value("cxFormat"):
            distribution["encodingFormat"] = cx_format
        if distribution:
            output["distribution"] = distribution

        author = {}
        if author_name := get_value("author", "rightsHolder", "owner", "bel:author", "Author"):
            author["name"] = author_name
        if author:
            output["author"] = author

        health_condition_list = []
        if disease := get_value("disease"):
            disease_text = strip_html_tags(disease)
            for item in disease_text.split(","):
                cleaned_name = item.strip()
                if cleaned_name:
                    health_condition_list.append({"name": cleaned_name})
        if properties_disease := get_value("diseases_id"):
            for disease_item in properties_disease:
                cleaned_name = strip_html_tags(disease_item).strip()
                if cleaned_name:
                    health_condition_list.append({"name": cleaned_name})
        if health_condition_list:
            output["healthCondition"] = health_condition_list

        species_list = []
        seen_species = set()
        for field in ["organism", "species", "idmapper.species", "species_common_name", "ORGANISM"]:
            value = get_value(field)
            add_species(value, seen_species, species_list)
        if species_list:
            output["species"] = species_list

        if reference := get_value("reference"):
            pmids, extra_url = process_reference(reference)
            if pmids:
                output["pmids"] = pmids
            if extra_url:
                output["citation"] = {
                    "url": extra_url,
                    "description": f"{output['name']} is found accessible at {extra_url}",
                }

        # Process other related fields (still added as isRelatedTo)
        is_related_to_list = []
        properties_is_related_to = {}
        if figure_title := get_value("figureTitle"):
            properties_is_related_to["name"] = figure_title
        if figure_link := get_value("figureLink"):
            properties_is_related_to["url"] = figure_link
        if properties_is_related_to:
            is_related_to_list.append(properties_is_related_to)
        if uri := get_value("uri"):
            is_related_to_list.append({"url": uri})
        if is_related_to_list:
            output["isRelatedTo"] = is_related_to_list

        if rights := get_value("rights"):
            output["license"] = rights
        elif properties_license := get_value("license"):
            output["license"] = properties_license
        elif properties_copyright := get_value("bel:copyright"):
            output["license"] = properties_copyright

        if doi := get_value("doi"):
            output["doi"] = doi

        keywords = []
        if labels := get_value("labels"):
            if isinstance(labels, list):
                keywords.extend(labels.split(","))
            else:
                keywords.append(labels)
        if network_type := get_value("network type"):
            if isinstance(network_type, list):
                keywords.extend(network_type)
            else:
                keywords.append(network_type)
        if network_type2 := get_value("networkType"):
            if isinstance(network_type2, list):
                keywords.extend(network_type2)
            else:
                keywords.append(network_type2)
        if properties_keywords := get_value("dc:type"):
            if isinstance(properties_keywords, list):
                keywords.extend(properties_keywords)
            else:
                keywords.append(properties_keywords)
        if keywords:
            output["keywords"] = keywords

        same_as_list = []
        if wiki_pathways_iri := get_value("wikipathwaysIRI"):
            same_as_list.append(wiki_pathways_iri)
        if properties_uri := get_value("URI"):
            same_as_list.append(properties_uri)
        if kegg_pathways := get_value("KEGG_PATHWAY_LINK"):
            same_as_list.append(kegg_pathways)
        if same_as_list:
            output["sameAs"] = same_as_list

        sd_publisher_list = []
        if data_source := get_value("Data source"):
            sd_publisher_list.append({"name": data_source})
        if tcga_data_source := get_value("TCGA Data Source"):
            sd_publisher_list.append({"name": tcga_data_source})
        if properties_source := get_value("source"):
            sd_publisher_list.append({"name": properties_source})
        if sd_publisher_list:
            output["sdPublisher"] = sd_publisher_list

        if date := get_value("dc:date"):
            date_formatted = date.split(" ")[0]
            output["date"] = date_formatted

        is_based_on_dict = {}
        if properties_data_source := get_value("dataSource"):
            is_based_on_dict["url"] = properties_data_source
        if properties_source := get_value("Source"):
            is_based_on_dict["url"] = properties_source
        if was_derived_from := get_value("prov:wasDerivedFrom"):
            if "http" in was_derived_from:
                is_based_on_dict["url"] = was_derived_from
            elif "ftp://" in was_derived_from:
                ftp_url = re.search(r"ftp://[^ ]+", was_derived_from).group()
                is_based_on_dict["url"] = ftp_url
            else:
                is_based_on_dict["name"] = was_derived_from
        if is_based_on_dict:
            output["isBasedOn"] = is_based_on_dict

        if properties_treatment := get_value("Treatment"):
            output["variableMeasured"] = {"name": properties_treatment}

        measurement_technique_list = []
        measurement_technique_identifiers = set()
        if properties_method := get_value("methods"):
            methods = re.split(r",(?![^()]*\))", properties_method)
            for method in methods:
                method = method.strip()
                if method in technique_lookup:
                    measurement_technique = technique_lookup[method]
                    identifier = (measurement_technique["name"], measurement_technique["url"])
                    if identifier not in measurement_technique_identifiers:
                        measurement_technique_list.append(measurement_technique)
                        measurement_technique_identifiers.add(identifier)
        if measurement_technique_list:
            output["measurementTechnique"] = measurement_technique_list

        # Default topicCategory remains unchanged.
        output["topicCategory"] = [
            {
                "url": "http://edamontology.org/topic_0602",
                "identifier": "topic_0602",
                "name": "Molecular interactions, pathways and networks",
                "inDefinedTermSet": "EDAM",
            }
        ]

        # Only yield the record if it comes from one of the valid network sets or has a citation/pmid.
        if output["identifier"] not in valid_network_ids and ("citation" not in output or "pmids" not in output):
            logger.warning(
                f"Skipping network {output.get('identifier')} because it is not in a valid network set and has no reference"
            )
            continue

        yield output


def parse():
    logger.info("Starting NDEx parsing")
    # Initialize the NDEx client
    anon_ndex = nc.Ndex2("http://public.ndexbio.org")
    anon_ndex.update_status()

    def fetch_networks(start, size):
        return anon_ndex.search_networks(start=start, size=size)

    status = anon_ndex.status
    networks_count = status.get("networkCount")
    users_count = status.get("userCount")
    groups_count = status.get("groupCount")
    logger.info(f"anon client: {networks_count} networks, {users_count} users, {groups_count} groups")

    start = 0
    size = 1000
    networks = anon_ndex.search_networks(start=start, size=size)
    total_networks = networks.get("numFound", 0)
    logger.info(f"Found {total_networks} networks")

    valid_network_ids = set(get_valid_network_ids())
    logger.info(f"Found {len(valid_network_ids)} valid network IDs")

    while networks.get("networks"):
        for network in process_networks(networks, valid_network_ids):
            yield network
        start += 1
        logger.info(f"Retrieved {start * size} networks")
        networks = fetch_networks(start, size)
