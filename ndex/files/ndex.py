import datetime
import logging
import re

import ndex2.client as nc

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

test_networks = [
    "89274295-1730-11e7-b39e-0ac135e8bacf",
    "e6bf9a50-b666-11ee-8a13-005056ae23aa",
    "792f0c5c-22b6-11ea-bb65-0ac135e8bacf",
    "b164f77e-499c-11e8-a4bf-0ac135e8bacf",
]


def strip_html_tags(text: str) -> str:
    """Return text with any HTML tags removed."""
    return re.sub(r"<[^>]*>", "", text)


def process_networks(networks):
    for network in networks.get("networks", []):
        properties_dict = {}
        if properties := network.get("properties"):
            properties_dict = {prop["predicateString"]: prop["value"] for prop in properties}

        # Helper function to get value from network or properties_dict
        def get_value(*keys):
            for key in keys:
                if value := network.get(key):
                    return value
                elif value := properties_dict.get(key):
                    return value
            return None

        output = {
            "includedInDataCatalog": {
                "@type": "Dataset",
                "name": "NDEx",
                "url": "https://www.ndexbio.org/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Dataset",
        }

        # Use the helper function to get external_id
        if external_id := get_value("externalId"):
            if external_id in test_networks:
                logger.warning(f"Skipping test network {external_id}")
                continue
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
            if visibility == "PUBLIC":
                output["conditionsOfAccess"] = "Open"
            else:
                output["conditionsOfAccess"] = "Restricted"
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
        if organism := get_value("organism"):
            if "human" in organism.lower():
                organism = "human"
            if organism not in seen_species:
                species_list.append({"name": organism})
                seen_species.add(organism)
        if species := get_value("species"):
            if "human" in species.lower():
                species = "human"
            if species not in seen_species:
                species_list.append({"name": species})
                seen_species.add(species)
        if idmapper_species := get_value("idmapper.species"):
            if "human" in idmapper_species.lower():
                idmapper_species = "human"
            if idmapper_species not in seen_species:
                species_list.append({"name": idmapper_species})
                seen_species.add(idmapper_species)
        if properties_species := get_value("species_common_name"):
            if "human" in properties_species.lower():
                properties_species = "human"
            if properties_species not in seen_species:
                species_list.append({"name": properties_species})
                seen_species.add(properties_species)
        if properties_organism := get_value("ORGANISM"):
            if "human" in properties_organism.lower():
                properties_organism = "human"
            if properties_organism not in seen_species:
                species_list.append({"name": properties_organism})
                seen_species.add(properties_organism)
        if species_list:
            output["species"] = species_list

        is_related_to_list = []
        if reference := get_value("reference"):
            is_related_to_list.append({"url": reference})
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

        citation = {}
        if pmcid := get_value("pmcid"):
            citation["pmcid"] = pmcid
        if paper_title := get_value("paperTitle"):
            citation["name"] = paper_title
        if paper_link := get_value("paperLink"):
            citation["url"] = paper_link
        if citation:
            output["citation"] = citation

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

        # Default topicCategory to "Molecular interactions, pathways and networks"
        # https://github.com/NIAID-Data-Ecosystem/nde-crawlers/issues/152#issuecomment-2597016173
        output["topicCategory"] = [
            {
                "url": "http://edamontology.org/topic_0602",
                "identifier": "topic_0602",
                "name": "Molecular interactions, pathways and networks",
                "inDefinedTermSet": "EDAM",
            }
        ]

        yield output


def parse():
    logger.info("Starting NDEx parsing")
    # Initialize the NDEx client
    anon_ndex = nc.Ndex2("http://public.ndexbio.org")
    anon_ndex.update_status()

    def fetch_networks(start, size):
        return anon_ndex.search_networks(start=start, size=size)

    # Fetch general status information
    status = anon_ndex.status
    networks_count = status.get("networkCount")
    users_count = status.get("userCount")
    groups_count = status.get("groupCount")
    logger.info(f"anon client: {networks_count} networks, {users_count} users, {groups_count} groups")

    # Initialize variables
    start = 0
    size = 1000

    networks = anon_ndex.search_networks(start=start, size=size)
    total_networks = networks.get("numFound", 0)
    logger.info(f"Found {total_networks} networks")

    while networks.get("networks"):
        for network in process_networks(networks):
            yield network
        start += 1
        logger.info(f"Retrieved {start * size} networks")
        networks = fetch_networks(start, size)
