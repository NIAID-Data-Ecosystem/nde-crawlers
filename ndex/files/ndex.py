import datetime
import logging

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


def process_networks(networks):
    for network in networks.get("networks", []):
        properties_dict = {}
        if properties := network.get("properties"):
            properties_dict = {prop["predicateString"]: prop["value"] for prop in properties}
        output = {
            "includedInDataCatalog": {
                "@type": "Dataset",
                "name": "NDEx",
                "url": "https://www.ndexbio.org/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Dataset",
        }
        if external_id := network.get("externalId"):
            if external_id in test_networks:
                logger.warning(f"Skipping test network {external_id}")
                continue
            output["identifier"] = external_id
            output["url"] = f"https://www.ndexbio.org/viewer/networks/{external_id}"
            output["_id"] = f"ndex_{external_id}"
        else:
            logger.warning("Network missing externalId")
            logger.warning(network)
            continue
        if network.get("isDeleted"):
            logger.warning(f"Network {external_id} is deleted")
            continue
        if name := network.get("name"):
            output["name"] = name
        if description := network.get("description"):
            output["description"] = description
        if creation_time := network.get("creationTime"):
            output["dateCreated"] = datetime.datetime.fromtimestamp(creation_time / 1000).isoformat()
        if modification_time := network.get("modificationTime"):
            output["dateModified"] = datetime.datetime.fromtimestamp(modification_time / 1000).isoformat()
        elif properties_modification_time := properties_dict.get("lastmodifieddate"):
            output["dateModified"] = properties_modification_time
        elif properties_modification_time := properties_dict.get("ndex:modificationTime"):
            output["dateModified"] = properties_modification_time
        if visibility := network.get("visibility"):
            if visibility == "PUBLIC":
                output["conditionsOfAccess"] = "Open"
            else:
                output["conditionsOfAccess"] = "Restricted"
        elif properties_visibility := properties_dict.get("visibility"):
            if properties_visibility is False:
                output["conditionsOfAccess"] = "Open"
        # if version := network.get("version"):
        #     output["version"] = version

        distribution = {}
        if cx_file_size := network.get("cxFileSize"):
            distribution["contentSize"] = cx_file_size
        if cx_format := network.get("cxFormat"):
            distribution["encodingFormat"] = cx_format
        if distribution:
            output["distribution"] = distribution

        author = {}
        if author_name := network.get("author"):
            author["name"] = author_name
        elif rights_holder := network.get("rightsHolder"):
            author["name"] = rights_holder
        elif owner := network.get("owner"):
            author["name"] = owner
        elif properties_author := properties_dict.get("bel:author"):
            author["name"] = properties_author
        elif properties_author_2 := properties_dict.get("Author"):
            author["name"] = properties_author_2
        if author:
            output["author"] = author

        health_condition_list = []
        if disease := network.get("disease"):
            health_condition_list.append({"name": disease})
        if properties_disease := properties_dict.get("diseases_id"):
            for disease in properties_disease:
                health_condition_list.append({"name": disease})
        if health_condition_list:
            output["healthCondition"] = health_condition_list

        species_list = []
        seen_species = set()
        if organism := network.get("organism"):
            if organism in seen_species:
                continue
            species_list.append({"name": organism})
            seen_species.add(organism)
        if species := network.get("species"):
            if species in seen_species:
                continue
            species_list.append({"name": species})
            seen_species.add(species)
        if idmapper_species := properties_dict.get("idmapper.species"):
            if idmapper_species in seen_species:
                continue
            species_list.append({"name": idmapper_species})
            seen_species.add(idmapper_species)
        if properties_species := properties_dict.get("species_common_name"):
            if properties_species in seen_species:
                continue
            species_list.append({"name": properties_species})
            seen_species.add(properties_species)
        if properties_organism := properties_dict.get("ORGANISM"):
            if properties_organism in seen_species:
                continue
            species_list.append({"name": properties_organism})
            seen_species.add(properties_organism)
        if species_list:
            output["species"] = species_list

        is_related_to_list = []
        if reference := network.get("reference"):
            is_related_to_list.append({"url": reference})
        properties_is_related_to = {}
        if figure_title := properties_dict.get("figureTitle"):
            properties_is_related_to["name"] = figure_title
        if figure_link := properties_dict.get("figureLink"):
            properties_is_related_to["url"] = figure_link
        if properties_is_related_to:
            is_related_to_list.append(properties_is_related_to)
        if uri := network.get("uri"):
            is_related_to_list.append({"url": uri})
        if is_related_to_list:
            output["isRelatedTo"] = is_related_to_list

        if rights := network.get("rights"):
            output["license"] = rights
        elif properties_license := properties_dict.get("license"):
            output["license"] = properties_license
        elif properties_copyright := properties_dict.get("bel:copyright"):
            output["license"] = properties_copyright

        if doi := network.get("doi"):
            output["doi"] = doi

        keywords = []
        if labels := network.get("labels"):
            keywords.extend(labels.split(","))
        if network_type := network.get("network type"):
            keywords.append(network_type)
        if network_type2 := network.get("networkType"):
            keywords.append(network_type2)
        if properties_keywords := properties_dict.get("dc:type"):
            keywords.append(properties_keywords)
        if keywords:
            output["keywords"] = keywords

        same_as_list = []
        if wiki_pathways_iri := network.get("wikipathwaysIRI"):
            same_as_list.append(wiki_pathways_iri)
        if properties_uri := properties_dict.get("URI"):
            same_as_list.append(properties_uri)
        if kegg_pathways := properties_dict.get("KEGG_PATHWAY_LINK"):
            same_as_list.append(kegg_pathways)
        if same_as_list:
            output["sameAs"] = same_as_list

        if was_derived_from := network.get("prov:wasDerivedFrom"):
            output["wasDerivedFrom"] = was_derived_from

        citation = {}
        if pmcid := properties_dict.get("pmcid"):
            citation["pmcid"] = pmcid
        if paper_title := properties_dict.get("paperTitle"):
            citation["name"] = paper_title
        if paper_link := properties_dict.get("paperLink"):
            citation["url"] = paper_link
        if citation:
            output["citation"] = citation

        sd_publisher_list = []
        if data_source := properties_dict.get("Data source"):
            sd_publisher_list.append({"name": data_source})
        if tcga_data_source := properties_dict.get("TCGA Data Source"):
            sd_publisher_list.append({"name": tcga_data_source})
        if properties_source := properties_dict.get("source"):
            sd_publisher_list.append({"name": properties_source})
        if sd_publisher_list:
            output["sdPublisher"] = sd_publisher_list

        if date := properties_dict.get("dc:date"):
            date_formatted = date.split(" ")[0]
            output["date"] = date_formatted

        is_based_on_list = []
        if properties_data_source := properties_dict.get("dataSource"):
            is_based_on_list.append({"url": properties_data_source})
        if properties_source := properties_dict.get("Source"):
            is_based_on_list.append({"url": properties_source})
        if is_based_on_list:
            output["isBasedOn"] = is_based_on_list

        if properties_treatment := properties_dict.get("Treatment"):
            output["variableMeasured"] = {"name": properties_treatment}

        if properties_method := properties_dict.get("methods"):
            if properties_method in technique_lookup:
                output["measurementTechnique"] = technique_lookup[properties_method]

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
