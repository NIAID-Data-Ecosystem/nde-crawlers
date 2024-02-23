# # what dde gives us https://www.ebi.ac.uk/ols/api/terms?iri=https://www.ontobee.org/ontology/NCBITaxon?iri=http://purl.obolibrary.org/obo/NCBITaxon_9554
# def get_uniprot_details(url):
#     logger.info(f"Getting details for {url}")
#     # get id from url

#     # Fetch details from the UniProt API
#     identifier = identifier.split("*")[-1]
#     # try:
#     species_info = requests.get(f"https://rest.uniprot.org/taxonomy/{identifier}")
#     species_info.raise_for_status()
#     species_info = species_info.json()
#     standard_dict = {
#         "identifier": identifier,
#         "inDefinedTermSet": "UniProt",
#         "url": f"https://www.uniprot.org/taxonomy/{identifier}",
#         "originalName": original_name,
#         "isCurated": True,
#         "curatedBy": {
#             "name": "PubTator",
#             "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
#             "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
#         },
#     }
#     if scientific_name := species_info.get("scientificName"):
#         standard_dict["name"] = scientific_name
#     else:
#         standard_dict["name"] = original_name

#     alternative_names = []
#     if common_name := species_info.get("commonName"):
#         standard_dict["commonName"] = common_name
#         alternative_names.append(common_name)

#         standard_dict["displayName"] = f"{common_name} | {scientific_name}"

#     if other_names := species_info.get("otherNames"):
#         alternative_names.extend(other_names)

#     if alternative_names:
#         standard_dict["alternateName"] = alternative_names

#     if lineage := species_info.get("lineage"):
#         standard_dict["classification"] = classify_as_host_or_agent(lineage)
#     else:
#         logger.warning(f"No lineage found for {identifier}")

#     return standard_dict

# {
#     "species": {
#         "http://purl.obolibrary.org/obo/NCBITaxon_": "UniProt"
#     },
#     "infectiousAgent": {
#         "http://purl.obolibrary.org/obo/NCBITaxon_": "UniProt"
#     },

# }

# import requests
# from biothings.utils.dataload import tab2dict

# https://docs.google.com/spreadsheets/d/107WVX39r_a6xBGZ_gCku0LBNRWWBmk9x7Dg53wj1SiI/edit#gid=0


# result = requests.get("https://docs.google.com/spreadsheets/d/107WVX39r_a6xBGZ_gCku0LBNRWWBmk9x7Dg53wj1SiI/export?format=csv")

import csv
from io import StringIO

import requests

response = requests.get(
    "https://docs.google.com/spreadsheets/d/107WVX39r_a6xBGZ_gCku0LBNRWWBmk9x7Dg53wj1SiI/export?format=csv"
).text


def process_csv_data(csv_content):
    # Prepare the CSV data for reading
    csv_file = StringIO(csv_content)
    csv_reader = csv.DictReader(csv_file)

    # Initialize a dictionary to hold the processed data
    properties = {}

    # Variable to keep track of the current property being processed
    current_property = None

    for row in csv_reader:
        if ".inDefinedTermSet" in row["Old button text"]:
            # Start of a new property section
            current_property = row["Old button text"].split(".")[0]
            properties[current_property] = {}
        elif row["base url (for mapping purposes)"] and current_property:
            if "http" not in row["base url (for mapping purposes)"]:
                continue
            # Row belongs to the current property
            base_url = row["base url (for mapping purposes)"]
            # The script assumes the 'InDefinedTermSet' field to hold the intended value
            value = row["InDefinedTermSet"]
            properties[current_property][base_url] = value

    return properties


def species_func(species, term_set):
    if term_set == "UniProt":
        species["name"] = "WE GOT THIS WORKING"


def get_term_set_helper(url, property_dict):
    for prop, term_set in property_dict.items():
        if prop in url:
            return term_set
    else:
        return None


def get_in_defined_term_set(doc, properties_dict):
    # http://purl.obolibrary.org/obo/NCBITaxon_9554
    # https://www.geeksforgeeks.org/python-store-function-as-dictionary-value/
    # nde_properties = {
    #     "species": species_func,
    #     "infectiousAgent": infectious_agent_func,
    #     "healthCondition": health_condition_func,
    #     "measurementTechnique": measurement_technique_func,
    #     "variableMeasured": variable_measured_func,
    #     "keywords": keywords_func,
    #     "topicCategory": topic_category_func,
    # }

    nde_properties = {
        "species": species_func,
    }

    for nde_property, func in nde_properties.items():
        if doc_property := doc.get(nde_property):
            if isinstance(doc_property, list):
                for item in doc_property:
                    if url := item.get("url"):
                        if term_set := get_term_set_helper(url, properties_dict[nde_property]):
                            func(item, term_set)
            else:
                if url := doc_property.get("url"):
                    if term_set := get_term_set_helper(url, properties_dict[nde_property]):
                        func(doc_property, term_set)
