import datetime
import logging
import re
import time

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def query_ols(iri):
    """Gets the name field of measurementTechnique, infectiousAgent, healthCondition (infectiousDisease), species, and variableMeasured in our nde schema

    ols api doc here: https://www.ebi.ac.uk/ols4/swagger-ui/index.html
    Returns the formatted dictionary {name: ####, url: ####} if an url was given or {name: ####}
    """

    url = "https://www.ebi.ac.uk/ols/api/terms?"
    pattern = re.compile("^https?://")

    if pattern.match(iri):
        params = {
            # isnt the best way but good enough to get first instance of http while ignoring https:
            # https://stackoverflow.com/questions/9760588/how-do-you-extract-a-url-from-a-string-using-python
            "iri": re.search(r"(?P<url>http?://[^\s]+)", iri).group("url")
        }

        request = requests.get(url, params).json()
        # no documentation on how many requests can be made
        time.sleep(0.5)

        lookup = {
            "name": request["_embedded"]["terms"][0]["label"],
            "url": iri,
            "inDefinedTermSet": request["_embedded"]["terms"][0]["ontology_prefix"],
            "curatedBy": {
                "name": "Data Discovery Engine",
                "url": "https://discovery.biothings.io/",
                "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
            },
            "isCurated": True,
        }
        if alternate_name := request["_embedded"]["terms"][0]["annotation"].get("alternative label"):
            lookup["alternateName"] = alternate_name
        elif alternate_name := request["_embedded"]["terms"][0].get("synonyms"):
            lookup["alternateName"] = alternate_name
        else:
            logger.info(f"Does not have an alternateName: {iri}")

        return lookup, True
    else:
        return iri, False


def format_query_ols(iri):
    """
    Formats the query_ols return value to fit our schema for measurementTechnique, infectiousAgent, infectiousDisease, and species
    Returns the formatted dictionary {name: ####, url: ####} if an url was given or {name: ####}
    """
    info, is_url = query_ols(iri)
    if is_url:
        return info
    else:
        return {"name": info}


def format_query_ky(iri):
    """
    Formats the query_ols return value to fit our schema for variableMeasured and keywords
    Returns the formatted keyword
    """
    info, is_url = query_ols(iri)
    if is_url:
        return info["name"]
    else:
        return info


def parse():
    # initial request to find total number of hits
    url = "https://discovery.biothings.io/api/dataset/query?q=_meta.guide:%22/guide/niaid/ComputationalTool%22%20OR%20_meta.guide:%22/guide/niaid%22%20OR%20_meta.guide:%22/guide/nde/ResourceCatalog%22%20OR%20_meta.guide:%22/guide/creid%22&sort=-_ts.last_updated"
    request = requests.get(url).json()
    # get the number of pages to paginate through
    total = request["total"]
    pages = (total - 1) // 1000
    count = 0
    logger.info("Starting requests...")
    # paginate through the requests
    for page in range(pages + 1):
        url = (
            "https://discovery.biothings.io/api/dataset/query?q=_meta.guide:%22/guide/niaid/ComputationalTool%22%20OR%20_meta.guide:%22/guide/niaid%22%20OR%20_meta.guide:%22/guide/nde/ResourceCatalog%22%20OR%20_meta.guide:%22/guide/creid%22&sort=-_ts.last_updated&size=1000&from="
            + str(page * 1000)
        )
        request = requests.get(url).json()
        for hit in request["hits"]:
            count += 1
            if count % 100 == 0:
                logger.info("Parsed %s records", count)

            # adjust @type value to fit our schema
            if nde_type := hit.pop("@type", None):
                nde_type = nde_type.split(":")[-1]
                if "Dataset" in nde_type:
                    nde_type = "Dataset"
                elif "ComputationalTool" in nde_type:
                    nde_type = "ComputationalTool"

                hit["@type"] = nde_type

            included_in_data_catalog = [
                {
                    "@type": "DataCatalog",
                    "name": "Data Discovery Engine",
                    "url": "https://discovery.biothings.io/",
                    "versionDate": datetime.date.today().isoformat(),
                }
            ]

            # add sourceOrganization (as determined by the DDE portal to which the record was submitted) - this can be determined based on the @context
            if hit.get("@context") and "nde" in hit.get("@context"):
                source_organization = [
                    {
                        "@type": "ResearchProject",
                        "name": "NIAID Data Ecosystem",
                        "url": "https://data.niaid.nih.gov",
                        "parentOrganization": "NIAID"
                    }
                ]

            # add creid to included in sourceOrganization
            if hit.get("@context") and "creid" in hit.get("@context"):
                source_organization = [
                    {
                        "@type": "ResearchProject",
                        "name": "NIAID CREID Network",
                        "description": "The Centers for Research in Emerging Infectious Diseases (CREID) Network is a coordinated group of emerging infectious disease research centers situated in regions around the globe where emerging and re-emerging infectious disease outbreaks are likely to occur.",
                        "alternateName": ["CREID","Centers for Research in Emerging Infectious Disease"],
                        "url": "https://creid-network.org/",
                        "parentOrganization": "NIAID"
                    },
                    {
                        "@type": "ResearchProject",
                        "name": "NIAID Data Ecosystem",
                        "url": "https://data.niaid.nih.gov",
                        "parentOrganization": "NIAID"
                    }
                ]

            # all of niaid systems biology is a subset of niaid data ecosystem but if nde is in the context then it is not part of niaid systems biology
            # if the above comment confuses the reader, think of it like the story of oedipus
            if (
                hit.get("@context")
                and "niaid" in hit.get("@context")
                and "nde" not in hit.get("@context")
                and "creid" not in hit.get("@context")
            ):
                source_organization = [
                    {
                        "@type": "ResearchProject",
                        "name": "NIAID Systems Biology",
                        "description": "The NIAID/Division of Microbiology and Infectious Diseases (DMID) Systems Biology Consortium for Infectious Diseases is a group of interdisciplinary scientists that bridge disparate scientific disciplines including microbiology, immunology, infectious diseases, microbiome, mathematics, physics, bioinformatics, computational biology, machine learning, statistical methods, and mathematical modeling.",
                        "alternateName": ["NIAID Systems Biology Consortium for Infectious Diseases", "NIAID SysBio"],
                        "url": "https://www.niaid.nih.gov/research/systems-biology-consortium",
                        "parentOrganization": "NIAID"
                    },
                    {
                        "@type": "ResearchProject",
                        "name": "NIAID Data Ecosystem",
                        "url": "https://data.niaid.nih.gov",
                        "parentOrganization": "NIAID"
                    }
                ]

            hit["includedInDataCatalog"] = included_in_data_catalog
            hit["sourceOrganization"] = source_organization

            if citations := hit.get("citation"):
                if not isinstance(citations, list):
                    citations = [citations]
                for citation in citations:
                    if citation.get("@type") == "ScholarlyArticle":
                        if pmid := citation.get("pmid"):
                            hit["pmids"] = (
                                hit.get("pmids") + "," + str(pmid).lstrip("0") if hit.get("pmids") else str(pmid)
                            )

                # Use list comprehension to filter out citations with a PMID
                hit["citation"] = [
                    citation
                    for citation in citations
                    if not (citation.get("pmid") and citation.get("@type") == "ScholarlyArticle")
                ]
                if not hit["citation"]:
                    del hit["citation"]

            # remove nonetypes from distribution
            if distribution := hit.pop("distribution", None):
                if isinstance(distribution, list):
                    hit["distribution"] = [{k: v for k, v in d.items() if v is not None} for d in distribution]
                else:
                    hit["distribution"] = {k: v for k, v in distribution.items() if v is not None}

            # rename our id value and creator to author
            if authors := hit.pop("creator", None):
                if type(authors) is list:
                    for author in authors:
                        if affiliation := author.get("affiliation"):
                            author["affiliation"] = {"name": affiliation}
                else:
                    if affiliation := authors.get("affiliation"):
                        authors["affiliation"] = {"name": affiliation}
                hit["author"] = authors

            # Adjust URL based on type
            if hit["@type"] == "Dataset":
                hit["url"] = "https://discovery.biothings.io/dataset/" + hit["_id"]

            hit["_id"] = "dde_" + hit["_id"].lower()

            # adjust date values
            if dates := hit.pop("_ts", None):
                hit["dateCreated"] = datetime.datetime.fromisoformat(dates["date_created"]).date().isoformat()
                hit["dateModified"] = datetime.datetime.fromisoformat(dates["last_updated"]).date().isoformat()

            # WE DONT NEED TO CHANGE THIS
            # adjust applicationSubCategory to fit our schema
            # if app_subs := hit.pop('applicationSubCategory', None):
            # hit['applicationSubCategory'] = []
            # for app_sub in app_subs:
            #     hit['applicationSubCategory'].append(
            #         {'name': app_sub.get('name')})

            # query the ols to get measurementTechnique, infectiousAgent, healthCondition (infectiousDisease), and species
            if mts := hit.pop("measurementTechnique", None):
                if type(mts) is list:
                    hit["measurementTechnique"] = []
                    for mt in mts:
                        hit["measurementTechnique"].append(format_query_ols(mt))
                else:
                    hit["measurementTechnique"] = format_query_ols(mts)

            if ias := hit.pop("infectiousAgent", None):
                if type(ias) is list:
                    hit["infectiousAgent"] = []
                    for ia in ias:
                        hit["infectiousAgent"].append(format_query_ols(ia))
                else:
                    hit["infectiousAgent"] = format_query_ols(ias)

            if hcs := hit.pop("healthCondition", None):
                hit["healthCondition"] = []
                if type(hcs) is list:
                    for hc in hcs:
                        hit["healthCondition"].append(format_query_ols(hc))
                else:
                    hit["healthCondition"].append(format_query_ols(hcs))

            # infectious disease is deprecated change to healthCondition
            if ids := hit.pop("infectiousDisease", None):
                if not hit.get("healthCondition"):
                    hit["healthCondition"] = []
                if type(ids) is list:
                    for i_d in ids:
                        hit["healthCondition"].append(format_query_ols(i_d))
                else:
                    hit["healthCondition"].append(format_query_ols(ids))

            if species := hit.pop("species", None):
                if type(species) is list:
                    hit["species"] = []
                    for a_species in species:
                        hit["species"].append(format_query_ols(a_species))
                else:
                    hit["species"] = format_query_ols(species)

            # TODO update variableMeasured when we have update pubtator helper
            if vms := hit.pop("variableMeasured", None):
                if type(vms) is list:
                    hit["variableMeasured"] = []
                    for vm in vms:
                        hit["variableMeasured"].append(format_query_ky(vm))
                else:
                    hit["variableMeasured"] = format_query_ky(vms)

            if keywords := hit.pop("keywords", None):
                if type(keywords) is list:
                    hit["keywords"] = []
                    for keyword in keywords:
                        hit["keywords"].append(format_query_ky(keyword))
                else:
                    hit["keywords"] = format_query_ky(keywords)

            # fix temporalCoverage
            if temporalCoverage := hit.pop("temporalCoverage", None):
                if isinstance(temporalCoverage, list):
                    hit["temporalCoverage"] = [{"temporalInterval": tc} for tc in temporalCoverage]
                else:
                    hit["temporalCoverage"] = {"temporalInterval": temporalCoverage}

            # fix topicCategory
            if topicCategory := hit.pop("topicCategory", None):
                hit["topicCategory"] = {"url": topicCategory}

            # remove unnecessary values
            hit.pop("_meta", None)
            hit.pop("_score", None)
            hit.pop("@context", None)
            hit.pop("version", None)

            yield hit

    if count == total:
        logger.info("Total number of documents parsed: %s", total)
    else:
        logger.warning(
            "Did not parse all the records \n" "Total number parsed: %s \n Total number of documents: %s", count, total
        )
