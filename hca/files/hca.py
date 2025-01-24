import datetime
import logging

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def retrieve_ids():
    """Get all project IDs from the HCA and convert to UUID format"""
    url = "https://service.azul.data.humancellatlas.org/index/catalogs"

    sources = requests.get(url).json()["catalogs"]
    for key in sources:
        sources = sources[key]["plugins"]["repository"]["sources"]
        break

    logger.info("Retrieving %s project IDs", len(sources))

    uuids = []
    for source in sources:
        id = source.split("_")[2]
        uuid = id[:8] + "-" + id[8:12] + "-" + id[12:16] + "-" + id[16:20] + "-" + id[20:]
        uuids.append(uuid)
    return uuids


def retrieve_project_metadata():
    """Get project metadata from UUIDS"""
    uuids = retrieve_ids()
    logger.info("Retrieving individual project metadata")
    count = 0
    meta = []
    for uuid in uuids:
        count += 1
        logger.info(f"Retrieving project {count} of {len(uuids)}")
        url = f"https://service.azul.data.humancellatlas.org/index/projects/{uuid}"
        r = requests.get(url)
        meta.append(r.json())
    return meta


def get_same_as(accession):
    """Get the url for a given accession"""
    if accession.startswith(("SRP", "ERP", "DRP", "SRR", "SRX")):
        return f"https://www.ebi.ac.uk/ena/data/view/{accession}"
    elif accession.startswith(("GSE", "GSM", "GDS")):
        return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    elif accession.startswith(("E-")):
        return f"https://www.ebi.ac.uk/arrayexpress/experiments/{accession}"
    elif accession.startswith(("PRJ")):
        return f"https://www.ncbi.nlm.nih.gov/bioproject/{accession}"
    elif accession.startswith(("EGA")):
        return f"https://ega-archive.org/studies/{accession}"
    elif accession.startswith(("PRJEB")):
        return f"https://www.ebi.ac.uk/ena/browser/view/{accession}"
    elif accession.startswith(("S-")):
        return f"https://www.ebi.ac.uk/biostudies/dsp/studies/{accession}"
    elif accession.startswith(("phs")):
        return f"https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={accession}"
    else:
        logger.info(f"No sameAs found for {accession}")
        return None


def parse():
    count = 0
    project_meta = retrieve_project_metadata()
    for metadata in project_meta:
        count += 1
        logger.info(f"Parsing project {count} of {len(project_meta)}")
        output = {
            "includedInDataCatalog": {
                "@type": "Dataset",
                "name": "Human Cell Atlas",
                "url": "https://data.humancellatlas.org/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "Dataset",
        }
        if entry_id := metadata.get("entryId"):
            url = f"https://data.humancellatlas.org/explore/projects/{entry_id}"
            output["url"] = url
            output["includedInDataCatalog"]["dataset"] = url
            output["_id"] = f"HCA_{entry_id}"
        else:
            logger.error("No entryId found for project. Skipping...")
            continue

        if projects := metadata.get("projects"):
            for project in projects:
                if project_title := project.get("projectTitle"):
                    output["name"] = project_title

                if project_description := project.get("projectDescription"):
                    output["description"] = project_description

                if contributors := project.get("contributors"):
                    authors = []
                    for contributor in contributors:
                        author_dict = {}
                        if institution := contributor.get("institution"):
                            author_dict["affiliation"] = {"name": institution}
                        if name := contributor.get("contactName"):
                            given_name = name.split(",")[0]
                            family_name = name.split(",")[-1]
                            author_dict["givenName"] = given_name
                            author_dict["familyName"] = family_name
                            author_dict["name"] = given_name + " " + family_name
                        if role := contributor.get("projectRole"):
                            author_dict["role"] = role
                        if laboratory := contributor.get("laboratory"):
                            if "affiliation" in author_dict:
                                author_dict["affiliation"]["name"] += ", " + laboratory
                            else:
                                author_dict["affiliation"] = {"name": laboratory}
                        if email := contributor.get("email"):
                            author_dict["email"] = email
                        if author_dict:
                            authors.append(author_dict)
                    if len(authors):
                        output["author"] = authors

                if publications := project.get("publications"):
                    citations = []
                    for publication in publications:
                        citation_dict = {}
                        if title := publication.get("publicationTitle"):
                            citation_dict["name"] = title
                        if url := publication.get("publicationUrl"):
                            citation_dict["url"] = url
                        if doi := publication.get("doi"):
                            citation_dict["doi"] = doi
                        if citation_dict:
                            citations.append(citation_dict)
                    if len(citations):
                        output["citation"] = citations

                if links := project.get("supplementaryLinks"):
                    other_links = []
                    github_links = []
                    for link in links:
                        if link is not None:
                            if "github" in link:
                                github_links.append(link)
                            else:
                                other_links.append(link)
                    if len(other_links):
                        output["mainEntityOfPage"] = other_links
                    if len(github_links):
                        is_based_on = []
                        for link in github_links:
                            is_based_on.append({"@type": "ComputationalTool", "codeRepository": link})
                        output["isBasedOn"] = is_based_on

                if accessions := project.get("accessions"):
                    for accession_dict in accessions:
                        if accession := accession_dict.get("accession"):
                            same_as = get_same_as(accession)
                            if same_as is not None:
                                output["sameAs"] = same_as

                if accessible := project.get("accessible"):
                    output["isAccessibleForFree"] = accessible

            if samples := metadata.get("samples"):
                keywords = []
                for sample in samples:
                    if eff_organs := sample.get("effectiveOrgan"):
                        for eff_organ in eff_organs:
                            if eff_organ is not None:
                                keywords.append(eff_organ)
                    if organs := sample.get("organ"):
                        for organ in organs:
                            if organ is not None:
                                keywords.append(organ)
                    if organ_parts := sample.get("organPart"):
                        for organ_part in organ_parts:
                            if organ_part is not None:
                                keywords.append(organ_part)
                    if diseases := sample.get("disease"):
                        health_conditions = []
                        for disease in diseases:
                            if disease == "normal":
                                health_conditions.append({"name": "healthy"})
                            elif disease is not None:
                                health_conditions.append({"name": disease})
                        if len(health_conditions):
                            output["healthCondition"] = health_conditions
                if len(keywords):
                    output["keywords"] = list(set(keywords))

            if donor_organisms := metadata.get("donorOrganisms"):
                for donor_organism in donor_organisms:
                    if species := donor_organism.get("genusSpecies"):
                        species_list = []
                        for single_species in species:
                            if species is not None:
                                species_list.append({"name": single_species})
                        if len(species_list):
                            output["species"] = species_list

            if dates := metadata.get("dates"):
                for date in dates:
                    if last_modified := date.get("lastModifiedDate"):
                        output["dateModified"] = datetime.datetime.strptime(
                            last_modified, "%Y-%m-%dT%H:%M:%S.%f%z"
                        ).strftime("%Y-%m-%d")
                    if submission_date := date.get("submissionDate"):
                        output["datePublished"] = datetime.datetime.strptime(
                            submission_date, "%Y-%m-%dT%H:%M:%S.%f%z"
                        ).strftime("%Y-%m-%d")
        yield output
