import datetime
import logging
import time

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def download_jsondocs():
    logger.info("Retrieving Metadata From API")

    biotoolsapiurl = "https://bio.tools/api/t"
    jsondoclist = []
    payloads = {"format": "json", "page": 1}

    while True:
        response = requests.get(biotoolsapiurl, params=payloads, timeout=20)
        response.raise_for_status()
        data = response.json()

        jsondoclist.extend(data["list"])
        logger.info("Retrieved %d records so far", len(jsondoclist))

        # Get the 'next' page value from the response and update 'payloads'
        biotoolsapipage = data.get("next")
        if not biotoolsapipage:
            break
        else:
            # Extract the next page number from the 'next' query parameter
            next_page = biotoolsapipage.split("=")[-1]
            payloads["page"] = next_page  # Update the page in the payloads

    logger.info("Finished retrieving all records. Total: %d", len(jsondoclist))
    return jsondoclist


def parse():
    count = 0
    logger.info("Parsing biotools metadata")

    jsondoclist = download_jsondocs()
    for tool in jsondoclist:
        count += 1
        logger.info(f"Parsing tool {count} of {len(jsondoclist)}")
        if count % 50 == 0:
            logger.info("Parsed %s tools", count)

        # Initialize the output structure
        output = {
            "@type": "ComputationalTool",
            "includedInDataCatalog": {
                "@type": "Dataset",
                "name": "Bio.tools",
                "url": "https://bio.tools/",
                "versionDate": datetime.date.today().isoformat(),
            },
        }

        # Direct mappings
        if tool_name := tool.get("name"):
            output["name"] = tool_name

        if description := tool.get("description"):
            output["description"] = description

        if homepage := tool.get("homepage"):
            output["mainEntityOfPage"] = homepage

        # Biotools ID needs to be split into identifier and URL
        if biotools_id := tool.get("biotoolsID"):
            output["identifier"] = biotools_id
            output["url"] = f"https://bio.tools/{biotools_id}"
            output["_id"] = f"biotools_{biotools_id}"
        else:
            logger.warning("Skipping tool without biotoolsID")
            continue

        # Conditional mappings for date fields with formatting
        if addition_date := tool.get("additionDate"):
            output["datePublished"] = addition_date.split("T")[0]
        if last_update := tool.get("lastUpdate"):
            output["dateModified"] = last_update.split("T")[0]

        # Topic and Function mappings
        if topics := tool.get("topic"):
            output["topicCategory"] = [
                {
                    "url": topic.get("uri"),
                    "identifier": topic.get("uri"),
                    "name": topic.get("term"),
                }
                for topic in topics
                if topic.get("uri") and topic.get("term")
            ]

        if functions := tool.get("function"):
            output["featureList"] = [
                {
                    "url": op.get("uri"),
                    "name": op.get("term"),
                }
                for func in functions
                for op in func.get("operation", [])
                if op.get("uri") and op.get("term")
            ]

        # Credit handling based on typeEntity and typeRole
        if credits := tool.get("credit"):
            funders = []
            authors = []
            contributors = []

            # Sets to track unique entries by name or URL to prevent duplicates
            unique_authors = set()
            unique_contributors = set()
            unique_funders = set()

            for credit in credits:
                if not credit.get("name"):
                    logger.warning(f"Skipping incomplete credit entry without a name: {credit}")
                    continue

                credit_entry = {}
                if credit.get("name"):
                    credit_entry["name"] = credit.get("name")
                if credit.get("email"):
                    credit_entry["email"] = credit.get("email")
                if credit.get("orcidid") or credit.get("url"):
                    credit_entry["url"] = credit.get("orcidid") or credit.get("url")

                identifier = credit.get("name") or credit.get("url")

                if credit.get("typeEntity") == "Person":
                    credit_entry["@type"] = "Person"
                elif credit.get("typeEntity") == "Organization":
                    credit_entry["@type"] = "Organization"

                type_roles = credit.get("typeRole") or ["Developer"]
                is_author = False

                if isinstance(type_roles, list):
                    for role in type_roles:
                        if role in ["Developer", "Primary contact"]:
                            # Prioritize authorship over contribution
                            if identifier not in unique_authors:
                                authors.append(credit_entry)
                                unique_authors.add(identifier)
                                is_author = True
                            break  # Stop further processing if added as an author
                    # Only add to contributors if not already added as an author
                    if not is_author:
                        for role in type_roles:
                            if role in [
                                "Maintainer",
                                "Provider",
                                "Contributor",
                                "Support",
                                "Documentor",
                            ]:
                                if identifier not in unique_contributors:
                                    contributors.append(credit_entry)
                                    unique_contributors.add(identifier)
                                break
                else:
                    if type_roles in ["Developer", "Primary contact"]:
                        if identifier not in unique_authors:
                            authors.append(credit_entry)
                            unique_authors.add(identifier)
                            is_author = True
                    elif not is_author and type_roles in [
                        "Maintainer",
                        "Provider",
                        "Contributor",
                        "Support",
                        "Documentor",
                    ]:
                        if identifier not in unique_contributors:
                            contributors.append(credit_entry)
                            unique_contributors.add(identifier)

                if credit.get("typeEntity") == "Funding agency":
                    if identifier not in unique_funders:
                        funders.append(credit_entry)
                        unique_funders.add(identifier)

            if authors:
                output["author"] = authors
            if contributors:
                output["contributor"] = contributors
            if funders:
                output["funding"] = {"funder": funders}

        # Publication handling based on type
        if publications := tool.get("publication"):
            for pub in publications:
                pub_entry = {}

                if doi := pub.get("doi"):
                    pub_entry["doi"] = doi

                if metadata := pub.get("metadata"):
                    if title := metadata.get("title"):
                        pub_entry["name"] = title

                    if date := metadata.get("date"):
                        pub_entry["datePublished"] = date.split("T")[0]

                    if journal := metadata.get("journal"):
                        pub_entry["journal"] = journal

                    authors = []
                    if metadata.get("authors"):
                        for author in metadata.get("authors", []):
                            if author_name := author.get("name"):
                                authors.append({"name": author_name})
                        if authors:
                            pub_entry["author"] = authors

                    if abstract := metadata.get("abstract"):
                        pub_entry["abstract"] = abstract

                if pmid := pub.get("pmid"):
                    pub_entry["pmid"] = pmid

                if pmcid := pub.get("pmcid"):
                    pub_entry["pmcid"] = pmcid

                if version := pub.get("version"):
                    pub_entry["version"] = version

                pub_type = pub.get("type")
                if isinstance(pub_type, list) and len(pub_type) == 1:
                    pub_type = pub_type[0]
                if pub_type == "Primary":
                    output.setdefault("citation", []).append(pub_entry)
                elif pub_type == "Method":
                    output.setdefault("isBasedOn", []).append(pub_entry)
                elif pub_type == "Usage":
                    output.setdefault("citedBy", []).append(pub_entry)
                elif pub_type == "Benchmarking study":
                    output.setdefault("isBasisFor", []).append(pub_entry)
                elif pub_type in ["Review", "Other"]:
                    output.setdefault("citedBy", []).append(pub_entry)
                else:
                    logger.warning("Unknown publication type: %s", pub_type)
                    if pub_entry:
                        logger.warning("Adding publication entry to 'citedBy' by default")
                        logger.warning(pub_entry)
                        output.setdefault("citedBy", []).append(pub_entry)

        # Relations
        if relations := tool.get("relation"):
            for relation in relations:
                relation_entry = relation.get("biotoolsID")
                relation_type = relation.get("type")
                if not relation_entry:
                    continue

                if relation_type in ["isNewVersionOf", "hasNewVersion"]:
                    output.setdefault("sameAs", []).append(relation_entry)
                elif relation_type == "uses":
                    output.setdefault("isBasedOn", []).append({"identifier": relation_entry})
                elif relation_type == "usedBy":
                    output.setdefault("isBasisFor", []).append({"identifier": relation_entry})
                elif relation_type == "includes":
                    output.setdefault("hasPart", []).append({"identifier": relation_entry})
                elif relation_type == "includedIn":
                    output.setdefault("isPartOf", []).append({"identifier": relation_entry})
                else:
                    output.setdefault("isRelatedTo", []).append({"identifier": relation_entry})
                    logger.warning("Unknown relation type: %s", relation_type)

        # Software-specific fields
        if version := tool.get("version"):
            output["softwareVersion"] = version

        if license := tool.get("license"):
            output["license"] = license

        if programming_language := tool.get("language"):
            output["programmingLanguage"] = programming_language

        if os := tool.get("operatingSystem"):
            output["operatingSystem"] = os

        if links := tool.get("link"):
            for link in links:
                link_type = link.get("type")
                link_url = link.get("url")
                if not link_url:
                    continue
                if link_type == "Software catalogue":
                    output.setdefault("sdPublisher", []).append({"url": link_url})
                elif link_type == "Repository":
                    output.setdefault("codeRepository", []).append({"url": link_url})
                else:
                    output.setdefault("isRelatedTo", []).append({"url": link_url})
                    logger.warning("Unknown link type: %s", link_type)

        # Download URLs
        if downloads := tool.get("download"):
            for download in downloads:
                download_url = download.get("url")
                download_type = download.get("type")

                if not download_url:
                    continue  # Skip if either URL or type is missing

                if download_type == "API specification":
                    output.setdefault("softwareHelp", []).append({"url": download_url})
                elif download_type == "Biological data":
                    output.setdefault("softwareAddOn", []).append({"url": download_url})
                elif download_type == "Binaries":
                    output.setdefault("downloadUrl", []).append({"name": download_url})
                elif download_type == "Command-line specificaiton":
                    output.setdefault("softwareHelp", []).append({"url": download_url})
                elif download_type == "Container file":
                    output.setdefault("downloadUrl", []).append({"name": download_url})
                elif download_type == "Icon":
                    output.setdefault("thumbnailUrl", []).append(download_url)
                elif download_type == "Screenshot":
                    output.setdefault("thumbnailUrl", []).append(download_url)
                elif download_type == "Source code":
                    output.setdefault("codeRepository", []).append(download_url)
                elif download_type == "Software package":
                    output.setdefault("downloadUrl", []).append({"name": download_url})
                elif download_type == "Test data":
                    output.setdefault("softwareAddOn", []).append({"url": download_url})
                elif download_type == "Test script":
                    output.setdefault("softwareAddOn", []).append({"url": download_url})
                elif download_type == "Tool wrapper (CWL)":
                    output["programmingLanguage"] = "CWL"
                elif download_type == "Tool wrapper (galaxy)":
                    output["isAvailableOnDevice"] = "Galaxy"
                elif download_type == "Tool wrapper (taverna)":
                    output["applicationSuite"] = "Taverna"
                elif download_type == "VM image":
                    output.setdefault("downloadUrl", []).append({"name": download_url})
                elif download_type == "Downloads page":
                    output.setdefault("downloadUrl", []).append({"name": download_url})
                elif download_type == "Other":
                    output.setdefault("downloadUrl", []).append({"name": download_url})
                else:
                    output.setdefault("downloadUrl", []).append({"name": download_url})
                    logger.warning("Unknown download type: %s", download_type)

        # Cost and accessibility
        if cost := tool.get("cost"):
            if cost.lower() == "free of charge":
                output["isAccessibleForFree"] = True
            elif cost.lower() == "commercial":
                output["isAccessibleForFree"] = False
            elif cost.lower() == "free of charge (with restrictions)":
                output["isAccessibleForFree"] = True
            else:
                logger.warning("Unknown cost value: %s", cost)

        if accessibility := tool.get("accessibility"):
            if accessibility.lower() == "open access":
                output["conditionsOfAccess"] = "Open"
            elif accessibility.lower() == "restricted access":
                output["conditionsOfAccess"] = "Closed"
            elif accessibility.lower() == "open access (with restrictions)":
                output["conditionsOfAccess"] = "Restricted"
            else:
                logger.warning("Unknown accessibility value: %s", accessibility)

        if documentation := tool.get("documentation"):
            software_help = []
            for doc in documentation:
                software_help_entry = {}
                if url := doc.get("url"):
                    software_help_entry["url"] = url
                if doc_type := doc.get("type"):
                    software_help_entry["name"] = doc_type
                if software_help_entry:
                    software_help.append(software_help_entry)
            if software_help:
                if "softwareHelp" in output:
                    output["softwareHelp"].extend(software_help)
                else:
                    output["softwareHelp"] = software_help

        if other_ids := tool.get("otherID"):
            for oid in other_ids:
                if oid.get("type") == "doi":
                    output["doi"] = oid.get("value")

        if tool_type := tool.get("toolType"):
            output["applicationCategory"] = tool_type

        if collection_id := tool.get("collectionID"):
            output["keywords"] = collection_id

        if maturity := tool.get("maturity"):
            output["creativeWorkStatus"] = maturity

        # Initialize lists to gather all inputs and outputs
        all_inputs = []
        all_outputs = []

        if functions := tool.get("function"):
            for func in functions:
                # Process inputs for each function
                if inputs := func.get("input"):
                    for i in inputs:
                        input_entry = {}

                        if url := i.get("data", {}).get("uri"):
                            input_entry["url"] = url

                        if name := i.get("data", {}).get("term"):
                            input_entry["name"] = name

                        if formats := i.get("format"):
                            encoding_format = {}
                            for f in formats:
                                if encoding_format_url := f.get("uri"):
                                    encoding_format["url"] = encoding_format_url
                                if encoding_format_name := f.get("term"):
                                    encoding_format["name"] = encoding_format_name
                            if encoding_format:
                                input_entry["encodingFormat"] = encoding_format

                        if input_entry:
                            all_inputs.append(input_entry)

                # Process outputs for each function
                if outputs := func.get("output"):
                    for o in outputs:
                        output_entry = {}

                        if url := o.get("data", {}).get("uri"):
                            output_entry["url"] = url

                        if name := o.get("data", {}).get("term"):
                            output_entry["name"] = name

                        if formats := o.get("format"):
                            encoding_format = {}
                            for f in formats:
                                if encoding_format_url := f.get("uri"):
                                    encoding_format["url"] = encoding_format_url
                                if encoding_format_name := f.get("term"):
                                    encoding_format["name"] = encoding_format_name
                            if encoding_format:
                                output_entry["encodingFormat"] = encoding_format

                        if output_entry:
                            all_outputs.append(output_entry)

        if all_inputs:
            output["input"] = all_inputs

        if all_outputs:
            output["output"] = all_outputs

        yield output

    logger.info("Finished parsing all tools.")
