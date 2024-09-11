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

        # Respect API rate limiting
        time.sleep(1)

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

            for credit in credits:
                credit_entry = {
                    "name": credit.get("name"),
                    "email": credit.get("email"),
                    "url": credit.get("orcidid") or credit.get("url"),
                }

                if credit.get("typeEntity") == "Person":
                    credit_entry["@type"] = "Person"
                else:
                    credit_entry["@type"] = "Organization"

                if credit.get("typeEntity") == "Funding agency":
                    funders.append(credit_entry)
                else:
                    if credit.get("typeRole") in ["Developer", "Primary contact"]:
                        authors.append(credit_entry)
                    elif credit.get("typeRole") in ["Maintainer", "Provider", "Contributor", "Support"]:
                        contributors.append(credit_entry)

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

                # Handle DOI
                if doi := pub.get("doi"):
                    pub_entry["doi"] = doi

                # Safely handle metadata and its fields
                if metadata := pub.get("metadata"):
                    # Handle title (name)
                    if title := metadata.get("title"):
                        pub_entry["name"] = title

                    # Handle datePublished, ensuring it's not None
                    if date := metadata.get("date"):
                        pub_entry["datePublished"] = date.split("T")[0]  # Only store YYYY-MM-DD format

                    # Handle journal
                    if journal := metadata.get("journal"):
                        pub_entry["journal"] = journal

                    # Handle authors
                    authors = []
                    if metadata.get("authors"):
                        for author in metadata.get("authors", []):
                            if author_name := author.get("name"):
                                authors.append({"name": author_name})
                        if authors:
                            pub_entry["author"] = authors

                    # Handle abstract
                    if abstract := metadata.get("abstract"):
                        pub_entry["abstract"] = abstract

                # Handle PMID
                if pmid := pub.get("pmid"):
                    pub_entry["pmid"] = pmid

                # Handle PMCID
                if pmcid := pub.get("pmcid"):
                    pub_entry["pmcid"] = pmcid

                # Handle version
                if version := pub.get("version"):
                    pub_entry["version"] = version

                # Conditional logic based on publication type
                pub_type = pub.get("type")
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

        # Relations
        if relations := tool.get("relation"):
            for relation in relations:
                relation_entry = {"identifier": relation.get("biotoolsID")}
                relation_type = relation.get("type")

                if relation_type == "isNewVersion":
                    output.setdefault("isRelatedTo", []).append(relation_entry)
                elif relation_type == "uses":
                    output.setdefault("isBasedOn", []).append(relation_entry)
                elif relation_type == "usedBy":
                    output.setdefault("isBasisFor", []).append(relation_entry)
                elif relation_type == "includes":
                    output.setdefault("hasPart", []).append(relation_entry)
                elif relation_type == "includedIn":
                    output.setdefault("isPartOf", []).append(relation_entry)
                else:
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

        # Download URLs
        if downloads := tool.get("download"):
            for download in downloads:
                download_entry = {"url": download.get("url")}
                if download.get("version"):  # Only add 'version' if it's not None
                    download_entry["version"] = download.get("version")
                output.setdefault("downloadURL", []).append(download_entry)

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
            for doc in documentation:
                software_help = []
                if url := doc.get("url"):
                    software_help.append({"url": url})
                if doc_type := doc.get("type"):
                    software_help.append({"name": doc_type})
                if software_help:
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

        # Iterate over each function in the list
        if functions := tool.get("function"):
            for func in functions:
                # Process inputs for each function
                if inputs := func.get("input"):
                    for i in inputs:
                        input_entry = {}

                        # Add URL if available
                        if url := i.get("data", {}).get("uri"):
                            input_entry["url"] = url

                        # Add name if available
                        if name := i.get("data", {}).get("term"):
                            input_entry["name"] = name

                        # Add encoding format if available
                        if formats := i.get("format"):
                            encoding_format = {}
                            for f in formats:
                                if encoding_format_url := f.get("uri"):
                                    encoding_format["url"] = encoding_format_url
                                if encoding_format_name := f.get("term"):
                                    encoding_format["name"] = encoding_format_name
                            if encoding_format:
                                input_entry["encodingFormat"] = encoding_format

                        # Append the input entry if any data was collected
                        if input_entry:
                            all_inputs.append(input_entry)

                # Process outputs for each function
                if outputs := func.get("output"):
                    for o in outputs:
                        output_entry = {}

                        # Add URL if available
                        if url := o.get("data", {}).get("uri"):
                            output_entry["url"] = url

                        # Add name if available
                        if name := o.get("data", {}).get("term"):
                            output_entry["name"] = name

                        # Add encoding format if available
                        if formats := o.get("format"):
                            encoding_format = {}
                            for f in formats:
                                if encoding_format_url := f.get("uri"):
                                    encoding_format["url"] = encoding_format_url
                                if encoding_format_name := f.get("term"):
                                    encoding_format["name"] = encoding_format_name
                            if encoding_format:
                                output_entry["encodingFormat"] = encoding_format

                        # Append the output entry if any data was collected
                        if output_entry:
                            all_outputs.append(output_entry)

        # Only add 'input' and 'output' to the final output if we processed any inputs or outputs
        if all_inputs:
            output["input"] = all_inputs

        if all_outputs:
            output["output"] = all_outputs

        yield output

    logger.info("Finished parsing all tools.")
