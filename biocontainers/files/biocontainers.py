import datetime
import logging

import requests
from requests.adapters import HTTPAdapter, Retry

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def parse():
    logger.info("Retrieving Metadata")

    offset = 0
    limit = 1000
    total = 0

    tool_metadata = []
    s = requests.Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    s.mount("https://", HTTPAdapter(max_retries=retries))

    while True:
        url = f"https://api.biocontainers.pro/ga4gh/trs/v2/tools?limit={limit}&offset={offset}"
        try:
            all_data = s.get(url)
            all_data.raise_for_status()
        except requests.HTTPError as e:
            logger.error(f"HTTP Error: {e}")
            raise Exception(f"HTTP Error: {e}")
        if all_data.status_code == 204:
            # 204 means we've successfully reached the end of the metadata through pagination
            break
        for obj in all_data.json():
            total += 1
            tool_metadata.append(obj)

        logger.info("Tools Retrieved: %s", total)

        offset += 1000

    logger.info("Retrieved %s Tools", total)
    logger.info("Parsing %s Tools", total)

    count = 0
    for metadata in tool_metadata:
        count += 1
        if count % 100 == 0:
            logger.info("Parsed %s Tools", count)

        output = {
            "includedInDataCatalog": {
                "@type": "ComputationalTool",
                "name": "BioContainers",
                "url": "https://biocontainers.pro/",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@type": "ComputationalTool",
        }

        if description := metadata.get("description"):
            output["description"] = description

        if name := metadata.get("name"):
            output["name"] = name
            output["url"] = f"https://biocontainers.pro/tools/{name}"
            output["_id"] = "biocontainers_" + name

            versions_response = requests.get(f"https://api.biocontainers.pro//ga4gh/trs/v2/tools/{name}/versions")

            # similars_response = requests.get(
            #     f"https://api.biocontainers.pro//ga4gh/trs/v2/tools/{name}/similars")
            # if similars_response.status_code != 200:
            #     print(name, similars_response.status_code, 'similars')

        if identifiers := metadata.get("identifiers"):
            same_as = []
            for identifier in identifiers:
                if identifier.startswith("biotools:"):
                    same_as.append(f'https://bio.tools/{identifier.split(":")[1]}')
                if identifier.startswith("PMID:"):
                    output["pmids"] = identifier.split(":")[1]
            output["sameAs"] = same_as

        if license := metadata.get("license"):
            output["license"] = license

        if pulls := metadata.get("pulls"):
            output["interactionStatistic"] = {
                "userInteractionCount": pulls,
                "interactionType": "Number of downloads since release",
            }

        if tool_tags := metadata.get("tool_tags"):
            output["keywords"] = tool_tags

        if tool_url := metadata.get("tool_url"):
            if "github.com" in tool_url:
                output["codeRepository"] = tool_url
            else:
                output["mainEntityOfPage"] = tool_url

        if toolclass := metadata.get("toolclass"):
            output["applicationCategory"] = toolclass["description"]

        dates = []
        software_requirements = []
        if versions_response.status_code == 200:
            for version in versions_response.json():
                if images := version.get("images"):
                    for image in images:
                        if updated := image.get("updated"):
                            dates.append(datetime.datetime.strptime(updated, "%Y-%m-%dT%H:%M:%SZ").strftime("%Y-%m-%d"))
                        if software_requirement := version.get("image_type"):
                            if software_requirement not in software_requirements:
                                software_requirements.append(software_requirement)

        if len(dates):
            output["datePublished"] = sorted(dates)[0]
            output["dateModified"] = sorted(dates)[-1]

        if len(software_requirements):
            output["availableOnDevice"] = software_requirements

        # similars = []
        # if similars_response.status_code == 200:
        #     for similar in similars_response.json():
        #         similar_dict = {}
        #         if name := similar.get('name'):
        #             similar_dict['name'] = name
        #             similar_dict['_id'] = 'biocontainers_' + name
        #             similar_dict['identifier'] = name
        #             similar_dict['url'] = f'https://biocontainers.pro/tools/{name}'
        #         if bool(similar_dict):
        #             similars.append(similar_dict)
        # if similars:
        #     output['isSimilarTo'] = similars

        yield output
    logger.info("Finished Parsing %s Tools", count)
