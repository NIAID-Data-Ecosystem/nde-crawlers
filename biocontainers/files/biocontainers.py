import datetime
import logging
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')

logger.info('Retrieving Metadata')
# metadata = 'https://api.biocontainers.pro/ga4gh/trs/v2/tools?all_fields_search=&sort_order=asc&sort_field=default&offset=0&limit=1000'
# versions https://api.biocontainers.pro//ga4gh/trs/v2/tools/samtools/versions
# similars https://api.biocontainers.pro//ga4gh/trs/v2/tools/samtools/similars

# total 10556
offset = 0
limit = 1000
total = 0
# save metadata for every tool/workflow
tool_metadata = []
while True:
    data_count = 0
    url = f'https://api.biocontainers.pro/ga4gh/trs/v2/tools?limit={limit}&offset={offset}'
    all_data = requests.get(url)
    if all_data.status_code != 200:
        break
    for obj in all_data.json():
        total += 1
        tool_metadata.append(obj)
        data_count += 1

    logger.info('Tools Retrieved: %s', total)

    offset += 1000
    break
logger.info('Parsing %s Tool', total)
for metadata in tool_metadata:
    start = datetime.datetime.now()
    output = {'includedInDataCatalog': {
        '@type': 'ComputationalTool',
        'name': 'BioContainers',
        'url': 'https://biocontainers.pro/',
        'versionDate': datetime.date.today().isoformat()
    },
        '@type': 'ComputationalTool',
    }

    if description := metadata.get('description'):
        output['description'] = description

    if name := metadata.get('name'):
        output['name'] = name
        output['url'] = f'https://biocontainers.pro/tools/{name}'
        output['_id'] = 'Biocontainers_' + name
        versions_response = requests.get(
            f'https://api.biocontainers.pro//ga4gh/trs/v2/tools/{name}/versions')
        similars_response = requests.get(
            f'https://api.biocontainers.pro//ga4gh/trs/v2/tools/{name}/similars')

    if identifiers := metadata.get('identifiers'):
        for identifier in identifiers:
            if identifier.startswith('biotools:'):
                output['mainEntityOfPage'] = f'https://bio.tools/{identifier.split(":")[1]}'
            if identifier.startswith('pmid:'):
                output['pmids'] = identifier.split(":")[1]

    if pulls := metadata.get('pulls'):
        output['interactionStatistic'] = {
            'userInteractionCount': pulls,
            'interactionType': 'Number of downloads since release'
        }

    if tool_tags := metadata.get('tool_tags'):
        output['keywords'] = tool_tags

    if tool_url := metadata.get('tool_url'):
        if 'github.com' in tool_url:
            output['codeRepository'] = tool_url
        else:
            output['mainEntityOfPage'] = tool_url

    if toolclass := metadata.get('toolclass'):
        output['applicationCategory'] = toolclass['description']

    dates = []
    for version in versions_response.json():
        if images := version.get('images'):
            for image in images:
                if updated := image.get('updated'):
                    dates.append(datetime.datetime.strptime(
                        updated, '%Y-%m-%dT%H:%M:%SZ').strftime('%Y-%m-%d'))
    if dates:
        output['datePublished'] = sorted(dates)[0]
        output['dateModified'] = sorted(dates)[-1]

    similiars = []
    for similar in similars_response.json():
        similar_dict = {}
        if name := similar.get('name'):
            similar_dict['name'] = name
            similar_dict['_id'] = 'Biocontainers_' + name
            similar_dict['identifier'] = name
            similar_dict['url'] = f'https://biocontainers.pro/tools/{name}'
        if bool(similar_dict):
            similiars.append(similar_dict)
    output['isSimilarTo'] = similiars
    end = datetime.datetime.now()
    print(end - start)
