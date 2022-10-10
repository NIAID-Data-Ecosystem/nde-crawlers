import datetime
import logging
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def get_ids():
    """Get the ids of the workflows from the workflowhub."""
    logger.info('Getting workflow ids')
    count = 0
    url = 'https://workflowhub.eu/ga4gh/trs/v2/tools'
    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.HTTPError as err:
        logger.error(err)
    result = []
    for workflow in response.json():
        count += 1
        if count % 50 == 0:
            logger.info('Retrieved %s ids.', count)
        result.append(workflow['id'])
    return result


def get_metadata(workflow_id):
    """Get the metadata of a workflow."""
    url = f'https://workflowhub.eu/workflows/{workflow_id}.json'
    response = requests.get(url)
    return response.json()['data']


def parse_metadata(metadata):
    """Parse the metadata of a workflow."""
    output = {'includedInDataCatalog': {
        '@type': 'ComputationalTool',
        'name': 'WorkflowHub',
        'url': 'https://workflowhub.eu/',
        'versionDate': datetime.date.today().isoformat()
    },
        '@type': 'ComputationalTool',
    }

    if type := metadata.get('type'):
        output['applicationCategory'] = type

    if workflow_id := metadata.get('id'):
        output['_id'] = 'workflowhub_workflow' + workflow_id

    if attributes := metadata.get('attributes'):
        if discussion_links := attributes.get('discussion_links'):
            for discussion_link in discussion_links:
                if url := discussion_link.get('url'):
                    output['discussionUrl'] = url

        if title := attributes.get('title'):
            output['name'] = title

        if license := attributes.get('license'):
            output['license'] = license

        if description := attributes.get('description'):
            output['description'] = description

        if latest_version := attributes.get('latest_version'):
            output['softwareVersion'] = str(latest_version)

        if tags := attributes.get('tags'):
            output['keywords'] = tags
            if 'covid-19' in tags:
                output['topicCategory'] = {
                    'name': 'COVID-19',
                    'url': 'https://workflowhub.eu/workflows?filter%5Btag%5D=covid-19',
                    'curatedBy': {
                        'name': 'WorkflowHub',
                        'url': 'https://workflowhub.eu/'
                    }
                }

        if versions := attributes.get('versions'):
            if isinstance(versions, list):
                for version in versions:
                    if version.get('version') == latest_version:
                        if remote := version.get('remote'):
                            output['codeRepository'] = remote

        if created := attributes.get('created_at'):
            output['dateCreated'] = datetime.datetime.strptime(
                created, '%Y-%m-%dT%H:%M:%S.%fZ').strftime('%Y-%m-%d')

        if updated := attributes.get('updated_at'):
            output['dateModified'] = datetime.datetime.strptime(
                updated, '%Y-%m-%dT%H:%M:%S.%fZ').strftime('%Y-%m-%d')

        if doi := attributes.get('doi'):
            output['doi'] = doi

        if creators := attributes.get('creators'):
            author_list = []
            if isinstance(creators, list):
                for creator in creators:
                    author_dict = {}
                    if profile := creator.get('profile'):
                        author_dict['url'] = f'https://workflowhub.eu/{profile}'
                    if family_name := creator.get('family_name'):
                        author_dict['familyName'] = family_name
                    if given_name := creator.get('given_name'):
                        author_dict['givenName'] = given_name
                    if family_name and given_name:
                        author_dict['name'] = f'{given_name} {family_name}'
                    if affiliation := creator.get('affiliation'):
                        author_dict['affiliation'] = {'name': affiliation}
                    if orcid := creator.get('orcid'):
                        author_dict['identifier'] = orcid
                    author_list.append(author_dict)

            # if other_creators := attributes.get('other_creators'):
            #     other_creators = other_creators.split(', ')
            #     other_creators = [creator.split(' and ')
            #                       for creator in other_creators]
            #     for creator in other_creators:
            #         author_list.append({'name': creator})

            if len(author_list):
                output['author'] = author_list

        if workflow_class := attributes.get('workflow_class'):
            if title := workflow_class.get('title'):
                output['programmingLanguage'] = title

        if edam_operations := attributes.get('edam_operations'):
            feature_list = []
            if isinstance(edam_operations, list):
                for operation in edam_operations:
                    if identifier := operation.get('identifier'):
                        feature_list.append(identifier)
                if len(feature_list):
                    output['featureList'] = feature_list

        if edam_topics := attributes.get('edam_topics'):
            topic_list = []
            if isinstance(edam_topics, list):
                for topic in edam_topics:
                    topic_dict = {}
                    if identifier := topic.get('identifier'):
                        topic_dict['url'] = identifier
                        topic_dict['identifier'] = identifier.split('/')[-1]
                    if label := topic.get('label'):
                        topic_dict['name'] = label
                    topic_list.append(topic_dict)
                if len(topic_list):
                    output['applicationSubCategory'] = topic_list

        if internals := attributes.get('internals'):
            if inputs := internals.get('inputs'):
                if isinstance(inputs, list):
                    input_list = []
                    for input in inputs:
                        input_dict = {}
                        if id := input.get('id'):
                            input_dict['identifier'] = id
                        if name := input.get('name'):
                            input_dict['name'] = name
                        if description := input.get('description'):
                            input_dict['description'] = description
                        input_list.append(input_dict)
                    if len(input_list):
                        output['input'] = input_list

            if outputs := internals.get('outputs'):
                output_list = []
                if isinstance(outputs, list):
                    for single_output in outputs:
                        output_dict = {}
                        if id := single_output.get('id'):
                            output_dict['identifier'] = id
                        if name := single_output.get('name'):
                            output_dict['name'] = name
                        if description := single_output.get('description'):
                            output_dict['description'] = description
                        if type := single_output.get('type'):
                            if isinstance(type, list):
                                output_dict['type'] = type[0]['type']
                            elif type != None and type != '':
                                output_dict['type'] = type
                        output_list.append(output_dict)
                    if len(output_list):
                        output['output'] = output_list

    if links := metadata.get('links'):
        if self_link := links.get('self'):
            output['url'] = f'https://workflowhub.eu{self_link}'

    return output


def parse():
    """Parse the workflows from the workflowhub."""
    metadata_list = []
    count = 0
    workflow_ids = get_ids()
    for workflow_id in workflow_ids:
        count += 1
        if count % 50 == 0:
            logger.info('Retrieved %s workflows.', count)
        metadata = get_metadata(workflow_id)
        metadata_list.append(metadata)

    logger.info('Parsing Metadata')
    count = 0
    for metadata in metadata_list:
        count += 1
        if count % 50 == 0:
            logger.info('Parsed %s workflows.', count)
        yield parse_metadata(metadata)
    logger.info('Finished parsing metadata. Total: %s workflows.', count)
