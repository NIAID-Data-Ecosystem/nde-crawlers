import datetime
import logging
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def get_categories():
    categories_url = 'https://dockstore.org/api/categories'
    categories = []
    response = requests.get(categories_url)
    for category in response.json():
        categories.append(category['name'])
    return categories


def filter_topics():
    url = 'https://dockstore.org/api/api/ga4gh/v2/extended/tools/entry/_search'
    topics = {}
    for category in get_categories():
        topics[category] = []
        data = {"size": 201, "_source": ["full_workflow_path"], "query": {"bool": {"filter": {"term": {
            "categories.name.keyword": f"{category}"}}, "must": [{"match": {"_index": "workflows"}}, {"match_all": {}}]}}}
        topic_response = requests.post(url, json=data)
        [topics[category].append(x['_source']['full_workflow_path'])
            for x in topic_response.json()['hits']['hits']]
    return topics


def parse():
    logger.info('Retrieving Metadata')

    # pagination using offset and limit, api using offset as page
    offset = 0
    limit = 100

    # save metadata for every tool/workflow
    tool_metadata = []
    while True:
        data_count = 0
        url = f"https://dockstore.org/api/api/ga4gh/v2/tools?limit={limit}&offset={offset}"
        all_data = requests.get(url)
        for obj in all_data.json():
            tool_metadata.append(obj)
            data_count += 1

        logger.info('Current Page: %s', offset)

        # condition to check if we've reached the last page
        if data_count != limit:
            break

        offset += 1

    logger.info('Retrieving Topic Categories')
    topic_dict = filter_topics()

    logger.info('Parsing Metadata')
    count = 0
    for metadata in tool_metadata:
        count += 1
        if count % 100 == 0:
            logger.info(f'Parsed {count} package metadata')

        output = {'includedInDataCatalog': {
            '@type': 'ComputationalTool',
            'name': 'Dockstore',
            'url': 'https://dockstore.org/',
            'versionDate': datetime.date.today().isoformat()
        },
            '@type': 'ComputationalTool',
        }

        if identifier := metadata.get('toolname'):
            output['identifier'] = identifier
            output['_id'] = 'dockstore_' + identifier.replace('/', '_')

        if doi := metadata.get('aliases'):
            output['doi'] = doi[0]

        author_list = []
        if author := metadata.get('author'):
            # if affiliation := metadata.get('organization'):
            #     authors = author.split(', ')
            #     for author in authors:
            #         if author != 'Unknown author' and author != 'docker':
            #             author_list.append({
            #                 'name': author, 'affiliation': {'name': affiliation}})
            # else:
            authors = author.split(', ')
            for author in authors:
                author_list.append({'name': author})
        if len(author_list):
            output['author'] = author_list

        if description := metadata.get('description'):
            output['description'] = description

        if url := metadata.get('id'):
            if url.startswith('#workflow/'):
                output['url'] = 'https://dockstore.org/workflows/' + url[10:]
                for key in topic_dict:
                    if url[10:] in topic_dict[key]:
                        output['topicCategory'] = {
                            'name': key,
                            'url': 'https://dockstore.org/search?categories.name.keyword=' + key,
                            'curatedBy': {
                                'name': 'Dockstore',
                                'url': 'https://dockstore.org/search?entryType=workflows&searchMode=files'
                            }
                        }
                if 'github.com' in url:
                    output['codeRepository'] = 'https://' + \
                        url[10:].split(':')[0]
            elif url.startswith('#service/'):
                output['url'] = 'https://dockstore.org/services/' + url[9:]
                if 'github.com' in url:
                    output['codeRepository'] = 'https://' + \
                        url[9:].split(':')[0]
            else:
                output['url'] = 'https://dockstore.org/containers/' + url
                if 'github' in url:
                    output['codeRepository'] = 'https://' + url.split(':')[0]

        if date_modified := metadata.get('meta_version'):
            if metadata.get('versions') == []:
                output['datePublished'] = datetime.datetime.strptime(
                    date_modified, '%Y-%m-%d %H:%M:%S.%f').strftime('%Y-%m-%d')
            else:
                output['dateModified'] = datetime.datetime.strptime(
                    date_modified, '%Y-%m-%d %H:%M:%S.%f').strftime('%Y-%m-%d')

        if toolclass := metadata.get('toolclass'):
            output['applicationCategory'] = toolclass['name']

        if name := metadata.get('toolname'):
            output['name'] = name

        if versions := metadata.get('versions'):
            if versions != []:
                oldest = versions[-1]
                if oldest['meta_version'] != 'Thu Jan 01 00:00:00 UTC 1970':
                    output['datePublished'] = datetime.datetime.strptime(
                        oldest['meta_version'], '%Y-%m-%d %H:%M:%S.%f').strftime('%Y-%m-%d')
                else:
                    output['datePublished'] = output['dateModified']
                    output.pop('dateModified')
                latest = versions[0]
                output['softwareVersion'] = latest['name']

                languages = []
                for version in versions:
                    [languages.append(
                        x) for x in version['descriptor_type'] if x not in languages]
                if languages != []:
                    for language in languages:
                        if language != 'CWL' and language != 'WDL':
                            languages.remove(language)
                    if languages != []:
                        output['programmingLanguage'] = languages
        yield output
