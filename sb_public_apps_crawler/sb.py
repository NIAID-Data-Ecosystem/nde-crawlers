import json
from pprint import pprint
import requests


def parse():

    url = "https://igor.sbgenomics.com/ns/amigo/api/v1/apps/explore?limit=20&offset=0&order_by=label"
    all_apps = requests.get(url)
    json_obj = json.loads(all_apps.text)
    obj_list = json_obj['data']
    public_ids = []
    for obj in obj_list:
        public_ids.append(obj['public_id'])

    all_app_meta_data = []
    for id in public_ids:
        app_meta_data = requests.get(
            "https://igor.sbgenomics.com/ns/brood/v1/raw/" + id)
        app_json = json.loads(app_meta_data.text)
        all_app_meta_data.append(app_json)

    # if no description use doc

    for data in all_app_meta_data:
        identifier = data.get('sbg:id')
        output = {
            "_id": "SB_Public_Apps_" + identifier,
            "includedInDataCatalog": {"name": "PublicApps@SevenBridges"},
            "@type": "Dataset",
            "url": "https://igor.sbgenomics.com/public/apps/" + identifier
        }
        if app_class := data.get('class'):
            output['applicationCategory'] = app_class
        if cwl_version := data.get('cwlVersion'):
            output['version'] = cwl_version
        if label := data.get('label'):
            output['name'] = label
        if description := data.get('description'):
            output['description'] = description
        if doc := data.get('doc'):
            output['description'] = doc

        # TODO value for input should be list of FormalParameter objects https://bioschemas.org/types/FormalParameter/1.0-RELEASE
        if input_object_list := data.get('inputs'):
            result_input_list = []
            for input_object in input_object_list:
                result_input_object = {}
                if input_id := input_object.get('id'):
                    result_input_object['identifier'] = input_id
                if name := input_object.get('label'):
                    result_input_object['name'] = name
                if description := input_object.get('description'):
                    result_input_object['description'] = description
                if doc := input_object.get('doc'):
                    result_input_object['description'] = doc
                required = input_object.get('required')
                if required is not None:
                    result_input_object['valueRequired'] = required
                # if file_type := input_object.get('sbg:fileTypes'):
                #     result_input_object['encodingFormat'] = file_type
                result_input_list.append(result_input_object)
            output['input'] = result_input_list
            output['original_input'] = input_object_list

            # TODO value for output should be list of FormalParameter objects https://bioschemas.org/types/FormalParameter/1.0-RELEASE
        if outputs := data.get('output'):
            output['bioschemas:outputs'] = outputs

        # commenting out for readablility
        # if steps := data.get('steps'):
        #     output['bioschemas:hasPart'] = steps

        if requirements := data.get('requirements'):
            result_list = []
            for obj in requirements:
                result_list.append(obj['class'])
            output['schema:softwareRequirements'] = result_list

        if image_url := data.get('sbg:image_url'):
            output['thumbnailUrl'] = image_url

        # TODO check applicationSuite exists, append version?
        if toolkit := data.get('sbg:toolkit'):
            output['applicationSuite'] = toolkit + \
                'version:' + data.get('sbg:toolkitVersion')

        if license := data.get('license'):
            output['license'] = license

        if links := data.get('sbg:links'):
            github_links = []
            for link_obj in links:
                if 'https://github.com' in link_obj['id']:
                    github_links.append(link_obj['id'])
            output['codeRepository'] = github_links

        if categories := data.get('sbg:categories'):
            output['applicationSubCategory'] = categories
        if revisions_info := data.get('sbg:revisionsInfo'):
            output['softwareVersion'] = revisions_info
        # if project_name := data.get('sbg:projectName'):
        if tool_author := data.get('sbg:toolAuthor'):
            output['creator'] = tool_author
        if app_version := data.get('appVersion'):
            output['version'] = app_version
        if modified_on := data.get('sbg:modifiedOn'):
            output['dateModified'] = modified_on
        if created_on := data.get('sbg:createdOn'):
            output['dateCreated'] = created_on
        if contributors := data.get('sbg:contributors'):
            output['contributor'] = contributors
        if publisher := data.get('sbg:publisher'):
            output['sdPublisher'] = publisher
        if workflow_language := data.get('sbg:workflowLanguage'):
            output['programmingLanguage'] = workflow_language

        yield output


json_list = []
for i in parse():
    json_list.append(i)

with open('json_data.json', 'w') as outfile:
    json.dump(json_list, outfile)
