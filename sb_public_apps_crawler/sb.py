import json
from pprint import pprint
import requests

url = "https://igor.sbgenomics.com/ns/amigo/api/v1/apps/explore?limit=20&offset=0&order_by=label"


def parse():

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
        if inputs := data.get('inputs'):
            output['bioschemas:input'] = inputs
        if outputs := data.get('outputs'):
            output['bioschemas:input'] = outputs
        if steps := data.get('steps'):
            output['bioschemas:hasPart'] = steps
        if requirements := data.get('requirements'):
            output['schema:softwareRequirements'] = requirements
        if image_url := data.get('sbg:image_url'):
            output['thumbnailUrl'] = image_url
        if toolkit := data.get('sbg:toolkit'):
            output['applicationSuite'] = toolkit
        if license := data.get('license'):
            output['license'] = license
        if links := data.get('sbg:links'):
            output['codeRepository'] = links
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
