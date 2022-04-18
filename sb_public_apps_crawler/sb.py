import json
from pprint import pprint
import requests


def parse():

    mime_type = {
        "txt": "text/plain", "fastq": "text/fastq", "fastq.gz": "text/fastq", "fasta": "application/x-fasta", "fasta.gz": "application/x-fasta", "tar": "application/x-tar", "tar.gz": "application/x-tar", "gtf": "application/x-gtf", "gtf.gz": "application/x-gtf", "html": "text/html", "sam": "application/x-sam", "bam": "application/x-bam", "zip": "application/zip", "vcf": "application/x-vcf", "vcf.gz": "application/x-vcf", "bed": "text/x-bed", "bed.gz": "text/x-bed", "hdf5": "application/x-hdf5"
    }

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

        # metadata has either doc or description, make sure we get both
        if description := data.get('description'):
            output['description'] = description
        if doc := data.get('doc'):
            output['description'] = doc

        # value for input should be list of FormalParameter objects https://bioschemas.org/types/FormalParameter/1.0-RELEASE
        # initial metadata is a list of objects
        if input_object_list := data.get('inputs'):
            result_input_list = []
            for input_object in input_object_list:
                # transform to a new object we create here
                formal_parameter_object = {}
                if name := input_object.get('label'):
                    formal_parameter_object['name'] = name

                if file_types := input_object.get('sbg:fileTypes'):
                    # incase of just one file_type, make an array so we don't iterate through the string, make lower case for mime_type dictionary keys
                    file_types = file_types.lower()
                    file_types = file_types.split(', ')
                # using dictionary mime_type defined above, find correct value transformations by iterating through keys
                    for file_type in file_types:
                        if file_type in mime_type.keys():
                            formal_parameter_object['encodingFormat'] = mime_type[file_type]

                if len(formal_parameter_object):
                    result_input_list.append(formal_parameter_object)

            if len(result_input_list):
                output['input'] = result_input_list
            # to compare transformed input
            # output['original_input'] = input_object_list

            # value for output should be list of FormalParameter objects https://bioschemas.org/types/FormalParameter/1.0-RELEASE
            # reference the input method above, same workflow
        if output_object_list := data.get('outputs'):
            result_output_list = []
            for output_object in output_object_list:
                result_output_object = {}
                if name := output_object.get('label'):
                    result_output_object['name'] = name
                if file_types := output_object.get('sbg:fileTypes'):
                    file_types = [file_types.lower()]
                    if ',' in file_types[0]:
                        file_types = file_types[0].split(', ')
                    for file_type in file_types:
                        if file_type in mime_type.keys():
                            result_output_object['encodingFormat'] = mime_type[file_type]

                if len(result_output_object):
                    result_output_list.append(result_output_object)

            if len(result_output_list):
                output['output'] = result_output_list

        # commenting out for readablility
        # if steps := data.get('steps'):
        #     output['bioschemas:hasPart'] = steps

        if requirements := data.get('requirements'):
            result_list = []
            for obj in requirements:
                result_list.append(obj['class'])
            output['softwareRequirements'] = result_list

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
        if project_name := data.get('sbg:projectName'):
            output['project'] = project_name
        if tool_author := data.get('sbg:toolAuthor'):
            output['creator'] = tool_author
        if app_version := data.get('appVersion'):
            output['version'] = app_version
        if modified_on := data.get('sbg:modifiedOn'):
            output['dateModified'] = modified_on
        if created_on := data.get('sbg:createdOn'):
            output['dateCreated'] = created_on
        if contributors := data.get('sbg:contributors'):
            output['contributor'] = contributors.split(', ')
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
