import datetime 
import logging
from typing import final  
import requests 
import time
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')



def record_generator():

    # Request API call that returns a list of all available data records 
    api_command = 'https://veupathdb.org/veupathdb/service/record-types/dataset/searches/AllDatasets/reports/standard?reportConfig={"attributes":["primary_key","organism_prefix","project_id","eupath_release","newcategory","summary","contact","wdk_weight","version","institution","build_number_introduced","pmids_download","release_policy","short_attribution","type","genecount"],"tables":["Publications","Contacts","GenomeHistory","DatasetHistory","Version","References","HyperLinks","GeneTypeCounts","TranscriptTypeCounts"],"attributeFormat":"text"}'
    request = requests.get(api_command)
    _record_dict = request.json()
    record_list = _record_dict['records']

    # paginate through records
    for _record_dict in record_list[:1]:
        record_dict = {} # initlialize the record dictionary

        # add the 
        record_dict = {
            'veupathdb_id': "veupathdb_"+_record_dict['id'][0]['value'],
            '@type': 'Dataset',
            'includedInDataCatalog.name': "VEuPathDB",
            'url': "https://veupathdb.org/veupathdb/app/record/dataset/",
            
        }
        record_dict['authors'] =  [{"author.name": hit['contact_name'], "author.affiliation": hit['affiliation'] } for hit in _record_dict['tables']['Contacts'] ] # SHOULD BE A LIST
        record_dict['name'] = _record_dict['displayName'] 
        record_dict['identifier'] = _record_dict['id'][0]['value']

        # attributes          
        record_dict['description'] = _record_dict['attributes']['summary']
        record_dict['measurementTechnique'] = _record_dict['attributes']['type']
        record_dict['dateUpdated'] = _record_dict['attributes']['version']
        record_dict['usageInfo'] = _record_dict['attributes']['release_policy']
        record_dict['sdPublisher.name'] = _record_dict['attributes']['project_id']
        record_dict['creditText'] = _record_dict['attributes']['short_attribution']
        #record_dict[''] = ## NEED TO NAME records.attributes.genecount

        # Publications Table -- have to use Jasons helper function
        #if len(_record_dict['tables']["Publications"]) > 1 : print(len(_record_dict['tables']["Publications"]))
        #record_dict['citation.pmid'] = _record_dict['tables']['Publications']['pmid']
        #record_dict['citation.url'] = _record_dict['tables']['Publications']['url']
        
        # Version Table
        #if len(_record_dict['tables']['Version']) > 1 : print(json.dumps(_record_dict['tables']['Version'], indent=4))
        record_dict['species'] = [hit['organism'] for hit in _record_dict['tables']['Version']]
        record_dict['datePublished'] = [hit['version'] for hit in _record_dict['tables']['Version']]

        # HyperLinks Table
        #if len(_record_dict['tables']['Version']) > 1 : print(json.dumps(_record_dict['tables']['Version'], indent=4))
        record_dict['distributions'] = [{"distribution.text": hit['text'], "distribution.url": hit['url']} for hit in _record_dict['tables']['HyperLinks']]

        # GeneTypeCounts Table
        record_dict['variableMeasured'] = {
            '@type': "PropertyValue",
            'value': [hit['gene_count'] for hit in _record_dict['tables']['GeneTypeCounts']], 
            'valueReference': [hit['gene_type'] for hit in _record_dict['tables']['GeneTypeCounts']]
            } 

        # TranscriptTypeCounts Table
        record_dict['TranscriptTypeCounts'] = _record_dict['tables']['TranscriptTypeCounts']

        # GenomeHistory Table
        record_dict['GenomeHistory'] = _record_dict['tables']['GenomeHistory']
        for rec in record_dict['GenomeHistory']:
            del rec['dataset_id']
        
        yield record_dict


generator = record_generator()
for record in generator:
    print(json.dumps(record, indent=4))