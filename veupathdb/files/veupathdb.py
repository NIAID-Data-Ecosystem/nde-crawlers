import logging
from platform import release
import requests 
import json
from datetime import datetime

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def record_generator():

    # Request API call that returns a list of all available data records 
    api_command = 'https://veupathdb.org/veupathdb/service/record-types/dataset/searches/AllDatasets/reports/standard?reportConfig={"attributes":["primary_key","organism_prefix","project_id","eupath_release","newcategory","summary","contact","wdk_weight","version","institution","build_number_introduced","pmids_download","release_policy","short_attribution","type","genecount"],"tables":["Publications","Contacts","GenomeHistory","DatasetHistory","Version","References","HyperLinks","GeneTypeCounts","TranscriptTypeCounts"],"attributeFormat":"text"}'
    request = requests.get(api_command)
    _record_dict = request.json()

    # paginate through records
    for _record_dict in list(_record_dict['records'])[:25]:

        # add the 
        _record_dict.update({
            'veupathdb_id': "veupathdb_"+_record_dict['id'][0]['value'],
            '@type': 'Dataset',
            'includedInDataCatalog.name': "VEuPathDB",
            'url': "https://veupathdb.org/veupathdb/app/record/dataset/",
            
        })
        _record_dict.pop('recordClassName')
        _record_dict.pop('tableErrors')
        _record_dict['name'] = _record_dict['displayName'] 
        _record_dict.pop('displayName') 
        _record_dict['identifier'] = _record_dict['id'][0]['value']
        _record_dict.pop('id')


        # tables.Contacts
        _record_dict['author']=[{'name': _dict.pop('contact_name'), "affiliation": _dict.pop("affiliation")} for _dict in _record_dict['tables']['Contacts']]
        _record_dict['tables'].pop('Contacts')

        # attributes          
        _record_dict['description'] = _record_dict['attributes'].pop('summary')
        _record_dict['measurementTechnique'] = _record_dict['attributes'].pop('type')
        _record_dict['dateUpdated'] = _record_dict['attributes'].pop('version')
        _record_dict['usageInfo'] = _record_dict['attributes'].pop('release_policy')
        _record_dict['sdPublisher.name'] = _record_dict['attributes'].pop('project_id')
        _record_dict['creditText'] = _record_dict['attributes'].pop('short_attribution')
        _record_dict.pop('attributes')
        

        # **** tables.GenomeHistory **** 
        release_dates = [hit['version'] for hit in _record_dict['tables']['Version']]

        # if multiple dates passed, keep the most recent date
        if len(release_dates) > 1:
            try:
                release_date = sorted(
                    release_dates,  # list of dates input 
                    key = lambda d: datetime.strptime(d, '%Y-%m-%d'),  
                    # convert each string into date
                    reverse=True  # for decreasing order 
                )[0]
            except:
                release_dates[0]
        elif not release_dates:
            release_date = None
        else:
            release_date = release_dates[0]
           

        _record_dict['dateUpdated'] = release_date
        _record_dict['tables'].pop('GenomeHistory')
    
        # **** tables.Publications ****
        # - have to use Jasons helper function
        for _dict in _record_dict['tables']['Publications']:
            print(_dict)
            #_dict.pop('dataset_id')
            #_dict.pop('pubmed_link')
            #_dict.pop('citation')
            
        # **** tables.Version ****
        dates = [hit['version'] for hit in _record_dict['tables']['Version']]

        # if multiple dates passed, keep the most recent date
        if len(dates) > 1: 
            try:
                recent_date = sorted(
                    dates,  # list of dates input 
                    key = lambda d: datetime.strptime(d, '%Y-%m-%d'),  
                    # convert each string into date
                    reverse=True  # for decreasing order 
                )[0]
            except:
                release_dates[0]
        elif not dates:
            recent_date = None
        else:
            recent_date = dates[0] 
            
        _record_dict['datePublished'] = recent_date
        _record_dict['species'] = [hit['organism'] for hit in _record_dict['tables']['Version']]
        _record_dict['tables'].pop('Version')

        # **** tables.HyperLinks ****
        _record_dict['distribution'] = [{"text": hit['text'], "url": hit['url']} for hit in _record_dict['tables']['HyperLinks']]
        _record_dict['tables'].pop('HyperLinks')
        

        # **** table.GeneTypeCounts ****
        values = [hit['gene_count'] for hit in _record_dict['tables']['GeneTypeCounts']]
        refs =[hit['gene_type'] for hit in _record_dict['tables']['GeneTypeCounts']]
        if not values:_value = None
        else:_value = values[0]
        if not refs: _ref = None
        else: _ref = refs[0]

        _record_dict['variableMeasured'] = {
            '@type': "PropertyValue",
            'value': _value , 
            'valueReference': _ref
            }
        _record_dict.pop('tables')

        yield _record_dict

rec_gen = record_generator()
for rec in rec_gen:
    print("\n-----\n\n[INFO] record: \n")
    print(json.dumps(rec, indent=4))
    print("\n-----\n")