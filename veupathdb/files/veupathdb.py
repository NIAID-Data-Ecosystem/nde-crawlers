import logging
import requests 
import datetime


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')
#debug_file='/Users/nacosta/Documents/NIAID-Projects/clinepidb/debugging/bad-dates.log'
#logging.basicConfig(filename=debug_file, level=logging.DEBUG)


def record_generator():
   
    # API call that returns a list of all available data records 
    api_command = 'https://veupathdb.org/veupathdb/service/record-types/dataset/searches/AllDatasets/reports/standard?reportConfig={"attributes":["primary_key","organism_prefix","project_id","eupath_release","newcategory","summary","contact","wdk_weight","version","institution","build_number_introduced","pmids_download","release_policy","short_attribution","type","genecount"],"tables":["Publications","Contacts","GenomeHistory","DatasetHistory","Version","References","HyperLinks","GeneTypeCounts","TranscriptTypeCounts"],"attributeFormat":"text"}'
    # send and retrieve request call
    request = requests.get(api_command)
    json_records = request.json()

    # paginate through records
    for _record_dict in json_records['records']:
        # add custom values to the record
        _record_dict.update({
            '_id': "veupathdb_"+_record_dict['id'][0]['value'],
            '@type': 'Dataset',
            'includedInDataCatalog': {'name': "VEuPathDB"},
            'url': "https://veupathdb.org/veupathdb/app/record/dataset/"+_record_dict['id'][0]['value']
            
        })
    
        _record_dict['name'] = _record_dict['displayName']  # set name to records.displayName
        _record_dict['identifier'] = _record_dict['id'][0]['value'] # set identifier to records.id.value
        
        # get pmid for, set as string for helper function
        pmids_list = [ _dict.pop('pmid') for _dict in _record_dict['tables']['Publications']]
        if pmids_list:
            if len(pmids_list) == 1: 
                 _record_dict["pmids"] = pmids_list[0]
            else:
                _record_dict['pmids'] = ','.join(pmids_list)
        # attributes  
        _record_dict['description'] = _record_dict['attributes'].pop('summary')
        _record_dict['measurementTechnique'] = {'name': _record_dict['attributes'].pop('type')}
        _record_dict['sdPublisher'] = {'name': _record_dict['attributes'].pop('project_id')}
        _record_dict['creditText'] = _record_dict['attributes'].pop('short_attribution')  

        if _record_dict['attributes']['release_policy']:
            _record_dict['conditionOfAccess'] = _record_dict['attributes'].pop('release_policy')

        
        if _record_dict['attributes']['version']:
            try:
                date_modified = _record_dict['attributes']['version']
                date_modified = datetime.datetime.strptime(date_modified, '%Y-%m-%d').date().isoformat()
                _record_dict['dateModified'] = date_modified
            except:
                logging.debug("[INFO] BAD DATE FROM _record_dict['attributes']['version']: %s"%_date)

        # tables.Contacts 
        _record_dict['author']=[{'name': _dict.pop('contact_name'), "affiliation": _dict.pop("affiliation")} for _dict in _record_dict['tables']['Contacts']]

        # tables.GenomeHistory
        release_dates = [hit['release_date'] for hit in _record_dict['tables']['GenomeHistory']]
        # if multiple dates passed, keep the most recent date
        if release_dates:
            try:   
                iso_list = [datetime.datetime.strptime(d, '%Y-%m-%d').date().isoformat() for d in release_dates]
                date_updated = sorted(iso_list)[0]
                _record_dict['dateUpdated'] = date_updated
            except:
                ...
                #logging.debug("[INFO] BAD DATE FROM _record_dict['tables']['Version']: %s"%dates)

        
        # tables.Version 
        published_dates = [hit['version'] for hit in _record_dict['tables']['Version']]
        if published_dates:
            try:
                _iso_list = [datetime.datetime.strptime(d, '%Y-%m-%d').date().isoformat() for d in published_dates]
                date_published = sorted(_iso_list)[0]
                _record_dict['datePublished'] = date_published
            except:
                logging.debug("[INFO] BAD DATE FROM _record_dict['tables']['Version']: %s"%dates)
                #print("----\n[INFO] BAD dates: ",dates)
            
        
        if _record_dict['tables']['Version']:
            _record_dict['species'] = {'name': [hit['organism'] for hit in _record_dict['tables']['Version']]} 

        # tables.HyperLinks
        if _record_dict['tables']['HyperLinks']:
            _record_dict['distribution'] = [{"name": hit['text'], "url": hit['url']} for hit in _record_dict['tables']['HyperLinks']]
        
        # table.GeneTypeCounts 
        gene_counts = [hit['gene_count'] for hit in _record_dict['tables']['GeneTypeCounts']]
        gene_refs =[hit['gene_type'] for hit in _record_dict['tables']['GeneTypeCounts']]
        if gene_refs:
            _record_dict['variableMeasured'] = gene_refs[0]
            """_record_dict['variableMeasured'] = {
            '@type': "PropertyValue",
            'identifier': gene_counts[0],
            'name': gene_refs[0]
            }"""

        # remove 
        _record_dict.pop('recordClassName')
        _record_dict.pop('tableErrors')
        _record_dict.pop('displayName') 
        _record_dict.pop('id')
        _record_dict.pop('attributes')
        _record_dict.pop('tables')

        yield _record_dict