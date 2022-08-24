import requests
import json
import xmltodict
import logging
import datetime

logger = logging.getLogger('clinepidb-logger')

def record_generator():

    logger.info('record generator function started')
    xml_url = "https://clinepidb.org/ce/sitemap-SitemapDatasets.xml"
    sitemap_request = requests.get(xml_url)
    dict_data = xmltodict.parse(sitemap_request.content)
    data_urls = [hit['loc'] for hit in dict_data["urlset"]["url"]]
    logger.info('%s datasets found'%len(data_urls))
    cond_acc_dict = {
        "public": "open",
        "protected": "restricted", 
        "controlled": "restricted",
        "prerelease": "restricted"
    }

    for url in data_urls:
        dataset_id=url.split("/")[-1]
        
        cookies = {
            'JSESSIONID': '34D7C5C6EDE1FBA34B0956E13896A19D',
        }

        headers = {
            'Accept': '*/*',
            'Accept-Language': 'en-US,en;q=0.9',
            'Cache-Control': 'no-cache',
            'Connection': 'keep-alive',
            'Origin': 'https://clinepidb.org',
            'Pragma': 'no-cache',
            'Referer': 'https://clinepidb.org/ce/app/workspace/analyses/%s/new/details'%dataset_id,
            'content-type': 'application/json'
            }

        data = '{"attributes":["dataset_id","eda_study_id","study_access", "disease","sex","Participant_Type","Sample_Type","Study_Type","Study_Design", \
                "HH_count","Variable_count","approved", "Part_count","percent_approved","total","total_number","Obser_count","Samp_count","Country","Years","AdditionalData","WHO", \
                "ProjectName","contact","institution","description","eupath_release","email","summary"], \
                "primaryKey":[{"name":"dataset_id","value":%s}], \
                "tables":["StudyCharacteristicTable","Publications","DownloadVersion","Contacts","AccessRequestStats","AccessRequest","References","HyperLinks"]}'%dataset_id

        response = requests.post('https://clinepidb.org/ce/service/record-types/dataset/records', headers=headers, cookies=cookies, data=data)
        record = response.json()

        
        # set initial available values
        record['identifier'] = record['id'][0].pop('value')
        record['name'] = record.pop('displayName')

        # add keys
        record.update(
            {
                '_id': "clinepidb_"+record['identifier'],
                '@type': "Dataset",
                'includedInDataCatalog': { 'name': 'ClinEpiDB'},
                'url': 'https://clinepidb.org/ce/app/record/dataset/'+record['identifier']
            }
        )
        try:
            # attributes 
            record['measurementTechnique'] = {'name': record['attributes'].pop('Study_Design')}
            record['description'] = record['attributes'].pop('description')
            record['description'] += record['attributes'].pop('summary')
            record['keywords'] = record['attributes'].pop('WHO')
            record['keywords'] += record['attributes'].pop('Participant_Type')
            record['keywords'] += record['attributes'].pop('Study_Type')
            record['spatialCoverage'] = {'name': record['attributes'].pop('Country')}

            ##version_date =  record['attributes'].pop('eupath_release').split(",")[0]
            #iso_version_date = datetime.datetime.strptime(version_date, "%Y-%m-%d").date().isoformat()
            #record['version'] = {str( record['attributes'].pop('eupath_release'))}
            
            cond_acc = cond_acc_dict[record['attributes']['study_access'].lower()]
            record['conditionOfAccess'] = cond_acc

            temp_cov = record['attributes'].pop('Years')
            record['temporalCoverage'] = {
                'temporalInterval': 
                    {
                        'startDate': (temp_cov.split(' ')[0]), 
                        'endDate': (temp_cov.split(' ')[-1])
                    }
            }
            if not temp_cov:
                record['temporalCoverage'] = {
                    'temporalInterval': 
                        {
                            'startDate':record['tables']['StudyCharacteristicTable'][0]['Years'][0],
                            'endDate': record['tables']['StudyCharacteristicTable'][0]['Years'][-1]
                        }
                }
                        
           
            # tables.StudyCharacteriticTable -- get these variables if the variables above are missing
            health_cond = record['attributes'].pop('disease')
            record['healthCondition'] =  { 'name': health_cond }
            if not health_cond: 
                record['healthCondition'] = {'name': record['tables']['StudyCharacteristicTable'][0].pop('disease')}

            # tables.Publications -- helper function to create citations 
            pmid_list = [hit.pop('pmid') for hit in record['tables']['Publications']]

            if pmid_list:
                if len(pmid_list) == 1:
                    record['pmids'] = pmid_list[0]
                else:
                    record['pmids'] =  ','.join(pmid_list)

            # tables.DowloadVersion
            date_dict = { "JAN": '01', 'FEB': '02', 'MAR': '03', 'APR': '04', 'MAY': '05',
                        "JUN": '06', 'JUL': '07', 'AUG': '08', 'SEP': '09', 'OCT': '10', 
                        'NOV': '11', 'DEC': '12'}
            distr_dic = {}
            recent_date=None
            for hit in record['tables']['DownloadVersion']:
                date = hit['release_date']
                month = date.split("-")[1]
                if month in date_dict.keys():
                    new_date = date.split("-")[0]+'-'+date_dict[month]+'-'+date.split("-")[2]
                    distr_dic[new_date] = hit
            recent_date = sorted([datetime.datetime.strptime(date, "%Y-%m-%d").date().isoformat() for date in distr_dic])[-1]
            record['distribution'] = {
              'dateModified': recent_date,
              'name': distr_dic[recent_date].pop('dataset_name')
             }

            # tables.Contacts -- multiple                
            record['author'] = [{'name': hit['contact_name'], 'affiliation': hit['affiliation']} for hit in record['tables']['Contacts']]

            # tables.HyperLinks
            hl_list = []
            for dict_ in record['tables']['HyperLinks']:
                url_ = dict_.pop('url')
                hl_url = dict_['hyper_link'].pop('url')
                if "NCT" in hl_url:
                    hl_list.append({
                        'identifier': hl_url.split("/")[-1],
                        'url': hl_url
                        })
                else:
                    hl_list.append({'url': hl_url})

            record['isBasedOn'] = hl_list

        except Exception as e:
            logger.error('identifier:' + record['identifier'] + ": " + str(e))
            
        # Remove excess data 
        record.pop('attributes')
        record.pop('tables')
        record.pop('id')
        record.pop('recordClassName')
        record.pop('tableErrors')
        
        yield record