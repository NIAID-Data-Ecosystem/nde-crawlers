import requests

def record_generator():

    xml_url = "https://clinepidb.org/ce/sitemap-SitemapDatasets.xml"
    sitemap_request = requests.get(xml_url)
    dict_data = xmltodict.parse(sitemap_request.content)
    data_urls = [hit['loc'] for hit in dict_data["urlset"]["url"]]

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
        'Referer': 'https://clinepidb.org/ce/app/workspace/analyses/%s/new/details'%ID_HERE,
        'Sec-Fetch-Dest': 'empty',
        'Sec-Fetch-Mode': 'cors',
        'Sec-Fetch-Site': 'same-origin',
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Safari/537.36',
        'content-type': 'application/json',
        'sec-ch-ua': '".Not/A)Brand";v="99", "Google Chrome";v="103", "Chromium";v="103"',
        'sec-ch-ua-mobile': '?0',
        'sec-ch-ua-platform': '"macOS"',
        'sec-gpc': '1',
        'x-client-retry-count': '0',
        'x-client-wdk-timestamp': '1656002337969',
    }

    data = '{"attributes":[],"primaryKey":[{"name":"dataset_id","value":"DS_010e5612b8"}],"tables":["StudyCharacteristicTable","Publications","DownloadVersion","Contacts","AccessRequestStats","AccessRequest","References","HyperLinks"]}'
    response = requests.post('https://clinepidb.org/ce/service/record-types/dataset/records', headers=headers, cookies=cookies, data=data)
    record = response.json()
    
    # ---- attributes 
    record['identifier'] = record['id'][0].pop('value')
    record['name'] = record.pop('displayName')
    record['measurementTechnique'] = record['attributes'].pop('Study_Design')
    record['description'] = record['attributes'].pop('description')
    record['description'] += record['attributes'].pop('summary')
    record['conditionOfAccess'] = record['attributes'].pop('study_access')
    record['keywords'] = record['attributes'].pop('WHO')
    record['keywords'] += record['attributes'].pop('Participant_Type')
    record['keywords'] += record['attributes'].pop('Study_Type')
    record['healthCondition'] = record['attributes'].pop('disease')
    record['temporalCoverage'] = record['attributes'].pop('Years')
    record['spatialCoverage'] = record['attributes'].pop('Country')
    record['version'] = record['attributes'].pop('eupath_release')

    # tables.StudyCharacteriticTable
    if not record['healthCondition']:
        record['healthCondition'] = record['tables']['StudyCharacteristicTable'][0]['disease']
    if not record['temporalCoverage']:
        record['temporalCoverage'] = record['tables']['StudyCharacteristicTable'][0]['Years']

    # tables.Publications
    record['citation'] = {
        'name': record['tables']['Publications'][0]['pubmed_link']['displayText'],
        'url': record['tables']['Publications'][0]['pubmed_link']['url'],
        'pmid' : record['tables']['Publications'][0]['pmid']
        }

    # tables.DowloadVersion
    record['distribution'] = {
        'dateModified': record['tables']['DownloadVersion'][0]['release_date'],
        'name': record['tables']['DownloadVersion'][0]['dataset_name']
        }

    # tables.Contacts
    record['author'] = {
        'name': record['tables']['Contacts'][0]['contact_name'], 
        'affiliation': record['tables']['Contacts'][0]['affiliation']
        }

    # tables.HyperLinks
    if "NCT" in record['tables']['HyperLinks'][0]['hyper_link']['url']:
    record['isBasedOn'] = {
        'identifier': record['tables']['HyperLinks'][0]['hyper_link']['url'].split("/")[-1],
        'url': record['tables']['HyperLinks'][0]['url']
            }
    else:
        record['isBasedOn'] = {
            'url': record['tables']['HyperLinks'][0]['url']
            }

    # Remove excess data 
    record.pop('attributes')
    record.pop('tables')
    record.pop('id')
    record.pop('recordClassName')
    record.pop('tableErrors')


    yield record