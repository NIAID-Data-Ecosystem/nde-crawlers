import requests


def get_record():
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
    #yield record