# get request and print the response of https://bio.tools/api/t with these headers Accept: application/json, text/plain, */*
# paginate with next and previous links

import time
import requests
import json


headers = {'Accept': 'application/json, text/plain, */*'}
response = requests.get('https://bio.tools/api/t', headers=headers)
all_data = []
data = response.json()
for meta in data['list']:
    all_data.append(meta)
page = data['next']
while page is not None:
    start = time.time()
    response = requests.get(f'https://bio.tools/api/t{page}', headers=headers)
    data = response.json()
    for meta in data['list']:
        # if meta['credit'] != []:
        #     for obj in meta['credit']:
        #         if obj['note'] is not None:
        #             print(obj['note'])
        # if meta['community'] is not None:
        #     print(meta['community'])
        meta['test'] = 'https://bio.tools/' + meta['biotoolsID']
        all_data.append(meta)

    if page == '?page=300':
        break
    page = data['next']
    print(page)
# send all data to json file
with open('data.json', 'w') as outfile:
    json.dump(all_data, outfile)
