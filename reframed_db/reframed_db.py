from datetime import datetime
import logging
import json
import requests


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')

assay_list_url = f"https://reframedb.org/api/assay_list"
response = requests.get(assay_list_url)
assays = json.loads(response.text)

assay_ids = []
for assay in assays:
    assay_ids.append(assay['assay_id'])

bad_ids = []
count = 0
for id in assay_ids:
    individual_assay_url = f"https://reframedb.org/api/assay_details?aid={id}"
    response = requests.get(individual_assay_url)
    if response.status_code != 200:
        bad_ids.append(id)
    else:
        count += 1
        metadata = json.loads(response.text)
        if metadata[0]['bibliography']:
            print(metadata)
    if count % 50 == 0:
        print(count)
print(bad_ids)
