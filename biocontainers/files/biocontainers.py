import datetime
import logging
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')

url = 'https://api.biocontainers.pro/ga4gh/trs/v2/tools?all_fields_search=&sort_order=asc&sort_field=default&offset=0&limit=1000'
logger.info('Retrieving Metadata')

# total 10556
offset = 0
limit = 1000
total = 0
# save metadata for every tool/workflow
tool_metadata = []
while True:
    data_count = 0
    url = f'https://api.biocontainers.pro/ga4gh/trs/v2/tools?limit={limit}&offset={offset}'
    all_data = requests.get(url)
    if all_data.status_code != 200:
        break
    for obj in all_data.json():
        total += 1
        if obj['contains'] != []:
            print(obj)
        tool_metadata.append(obj)
        data_count += 1

    logger.info('Current Page: %s', offset)

    # condition to check if we've reached the last page
    # if data_count != limit:
    #     break

    offset += 1000
print(total)
