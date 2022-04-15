import json
from pprint import pprint
import requests

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


# with open('json_data.json', 'a') as outfile:
#     json.dump(all_app_meta_data, outfile)


# for data in all_app_meta_data:
