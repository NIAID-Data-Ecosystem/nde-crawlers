import requests

url = 'https://service.azul.data.humancellatlas.org/index/catalogs'
sources = requests.get(url).json()[
    'catalogs']['dcp20']['plugins']['repository']['sources']
uuids = []
for source in sources:
    id = source.split('_')[2]
    uuid = id[:8] + '-' + id[8:12] + '-' + \
        id[12:16] + '-' + id[16:20] + '-' + id[20:]
    uuids.append(uuid)

meta = []
for uuid in uuids:
    url = f'https://service.azul.data.humancellatlas.org/index/projects/{uuid}'
    r = requests.get(url)
    print(r.status_code)
    meta.append(r.json())

for data in meta:
    projects = data['projects']
    for project in projects:
        if project['contributedAnalyses'] != []:
            print(project['contributedAnalyses'])
