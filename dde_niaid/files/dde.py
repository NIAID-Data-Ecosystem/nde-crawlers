import requests
import datetime
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def parse():
    # initial request to find total number of hits
    url = "https://discovery.biothings.io/api/dataset/query?q=_meta.guide:%22/guide/niaid%22&sort=-_ts.last_updated"
    request = requests.get(url).json()
    # get the number of pages to paginate through
    total = request['total']
    pages = (total-1)//1000
    count = 0

    # paginate through the requests
    for page in range(pages+1):
        url = "https://discovery.biothings.io/api/dataset/query?q=_meta.guide:%22/guide/niaid%22&sort=-_ts.last_updated&size=1000&from=" + str(page*1000)
        request = requests.get(url).json()
        for hit in request['hits']:
            count += 1
            if count % 1000 == 0:
                logger.info("Parsed %s records", count)

            # add included in data catalog
            hit['includedInDataCatalog'] = {
                '@type': 'Dataset',
                'name': 'Data Discovery Engine',
                'url': 'https://discovery.biothings.io/',
                'versionDate': datetime.date.today().isoformat()
                }

            # remove unnecessary values
            hit.pop('_meta', None)
            hit.pop('_score', None)

            # adjust date values
            if dates := hit.pop('_ts', None):
                hit['dateCreated'] = datetime.datetime.fromisoformat(dates['date_created']).date().isoformat()
                hit['dateModified'] = datetime.datetime.fromisoformat(dates['last_updated']).date().isoformat()
            
            # adjust @type value to fit our schema
            if nde_type := hit.pop('@type', None):
                nde_type = nde_type.split(":")[-1]
                hit['@type'] = nde_type

            # rename values
            if author := hit.pop('creator', None):
                hit['author'] = author
            hit['_id'] = "DDE_" + hit['_id']

            yield hit

    if count == total:
        logger.info("Total number of documents parsed: %s", total)
    else:
        logger.warning("Did not parse all the records \n"
                       "Total number parsed: %s \n Total number of documents: %s", count, total)













