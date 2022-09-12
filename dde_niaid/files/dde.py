import requests
import datetime
import re
import time
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def query_ols(iri):
    """ Gets the name field of measurementTechnique, infectiousAgent, infectiousDisease, and species in our nde schema

        ols api doc here: https://www.ebi.ac.uk/ols/docs/api
        Returns the formatted dictionary {name: ####, url: ####} if an url was given or {name: ####}
    """

    url = "https://www.ebi.ac.uk/ols/api/terms?"
    pattern = re.compile("^https?://")

    if pattern.match(iri):
        params = {
            # isnt the best way but good enough to get http: https://stackoverflow.com/questions/9760588/how-do-you-extract-a-url-from-a-string-using-python
            "iri": re.search("(?P<url>http?://[^\s]+)", iri).group("url")
        }

        request = requests.get(url, params).json()
        # no documentation on how many requests can be made
        time.sleep(0.5)
        return {'name': request['_embedded']['terms'][0]['label'], 'url': iri}
    else:
        return {'name': iri}


def parse():
    # initial request to find total number of hits
    url = "https://discovery.biothings.io/api/dataset/query?q=_meta.guide:%22/guide/niaid/ComputationalTool%22%20OR%20_meta.guide:%22/guide/niaid%22&sort=-_ts.last_updated"
    request = requests.get(url).json()
    # get the number of pages to paginate through
    total = request['total']
    pages = (total - 1) // 1000
    count = 0

    # paginate through the requests
    for page in range(pages + 1):
        url = "https://discovery.biothings.io/api/dataset/query?q=_meta.guide:%22/guide/niaid/ComputationalTool%22%20OR%20_meta.guide:%22/guide/niaid%22&sort=-_ts.last_updated&size=1000&from=" \
              + str(page * 1000)
        request = requests.get(url).json()
        for hit in request['hits']:
            count += 1
            if count % 100 == 0:
                logger.info("Parsed %s records", count)

            # add included in data catalog
            hit['includedInDataCatalog'] = {
                '@type': 'Dataset',
                'name': 'Data Discovery Engine',
                'url': 'https://discovery.biothings.io/',
                'versionDate': datetime.date.today().isoformat()
            }

            # rename our id value and creator to author
            if authors := hit.pop('creator', None):
                if type(authors) is list:
                    for author in authors:
                        if affiliation := author.get('affiliation'):
                            author['affiliation'] = {'name': affiliation}
                else:
                    if affiliation := authors.get('affiliation'):
                        authors['affiliation'] = {'name': affiliation}
                hit['author'] = authors
            hit['_id'] = "DDE_" + hit['_id']

            # adjust date values
            if dates := hit.pop('_ts', None):
                hit['dateCreated'] = datetime.datetime.fromisoformat(
                    dates['date_created']).date().isoformat()
                hit['dateModified'] = datetime.datetime.fromisoformat(
                    dates['last_updated']).date().isoformat()

            # adjust applicationSubCategory to fit our schema
            if app_subs := hit.pop('applicationSubCategory', None):
                hit['applicationSubCategory'] = []
                for app_sub in app_subs:
                    hit['applicationSubCategory'].append(
                        {'name': app_sub.get('name')})

            # adjust @type value to fit our schema
            if nde_type := hit.pop('@type', None):
                nde_type = nde_type.split(":")[-1]
                if "Dataset" in nde_type:
                    hit['@type'] = "Dataset"
                else:
                    hit['@type'] = nde_type

            # query the ols to get measurementTechnique, infectiousAgent, healthCondition (infectiousDisease), and species
            if mts := hit.pop('measurementTechnique', None):
                if type(mts) is list:
                    hit['measurementTechnique'] = []
                    for mt in mts:
                        hit['measurementTechnique'].append(query_ols(mt))
                else:
                    hit['measurementTechnique'] = query_ols(mts)

            if ias := hit.pop('infectiousAgent', None):
                if type(ias) is list:
                    hit['infectiousAgent'] = []
                    for ia in ias:
                        hit['infectiousAgent'].append(query_ols(ia))
                else:
                    hit['infectiousAgent'] = query_ols(ias)

            if hcs := hit.pop('healthCondition', None):
                hit['healthCondition'] = []
                if type(hcs) is list:
                    for hc in hcs:
                        hit['healthCondition'].append(query_ols(hc))
                else:
                    hit['healthCondition'].append(query_ols(hcs))

            # infectious disease is deprecated change to healthCondition
            if ids := hit.pop('infectiousDisease', None):
                if not hit.get('healthCondition'):
                    hit['healthCondition'] = []
                if type(ids) is list:
                    for i_d in ids:
                        hit['healthCondition'].append(query_ols(i_d))
                else:
                    hit['healthCondition'].append(query_ols(ids))

            if species := hit.pop('species', None):
                if type(species) is list:
                    hit['species'] = []
                    for a_species in species:
                        hit['species'].append(query_ols(a_species))
                else:
                    hit['species'] = query_ols(species)

            # remove unnecessary values
            hit.pop('_meta', None)
            hit.pop('_score', None)
            hit.pop('@context', None)

            yield hit

    if count == total:
        logger.info("Total number of documents parsed: %s", total)
    else:
        logger.warning("Did not parse all the records \n"
                       "Total number parsed: %s \n Total number of documents: %s", count, total)
