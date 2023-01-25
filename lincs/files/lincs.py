import requests
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('nde-logger')

class LINCS():
  def parser(self):
    lincsportal_url = "https://lincsportal.ccs.miami.edu/dcic/api/fetchdata?limit=1000&searchTerm=*"
    lincsportal_res = requests.get(lincsportal_url)
    logger.info(f"Started extraction of LINCS database: {lincsportal_url}.")
    logger.info(f"Request status: {lincsportal_res.status_code}")
    lincsportal_data = lincsportal_res.json()

    doc_ct = 0
    success_ct = 0

    for document in lincsportal_data['results']["documents"]:
        # modify doc
        doc_ct += 1

        document['@type'] = "Dataset"
        document['includedInDataCatalog'] = {'name': 'LINCS'}

        if "centerdatasetid" in document:
          document["url"] = document.pop("centerdatasetid")
        #document["description"] = []
        if "assayoverview" in document:
          document["description"] = document.pop('assayoverview')
        if "centerurl" in document:
          document["author"] = {
            "name": document.pop("principalinvestigator"),
            "url": document.pop("centerurl"),
            "affiliation": [{'name': document.pop("centerfullname")}],
            } 
        if "funding" in document:
          document["funding"] = {'identfier': document.pop("funding")}
        if 'datemodified' in document:
          document["dateUpdated"] = document.pop('datemodified')
        if 'screeninglabinvestigator' in document:
          document['author']=[{"name":document.pop("screeninglabinvestigator")}]
        if "datasetname" in document:
          document["name"] = document.pop("datasetname")
        if "datasetid" in document:
          document["identifier"] = document.pop("datasetid")
          document["_id"] = document["identifier"]
        if "datereleased" in document:
          document["datePublished"] = document.pop("datereleased")
        if "assayname" in document:
          measure_tech = document.pop("assayname")
          if 'assayformat' in document:
            measure_tech.append(document.pop("assayformat"))
          document["measurementTechnique"] = [ {'description': measure_tech} ]
        document["keywords"] = []
        if 'assaydesignmethod' in document: 
          document["keywords"] = document.pop("assaydesignmethod")
        if 'physicaldetection' in document:
          document["variableMeasured"] = document.pop("physicaldetection")
        if 'biologicalprocess' in document:
          document["keywords"] = document["keywords"] + document.pop("biologicalprocess")
        if 'technologies' in document:
          document["keywords"].append(document.pop("technologies"))
        if 'biologicalbucket' in document:
          document["keywords"].append(document.pop("biologicalbucket"))
        if 'endpointcategorization' in document:
          document["keywords"].append(document.pop("endpointcategorization"))
        if 'protein' in document:
          for x in document['protein']:
            document['keywords'].append(x)
          document.pop('protein')
        if 'size' in document:
          document["distribution"] = [{"contentsize": document.pop("size")}]
        if "tool" in document and "toollink" in document:
          document["isBasedOn"] = [{"name": document.pop("tool"), "url": document.pop("toollink")}]
        elif "tool" in document:
          document["isBasedOn"] = [{"name": document.pop("tool")}]
        elif "toollink" in document:
          document["isBasedOn"] = [{"url": document.pop("toollink")}]
        if 'datasetgroup' in document:
          document['isRelatedTo'] = [{'url': document.pop('datasetgroup')}]
        if 'cellline' in document:
          for x in document['cellline']:
            document['keywords'].append(x)

        
        if 'concentrations' in document: document.pop("concentrations")
        if 'timepoints' in document: document.pop("timepoints")
        if 'expentimentalcomments' in document: document.pop("expentimentalcomments")
        if 'centerletter' in document: document.pop("centerletter")
        if 'id' in document: document.pop("id")
        if 'path' in document: document.pop("path")
        if 'protocol' in document: document.pop("protocol")
        if 'centername' in document: document.pop("centername")
        if 'versions' in document: document.pop("versions")
        if 'latestversions' in document: document.pop("latestversions")
        if 'datalevels' in document: document.pop("datalevels")
        if 'projectname' in document: document.pop('projectname')
        if 'ldplink' in document: document.pop("ldplink")
        if 'statsvalues' in document: document.pop('statsvalues')
        if 'statsfields' in document: document.pop('statsfields')
        if 'counts' in document: document.pop('counts')
        if 'levelspath' in document: document.pop('levelspath')
        if 'datasetlevels' in document: document.pop('datasetlevels')
        if '_version_' in document: document.pop('_version_')
        if 'endpoints' in document: document.pop('endpoints')
        if 'smlincsidentifier' in document: document.pop('smlincsidentifier')
        if 'pipeline' in document: document.pop('pipeline')

        yield document
        success_ct += 1
    logger.info(f"{doc_ct} documents found from https://lincsportal.ccs.miami.edu/, and {success_ct} documents successfully parsed from that set.")