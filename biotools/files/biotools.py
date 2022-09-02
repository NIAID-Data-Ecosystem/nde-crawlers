import pandas as pd
import requests
import math
import time
import copy
import logging

from datetime import datetime

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def cleanNullTerms(d):
    clean = {}
    for k, v in d.items():
        if v is not None and len(v) > 0:
            clean[k] = v
    return clean


def get_pages(biotoolsapiurl):
    biotoolsapiurl
    # payloads = {'format': 'json', 'domain': 'covid-19', 'page': 'i'}
    payloads = {'format': 'json', 'page': 'i'}
    r = requests.get(biotoolsapiurl, params=payloads, timeout=20).json()
    count = r['count']
    list_num = len(r['list'])
    total_pages = math.ceil(count/list_num)
    return total_pages


def get_dict_key(single_dict):
    keylist = list(single_dict.keys())
    only_key = keylist[0]
    return only_key


def parse_defined_terms(definedtermlist):
    try:
        termdf = pd.DataFrame(definedtermlist)
        termdf.rename(columns={'uri': 'url', 'term': 'name'}, inplace=True)
        termdf['inDefinedTermSet'] = 'http://edamontology.org/'
        termdf['@id'] = termdf['url'].str.replace(
            'http://edamontology.org/', '', regex=True)
        # termdf['@type'] = 'schema:DefinedTerm'
        termdf['termCode'] = termdf['@id']
        cleandf = termdf[['@type', '@id',
                          'inDefinedTermSet', 'termCode', 'url', 'name']]
        termjson = cleandf.to_dict(orient="records")
    except:
        termjson = -1
    return termjson


def parse_parameters(parameterdict):
    datalist = []
    formatlist = []
    for eachparameter in parameterdict:
        content_types = list(eachparameter.keys())
        if 'data' in content_types:
            if isinstance(eachparameter['data'], dict) == True:
                datadictlist = []
                datadictlist.append(eachparameter['data'])
            elif isinstance(eachparameter['data'], list) == True:
                datadictlist = eachparameter['data']
            datavalue = parse_defined_terms(datadictlist)
            if datavalue != -1:
                datalist.extend(datavalue)
        if 'format' in content_types:
            if isinstance(eachparameter['format'], dict) == True:
                formatdictlist = []
                formatdictlist.append(eachparameter['format'])
            elif isinstance(eachparameter['format'], list) == True:
                formatdictlist = eachparameter['format']
            formatvalue = parse_defined_terms(formatdictlist)
            if formatvalue != -1:
                formatlist.extend(formatvalue)
    if len(datalist) > 0 or len(formatlist) > 0:
        parameterjson = {'@type': 'bioschemastypes:FormalParameter',
                         'encodingFormat': datalist,
                         'defaultValue': formatlist}
    else:
        parameterjson = -1
    return parameterjson


class cleandoc:
    now = datetime.now()
    basejson = {
        "@type": "ComputationalTool",  # gettype
        # "@context": {
        #     "schema": "http://schema.org/",
        #     "outbreak": "https://discovery.biothings.io/view/outbreak/",
        #     "dct": "http://purl.org/dc/terms/"
        # },
        "@context": 'https://schema.org/',
        "includedInDataCatalog": {
            "@type": "ComputationalTool",  # gettype
            "name": "bio.tools",  # check period
            "url": "https://bio.tools/",
            'versionDate': datetime.date.today().isoformat()
        }
    }

    def add_basic_info(biotooljsonhit):
        cleanjson = cleandoc.basejson
        cleanjson['name'] = biotooljsonhit['name']
        cleanjson['_id'] = 'biotools_' + biotooljsonhit['name']
        cleanjson['description'] = biotooljsonhit['description']
        cleanjson['identifier'] = biotooljsonhit['biotoolsID']
        cleanjson['url'] = 'https://bio.tools/' + biotooljsonhit['biotoolsID']
        # cleanjson['url'] = biotooljsonhit['homepage']
        # if homepage is github link, add github repo link
        if 'github' in biotooljsonhit['homepage']:
            cleanjson['codeRepository'] = biotooljsonhit['homepage']
        else:
            cleanjson['contentUrl'] = biotooljsonhit['homepage']
        cleanjson['softwareVersion'] = biotooljsonhit['version']
        cleanjson['applicationCategory'] = biotooljsonhit['toolType']
        cleanjson['license'] = biotooljsonhit['license']
        cleanjson['programmingLanguage'] = biotooljsonhit['language']
        if biotooljsonhit['accessibility'] != None:
            cleanjson['conditionsOfAccess'] = biotooljsonhit['accessibility']
        elif biotooljsonhit['cost'] != None:
            cleanjson['conditionsOfAccess'] = biotooljsonhit['cost']
        try:
            cleanjson['dateModified'] = biotooljsonhit['lastUpdate'].split('T')[
                0]
        except:
            cleanjson['dateModified'] = biotooljsonhit['lastUpdate']
        try:
            cleanjson['dateCreated'] = biotooljsonhit['additionDate'].split('T')[
                0]
        except:
            cleanjson['dateCreated'] = biotooljsonhit['additionDate']
        return cleanjson

    def add_app_sub_cat(cleanjson, biotooljsonhit):
        try:
            alltopics = biotooljsonhit['topic']
            topicjson = parse_defined_terms(alltopics)
            if topicjson != -1:
                cleanjson['applicationSubCategory'] = {'name': topicjson}
        except:
            pass
        return cleanjson

    def add_features(cleanjson, biotooljsonhit):
        available_functions = biotooljsonhit['function']
        if len(available_functions) > 0:
            operationlist = []
            inputlist = []
            outputlist = []
        for eachfunction in available_functions:
            # parse operations
            operations = eachfunction['operation']
            features = parse_defined_terms(operations)
            if features != -1:
                operationlist.extend(features)
            # parse inputs
            inparameterdict = eachfunction['input']
            inparameters = parse_parameters(inparameterdict)
            if inparameters != -1:
                inputlist.append(inparameters)
            # parse outputs
            outparameterdict = eachfunction['output']
            outparameters = parse_parameters(outparameterdict)
            if outparameters != -1:
                outputlist.append(outparameters)
                # Note, additional content available from bio.tools in this section includes 'cmd', and 'notes'
            if len(operationlist) > 0:
                cleanjson['featureList'] = operationlist
            if len(inputlist) > 0:
                cleanjson['input'] = inputlist
            if len(outputlist) > 0:
                cleanjson['output'] = outputlist
        return cleanjson

    # get the downloadUrl and codeRepository from either link or download
    def add_links(cleanjson, biotooljsonhit):
        # pool links
        biotoollinks = []
        if isinstance(biotooljsonhit['link'], list) == True:
            for eachdict in biotooljsonhit['link']:
                biotoollinks.append(eachdict)
        if isinstance(biotooljsonhit['download'], list) == True:
            for eachdict in biotooljsonhit['link']:
                biotoollinks.append(eachdict)
        if isinstance(biotooljsonhit['link'], dict) == True:
            biotoollinks.append(biotooljsonhit['link'])
        if isinstance(biotooljsonhit['download'], dict) == True:
            biotoollinks.append(biotooljsonhit['download'])

        # Parse links
        if len(biotoollinks) > 0:
            codeRepository = []
            discussionUrl = []
            downloadUrl = []
            for eachitem in biotoollinks:
                if isinstance(eachitem, dict) == True and eachitem['type'] != None:
                    if ('Repository' or 'repository') in eachitem['type']:
                        codeRepository.append(eachitem['url'])
                    if ('Issue' or 'issue') in eachitem['type']:
                        discussionUrl.append(eachitem['url'])
                    if ('file' in eachitem['type']):
                        downloadUrl.append(eachitem['url'])
                    if ('note' in eachitem.keys()) and (eachitem['note'] != None) and ('image' in eachitem['note']):
                        downloadUrl.append(eachitem['url'])
            if len(codeRepository) > 0:
                cleanjson['codeRepository'] = list(set(codeRepository))
            if len(discussionUrl) > 0:
                cleanjson['discussionUrl'] = list(set(discussionUrl))
            if len(downloadUrl) > 0:
                cleanjson['downloadUrl'] = list(set(downloadUrl))
        return cleanjson

    def add_softwarehelp(cleanjson, biotooljsonhit):
        if len(biotooljsonhit['documentation']) > 0:
            if isinstance(biotooljsonhit['documentation'], list) == True:
                cleanjson['softwareHelp'] = [x['url']
                                             for x in biotooljsonhit['documentation']]
            if isinstance(biotooljsonhit['documentation'], dict) == True:
                cleanjson['softwareHelp'] = [biotooljsonhit['documentation']]
        return cleanjson

    def add_author(cleanjson, biotooljsonhit):
        authorlist = []
        for eachhit in biotooljsonhit['credit']:
            newdict = copy.deepcopy(eachhit)
            if newdict['typeEntity'] == 'Person':
                newdict['@type'] = 'Person'
            elif newdict['typeEntity'] == None:
                if newdict['orcidid'] != None:
                    newdict['@type'] = 'Person'
                # else:
                #     newdict['@type'] = 'dct:Agent'
            else:
                newdict['@type'] = 'Organization'
            newdict['identifier'] = newdict.pop('orcidid', None)
            newdict['role'] = newdict.pop('typeRole', None)
            newdict.pop('gridid', None)
            newdict.pop('rorid', None)
            newdict.pop('fundrefid', None)
            newdict.pop('note', None)
            newdict.pop('typeEntity', None)
            for key, value in dict(newdict).items():
                if (value is None) or (len(value) == 0):
                    del newdict[key]
            authorlist.append(newdict)
        cleanjson['author'] = authorlist
        return cleanjson

    def add_citations(cleanjson, biotooljsonhit):
        isBasisFor = []
        for eachpub in biotooljsonhit['publication']:
            tmppub = copy.deepcopy(eachpub)
            tmppub['@type'] = 'Publication'
            try:
                tmppub['name'] = tmppub['metadata']['title']
            except:
                tmppub['name'] = None
            tmppub.pop('type')
            tmppub.pop('version')
            tmppub.pop('note')
            tmppub.pop('metadata')
            for key, value in dict(tmppub).items():
                if (value is None) or (len(value) == 0):
                    del tmppub[key]
            isBasisFor.append(tmppub)
        cleanjson['isBasisFor'] = isBasisFor

        return cleanjson


def download_jsondocs():
    jsondoclist = []
    biotoolsapiurl = 'https://bio.tools/api/t'
    payloads = {'format': 'json', 'page': '1'}
    # payloads = {'format': 'json', 'domain': 'covid-19', 'page': '1'}
    #payloads = {'format': 'json', 'q': 'COVID-19','page':'1'}
    r = requests.get(biotoolsapiurl, params=payloads, timeout=5).json()
    count = r['count']
    list_num = len(r['list'])
    total_pages = math.ceil(count/list_num)
    i = 1
    while i < total_pages+1:
        logger.info(i)
        if i == 10:
            break
        payloads = {'format': 'json', 'page': i}
        r = requests.get(biotoolsapiurl, params=payloads, timeout=20).json()
        time.sleep(1)
        jsondoclist.extend(r['list'])
        i = i+1
    return jsondoclist


def transform_json(jsondoclist):
    for i in range(len(jsondoclist)):
        biotooljsonhit = jsondoclist[i]
        cleanjson = cleandoc.add_basic_info(biotooljsonhit)
        cleanjson = cleandoc.add_app_sub_cat(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_features(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_links(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_softwarehelp(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_author(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_citations(cleanjson, biotooljsonhit)
        yield cleanNullTerms(cleanjson)


def parse():
    jsondoclist = download_jsondocs()
    doclist = transform_json(jsondoclist)
    yield from doclist
