import copy
import datetime
import logging
import math
import time

import pandas as pd
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def cleanNullTerms(d):
    clean = {}
    for k, v in d.items():
        if not isinstance(v, bool):
            if v is not None and len(v) > 0:
                clean[k] = v
    return clean


def get_pages(biotoolsapiurl):
    biotoolsapiurl
    # payloads = {'format': 'json', 'domain': 'covid-19', 'page': 'i'}
    payloads = {"format": "json", "page": "i"}
    r = requests.get(biotoolsapiurl, params=payloads, timeout=20).json()
    count = r["count"]
    list_num = len(r["list"])
    total_pages = math.ceil(count / list_num)
    return total_pages


def get_dict_key(single_dict):
    keylist = list(single_dict.keys())
    only_key = keylist[0]
    return only_key


def parse_defined_terms(definedtermlist):
    try:
        termdf = pd.DataFrame(definedtermlist)
        termdf.rename(columns={"uri": "url", "term": "name"}, inplace=True)
        termdf["inDefinedTermSet"] = "http://edamontology.org/"
        termdf["@id"] = termdf["url"].str.replace("http://edamontology.org/", "", regex=True)
        # termdf['@type'] = 'schema:DefinedTerm'
        termdf["termCode"] = termdf["@id"]
        cleandf = termdf[["@type", "@id", "inDefinedTermSet", "termCode", "url", "name"]]
        termjson = cleandf.to_dict(orient="records")
    except Exception:
        termjson = -1
    return termjson


def parse_parameters(parameterdict):
    datalist = []
    formatlist = []
    for eachparameter in parameterdict:
        content_types = list(eachparameter.keys())
        if "data" in content_types:
            if isinstance(eachparameter["data"], dict) == True:
                datadictlist = []
                datadictlist.append(eachparameter["data"])
            elif isinstance(eachparameter["data"], list) == True:
                datadictlist = eachparameter["data"]
            datavalue = parse_defined_terms(datadictlist)
            if datavalue != -1:
                datalist.extend(datavalue)
        if "format" in content_types:
            if isinstance(eachparameter["format"], dict) == True:
                formatdictlist = []
                formatdictlist.append(eachparameter["format"])
            elif isinstance(eachparameter["format"], list) == True:
                formatdictlist = eachparameter["format"]
            formatvalue = parse_defined_terms(formatdictlist)
            if formatvalue != -1:
                formatlist.extend(formatvalue)
    if len(datalist) > 0 or len(formatlist) > 0:
        parameterjson = {
            "@type": "bioschemastypes:FormalParameter",
            "encodingFormat": datalist,
            "defaultValue": formatlist,
        }
    else:
        parameterjson = -1
    return parameterjson


class cleandoc:
    basejson = {
        "@type": "ComputationalTool",
        # "@context": {
        #     "schema": "http://schema.org/",
        #     "outbreak": "https://discovery.biothings.io/view/outbreak/",
        #     "dct": "http://purl.org/dc/terms/"
        # },
        "@context": "https://schema.org/",
        "includedInDataCatalog": {
            "@type": "ComputationalTool",
            "name": "biotools",
            "url": "https://bio.tools/",
            "versionDate": datetime.date.today().isoformat(),
        },
    }

    def add_basic_info(biotooljsonhit, domain_dict):
        cleanjson = cleandoc.basejson
        cleanjson["name"] = biotooljsonhit["name"]
        cleanjson["_id"] = "biotools_" + biotooljsonhit["biotoolsID"]

        for key in domain_dict:
            if biotooljsonhit["biotoolsID"] in domain_dict[key]["ids"]:
                cleanjson["topicCategory"] = {
                    "name": domain_dict[key]["title"],
                    "url": "https://bio.tools/t?domain=" + domain_dict[key]["domain"],
                    "curatedBy": {"name": "Biotools", "url": "https://bio.tools/domains"},
                }
                if domain_dict[key]["description"] is not None:
                    cleanjson["topicCategory"]["description"] = domain_dict[key]["description"]

        cleanjson["description"] = biotooljsonhit["description"]
        cleanjson["identifier"] = biotooljsonhit["biotoolsID"]
        cleanjson["url"] = "https://bio.tools/" + biotooljsonhit["biotoolsID"]
        # cleanjson['url'] = biotooljsonhit['homepage']
        # if homepage is github link, add github repo link
        if "github" in biotooljsonhit["homepage"]:
            cleanjson["codeRepository"] = biotooljsonhit["homepage"]
        else:
            cleanjson["mainEntityOfPage"] = biotooljsonhit["homepage"]
        cleanjson["softwareVersion"] = biotooljsonhit["version"]
        cleanjson["applicationCategory"] = biotooljsonhit["toolType"]
        cleanjson["license"] = biotooljsonhit["license"]
        cleanjson["programmingLanguage"] = biotooljsonhit["language"]

        if biotooljsonhit["accessibility"] == "Open access":
            cleanjson["conditionsOfAccess"] = "Open"
        elif biotooljsonhit["accessibility"] == "Open access(with restrictions)":
            cleanjson["conditionsOfAccess"] = "Restricted"
        elif biotooljsonhit["accessibility"] == "Restricted access":
            cleanjson["conditionsOfAccess"] = "Closed"

        if biotooljsonhit["cost"] == "Free of charge":
            cleanjson["isAccessibleForFree"] = True
        elif biotooljsonhit["cost"] == "Free of charge (with restrictions)":
            cleanjson["isAccessibleForFree"] = True
        elif biotooljsonhit["cost"] == "Commercial":
            cleanjson["isAccessibleForFree"] = False

        try:
            cleanjson["dateModified"] = biotooljsonhit["lastUpdate"].split("T")[0]
        except Exception:
            cleanjson["dateModified"] = biotooljsonhit["lastUpdate"]
        try:
            cleanjson["dateCreated"] = biotooljsonhit["additionDate"].split("T")[0]
        except Exception:
            cleanjson["dateCreated"] = biotooljsonhit["additionDate"]

        # temp fix for outbreak.info
        if cb_outbreak := cleanjson.get("includedInDataCatalog"):
            cleanjson["curatedBy"] = cb_outbreak

        return cleanjson

    def add_app_sub_cat(cleanjson, biotooljsonhit):
        try:
            alltopics = biotooljsonhit["topic"]
            # topicjson = parse_defined_terms(alltopics)
            # if topicjson != -1:
            topics = [x["term"] for x in alltopics]
            cleanjson["keywords"] = topics
        except Exception:
            pass
        return cleanjson

    def add_features(cleanjson, biotooljsonhit):
        available_functions = biotooljsonhit["function"]
        if len(available_functions) > 0:
            operationlist = []
            inputlist = []
            outputlist = []
        for eachfunction in available_functions:
            # parse operations
            operations = eachfunction["operation"]
            features = parse_defined_terms(operations)
            if features != -1:
                operationlist.extend(features)
            # parse inputs
            inparameterdict = eachfunction["input"]
            inparameters = parse_parameters(inparameterdict)
            if inparameters != -1:
                inputlist.append(inparameters)
            # parse outputs
            outparameterdict = eachfunction["output"]
            outparameters = parse_parameters(outparameterdict)
            if outparameters != -1:
                outputlist.append(outparameters)
                # Note, additional content available from bio.tools in this section includes 'cmd', and 'notes'
            if len(operationlist) > 0:
                cleanjson["featureList"] = operationlist
            if len(inputlist) > 0:
                cleanjson["input"] = inputlist
            if len(outputlist) > 0:
                cleanjson["output"] = outputlist
        return cleanjson

    # get the downloadUrl and codeRepository from either link or download
    def add_links(cleanjson, biotooljsonhit):
        # pool links
        biotoollinks = []
        if isinstance(biotooljsonhit["link"], list) == True:
            for eachdict in biotooljsonhit["link"]:
                biotoollinks.append(eachdict)
        if isinstance(biotooljsonhit["download"], list) == True:
            for eachdict in biotooljsonhit["download"]:
                biotoollinks.append(eachdict)
        if isinstance(biotooljsonhit["link"], dict) == True:
            biotoollinks.append(biotooljsonhit["link"])
        if isinstance(biotooljsonhit["download"], dict) == True:
            biotoollinks.append(biotooljsonhit["download"])

        # Parse links
        if len(biotoollinks) > 0:
            codeRepository = []
            discussionUrl = []
            downloadUrl = []
            for eachitem in biotoollinks:
                if isinstance(eachitem, dict) == True and eachitem["type"] != None:
                    if ("Repository" or "repository") in eachitem["type"]:
                        codeRepository.append(eachitem["url"])
                    if ("Issue" or "issue") in eachitem["type"]:
                        discussionUrl.append(eachitem["url"])
                    if "file" in eachitem["type"]:
                        downloadUrl.append({"name": eachitem["url"]})
                    if ("note" in eachitem.keys()) and (eachitem["note"] != None) and ("image" in eachitem["note"]):
                        downloadUrl.append({"name": eachitem["url"]})
            if len(codeRepository) > 0:
                cleanjson["codeRepository"] = list(set(codeRepository))
            if len(discussionUrl) > 0:
                cleanjson["discussionUrl"] = list(set(discussionUrl))
            if len(downloadUrl) > 0:
                cleanjson["downloadUrl"] = downloadUrl
        return cleanjson

    def add_softwarehelp(cleanjson, biotooljsonhit):
        if len(biotooljsonhit["documentation"]) > 0:
            if isinstance(biotooljsonhit["documentation"], list) == True:
                urls = [{"url": x["url"]} for x in biotooljsonhit["documentation"]]
                cleanjson["softwareHelp"] = urls
            if isinstance(biotooljsonhit["documentation"], dict) == True:
                urls = [{"url": biotooljsonhit["documentation"]}]
                cleanjson["softwareHelp"] = urls
        return cleanjson

    def add_author(cleanjson, biotooljsonhit):
        authorlist = []
        for eachhit in biotooljsonhit["credit"]:
            newdict = copy.deepcopy(eachhit)
            if newdict["typeEntity"] == "Person":
                newdict["@type"] = "Person"
            elif newdict["typeEntity"] == None:
                if newdict["orcidid"] != None:
                    newdict["@type"] = "Person"
                # else:
                #     newdict['@type'] = 'dct:Agent'
            else:
                newdict["@type"] = "Organization"
            newdict["identifier"] = newdict.pop("orcidid", None)
            newdict["role"] = newdict.pop("typeRole", None)
            newdict.pop("gridid", None)
            newdict.pop("rorid", None)
            newdict.pop("fundrefid", None)
            newdict.pop("note", None)
            newdict.pop("typeEntity", None)
            for key, value in dict(newdict).items():
                if (value is None) or (len(value) == 0):
                    del newdict[key]
            authorlist.append(newdict)
        cleanjson["author"] = authorlist
        return cleanjson

    def add_citations(cleanjson, biotooljsonhit):
        citation = []
        for eachpub in biotooljsonhit["publication"]:
            tmppub = copy.deepcopy(eachpub)
            tmppub["@type"] = "Publication"
            try:
                tmppub["name"] = tmppub["metadata"]["title"]
            except Exception:
                tmppub["name"] = None
            tmppub.pop("type")
            tmppub.pop("version")
            tmppub.pop("note")
            tmppub.pop("metadata")
            for key, value in dict(tmppub).items():
                if (value is None) or (len(value) == 0):
                    del tmppub[key]
            citation.append(tmppub)
        cleanjson["citation"] = citation

        return cleanjson


def download_jsondocs():
    logger.info("Retrieving Metadata From API")
    jsondoclist = []
    biotoolsapiurl = "https://bio.tools/api/t"
    payloads = {"format": "json", "page": "1"}
    # payloads = {'format': 'json', 'domain': 'covid-19', 'page': '1'}
    # payloads = {'format': 'json', 'q': 'COVID-19','page':'1'}
    r = requests.get(biotoolsapiurl, params=payloads, timeout=5).json()
    count = r["count"]
    list_num = len(r["list"])
    total_pages = math.ceil(count / list_num)
    i = 1
    while i < total_pages + 1:
        if i % 100 == 0:
            logger.info("Retrieved %s of %s pages" % (i, total_pages))
        if i % 400 == 0:
            break
        payloads = {"format": "json", "page": i}
        r = requests.get(biotoolsapiurl, params=payloads, timeout=20).json()
        time.sleep(1)
        jsondoclist.extend(r["list"])
        i = i + 1
    return jsondoclist


def transform_json(jsondoclist, domain_dict):
    logger.info("Parsing Metadata")
    for i in range(len(jsondoclist)):
        if i % 1000 == 0 and i != 0:
            logger.info("Parsed %s of %s" % (i, len(jsondoclist)))
        biotooljsonhit = jsondoclist[i]
        cleanjson = cleandoc.add_basic_info(biotooljsonhit, domain_dict)
        cleanjson = cleandoc.add_app_sub_cat(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_features(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_links(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_softwarehelp(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_author(cleanjson, biotooljsonhit)
        cleanjson = cleandoc.add_citations(cleanjson, biotooljsonhit)
        yield cleanNullTerms(cleanjson)
    logger.info("Finished Parsing. Total records: %s" % len(jsondoclist))


def get_domain_list():
    domain_dict = {}
    biotools_domains = "https://bio.tools/api/d/all"
    r = requests.get(biotools_domains, timeout=5).json()
    for domain in r["data"]:
        domain_dict[domain["title"]] = {
            "ids": [x["biotoolsID"] for x in domain["resources"]],
            "title": domain["title"],
            "description": domain["description"],
            "domain": domain["domain"],
        }
    return domain_dict


def parse():
    jsondoclist = download_jsondocs()
    domain_dict = get_domain_list()
    doclist = transform_json(jsondoclist, domain_dict)
    yield from doclist
