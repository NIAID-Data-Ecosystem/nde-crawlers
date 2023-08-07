import datetime
import logging

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


class LINCS:
    def parser(self):
        lincsportal_url = "https://lincsportal.ccs.miami.edu/dcic/api/fetchdata?limit=1000&searchTerm=*"
        lincsportal_res = requests.get(lincsportal_url)
        logger.info(f"Started extraction of LINCS database: {lincsportal_url}.")
        logger.info(f"Request status: {lincsportal_res.status_code}")
        lincsportal_data = lincsportal_res.json()

        doc_ct = 0
        success_ct = 0

        for document in lincsportal_data["results"]["documents"]:
            doc_ct += 1

            document["@type"] = "Dataset"
            document["includedInDataCatalog"] = {
                "name": "LINCS",
                "url": "https://lincsportal.ccs.miami.edu/datasets/",
                "@type": "Dataset",
                "versionDate": datetime.datetime.today().strftime("%Y-%m-%d"),
            }
            document["url"] = f'https://lincsportal.ccs.miami.edu/datasets/view/{document["datasetid"]}'

            if "assayoverview" in document:
                document["description"] = document.pop("assayoverview")

            if "centerurl" in document:
                if ";" in document["principalinvestigator"]:
                    authors = [author.strip() for author in document["principalinvestigator"].split(";")]
                else:
                    authors = [author.strip() for author in document["principalinvestigator"].split(",")]

                if len(authors) > 1:
                    author_list = []
                    for author in authors:
                        author_list.append(
                            {
                                "name": author,
                                "url": document["centerurl"],
                                "affiliation": {"name": document["centerfullname"]},
                            }
                        )
                    document["author"] = author_list
                else:
                    document["author"] = {
                        "name": authors[0],
                        "url": document["centerurl"],
                        "affiliation": {"name": document["centerfullname"]},
                    }
                document.pop("principalinvestigator")
                document.pop("centerurl")
                document.pop("centerfullname")

            if "funding" in document:
                document["funding"] = {"identifier": document.pop("funding")}

            if "datemodified" in document:
                document["dateModified"] = document.pop("datemodified")

            if "screeninglabinvestigator" in document:
                if ";" in document["screeninglabinvestigator"]:
                    authors = [author.strip() for author in document["screeninglabinvestigator"].split(";")]
                else:
                    authors = [author.strip() for author in document["screeninglabinvestigator"].split(",")]

                screening_authors_list = []
                for author in authors:
                    screening_authors_list.append({"name": author})

                # Check if document["author"] exists and it's a list
                if "author" in document and isinstance(document["author"], list):
                    existing_author_names = [
                        author["name"] for author in document["author"]
                    ]  # Get the existing author names
                    for author in screening_authors_list:
                        if author["name"] not in existing_author_names:  # check if the author is not already present
                            document["author"].append(author)
                # Check if document["author"] exists and it's a dict (single author)
                elif "author" in document and isinstance(document["author"], dict):
                    if document["author"]["name"] not in [
                        author["name"] for author in screening_authors_list
                    ]:  # if the existing author is not in the new authors list
                        screening_authors_list.insert(0, document["author"])
                    document["author"] = screening_authors_list
                else:
                    document["author"] = screening_authors_list

                document.pop("screeninglabinvestigator")

            if "datasetname" in document:
                document["name"] = document.pop("datasetname")

            if "datasetid" in document:
                document["identifier"] = document.pop("datasetid")
                document["_id"] = "LINCS_" + document["identifier"]

            if "datereleased" in document:
                document["datePublished"] = document.pop("datereleased")

            if "assayname" in document:
                assayname = document.pop("assayname")
                document["measurementTechnique"] = {"name": ",".join(assayname)}
                if "assayformat" in document:
                    document["measurementTechnique"]["description"] = document.pop("assayformat")

            if "size" in document:
                if len(set(document["size"])) > 1:
                    size_list = []
                    for size in document["size"]:
                        size_list.append({"contentSize": size})
                    document["distribution"] = size_list
                    document.pop("size")
                else:
                    document["distribution"] = {"contentSize": "".join(document.pop("size"))}

            if "physicaldetection" in document:
                document["variableMeasured"] = document.pop("physicaldetection")

            document["keywords"] = []
            if "assaydesignmethod" in document:
                document["keywords"] = document.pop("assaydesignmethod")
            if "biologicalprocess" in document:
                document["keywords"] = document["keywords"] + document.pop("biologicalprocess")
            if "technologies" in document:
                document["keywords"].append(document.pop("technologies"))
            if "biologicalbucket" in document:
                document["keywords"].append(document.pop("biologicalbucket"))
            if "endpointcategorization" in document:
                document["keywords"].append(document.pop("endpointcategorization"))
            if "protein" in document:
                for x in document["protein"]:
                    document["keywords"].append(x)
                document.pop("protein")

            isBasedOn_list = []
            if "tool" in document and "toollink" in document:
                for index, tool in enumerate(document["tool"]):
                    temp_dict = {"name": tool, "url": document["toollink"][index]}
                    isBasedOn_list.append(temp_dict)
                document.pop("tool")
                document.pop("toollink")
            if "protocol" in document:
                isBasedOn_list.append({"url": document.pop("protocol"), "name": f"{assayname[0]} protocol"})
            if len(isBasedOn_list) > 0:
                document[
                    "isBasedOn"
                ] = isBasedOn_list  # [{"name": document.pop("tool"), "url": document.pop("toollink")}]
            elif "tool" in document:
                document.pop("tool")
            elif "toollink" in document:
                document.pop("toollink")

            if "datasetgroup" in document:
                rt_id = document.pop("datasetgroup")
                document["isRelatedTo"] = {
                    "_id": "LINCS_" + rt_id,
                    "identifier": rt_id,
                    "name": rt_id,
                    "url": f"https://lincsportal.ccs.miami.edu/datasets/view/{rt_id}",
                }

            if "cellline" in document:
                for x in document["cellline"]:
                    document["keywords"].append(x)
                document.pop("cellline")
            document["keywords"] = ",".join(document["keywords"])

            if "concentrations" in document:
                document.pop("concentrations")
            if "timepoints" in document:
                document.pop("timepoints")
            if "expentimentalcomments" in document:
                document.pop("expentimentalcomments")
            if "centerletter" in document:
                document.pop("centerletter")
            if "id" in document:
                document.pop("id")
            if "path" in document:
                document.pop("path")
            if "centername" in document:
                document.pop("centername")
            if "versions" in document:
                document.pop("versions")
            if "latestversions" in document:
                document.pop("latestversions")
            if "datalevels" in document:
                document.pop("datalevels")
            if "projectname" in document:
                document.pop("projectname")
            if "ldplink" in document:
                document.pop("ldplink")
            if "statsvalues" in document:
                document.pop("statsvalues")
            if "statsfields" in document:
                document.pop("statsfields")
            if "counts" in document:
                document.pop("counts")
            if "levelspath" in document:
                document.pop("levelspath")
            if "datasetlevels" in document:
                document.pop("datasetlevels")
            if "_version_" in document:
                document.pop("_version_")
            if "endpoints" in document:
                document.pop("endpoints")
            if "smlincsidentifier" in document:
                document.pop("smlincsidentifier")
            if "pipeline" in document:
                document.pop("pipeline")
            if "phosphoprotein" in document:
                document.pop("phosphoprotein")
            if "dockerized_container" in document:
                document.pop("dockerized_container")
            if "iPSC" in document:
                document.pop("iPSC")
            if "primarycell" in document:
                document.pop("primarycell")
            if "differentiatediPSC" in document:
                document.pop("differentiatediPSC")
            if "antibody" in document:
                document.pop("antibody")
            if "centerdatasetid" in document:
                document.pop("centerdatasetid")

            if document["_id"]:
                yield document
                success_ct += 1
        logger.info(
            f"{doc_ct} documents found from https://lincsportal.ccs.miami.edu/, and {success_ct} documents successfully parsed from that set."
        )
