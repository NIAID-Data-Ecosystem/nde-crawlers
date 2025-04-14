import datetime
import logging

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")


class LINCS:
    def parse_authors(self, authors_string):
        # If there is only one author
        if "," not in authors_string:
            # Check for initials without periods and add periods if needed
            authors_string = " ".join([name + "." if len(name) == 1 else name for name in authors_string.split()])
            return [authors_string]
        else:
            # Handle multiple authors, first replace semicolons with commas if needed
            if ";" in authors_string:
                authors_string = authors_string.replace(";", ",")
            potential_authors = authors_string.split(",")
            authors = []
            i = 0
            while i < len(potential_authors):
                potential_author = potential_authors[i].strip()
                if " " in potential_author or i + 1 == len(potential_authors):
                    authors.append(potential_author)
                else:
                    author = potential_author + ", " + potential_authors[i + 1].strip()
                    authors.append(author)
                    i += 1
                i += 1
            authors = [
                " ".join([name + "." if len(name) == 1 else name for name in author.split()]) for author in authors
            ]
            return list(set(authors))

    def parser(self):
        lincsportal_url = "https://lincsportal.ccs.miami.edu/dcic/api/fetchdata?limit=1000&searchTerm=*"
        lincsportal_res = requests.get(lincsportal_url, verify=False)
        logger.info(f"Started extraction of LINCS database: {lincsportal_url}.")
        logger.info(f"Request status: {lincsportal_res.status_code}")
        lincsportal_data = lincsportal_res.json()

        docs = []  # Collect all processed documents here

        for document in lincsportal_data["results"]["documents"]:
            document["@type"] = "Dataset"
            document["includedInDataCatalog"] = {
                "name": "LINCS",
                "url": "https://lincsportal.ccs.miami.edu/datasets/",
                "@type": "DataCatalog",
                "versionDate": datetime.datetime.today().strftime("%Y-%m-%d"),
            }
            url = f'https://lincsportal.ccs.miami.edu/datasets/view/{document["datasetid"]}'
            document["url"] = url
            document["includedInDataCatalog"]["dataset"] = url

            if "assayoverview" in document:
                document["description"] = document.pop("assayoverview")

            author_list = []
            if "centerurl" in document:
                if "principalinvestigator" in document:
                    authors = self.parse_authors(document["principalinvestigator"])
                    for author in authors:
                        author_list.append(
                            {
                                "name": author,
                                "affiliation": {"name": document["centerfullname"]},
                                "url": document["centerurl"],
                            }
                        )
                document.pop("principalinvestigator")
                document.pop("centerurl")
                document.pop("centerfullname")

            if "funding" in document:
                document["funding"] = {"identifier": document.pop("funding")}

            if "datemodified" in document:
                try:
                    date = datetime.datetime.strptime(document["datemodified"], "%Y-%m-%d")
                    document["dateModified"] = date.strftime("%Y-%m-%d")
                except ValueError:
                    logger.warning("Invalid date format in datemodified: " + document["datemodified"])
                document.pop("datemodified")

            if "screeninglabinvestigator" in document:
                authors = self.parse_authors(document["screeninglabinvestigator"])
                for author in authors:
                    author_list.append({"name": author})
                document.pop("screeninglabinvestigator")

            # Remove duplicate authors
            no_dupes = []
            name_dict = {}
            for author_obj in author_list:
                if author_obj["name"] not in name_dict:
                    no_dupes.append(author_obj)
                    name_dict[author_obj["name"]] = author_obj
            document["author"] = no_dupes

            if "datasetname" in document:
                document["name"] = document.pop("datasetname")

            if "datasetid" in document:
                document["identifier"] = document.pop("datasetid")
                document["_id"] = document["identifier"]

            if "datereleased" in document:
                try:
                    date = datetime.datetime.strptime(document["datereleased"], "%Y-%m-%d")
                    document["datePublished"] = date.strftime("%Y-%m-%d")
                except ValueError:
                    logger.warning("Invalid date format in datereleased: " + document["datereleased"])
                document.pop("datereleased")

            if "assayname" in document:
                assay_data = document.pop("assayname")
                if isinstance(assay_data, str):
                    assay_data = [name.strip() for name in assay_data.split(",")]
                else:
                    assay_data = [name.strip() for name in assay_data]
                unique_assays = list(dict.fromkeys(assay_data))
                document["measurementTechnique"] = [{"name": assay} for assay in unique_assays]

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
                document["variableMeasured"] = {"name": document.pop("physicaldetection")}

            keywords_set = set()
            if "assaydesignmethod" in document:
                method = document.pop("assaydesignmethod")
                if isinstance(method, str):
                    keywords_set.update([k.strip() for k in method.split(",")])
                else:
                    keywords_set.update(method)
            if "biologicalprocess" in document:
                process = document.pop("biologicalprocess")
                if isinstance(process, str):
                    keywords_set.update([k.strip() for k in process.split(",")])
                else:
                    keywords_set.update(process)
            if "technologies" in document:
                keywords_set.add(document.pop("technologies"))
            if "biologicalbucket" in document:
                keywords_set.add(document.pop("biologicalbucket"))
            if "endpointcategorization" in document:
                keywords_set.add(document.pop("endpointcategorization"))
            if "protein" in document:
                for x in document["protein"]:
                    keywords_set.add(x)
                document.pop("protein")
            if "cellline" in document:
                for x in document["cellline"]:
                    keywords_set.add(x)
                document.pop("cellline")
            if keywords_set:
                document["keywords"] = list(keywords_set)

            # Leave isBasedOn and other fields unchanged...
            if "tool" in document and "toollink" in document:
                isBasedOn_list = []
                for index, tool in enumerate(document["tool"]):
                    isBasedOn_list.append({"name": tool, "url": document["toollink"][index]})
                document["isBasedOn"] = isBasedOn_list
                document.pop("tool")
                document.pop("toollink")
            if "protocol" in document:
                document.setdefault("isBasedOn", []).append({"url": document.pop("protocol"), "name": "protocol"})

            # Only change isRelatedTo: initially set it based on datasetgroup, then update later.
            if "datasetgroup" in document:
                rt_id = document.pop("datasetgroup")
                document["isRelatedTo"] = {
                    "_id": rt_id.lower(),
                    "identifier": rt_id,
                    "name": rt_id,  # temporary; will be updated below
                    "url": f"https://lincsportal.ccs.miami.edu/datasets/view/{rt_id}",
                }

            # Remove unused fields
            for key in [
                "concentrations",
                "timepoints",
                "expentimentalcomments",
                "centerletter",
                "id",
                "path",
                "centername",
                "versions",
                "latestversions",
                "datalevels",
                "projectname",
                "ldplink",
                "statsvalues",
                "statsfields",
                "counts",
                "levelspath",
                "datasetlevels",
                "_version_",
                "endpoints",
                "smlincsidentifier",
                "pipeline",
                "phosphoprotein",
                "dockerized_container",
                "iPSC",
                "primarycell",
                "differentiatediPSC",
                "antibody",
                "centerdatasetid",
            ]:
                document.pop(key, None)

            docs.append(document)

        # Build mapping from dataset identifier to proper dataset name.
        mapping = {doc["identifier"]: doc.get("name", doc["identifier"]) for doc in docs}

        # Update isRelatedTo field using the mapping.
        for doc in docs:
            if "isRelatedTo" in doc:
                rt_id = doc["isRelatedTo"]["identifier"]
                proper_name = mapping.get(rt_id, rt_id)
                doc["isRelatedTo"]["name"] = proper_name

        for doc in docs:
            yield doc
