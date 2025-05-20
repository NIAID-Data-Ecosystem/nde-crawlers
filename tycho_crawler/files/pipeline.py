# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface
# mapping https://docs.google.com/spreadsheets/d/19hQf4sQ6ZLwvt8ADYyvooYrqmQCsxh386kYnycX3xmg/edit#gid=0

import datetime
import json
import logging

import dateutil
from biothings.utils.dataload import tab2dict

logger = logging.getLogger("nde-logger")


__all__ = [
    "TychoItemProcessorPipeline",
]


class TychoItemProcessorPipeline:

    # map isAbout.value to variableMeasured field
    variable_measured_map = {
        "case": {"name": "Case", "url": "http://purl.obolibrary.org/obo/NCIT_C49152"},
        "cumulative incidence": {
            "name": "Cumulative incidence",
            "url": "http://purl.obolibrary.org/obo/EPO_0000061",
        },
        "death": {"name": "Dead", "url": "http://purl.obolibrary.org/obo/NCIT_C28554"},
        "dead": {"name": "Dead", "url": "http://purl.obolibrary.org/obo/NCIT_C28554"},
        "complete recovery": {
            "name": "Complete recovery",
            "url": "http://purl.obolibrary.org/obo/NCIT_C82467",
        },
    }

    def process_item(self, item: dict, spider):
        metadata = item.get("metadata")
        identifier = metadata["identifier"]["identifier"].split("/")[-1]
        _id = "tycho_" + identifier.casefold()
        if metadata.get("access"):
            if not isinstance(metadata["access"], list):
                metadata["access"] = [metadata["access"]]
            url = metadata["access"][0]["landingPage"]
        else:
            url = "https://www.tycho.pitt.edu/dataset/" + identifier
        output = {
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": _id,
            "identifier": identifier,
            "sameAs": item.get("zenodo_url"),
            "url": url,
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "Project Tycho",
                "url": "https://www.tycho.pitt.edu/",
                "versionDate": datetime.date.today().isoformat(),
                "dataset": url,
            },
            "conditionsOfAccess": "Closed",
            "species": {
                "@type": "DefinedTerm",
                "alternateName": [
                    "Human",
                    "Homo sapines",
                    "Home sapiens",
                    "Homo sapeins",
                    "human",
                    "Homo sapience",
                    "Homo sapiense",
                    "Homo sampiens",
                    "Homo spaiens",
                    "Homo sapiens (SIRT6)",
                    "Homo sapien",
                    "Homo sapients",
                    "Homo spiens",
                    "Homo sapiens (PARIS)",
                    "Humo sapiens",
                    "Homo sapian",
                    "Homo sapians",
                    "Homo sapiens Linnaeus, 1758",
                ],
                "classification": "host",
                "commonName": "Human",
                "curatedBy": {
                    "name": "Project Tycho",
                    "url": "https://www.tycho.pitt.edu/",
                    "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
                },
                "displayName": "Human | Homo sapiens",
                "identifier": "9606",
                "inDefinedTermSet": "UniProt",
                "isCurated": True,
                "name": "Homo sapiens",
                "originalName": "Homo sapiens",
                "url": "https://www.uniprot.org/taxonomy/9606",
            },
            "topicCategory": [
                {
                    "name": "Public health and epidemiology",
                    "url": "http://edamontology.org/topic_3305",
                    "identifier": "topic_3305",
                    "inDefinedTermSet": "EDAM",
                    "@type": "DefinedTerm",
                },
                {
                    "name": "Infectious disease",
                    "url": "http://edamontology.org/topic_3324",
                    "identifier": "topic_3324",
                    "inDefinedTermSet": "EDAM",
                    "@type": "DefinedTerm",
                },
            ],
        }

        if metadata.get("identifier") and metadata["identifier"].get("identifier"):
            output["identifier"] = metadata["identifier"]["identifier"]
            output["doi"] = output["identifier"]

        if name := metadata.get("title"):
            output["name"] = name

        if description := metadata.get("description"):
            output["description"] = description

        if dates := metadata.get("dates"):
            if not isinstance(metadata["dates"], list):
                dates = [metadata["dates"]]
            for date in dates:
                if date.get("date") and date.get("type") and date["type"].get("value"):
                    if date["type"]["value"].casefold() == "available":
                        output["datePublished"] = date["date"]
                    # TODO check if this is correct in dengue
                    elif date["type"]["value"].casefold() == "intervalequals":
                        output["temporalCoverage"] = [
                            {
                                "temporalInterval": {
                                    "startDate": date["date"].split("-")[0],
                                    "endDate": date["date"].split("-")[1],
                                }
                            }
                        ]
                    elif date["type"]["value"].casefold() == "created":
                        output["dateCreated"] = date["date"]

        if stored := metadata.get("storedIn"):
            if dates := stored.get("dates"):
                for date in dates:
                    if date.get("date") and date.get("type") and date["type"].get("value"):
                        if date["type"]["value"].casefold() == "available":
                            output["datePublished"] = date["date"]
                        elif date["type"]["value"].casefold() == "intervalequals":
                            output["temporalCoverage"] = {"temporalInterval": {"startDate": date["date"]}}
                        elif date["type"]["value"].casefold() == "created":
                            output["dateCreated"] = date["date"]

            if licenses := stored.get("licenses"):
                if not isinstance(metadata["storedIn"]["licenses"], list):
                    licenses = [metadata["storedIn"]["licenses"]]
                output["license"] = [
                    license["identifier"]["identifier"]
                    for license in licenses
                    if license.get("identifier") and license["identifier"].get("identifier")
                ]
                if len(output["license"]) == 1:
                    output["license"] = output["license"][0]

        s_cov = []
        if s_coverage := metadata.get("spatialCoverage"):
            if not isinstance(metadata["spatialCoverage"], list):
                s_coverage = [metadata["spatialCoverage"]]
            for coverage in s_coverage:
                c = {"@type": "Country"}
                if coverage.get("identifier") and coverage["identifier"].get("identifier"):
                    c["identifier"] = coverage["identifier"]["identifier"]
                if coverage.get("name"):
                    c["name"] = coverage["name"]
                s_cov.append(c)
        if s_cov:
            output["spatialCoverage"] = s_cov

        vm = {
            "name": "Count of disease cases",
            "url": "http://purl.obolibrary.org/obo/APOLLO_SV_00000497",
            "curatedBy": {"name": "NIAID Data Ecosystem"},
            "isCurated": True,
        }
        vm_set = set()
        vm_set.add(json.dumps(vm, sort_keys=True))
        mt_list = []
        if types := metadata.get("types"):
            if not isinstance(metadata["types"], list):
                types = [metadata["types"]]
            for t in types:
                vm = {}
                mt = {"@type": "DefinedTerm"}
                if information := t.get("information"):
                    if name := information.get("value"):
                        vm["name"] = name
                    if url := information.get("valueIRI"):
                        vm["identifier"] = url
                    if vm.get("name") or vm.get("identifier"):
                        vm_set.add(json.dumps(vm, sort_keys=True))
                if method := t.get("method"):
                    if name := method.get("value"):
                        if "surveillance" in name.casefold():
                            mt["name"] = "Surveillance"
                            mt["url"] = "http://purl.obolibrary.org/obo/NCIT_C15719"
                        elif "notification" in name.casefold():
                            mt["name"] = "Notification"
                            mt["url"] = "http://purl.obolibrary.org/obo/NCIT_C25297"
                        else:
                            logger.info("Unknown measurement technique type: %s", name)

                        if identifier := method.get("valueIRI"):
                            mt["identifier"] = identifier
                        if mt.get("name") or mt.get("identifier"):
                            mt_list.append(mt)

        if vm_set:
            output["variableMeasured"] = [json.loads(vm) for vm in vm_set]
        if mt_list:
            output["measurementTechnique"] = mt_list

        if keywords := metadata.get("refinement"):
            output["keywords"] = keywords

        if distributions := metadata.get("distributions"):
            d_list = []
            if not isinstance(metadata["distributions"], list):
                distributions = [metadata["distributions"]]

            for distribution in distributions:
                d = {}
                if distribution.get("access") and distribution["access"].get("landingPage"):
                    d["contentUrl"] = distribution["access"]["landingPage"]
                if formats := distribution.get("formats"):
                    if len(formats) == 1:
                        d["encodingFormat"] = formats[0]
                    else:
                        d["encodingFormat"] = formats
                if d:
                    d_list.append(d)

            if d_list:
                output["distribution"] = d_list

        a_list = []

        if authors := metadata.get("creators"):
            if not isinstance(metadata["creators"], list):
                authors = [metadata["creators"]]
            for author in authors:
                a = {}
                if author_id := author.get("identifier"):
                    if author_id.get("identifierSource") and author_id["identifierSource"] == "URL":
                        a["@type"] = "Organization"
                        a["url"] = author_id["identifier"]
                        a["name"] = author.get("name")
                    elif (
                        author_id.get("identifierSource") and "fairsharing" in author_id["identifierSource"].casefold()
                    ):
                        a["identifier"] = "https://orcid.org/" + author_id["identifier"]
                    elif author_id.get("identifier"):
                        a["identifier"] = author_id["identifier"]
                if name := author.get("fullName"):
                    a["name"] = name
                if first_name := author.get("firstName"):
                    a["givenName"] = first_name
                if last_name := author.get("lastName"):
                    a["familyName"] = last_name

                aff_list = []
                if affiliations := author.get("affiliation"):
                    if not isinstance(author["affiliation"], list):
                        affiliations = [author["affiliation"]]
                    for affiliation in affiliations:
                        aff = {}
                        if affiliation.get("name"):
                            aff["name"] = affiliation["name"]
                        if affiliation.get("identifier") and affiliation["identifier"].get("identifier"):
                            aff["url"] = affiliation["identifier"]["identifier"]
                        if aff:
                            aff_list.append(aff)
                if aff_list:
                    a["affiliation"] = aff_list
                if a:
                    a_list.append(a)

        f_list = []
        if acknowledges := metadata.get("acknowledges"):
            if not isinstance(metadata["acknowledges"], list):
                acknowledges = [metadata["acknowledges"]]
            for acknowledge in acknowledges:
                if funders := acknowledge.get("funders"):
                    if not isinstance(acknowledge["funders"], list):
                        funders = [acknowledge["funders"]]
                    for funder in funders:
                        f = {}
                        if funder.get("identifier") and funder["identifier"].get("identifier"):
                            f["identifier"] = funder["identifier"]["identifier"]
                        if funder.get("name"):
                            f["name"] = funder["name"]
                        if funder.get("abbreviation"):
                            f["alternateName"] = funder["abbreviation"]
                        if f:
                            f_list.append(f)

            if f_list:
                output["funding"] = {"funder": f_list}

                if authors := acknowledge.get("awardees"):
                    if not isinstance(acknowledge["awardees"], list):
                        authors = [acknowledge["awardees"]]
                    for author in authors:
                        a = {}
                        if name := author.get("fullName"):
                            a["name"] = name
                        if first_name := author.get("firstName"):
                            a["givenName"] = first_name
                        if last_name := author.get("lastName"):
                            a["familyName"] = last_name
                        if affiliations := author.get("affiliations"):
                            if not isinstance(author["affiliations"], list):
                                affiliations = [author["affiliations"]]
                            aff_list = []
                            for affiliation in affiliations:
                                aff = {}
                                if affiliation.get("identifier") and affiliation["identifier"].get("identifier"):
                                    aff["url"] = affiliation["identifier"]["identifier"]
                                if aff:
                                    aff_list.append(aff)
                            if aff_list:
                                a["affiliation"] = aff_list
                        if a:
                            a_list.append(a)
        if a_list:
            output["author"] = a_list

        is_based_on_list = []
        if relatedIdentifiers := metadata.get("relatedIdentifiers"):
            if not isinstance(metadata["relatedIdentifiers"], list):
                relatedIdentifiers = [metadata["relatedIdentifiers"]]
            for relatedIdentifier in relatedIdentifiers:
                is_based_on = {}
                if relatedIdentifier.get("relationType") and relatedIdentifier["relationType"].get("value"):
                    if relatedIdentifier["relationType"]["value"].casefold() == "isderivedfrom":
                        if relatedIdentifier.get("identifier") and relatedIdentifier["identifier"].get("identifier"):
                            is_based_on["url"] = relatedIdentifier["identifier"]["identifier"]
                            is_based_on_list.append(is_based_on)
        if is_based_on_list:
            output["isBasedOn"] = is_based_on_list

        credit_text = ""
        if acknowledgements := metadata.get("acknowledgements"):
            if not isinstance(metadata["acknowledgements"], list):
                acknowledgements = [metadata["acknowledgements"]]
            for acknowledgement in acknowledgements:
                if acknowledgement.get("name"):
                    credit_text += acknowledgement["name"] + ", "

        if title := metadata.get("title"):
            credit_text += title + "(version 2.0, April 1, 2018): Project Tycho data release, "

        credit_text += "DOI: " + metadata["identifier"]["identifier"]
        output["creditText"] = credit_text

        ia_list = []
        hc_list = []
        if about := metadata.get("isAbout"):
            if not isinstance(metadata["isAbout"], list):
                about = [metadata["isAbout"]]
            for a in about:
                if a.get("name") and a.get("identifier"):
                    url = a["identifier"]["identifier"]
                    if "http://purl.bioontology.org/ontology/NCBITAXON" in url:
                        ia = {
                            "@type": "DefinedTerm",
                            "name": a["name"],
                            "url": f"{url.replace('http://purl.bioontology.org/ontology/NCBITAXON/', 'https://www.uniprot.org/taxonomy/')}",
                            "isCurated": True,
                            "curatedBy": {
                                "name": "Project Tycho",
                                "url": "https://www.tycho.pitt.edu/",
                                "dateModified": datetime.datetime.now().strftime("%Y-%m-%d"),
                            },
                            "inDefinedTermSet": "UniProt",
                        }
                        ia_list.append(ia)
                    elif "http://purl.bioontology.org/ontology/SNOMEDCT" in url:
                        d = tab2dict("./healthConditions.csv", cols=(1, 2, 3, 4), key=0, sep=",")
                        hc = {}
                        if hc_info := d.get(url):
                            hc = {
                                "originalName": hc_info[0],
                                "name": hc_info[1],
                                "url": hc_info[2],
                                "isCurated": True,
                            }
                        else:
                            hc = {
                                "name": a["name"],
                                "url": url,
                            }
                        hc_list.append(hc)

                if a.get("value") and a.get("valueIRI"):

                    # variable measured
                    if a["value"].casefold() in self.variable_measured_map:
                        vm = self.variable_measured_map[a["value"].casefold()]
                        if output.get("variableMeasured"):
                            output["variableMeasured"].append(vm)
                        else:
                            output["variableMeasured"] = [vm]

                    # spatial coverage
                    if "https://www.iso.org/obp/ui" in a["valueIRI"]:
                        sp = {
                            "@type": "AdministrativeArea",
                            "identifier": a["valueIRI"],
                            "name": a["value"],
                            "country": item.get("country"),
                            "locationType": "other",
                        }
                        if output.get("spatialCoverage"):
                            output["spatialCoverage"].append(sp)
                        else:
                            output["spatialCoverage"] = [sp]

                    # temporal coverage
                    if "http://reference.data.gov.uk/id/gregorian-interval/" in a["valueIRI"]:
                        dates = a["value"]
                        parts = dates.split("-")
                        start_date = "-".join(parts[:3])
                        end_date = "-".join(parts[3:])
                        start_date = dateutil.parser.parse(start_date, ignoretz=True).date()
                        end_date = dateutil.parser.parse(end_date, ignoretz=True).date()
                        tc = {
                            "temporalInterval": {
                                "startDate": start_date.isoformat(),
                                "endDate": end_date.isoformat(),
                                "duration": str(end_date - start_date),
                                "name": "count period",
                                "temporalType": "other",
                            }
                        }
                        if output.get("temporalCoverage"):
                            output["temporalCoverage"].append(tc)
                        else:
                            output["temporalCoverage"] = [tc]

        if ia_list:
            output["infectiousAgent"] = ia_list
        if hc_list:
            output["healthCondition"] = hc_list

        return output
