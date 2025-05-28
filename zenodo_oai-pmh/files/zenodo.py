# to query one document: https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:5680920&metadataPrefix=oai_datacite
# to query all documents: https://zenodo.org/oai2d?verb=ListRecords&metadataPrefix=oai_datacite
# import re
import datetime
import json
import logging
from xml.etree import ElementTree

import dateutil
from sql_database import OAIDatabase

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


class Zenodo(OAIDatabase):
    # override variables
    SQL_DB = "zenodo.db"
    EXPIRE = datetime.timedelta(days=180)

    INTERVAL = {"days": 60}

    HOST = "https://zenodo.org/oai2d"
    METADATA_PREFIX = "oai_datacite"
    SLEEP_COUNT = 1000
    SLEEP_LENGTH = 0.5

    def __init__(self, sql_db=None):
        super().__init__(sql_db)

    def record_data(self, record):
        """Retrives the raw data using a sickle request and formats so dump can store it into the cache table
        Returns:
            A tuple (_id, data formatted as a json string)
        """
        # in each doc we want record.identifier and record stored
        doc = {
            "header": dict(record.header),
            "metadata": record.metadata,
            "xml": ElementTree.tostring(record.xml, encoding="unicode"),
        }

        return record.header.identifier, json.dumps(doc)

    def get_version_id(self, output, related_ids, url):
        version_id = []
        for related_id in related_ids:
            keys = related_id.keys()
            if sorted(keys) == sorted(["relatedIdentifierType", "relationType"]):
                if (
                    related_id.get("relatedIdentifierType") == "DOI"
                    and related_id.get("relationType") == "IsVersionOf"
                    and related_id.text
                    and "zenodo" in related_id.text
                ):
                    version_id.append(related_id.text)

        # It looks like zenodo oai has fixed their schema where they now specify IsVersionOf instead of a vague
        # IsPartOf for everything so we dont need to query the datacite anymore to find the original version
        # if version_id:
        #     if len(version_id) > 1:
        #         logger.info(
        #             "There is more than one version recordID in %s. Querying datacite to find correct version"
        #             % output["_id"]
        #         )
        #         # reset version_id and use datacite query to find correct version
        #         version_id = []
        #         record = sickle.GetRecord(identifier=f"oai:zenodo.org:{identifier}", metadataPrefix="datacite")
        #         versions = record.xml
        #         versions = versions.findall(".//{http://datacite.org/schema/kernel-4}relatedIdentifier")
        #         time.sleep(0.5)
        #         for version in versions:
        #             if (
        #                 version.get("relatedIdentifierType") == "DOI"
        #                 and version.get("relationType") == "IsVersionOf"
        #                 and "zenodo" in version.text
        #             ):
        #                 version_id.append(version.text)

        if version_id:
            version_id = list(set(version_id))
            if len(version_id) > 1:
                version_id = [min(version_id, key=lambda doi: int(doi.split(".")[-1]))]

            # sanity check
            assert len(version_id) <= 1, "There is more than one version per recordID: %s. Versions: %s" % (
                output["_id"],
                version_id,
            )
            output["versionId"] = version_id[0]
            output["sameAs"] = [url]
            output["doi"] = version_id[0]
            output["identifier"] = version_id[0].rsplit("/", 1)[-1]
            output["_id"] = "ZENODO_" + output["identifier"].rsplit(".", 1)[-1]
            output["url"] = "https://zenodo.org/record/" + output["identifier"].rsplit(".", 1)[-1]

        return output

    def parse_xml(self, xml, output, gen_type, missing_types, url, identifier):

        # use xml to query doi
        root = ElementTree.fromstring(xml)

        # used for testing to print out xml tags
        # for element in root.iter():
        #     print("%s - %s" % (element.tag, element.text))

        doi = root.find(".//{http://datacite.org/schema/kernel-4}identifier[@identifierType='DOI']")
        if doi is not None:
            output["doi"] = doi.text

        # use xml to get the type
        zenodo_type = root.find(".//{http://datacite.org/schema/kernel-4}resourceType[@resourceTypeGeneral]").get(
            "resourceTypeGeneral"
        )
        zenodo_type2 = root.find(".//{http://datacite.org/schema/kernel-4}resourceType[@resourceTypeGeneral]").text

        # format the types to be case insensitive and query the dictionary for the transformation
        zenodo_type = zenodo_type.lower()
        if zenodo_type2:
            zenodo_type2 = zenodo_type2.lower()
            zenodo_type2 = zenodo_type2.replace(" ", "")

        if zenodo_type is not None:
            if zenodo_type in gen_type.keys():
                output["@type"] = gen_type[zenodo_type]
            elif zenodo_type2 in gen_type.keys():
                output["@type"] = gen_type[zenodo_type2]
            else:
                missing_types[(zenodo_type, zenodo_type2)] = (zenodo_type, zenodo_type2)

        # xml to find all the authors and format them
        # we cannot use metadata to format due to creator and affiliation being in separate lists
        creators = root.findall(".//{http://datacite.org/schema/kernel-4}creator")
        for creator in creators:
            author = {}
            name = creator.find("./{http://datacite.org/schema/kernel-4}creatorName")
            affiliation = creator.find("./{http://datacite.org/schema/kernel-4}affiliation")
            orcid_id = creator.find(
                "./{http://datacite.org/schema/kernel-4}nameIdentifier[@nameIdentifierScheme='ORCID']"
            )
            if name is not None:
                author["name"] = name.text
            if affiliation is not None:
                # elasticsearch cannot index strings greater than 32766. Someone mistakenly put a WHOLE ARTICLE into the affiliation
                # https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:4675716&metadataPrefix=oai_datacite
                if len(affiliation.text) < 30000:
                    author["affiliation"] = {"name": affiliation.text}
            if orcid_id is not None:
                author["identifier"] = orcid_id.text
            if author:
                output["author"].append(author)

        # use xml to find the license and conditionOfAccess
        rights = root.findall(".//{http://datacite.org/schema/kernel-4}rights")
        for right in rights:
            if right.text and "access" in right.text.lower():
                enum = ["Open", "Restricted", "Closed", "Embargoed"]
                if right.text.split(" ", 1)[0] not in enum:
                    logger.error("Conditions of Access not valid: %s", right.text.split(" ", 1)[0])
                else:
                    output["conditionsOfAccess"] = right.text.split(" ", 1)[0]
            else:
                output["license"] = right.get("rightsURI")

        # use xml to find citedBy field
        related_ids = root.findall(".//{http://datacite.org/schema/kernel-4}relatedIdentifier")
        cited_by = []
        for related_id in related_ids:
            if related_id.get("relationType") == "IsCitedBy" and related_id.get("relatedIdentifierType") == "URL":
                cited_by.append({"url": related_id.text})
        if cited_by:
            output["citedBy"] = cited_by

        # use xml to find versioning
        related_ids = root.findall(".//{http://datacite.org/schema/kernel-4}relatedIdentifier")
        output = self.get_version_id(output, related_ids, url)

        # use xml to find funding
        contributors = root.findall(".//{http://datacite.org/schema/kernel-4}contributor[@contributorType='Funder']")
        fundings = []
        for contributor in contributors:
            name = contributor.find("./{http://datacite.org/schema/kernel-4}contributorName")
            if name is not None:
                fundings.append({"funder": {"name": name.text}})

        funders = root.findall(".//{http://datacite.org/schema/kernel-4}fundingReference")

        for funder in funders:
            funding = {}
            funder_identifier = funder.find(".//{http://datacite.org/schema/kernel-4}funderIdentifier")
            funder_name = funder.find(".//{http://datacite.org/schema/kernel-4}funderName")
            funding_name = funder.find(".//{http://datacite.org/schema/kernel-4}awardTitle")
            funding_identifier = funder.find(".//{http://datacite.org/schema/kernel-4}awardNumber")
            if funder_identifier is not None or funder_name is not None:
                funding = {"funder": {}}
                if funder_identifier is not None:
                    funding["funder"]["identifier"] = funder_identifier.text
                if funder_name is not None:
                    funding["funder"]["name"] = funder_name.text
            if funding_name is not None:
                funding["name"] = funding_name.text
            if funding_identifier is not None:
                funding["identifier"] = funding_identifier.text
            if funding:
                funding["@type"] = "MonetaryGrant"
                fundings.append(funding)
        if fundings:
            output["funding"] = fundings

        """ TODO try to get the codeRepository field. Not confident how to get it.
        "if this resource is a 'software' @type and has a relatedIdentifier with the isSupplementTo field".
        Examples:
        * 	https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:5680920&metadataPrefix=oai_datacite
                <relatedIdentifier relatedIdentifierType="URL" relationType="IsDerivedFrom" >https://gitlab.com/astron-idg/idg-fpga/</relatedIdentifier>
        * 	https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:1135290&metadataPrefix=oai_datacite
                <relatedIdentifier relatedIdentifierType="URL" relationType="IsSupplementTo" >https://github.com/ljcohen/planets/tree/v0.1</relatedIdentifier>
        """
        return output, missing_types

    def parse(self, records):
        """Transforms/pipeline data to the nde schema before writing the information into the ndson file"""
        # dictionary to convert general type to @type
        # https://docs.google.com/spreadsheets/d/1DOwMjvFL3CGPkdoaFCveKNb_pPVpboEgUjRTT_3eAW0/edit#gid=0
        gen_type = {
            "annotationcollection": "Collection",
            "article": "ScholarlyArticle",
            "audiovisual": "MediaObject",
            "book": "Book",
            "bookchapter": "Chapter",
            "collection": "Collection",
            "conferencepaper": "ScholarlyArticle",
            "datamanagementplan": "CreativeWork",
            "dataset": "Dataset",
            "deliverable": "CreativeWork",
            "diagram": "ImageObject",
            "drawing": "Drawing",
            "figure": "ImageObject",
            "image": "ImageObject",
            "interactiveresource": "CreativeWork",
            "journalarticle": "ScholarlyArticle",
            "lesson": "LearningResource",
            "other": "CreativeWork",
            "outputmanagementplan": "CreativeWork",
            "patent": "CreativeWork",
            "peerreview": "Review",
            "photo": "Photograph",
            "physicalobject": "Thing",
            "plot": "ImageObject",
            "poster": "Poster",
            "preprint": "ScholarlyArticle",
            "presentation": "PresentationDigitalDocument",
            "projectdeliverable": "CreativeWork",
            "projectmilestone": "CreativeWork",
            "proposal": "CreativeWork",
            "publication": "ScholarlyArticle",
            "report": "Report",
            "section": "CreativeWork",
            "software": "ComputationalTool",
            "softwaredocumentation": "TechArticle",
            "taxonomictreatment": "ScholarlyArticle",
            "technicalnote": "TechArticle",
            "thesis": "ScholarlyArticle",
            "video": "VideoObject",
            "workflow": "CreativeWork",
            "workingpaper": "ScholarlyArticle",
        }
        # dictionary log the types that cannot be converted
        missing_types = {}

        for record in records:
            data = json.loads(record[1])
            identifier = record[0]

            # format the identifier for _id, and identifier
            identifier = identifier.rsplit(":", 1)[-1]
            url = "https://zenodo.org/records/" + identifier
            distribution_base_url = "https://zenodo.org/api/records/"
            date_modified = (
                datetime.datetime.fromisoformat(data["header"]["datestamp"][:-1])
                .astimezone(datetime.timezone.utc)
                .date()
                .isoformat()
            )
            # use as much of the metadata variable or header variable to format the transformation
            output = {
                "@context": "https://schema.org/",
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "Zenodo",
                    "url": "https://zenodo.org/",
                    "versionDate": datetime.date.today().isoformat(),
                    "archivedAt": url,
                },
                "_id": "ZENODO_" + identifier,
                "name": data["metadata"].get("title")[0],
                "author": [],
                "identifier": "zenodo." + identifier,
                "dateModified": date_modified,
                "url": url,
                "distribution": [
                    {"contentUrl": distribution_base_url + identifier + "/files-archive", "dateModified": date_modified}
                ],
            }

            if description := data["metadata"].get("description"):
                output["description"] = description[0]

            if date_published := data["metadata"].get("date"):
                try:
                    output["datePublished"] = dateutil.parser.parse(date_published[0], ignoretz=True).isoformat()
                except (dateutil.parser._parser.ParserError, ValueError):
                    try:
                        logger.info("Could not parse date: %s" % date_published[0])
                        output["datePublished"] = dateutil.parser.parse(
                            date_published[0].split("/")[0].replace(" ", "").replace("_", "-"), ignoretz=True
                        ).isoformat()
                    except Exception:
                        # TODO REMOVE THIS LINE VERY SPECIFIC EXCEPTION WHERE USER PUT WRONG DATE
                        if date_published[0] == "2024006-24":
                            output["datePublished"] = dateutil.parser.parse("2024-06-24", ignoretz=True).isoformat()

            if language := data["metadata"].get("language"):
                output["inLanguage"] = {"name": language[0]}

            # zenodo uses different delimiters for keywords. See examples: oai:zenodo.org:1188946, oai:zenodo.org:1204780
            # there may be more delimiters than just '; ' and ', '
            # since zenodo does process the delimiters, for now we will not either example: oai:zenodo.org:1204780
            if keywords := data["metadata"].get("subject"):
                output["keywords"] = []
                for keyword in keywords:
                    # keyword = re.split('; |, ', keyword)
                    # for k in keyword:
                    output["keywords"].append(keyword)

            if sdpublishers := data["metadata"].get("relatedIdentifier"):
                sdps = []
                # some values may be null: https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:6670&metadataPrefix=oai_datacite
                for sdpublisher in sdpublishers:
                    if sdpublisher and "https://zenodo.org/communities/" in sdpublisher:
                        sdp = {"name": sdpublisher.rsplit("/", 1)[-1], "url": sdpublisher}
                        sdps.append(sdp)
                if sdps:
                    output["sdPublisher"] = sdps

            output, missing_types = self.parse_xml(data["xml"], output, gen_type, missing_types, url, identifier)

            # every doc has to have a type
            if output.get("@type"):
                yield output

        # output the missing transformations for @type
        if len(missing_types.keys()) > 0:
            logger.warning("Missing type transformation: {}".format(str(missing_types.keys())))
