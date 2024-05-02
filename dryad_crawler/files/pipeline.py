# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface

import datetime
import logging
import re

logger = logging.getLogger("nde-logger")


__all__ = [
    "DryadItemProcessorPipeline",
]


class DryadItemProcessorPipeline:
    def format_date(self, date_string):
        """Formats python and military datetime into an isoformat date

        Args:
            date_string: a date string
        Returns: An isoformatted date if there is a datetime if not then return None
        """
        if re.match(r"\d+-\d+-\d+", date_string):
            if "T" in date_string:
                date = datetime.datetime.fromisoformat(date_string.split("T")[0]).date().isoformat()
            else:
                date = datetime.datetime.fromisoformat(date_string).date().isoformat()
            return date
        return None

    def process_item(self, item: dict, spider):
        last_part = item["@id"].split("/")[-1]
        if "dryad." in last_part:
            last_part = last_part.replace("dryad.", "")
        output = {
            "@context": item.pop("@context"),
            "@type": item.pop("@type"),
            "url": item.pop("url"),
            "_id": "DRYAD_" + last_part,
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "Dryad Digital Repository",
                "url": "https://datadryad.org",
                "versionDate": datetime.date.today().isoformat(),
            },
        }

        if name := item.pop("name", None):
            output["name"] = name

        if description := item.pop("abstract", None):
            if methods := item.pop("methods", None):
                description = description + "\nMethods\n" + methods
            output["description"] = description

        if has_part := item.pop("usageNotes", None):
            output["hasPart"] = {"name": has_part}

        if identifier := item.pop("identifier", None):
            output["identifier"] = identifier
            if "doi" in identifier:
                output["doi"] = identifier

        if keywords := item.pop("keywords", None):
            output["keywords"] = keywords

        # we want to convert author.sameAs to author.url. If orcid url then convert to author.identifier
        if authors := item.pop("creator", None):
            author_list = []
            # checks if there is a list of authors
            if not isinstance(authors, list):
                if same_as := authors.pop("sameAs", None):
                    if "orcid" in same_as:
                        authors["identifier"] = same_as
                    else:
                        authors["url"] = same_as
                output["author"] = authors
            else:
                for author in authors:
                    if same_as := author.pop("sameAs", None):
                        if "orcid" in same_as:
                            author["identifier"] = same_as
                        else:
                            author["url"] = same_as
                    author_list.append(author)
                output["author"] = author_list

        # change distribution
        if distribution := item.pop("distribution", None):
            if distribution.get("hasPart"):
                raise AssertionError("Distribution contains hasPart")
            elif content_url := item.pop("contentUrl", None):
                distribution["hasPart"] = content_url
            else:
                pass
            output["distribution"] = distribution

        # There are 2 different cases and 4 different time format
        """
            Case 1: No list (contains one of the formats)
            Case 2: List (we assume it contains format 1 and  at least one format 2-4)
            Format 1: Y (2019)
            Format 2: Y-M-D (2018-12-13)
            Format 3: Python datetime (2015-10-16T00:00:00+00:00)
            Format 4: Military datetime (2015-10-16T00:00:00Z)
        """
        if temp_cov := item.pop("temporalCoverage", None):
            if isinstance(temp_cov, list):
                for tc in temp_cov:
                    # we assume there will be a format 2-4 in a list
                    if date := self.format_date(tc):
                        output["temporalCoverage"] = {"temporalInterval": {"endDate": date}}
            else:
                # we cannot assume format 2-4 is in the no list
                if date := self.format_date(temp_cov):
                    pass
                else:
                    date = datetime.datetime.strptime(temp_cov, "%Y").date().isoformat()
                output["temporalCoverage"] = {"temporalInterval": {"endDate": date}}

            if not output.get("temporalCoverage"):
                logger.info("No temporal coverage: " + " ".join(temp_cov))
                logger.info("ID is %s", output["_id"])

        if spatial_covs := item.pop("spatialCoverage", None):
            if isinstance(spatial_covs, list):
                sc = []
                for spatial_cov in spatial_covs:
                    if isinstance(spatial_cov, str):
                        sc.append({"name": spatial_cov})
                output["spatialCoverage"] = sc
            else:
                if isinstance(spatial_covs, str):
                    output["spatialCoverage"] = {"name": spatial_covs}

        if citation := item.pop("citation", None):
            output["citation"] = {"url": citation}

        if license_obj := item.pop("license", None):
            output["license"] = license_obj

        if date_pub := item.pop("publicationDate", None):
            date_pub = datetime.datetime.fromisoformat(date_pub).date().isoformat()
            output["datePublished"] = date_pub

        if is_accessible := item.pop("isAccessibleForFree", None):
            output["isAccessibleForFree"] = is_accessible

        if coa := item.pop("visibility", None):
            if coa.casefold() == "public":
                output["conditionsOfAccess"] = "Open"
            else:
                output["conditionsOfAccess"] = coa

        if funding := item.pop("funders", None):
            fds = []
            for funder in funding:
                fd = {}
                if name := funder.get("organization"):
                    fd["funder"] = {"name": name}
                if add_type := funder.get("identifierType"):
                    fd["additionalType"] = add_type
                if url := funder.get("identifier"):
                    fd["url"] = url
                if identifier := funder.get("awardNumber"):
                    fd["identifier"] = identifier
                fds.append(fd)
            output["funding"] = fds

        # if date_pub := item.pop('datePublished'):
        #     date_pub = datetime.datetime.strptime(date_pub.split(': ')[1], '%B %d, %Y').date().isoformat()
        #     output['datePublished'] = date_pub

        # remove unneeded values
        item.pop("provider", None)
        item.pop("publisher", None)
        item.pop("version", None)
        item.pop("links", None)
        item.pop("identifier", None)
        item.pop("id", None)
        item.pop("storageSize", None)
        item.pop("relatedPublicationISSN", None)
        # made using api
        item.pop("description", None)
        # provided in json+ld as name
        item.pop("title", None)
        # provided in json+ld as creators
        item.pop("authors", None)
        item.pop("relatedWorks", None)
        item.pop("versionNumber", None)
        item.pop("versionStatus", None)
        item.pop("curationStatus", None)
        item.pop("versionChanges", None)
        item.pop("lastModificationDate", None)
        item.pop("sharingLink", None)
        item.pop("userId", None)
        item.pop("changedFields", None)
        item.pop("license", None)
        # already in spatialCoverage
        item.pop("locations", None)
        item.pop("_links", None)
        item.pop("fieldOfScience", None)
        callback_url = item.pop("callback_url", None)

        if item:
            logger.warning("Haven't parsed all keys in dryad_crawler: " + ", ".join(item.keys()))
            logger.warning("Callback URL: %s", callback_url)

        return output
