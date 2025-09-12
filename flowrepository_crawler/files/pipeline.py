# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface
# mapping https://docs.google.com/spreadsheets/d/19hQf4sQ6ZLwvt8ADYyvooYrqmQCsxh386kYnycX3xmg/edit#gid=0

import datetime
import logging
import re

logger = logging.getLogger("nde-logger")


__all__ = [
    "FlowRepositoryItemProcessorPipeline",
]


class FlowRepositoryItemProcessorPipeline:
    def match(self, authors, author_name):
        match = False
        for name in authors:
            if name["name"] == author_name:
                match = True
        if not match:
            author = {"name": author_name}
            authors.append(author)

        return authors

    def process_item(self, item: dict, spider):
        # pop url first as it is the only key that doesnt contain a list
        url = item.pop("url")
        # convert keys that are not organizations or manuscripts
        for key, value in item.items():
            if key != "Organizations:" and key != "Manuscripts:":
                item[key] = " ".join(value)
            else:
                pass

        identifier = item.pop("Repository ID:")
        # logger.info(f"Url: {url}, ID: {identifier}")
        assert isinstance(item.get("Organizations:"), list)
        assert isinstance(item.get("Manuscripts:"), list)

        output = {
            "@type": "Dataset",
            "url": url,
            "_id": "FLOW_" + identifier,
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "Flow Repository",
                "url": "http://flowrepository.org/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
            "measurementTechnique": {
                "name": "flow cytometry method",
                "url": "http://purl.obolibrary.org/obo/MMO_0000617",
                "inDefinedTermSet": "MMO",
                "curatedBy": {"name": "Flow Repository", "url": "http://flowrepository.org/", "dateModified": datetime.date.today().isoformat()},
                "isCurated": True
            },
            "conditionsOfAccess": "Open",
            "license": "http://flowrepository.org/terms_of_service"

        }
        """ These are the keys:
            ['Repository ID:', 'Experiment name:', 'MIFlowCyt score:', 'Primary researcher:', 'PI/manager:',
            'Uploaded by:', 'Experiment dates:', 'Dataset uploaded:', 'Last updated:', 'Keywords:', 'Manuscripts:',
            'Organizations:', 'Purpose:', 'Conclusion:', 'Comments:', 'Funding:', 'Quality control:']
            Organizations and manuscripts can be more than one.
        """

        output["identifier"] = identifier
        if name := item.pop("Experiment name:", None):
            output["name"] = name

        # get the authors
        authors = []

        if author_name := item.pop("Primary researcher:", None):
            author = {"name": author_name}
            # if there is a Primary research then attach aff names
            if aff_names := item.pop("Organizations:", None):
                affiliations = []
                for aff_name in aff_names:
                    affiliation = {"name": aff_name}
                    affiliations.append(affiliation)
                author["affiliation"] = affiliations
            authors.append(author)
        else:
            # No primary researcher so we use author.name
            if author_names := item.pop("Organizations", None):
                for author_name in author_names:
                    author = {"name": author_name}
                    authors.append(author)

        if author_name := item.pop("PI/manager:", None):
            authors = self.match(authors, author_name)

        if author_name := item.pop("Uploaded by:", None):
            authors = self.match(authors, author_name)

        if authors:
            output["author"] = authors

        if temp_dates := item.pop("Experiment dates:", None):
            # All dates should be in this date pattern
            if dates := re.findall(r"\d+-\d+-\d+", temp_dates):
                dates = [datetime.datetime.fromisoformat(date_string).date().isoformat() for date_string in dates]
            else:
                logger.error("Could not parse dates. URL: %s", output["url"])
            # First date entry could be start or end date
            date = datetime.datetime.strptime(dates[0], "%Y-%m-%d").date().isoformat()
            assert len(dates) <= 2, f"There are more than 2 dates: {dates}. URL: {output['url']}"
            # check if there is more than one date
            if len(dates) > 1:
                end_date = datetime.datetime.strptime(dates[1], "%Y-%m-%d").date().isoformat()
                output["temporalCoverage"] = {"@type": "TemporalInterval", "temporalType": "study date", "startDate": date, "endDate": end_date}
            elif re.findall(r"\d+-\d+-\d+.*-", temp_dates):
                output["temporalCoverage"] = {"@type": "TemporalInterval", "temporalType": "study date", "startDate": date}
            else:
                output["temporalCoverage"] = {"@type": "TemporalInterval", "temporalType": "study date", "endDate": date}
                # I have not found an example of an end date yet.
                logger.info("There is only end date. URL: %s", output["url"])

        if date_pub := item.pop("Dataset uploaded:", None):
            date_iso = ""
            # This is to cover weird exception cases like http://flowrepository.org/id/FR-FCM-Z645.
            # Long year length
            if len(date_pub.split()[1]) == 4:
                date_iso = datetime.datetime.strptime(date_pub, "%b %Y").date().isoformat()
                output["datePublished"] = date_iso
            # Short Year length
            elif len(date_pub.split()[1]) == 2:
                date_iso = datetime.datetime.strptime(date_pub, "%b %y").date().isoformat()
            else:
                logger.warning("Date Uploaded has not been parsed. URL: %s", output["url"])
            if date_iso:
                output["datePublished"] = date_iso

        if date_mod := item.pop("Last updated:", None):
            date_iso = ""
            # This is to cover weird exception cases. Does not cover this execption case http://flowrepository.org/id/FR-FCM-Z3L9 "11:45AM"
            # Long year length
            if len(date_mod.split()[1]) == 4:
                date_iso = datetime.datetime.strptime(date_mod, "%b %Y").date().isoformat()
                output["dateModified"] = date_iso
            # Short Year length
            elif len(date_pub.split()[1]) == 2:
                date_iso = datetime.datetime.strptime(date_pub, "%b %y").date().isoformat()
            else:
                logger.warning("Date Uploaded has not been parsed. URL: %s", output["url"])
            if date_iso:
                output["dateModified"] = date_iso

            output["distribution"] = {
                "@type": "DataDownload",
                "name": "FCS files",
                "contentUrl": "http://flowrepository.org/experiments/510/download_ziped_files",
            }
            if output.get("dateModified"):
                output["distribution"]["dateModified"] = output["dateModified"]

        if keywords := item.pop("Keywords:", None):
            keywords = keywords.strip("[]").split("] [")
            keywords = [keyword for keyword in keywords if keyword.casefold() != "none"]
            if keywords:  # Only add to output if there are valid keywords
                output["keywords"] = keywords

        # split manuscripts into two comma separated strings: pmids and pmcs
        if manuscripts := item.pop("Manuscripts:", None):
            manuscripts = [manuscript.strip("[]") for manuscript in manuscripts]
            pmcs = []
            pmids = []
            for manuscript in manuscripts:
                if "PMC" in manuscript:
                    pmcs.append(manuscript)
                else:
                    pmids.append(manuscript)

            pmcs = ",".join([*set(pmcs)])
            pmids = ",".join([*set(pmids)])

            if pmcs:
                output["pmcs"] = pmcs
            if pmids:
                output["pmids"] = pmids

        # format the description
        description = ""
        # Do not include these in the description
        none_list = ["none", "n/a", "see linked paper", "not disclosed"]
        if purpose := item.pop("Purpose:", None):
            if purpose.casefold() not in none_list:
                description = description + purpose

        if conclusion := item.pop("Conclusion:", None):
            if conclusion.casefold() not in none_list:
                description = description + "\nConclusion: \n" + conclusion

        if notes := item.pop("Comments:", None):
            if notes.casefold() not in none_list:
                description = description + "\n Notes: \n" + notes

        if qc := item.pop("Quality control:", None):
            if qc.casefold() not in none_list:
                description = description + " " + qc

        if description:
            output["description"] = description

        # TODO Funding is under one string and there's no clear delimiter. We will keep it as if its only one funder for now.
        if funding := item.pop("Funding:", None):
            if funding.casefold() not in none_list:
                output["funding"] = {"description": funding}

        # pop unneeded values
        item.pop("MIFlowCyt score:", None)

        if item:
            logger.warning("Haven't parsed all keys in flowrepo_crawler: " + ", ".join(item.keys()))
            logger.warning("Callback URL: %s", output.get("url"))

        return output
