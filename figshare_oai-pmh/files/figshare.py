import datetime
import json

# import time
import logging
import re

# from xml.etree import ElementTree
from oai_helper import oai_helper

# from sickle import Sickle
from sql_database import NDEDatabase

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


class Figshare(NDEDatabase):
    SQL_DB = "figshare_filter.db"
    EXPIRE = datetime.timedelta(days=90)

    # sickle = Sickle('https://api.figshare.com/v2/oai',
    #                 max_retries=10, default_retry_after=20)

    def load_cache(self):
        # """Retrives the raw data using a sickle request and formats so dump can store it into the cache table
        # Returns:
        #     A tuple (_id, data formatted as a json string)
        # """
        # records = self.sickle.ListRecords(
        #     metadataPrefix='uketd_dc', ignore_deleted=True)

        # while True:
        #     try:
        #         # get the next item
        #         record = records.next()
        #         count += 1
        #         if count % 100 == 0:
        #             time.sleep(1)
        #             logger.info("Loading cache. Loaded %s records", count)

        #         if count % 1000 == 0:
        #             raise StopIteration

        #         # in each doc we want record.identifier and record stored
        #         doc = {'header': dict(record.header), 'metadata': record.metadata,
        #                'xml': ElementTree.tostring(record.xml, encoding='unicode')}

        #         yield (record.header.identifier, json.dumps(doc))

        #     except StopIteration:
        #         logger.info("Finished Loading. Total Records: %s", count)
        #         # if StopIteration is raised, break from loop
        #         break

        records = oai_helper()
        count = 0
        for record in records:
            count += 1
            yield record
        logger.info("Finished Loading. Total Records: %s", count)

    def parse(self, records):
        # used to test single record
        # record = sickle.GetRecord(
        #     identifier='oai:figshare.com:article/5849037', metadataPrefix='uketd_dc')
        # TODO parse time
        count = 0
        for record in records:
            data = json.loads(record[1])
            metadata = data["metadata"]

            # testing for missing properties
            # for key in metadata:
            #     if key not in properties:
            #         missing.append(key)

            count += 1
            if count % 1000 == 0:
                logger.info("Parsed %s records", count)
                # logging missing properties
                # if len(missing):
                #     logger.info(f'Missing {missing}')

            output = {
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "Figshare",
                    "versionDate": datetime.date.today().isoformat(),
                    "url": "https://figshare.com",
                },
            }
            if title := metadata.get("title"):
                output["name"] = title[0]

            if creators := metadata.get("creator"):
                creator_list = []
                for creator in creators:
                    creator_list.append({"name": creator})
                output["author"] = creator_list

            # TODO
            # checking if covid related article for outbreak api
            #     for keyword in subject:
            #         if keyword.lower() in covid_keywords:
            #             output['outbreakapi'] = True
            #         else:
            #             output['outbreakapi'] = False
            if subject := metadata.get("subject"):
                output["keywords"] = subject

            if description := metadata.get("description"):
                output["description"] = description[0]
            if description == None:
                if abstract := metadata.get("abstract"):
                    output["description"] = abstract

            if date := metadata.get("date"):
                output["dateModified"] = datetime.datetime.strptime(date[0], "%Y-%m-%dT%H:%M:%S%z").strftime("%Y-%m-%d")

            # sdPublsher = publisher > instituion/department > figshare
            if publisher := metadata.get("publisher"):
                output["sdPublisher"] = {"name": publisher[0]}
            if publisher == None:
                institution = metadata.get("institution")
                department = metadata.get("department")
                if institution and department:
                    if institution[0] != department[0]:
                        output["sdPublisher"] = {"name": institution[0] + "/" + department[0]}
                    else:
                        output["sdPublisher"] = {"name": institution[0]}
                else:
                    output["sdPublisher"] = {"name": "figshare"}

            if type := metadata.get("type"):
                if len(type) > 1:
                    if type[0] == type[1]:
                        if type[0] == "Software":
                            output["@type"] = "ComputationalTool"
                        else:
                            output["@type"] = type[0]
                    else:
                        output["@type"] = "Collection"
                        output["hasPart"] = {"@type": type}
                else:
                    output["@type"] = "Collection"
                    output["hasPart"] = {"@type": type[0]}

            if identifier := metadata.get("identifier"):
                for el in identifier:
                    # [None, 'https://ndownloader.figshare.com/files/35101450']
                    if el is not None:
                        if "ndownloader" in el:
                            output["distribution"] = {"url": el}
                        elif "10." in el:
                            output["doi"] = el
            if identifier == None:
                if reference := metadata.get("isReferencedBy"):
                    output["doi"] = reference[0]

            if language := metadata.get("language"):
                output["language"] = language[0]
            if relation := metadata.get("relation"):
                url = relation[0]
                output["url"] = url
                output["includedInDataCatalog"]["archivedAt"] = url
                output["_id"] = "Figshare_" + relation[0].split("/")[-1]
                output["identifier"] = relation[0].split("/")[-1]
            if license := metadata.get("license"):
                output["license"] = license[0]
            if issued := metadata.get("issued"):
                output["datePublished"] = issued[0]

            # Funding extraction:
            # Figshare's uketd_dc sometimes places the actual grant identifier in the <sponsor> tag
            # while <grantnumber> can be missing or a placeholder like "0". We attempt to:
            #  1. Use grantnumber as identifier when it looks valid.
            #  2. Detect NIH-style (or similar) grant identifiers inside sponsor when grantnumber is invalid.
            #  3. Support multiple sponsors/grantnumbers, returning a list if >1.
            #  4. Only include fields that have content.
            if sponsors := metadata.get("sponsor"):
                grantnumbers = metadata.get("grantnumber") or []

                def is_invalid(val: str) -> bool:
                    return val is None or val.strip().lower() in {"0", "na", "n/a", "", "none", "null"}

                # Regex: two letters followed by 6 digits (basic IC code + project number pattern)
                grant_pattern = re.compile(r"\b[A-Za-z]{2}\d{6}\b")

                fundings = []
                for idx, sponsor_val in enumerate(sponsors):
                    sponsor_val = sponsor_val.strip() if sponsor_val else sponsor_val
                    grant_val = grantnumbers[idx] if idx < len(grantnumbers) else None
                    grant_val = grant_val.strip() if isinstance(grant_val, str) else grant_val

                    funding_record = {}

                    # Case 1: valid grantnumber provided
                    if grant_val and not is_invalid(grant_val):
                        # sponsor is treated as funder name
                        if sponsor_val:
                            funding_record["funder"] = {"name": sponsor_val}
                        funding_record["identifier"] = grant_val
                    else:
                        # grantnumber invalid. Check if sponsor itself contains a grant id.
                        identifier_match = None
                        if sponsor_val:
                            # direct entire sponsor value matches the pattern
                            if grant_pattern.fullmatch(sponsor_val):
                                identifier_match = sponsor_val
                            else:
                                # search inside sponsor string
                                inner = grant_pattern.search(sponsor_val)
                                if inner:
                                    identifier_match = inner.group(0)
                        if identifier_match:
                            # Use detected identifier. If sponsor contained extra text, keep as funder name too.
                            funding_record["identifier"] = identifier_match
                            # If sponsor has more than just the identifier, store original as name
                            if sponsor_val and sponsor_val != identifier_match:
                                funding_record["funder"] = {"name": sponsor_val}
                        else:
                            # No identifier; just treat sponsor as funder name
                            if sponsor_val:
                                funding_record["funder"] = {"name": sponsor_val}

                    # Only append if we captured something meaningful
                    if funding_record:
                        fundings.append(funding_record)

                if fundings:
                    # Flatten to single dict if only one funding entry
                    output["funding"] = fundings[0] if len(fundings) == 1 else fundings

            yield output

    def update_cache(self):
        """If cache is not expired get the new records to add to the cache since last_updated"""

        last_updated = self.retreive_last_updated()

        # get all the records since last_updated to add into current cache
        # logger.info("Updating cache from %s", last_updated)
        # records = self.sickle.ListRecords(**{
        #     'metadataPrefix': 'uketd_dc',
        #     'ignore_deleted': True,
        #     'from': last_updated
        # }
        # )
        records = oai_helper(last_updated)
        count = 0
        for record in records:
            count += 1
            yield record
        logger.info("Finished updating cache. Total new records: %s", count)

        # # Very similar to load_cache()
        # count = 0
        # while True:
        #     try:
        #         # get the next item
        #         record = records.next()
        #         count += 1
        #         if count % 100 == 0:
        #             time.sleep(1)
        #             logger.info(
        #                 "Updating cache. %s new updated records", count)

        #         # in each doc we want record.identifier and record stored
        #         doc = {'header': dict(record.header), 'metadata': record.metadata,
        #                'xml': ElementTree.tostring(record.xml, encoding='unicode')}

        #         yield (record.header.identifier, json.dumps(doc))

        #     except StopIteration:
        #         logger.info(
        #             "Finished updating cache. Total new records: %s", count)
        #         # if StopIteration is raised, break from loop
        #         break

        self.insert_last_updated()
