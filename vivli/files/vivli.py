import datetime
import logging
import time

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

# for debugging
# import http.client
# http.client.HTTPConnection.debuglevel = 1

# # You must initialize logging, otherwise you'll not see debug output.
# logging.basicConfig()
# logging.getLogger().setLevel(logging.DEBUG)
# requests_log = logging.getLogger("requests.packages.urllib3")
# requests_log.setLevel(logging.DEBUG)
# requests_log.propagate = True

total = 0


def get_ids():
    # we need to create a session to keep-alive the http connection
    # vivli throttles requests (76 seconds!) when starting a new HTTPs connection
    # https://stackoverflow.com/questions/45783655/first-https-request-takes-much-more-time-than-the-rest
    # https://stackoverflow.com/questions/10115126/python-requests-close-http-connection

    with requests.Session() as session:
        url = (
            "https://vivli-prod-cus-srch.search.windows.net/indexes/studies/docs?api-version=2016-09-01&"
            "api-key=C8237BFE70B9CC48489DC7DD84D88379&search=*&$orderby=nctId%20desc,%20sponsorProtocolId%20desc&$count=true"
        )

        request = session.get(url).json()
        global total
        total = request.get("@odata.count")

        while True:
            # for debugging
            # start = time.time()
            # request = session.get(url).json()
            # end = time.time()
            # print('Time: {}. url: {}'.format(end - start, url))

            request = session.get(url).json()
            time.sleep(0.5)
            for value in request["value"]:
                yield value["id"]
            if url := request.get("@odata.nextLink"):
                continue
            else:
                break


def parse():
    # mapping https://docs.google.com/spreadsheets/d/1c7Iv2PrZPfs_QXIu5IWQfLKp0hc0HGK0Kl5Vqd9Dcpw/edit#gid=0
    ids = get_ids()
    with requests.Session() as session:
        logger.info("Getting IDS...")
        count = 0
        for _id in ids:
            url = "https://prod-api.vivli.org/api/studies/" + _id + "/metadata"
            request = session.get(url).json()
            # pprint(request)
            time.sleep(0.5)

            count += 1
            if count % 100 == 0:
                logger.info("Number of Records Parsed: %s", count)

            output = {
                "includedInDataCatalog": {
                    "@type": "Dataset",
                    "name": "Vivli",
                    "url": "https://vivli.org/",
                    "versionDate": datetime.date.today().isoformat(),
                },
                "@context": "http://schema.org/",
                "@type": "Dataset",
                "identifier": [],
                "_id": "VIVLI_" + request.get("id"),
            }

            if identifier := request.get("nctId"):
                output["identifier"].append(identifier)
            if identifiers := request.get("secondaryIds"):
                output["identifier"] += identifiers

            if sdPublishers := request["registryInfo"]:
                output["sdPublisher"] = []
                for sdPublisher in sdPublishers:
                    sd = {}
                    if name := sdPublisher["registryName"]:
                        sd["name"] = name
                    if identifier := sdPublisher["registryId"]:
                        sd["identifier"] = identifier
                    output["sdPublisher"].append(sd)

            # DOES NOT WORK
            # if author := request.get("principalInvestigator"):
            #     au = {"author": {}}
            #     if given_name := author.get("firstName"):
            #         au["author"]["givenName"] = given_name
            #     if family_name := author.get("lastName"):
            #         au["author"]["familyName"] = family_name
            #     if identifier := author.get("orcidId"):
            #         au["author"]["identifier"] = identifier
            #     if au["author"]:
            #         output["author"] = au

            if author := request.get("orgName"):
                output["author"] = {"name": author}

            if name := request.get("studyTitle"):
                output["name"] = "Dataset from " + name

            funding = []
            if funder := request.get("leadSponsor"):
                if name := funder.get("agency"):
                    funder_name = {"funder": {"name": name}}
                    funding.append(funder_name)

            if funders := request.get("collaborators"):
                for funder in funders:
                    if name := funder.get("agency"):
                        funder_name = {"funder": {"name": name}}
                        funding.append(funder_name)
            if funding:
                output["funding"] = funding

            temporal_coverage = {"temporalInterval": {}}
            if start_date := request.get("studyStartDate"):
                temporal_coverage["temporalInterval"]["startDate"] = (
                    datetime.datetime.fromisoformat(start_date).date().isoformat()
                )

            if end_date := request.get("actualStudyCompletionDate"):
                temporal_coverage["temporalInterval"]["endDate"] = (
                    datetime.datetime.fromisoformat(end_date).date().isoformat()
                )

            if temporal_coverage["temporalInterval"]:
                output["temporalCoverage"] = temporal_coverage

            if locations := request.get("locationsOfStudySites"):
                output["spatialCoverage"] = []
                for location in locations:
                    sc = {}
                    if name := location.get("name"):
                        sc["name"] = name
                    if identifier := location.get("code"):
                        sc["identifier"] = identifier
                    if sc:
                        output["spatialCoverage"].append(sc)

            keywords = []
            if keyword := request.get("phase"):
                keywords.append(keyword)

            if keyword := request.get("studyType"):
                keywords.append(keyword)

            # This should be equivalent to intervention names
            # if arms := request.get('arms'):
            #     for arm in arms:
            #         if interventions := arm.get('interventions'):
            #             for intervention in interventions:
            #                 if keyword := intervention.get('interventionName'):
            #                     keywords.append(keyword)
            if intervention_names := request.get("interventionNames"):
                keywords = keywords + intervention_names

            if pop_items := request.get("populationVocabularyItems"):
                for pop_item in pop_items:
                    if keyword := pop_item.get("term"):
                        keywords.append(keyword)

            if int_items := request.get("interventionVocabularyItems"):
                for int_item in int_items:
                    if keyword := int_item.get("term"):
                        keywords.append(keyword)

            if keywords:
                output["keywords"] = list(set(keywords))

            if conditions := request.get("conditions"):
                output["healthCondition"] = []
                for condition in conditions:
                    output["healthCondition"].append({"name": condition})

            vm = []
            if variable_measured := request.get("outcomeNames"):
                vm.append({"name": variable_measured})

            if outcomes := request.get("outcomes"):
                for outcome in outcomes:
                    if variable_measured := outcome.get("specificMeasurement"):
                        vm.append({"name": variable_measured})

            if vm:
                output["variableMeasured"] = vm

            if doi := request.get("digitalObjectId"):
                output["doi"] = doi
            elif doi := request.get("studyMetadataDoi"):
                output["doi"] = doi
            else:
                pass

            if url := output.get("doi"):
                output["url"] = doi

            if description := request.get("description"):
                output["description"] = description
            elif description := request.get("extractedBriefSummary"):
                output["description"] = description
            else:
                pass

            if date_created := request.get("draftCreatedDate"):
                output["dateCreated"] = datetime.datetime.fromisoformat(date_created.split(".")[0]).date().isoformat()
            elif date_created := request.get("submittedDate"):
                output["dateCreated"] = datetime.datetime.fromisoformat(date_created.split(".")[0]).date().isoformat()
            else:
                pass

            if date_published := request.get("postedDate"):
                output["datePublished"] = (
                    datetime.datetime.fromisoformat(date_published.split(".")[0]).date().isoformat()
                )

            if grant := request.get("grant"):
                logger.info("Grant Exists: %s. ID is: %s", grant, output["_id"])

            if funding := request.get("funder"):
                logger.info("Funding Exists: %s. ID is: %s", funding, output["_id"])

            if date_modified := request.get("updatedDate"):
                output["dateModified"] = datetime.datetime.fromisoformat(date_modified.split(".")[0]).date().isoformat()

            yield output

        if count == total:
            logger.info("Total number of documents parsed: %s", count)
        else:
            logger.warning(
                "Did not parse all the records \n" "Total number parsed: %s \n Total number of documents: %s",
                count,
                total,
            )
