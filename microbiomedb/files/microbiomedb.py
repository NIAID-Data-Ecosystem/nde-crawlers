import datetime
import json
import logging
import re

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


def parse():
    # mapping https://docs.google.com/spreadsheets/d/1uz1puYgLiiud2zuncXv43S1BkQvDNAeyl8pRiXgT14Y/edit?usp=sharing

    url = (
        "https://microbiomedb.org/mbio/service/record-types/dataset/searches/AllDatasets/reports/standard?"
        'reportConfig={"attributes":[acknowledgement,pmids,study_categories,description,type,institution,card_headline,'
        "contact,project_availability,study_access,email,is_species_scope,summary,card_points,dataset_id,build_number_introduced,"
        "short_display_name,display_name,dataset_name,primary_key,version,pmids_download,bulk_download_url,eda_study_id,"
        'short_attribution,is_public,eupath_release,caveat],"tables":[Contacts,DownloadVersion,Publications,'
        'AssociatedDatasets,HyperLinks],"attributeFormat":"text"}'
    )

    logger.info("Making request...")
    request = requests.get(url)
    logger.info("Request made. HTTP STATUS: %d", request.status_code)
    records = request.json()
    logger.info("Parsing records...")
    for count, record in enumerate(records["records"], start=1):
        if count % 1000 == 0:
            logger.info("Parsed %d records", count)

        output = {
            "includedInDataCatalog": {
                "@type": "Dataset",
                "name": "MicrobiomeDB",
                "url": "https://microbiomedb.org/mbio/app",
                "versionDate": datetime.date.today().isoformat(),
            },
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": "MICROBIOME_" + record.get("displayName").replace(" ", "_").replace("(", "-").replace(")", "-"),
        }

        if record_id := record["id"]:
            if len(record_id) > 1:
                logger.warning("More than one record id found. Using first record id.")
            record_id = record_id[0].get("value")
            output["url"] = f"https://microbiomedb.org/mbio/app/record/dataset/{record_id}"
            output["includedInDataCatalog"]["dataset"] = output["url"]

        if record_name := record.get("displayName"):
            output["name"] = record_name

        # set up keywords list
        k = []
        if keywords := record["attributes"].get("study_categories"):
            keywords = json.loads(keywords)
            k = keywords + k

        # type should only be a string
        if keyword := record["attributes"].get("type"):
            k.append(keyword)

        if k:
            output["keywords"] = k

        output_description = []

        if description := record["attributes"].get("description"):
            if isinstance(description, str):
                output_description.append(description)
            else:
                output_description.extend(description)

        if description := record["attributes"].get("card_points"):
            description = json.loads(description)
            if isinstance(description, str):
                output_description.append(description)
            else:
                output_description.extend(description)

        # Join all descriptions into one string
        if output_description:
            output["description"] = "\n".join(output_description)

        # set up contributor field
        contributor = {}
        if contributor_aff := record["attributes"].get("institution"):
            contributor["affiliation"] = contributor_aff

        if contributor_name := record["attributes"].get("name"):
            contributor["name"] = contributor_name

        if contributor_email := record["attributes"].get("email"):
            contributor["email"] = contributor_email

        if contributor:
            output["contributor"] = contributor

        if coa := record["attributes"].get("study_access"):
            if coa.casefold() == "public":
                coa = "Open"
            output["conditionsOfAccess"] = coa

        if abstract := record["attributes"].get("summary"):
            output["abstract"] = abstract

        if alternate_name := record["attributes"].get("dataset_name"):
            output["alternateName"] = alternate_name

        if pmids := record["attributes"].get("pmids_download"):
            output["pmids"] = pmids

        # if distribution_content_url := record['attributes'].get('bulk_url_download'):
        #     output['distribution'] = {'contentUrl': distribution_content_url}

        # set up author field
        if authors := record["tables"].get("Contacts"):
            a_list = []
            for author in authors:
                a = {}
                if name := author.get("contact_name"):
                    a["name"] = name
                if affiliation := author.get("affiliation"):
                    a["affiliation"] = {"name": affiliation}
                if a:
                    a_list.append(a)
            if a_list:
                output["author"] = a_list

        # set up distribution field
        if distributions := record["tables"].get("DownloadVersion"):
            d_list = []
            for distribution in distributions:
                d = {}
                if date := distribution.get("release_date"):
                    date = datetime.datetime.strptime(date, "%Y-%b-%d").date().isoformat()
                    d["datePublished"] = date
                if identifier := distribution.get("dataset_id"):
                    d["identifier"] = identifier
                if name := distribution.get("dataset_name"):
                    d["name"] = name
                if version_number_link := distribution.get("version_number_link"):
                    match = re.search(r'"downloadUrl":\s*"([^"]+)"', version_number_link)
                    if match:
                        contentUrl = match.group(1)
                        if contentUrl:
                            d["contentUrl"] = "https://microbiomedb.org" + contentUrl
                if d:
                    d_list.append(d)
            if d_list:
                output["distribution"] = d_list

        # get both isPartOf and isRelatedTo
        if ads := record["tables"].get("AssociatedDatasets"):
            # isPartOf
            ipo_list = []
            added_ipo = set()
            # isRelatedTo
            irt_list = []
            for ad in ads:
                ipo = {}
                if name := ad.get("project"):
                    ipo["name"] = name
                if alternate_name := ad.get("webapp_name"):
                    ipo["alternateName"] = alternate_name

                ipo_key = (ipo.get("name"), ipo.get("alternateName"))
                if ipo and ipo_key not in added_ipo:
                    ipo_list.append(ipo)
                    added_ipo.add(ipo_key)

                irt = {}
                if cite_url := ad.get("pubmed link"):
                    irt["citation"] = {"url": cite_url}

                if (cite_pmid := ad.get("pmid")) and cite_url:
                    irt["citation"]["pmid"] = cite_pmid
                elif cite_pmid := ad.get("pmid"):
                    irt["citation"] = {"pmid": cite_pmid}

                if url := ad.get("associated dataset"):
                    irt["url"] = url
                if name := ad.get("display_name"):
                    irt["name"] = name
                if identifier := ad.get("associated_dataset_id"):
                    irt["identifier"] = identifier
                if irt:
                    irt_list.append(irt)
            if ipo_list:
                output["isPartOf"] = ipo_list
            if irt_list:
                output["isRelatedTo"] = irt_list

        # grabs urls
        # Traditionally we've been using url as the url where the source contains the dataset.
        # Should only be one url
        if links := record["tables"].get("HyperLinks"):
            link = links[0]
            sd_publisher = {}
            if url := link.get("url"):
                sd_publisher["url"] = url
            if dataset_id := link.get("dataset_id"):
                sd_publisher["identifier"] = dataset_id
            if sd_publisher:
                sd_publisher["name"] = "ClinEpiDB"
                output["sdPublisher"] = sd_publisher
            if len(links) > 1:
                logger.warning("More than one link found. Using first link.")

        yield output

    logger.info("Finished Parsing. Total records: %d", count)
