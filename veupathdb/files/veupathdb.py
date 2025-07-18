import datetime
import logging

import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")


def record_generator():
    # API call that returns a list of all available data records
    api_command = 'https://veupathdb.org/veupathdb/service/record-types/dataset/searches/AllDatasets/reports/standard?reportConfig={"attributes":["primary_key","organism_prefix","project_id","eupath_release","newcategory","summary","contact","wdk_weight","version","institution","build_number_introduced","pmids_download","release_policy","short_attribution","type","genecount"],"tables":["Publications","Contacts","GenomeHistory","DatasetHistory","Version","References","HyperLinks","GeneTypeCounts","TranscriptTypeCounts"],"attributeFormat":"text"}'
    # send and retrieve request call
    request = requests.get(api_command)
    json_records = request.json()
    logging.info("[INFO] processing %s records...." % len(json_records["records"]))
    # paginate through records
    for _record_dict in json_records["records"]:
        # add custom values to the record
        _record_dict.update(
            {
                "_id": _record_dict["id"][0]["value"],
                "@type": "Dataset",
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "VEuPathDB",
                    "url": "https://veupathdb.org/veupathdb/app",
                    "versionDate": datetime.date.today().isoformat(),
                    "archivedAt": "https://veupathdb.org/veupathdb/app/record/dataset/"
                    + _record_dict["id"][0]["value"],
                },
                "url": "https://veupathdb.org/veupathdb/app/record/dataset/" + _record_dict["id"][0]["value"],
            }
        )

        # set name to records.displayName
        _record_dict["name"] = _record_dict["displayName"]
        # set identifier to records.id.value
        _record_dict["identifier"] = _record_dict["id"][0]["value"]

        # get pmid for, set as string for helper function
        pmids_list = [_dict.pop("pmid") for _dict in _record_dict["tables"]["Publications"]]
        if pmids_list:
            if len(pmids_list) == 1:
                _record_dict["pmids"] = pmids_list[0]
            else:
                _record_dict["pmids"] = ",".join(pmids_list)

        # attributes
        _record_dict["description"] = _record_dict["attributes"].pop("summary")
        _record_dict["measurementTechnique"] = {"name": _record_dict["attributes"].pop("type")}
        _record_dict["sdPublisher"] = {"name": _record_dict["attributes"].pop("project_id")}
        _record_dict["creditText"] = _record_dict["attributes"].pop("short_attribution")

        if _record_dict["attributes"]["release_policy"]:
            _record_dict["conditionOfAccess"] = _record_dict["attributes"].pop("release_policy")

        # tablexs.Contacts
        _record_dict["author"] = [
            {"name": _dict.pop("contact_name"), "affiliation": {"name": str(_dict.pop("affiliation"))}}
            for _dict in _record_dict["tables"]["Contacts"]
        ]

        if _record_dict["attributes"]["version"]:
            try:
                version_date = _record_dict["attributes"]["version"]
                date_updated = datetime.datetime.strptime(version_date, "%Y-%m-%d").date().isoformat()
            except Exception:
                date_updated = None
                logging.debug(
                    "[INFO] BAD DATE FROM _record_dict['attributes']['version']: %s"
                    % _record_dict["attributes"]["version"]
                )

        # tables.GenomeHistory
        release_dates = [hit["release_date"] for hit in _record_dict["tables"]["GenomeHistory"]]
        # if multiple dates passed, keep the most recent date
        if release_dates:
            try:
                if date_updated:
                    iso_list = [
                        datetime.datetime.strptime(d, "%Y-%m-%d").date().isoformat() for d in release_dates
                    ].append(date_updated)
                else:
                    iso_list = [datetime.datetime.strptime(d, "%Y-%m-%d").date().isoformat() for d in release_dates]
                date_updated = sorted(iso_list)[-1]
            except Exception:
                logging.debug("[INFO] BAD DATE FROM _record_dict['tables']['Version']: %s" % release_dates)

        if date_updated:
            _record_dict["dateModified"] = date_updated

        # tables.Version
        published_dates = [hit["version"] for hit in _record_dict["tables"]["Version"]]
        if published_dates:
            try:
                _iso_list = [datetime.datetime.strptime(d, "%Y-%m-%d").date().isoformat() for d in published_dates]
                date_published = sorted(_iso_list)[0]
                _record_dict["datePublished"] = date_published
            except Exception:
                logging.debug("[INFO] BAD DATE FROM _record_dict['tables']['Version']: %s" % published_dates)

        if _record_dict["tables"]["Version"]:
            _record_dict["species"] = [{"name": hit["organism"]} for hit in _record_dict["tables"]["Version"]]

        # tables.HyperLinks
        if _record_dict["tables"]["HyperLinks"]:
            _record_dict["distribution"] = [
                {"name": hit["text"], "url": hit["url"]} for hit in _record_dict["tables"]["HyperLinks"]
            ]

        # table.GeneTypeCounts
        # TODO DO WE NEED THIS?
        # gene_counts = [hit["gene_count"] for hit in _record_dict["tables"]["GeneTypeCounts"]]
        gene_refs = [hit["gene_type"] for hit in _record_dict["tables"]["GeneTypeCounts"]]
        if gene_refs:
            _record_dict["variableMeasured"] = {"name": gene_refs[0]}

        # set conditionsOfAccess and isAccessibleForFree
        _record_dict["conditionsOfAccess"] = "Closed"
        _record_dict["isAccessibleForFree"] = True

        # remove
        _record_dict.pop("recordClassName")
        _record_dict.pop("tableErrors")
        _record_dict.pop("displayName")
        _record_dict.pop("id")
        _record_dict.pop("attributes")
        _record_dict.pop("tables")

        yield _record_dict
