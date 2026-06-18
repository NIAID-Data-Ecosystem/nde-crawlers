"""Standalone Zenodo parser that reads the full XML export tarball.
to query one document: https://zenodo.org/oai2d?verb=GetRecord&identifier=oai:zenodo.org:5680920&metadataPrefix=oai_datacite
to query all documents: https://zenodo.org/oai2d?verb=ListRecords&metadataPrefix=oai_datacite

Instead of harvesting OAI-PMH with sickle, this streams Zenodo's bulk export
(https://zenodo.org/api/exporter/records-xml.tar.gz). Each member is one
<oai_datacite> file whose <payload><resource> block is DataCite kernel-4 XML.
We transform each record into the NDE schema and yield it.

Header fields that the OAI feed would have carried are recovered from the XML:`
  * record id  -> the member filename (e.g. 8435696.xml -> 8435696)
  * dateModified -> DataCite <date dateType="Updated"> (falls back to "Issued")
"""

import datetime
import itertools
import logging
import os
import tarfile
import xml.etree.ElementTree as ET

import dateutil.parser
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

TARBALL_URL = "https://zenodo.org/api/exporter/records-xml.tar.gz"

# DataCite kernel-4 namespace, used as an ElementTree tag prefix.
NS = "{http://datacite.org/schema/kernel-4}"

# schema.org types for (resourceTypeGeneral, resourceType) pairs. Built at runtime from Zenodo's resource type vocabulary.
def _build_type_mapping():
    mapping = {}
    request = requests.get("https://zenodo.org/api/vocabularies/resourcetypes?page=1").json()
    total = request["hits"]["total"]
    for page in range(1, (total // 100) + 2):
        request = requests.get(f"https://zenodo.org/api/vocabularies/resourcetypes?page={page}&size=100").json()
        for hit in request["hits"]["hits"]:
            prop = hit["props"]
            general = (prop["datacite_general"] or "").lower().replace(" ", "")
            specific = (prop["datacite_type"] or "").lower().replace(" ", "")
            mapping[(general, specific)] = prop["schema.org"].split("/")[-1]
    return mapping


def _find_date(root, *date_types):
    """Return the first <date> element matching the given dateTypes, else the
    first <date>. Uses explicit None checks: an Element with no children is
    falsy, so `find(...) or find(...)` would skip a real match."""
    for date_type in date_types:
        el = root.find(f".//{NS}date[@dateType='{date_type}']")
        if el is not None:
            return el
    return root.find(f".//{NS}date")


def _date_modified(root):
    """Return an ISO date string from DataCite <dates>, preferring Updated."""
    el = _find_date(root, "Updated", "Issued")
    if el is None or not el.text:
        return None
    try:
        return dateutil.parser.parse(el.text, ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, ValueError):
        return None


def _date_published(root):
    """Return an ISO datetime string from the Issued date (else the first date)."""
    el = _find_date(root, "Issued")
    if el is None or not el.text:
        return None
    try:
        return dateutil.parser.parse(el.text, ignoretz=True).isoformat()
    except (dateutil.parser.ParserError, ValueError):
        try:
            cleaned = el.text.split("/")[0].replace(" ", "").replace("_", "-")
            return dateutil.parser.parse(cleaned, ignoretz=True).isoformat()
        except Exception:
            logger.info("Could not parse date: %s", el.text)
            return None


def _apply_version(output, root, url):
    """If the record IsVersionOf a Zenodo DOI, point the doc at the concept (original) version."""
    version_ids = []
    for rel in root.findall(f".//{NS}relatedIdentifier"):
        # A valid version identifier carries only relatedIdentifierType and
        # relationType; any extra attribute (e.g. resourceTypeGeneral) marks a
        # malformed entry whose text is not a clean DOI, so skip it.
        if set(rel.attrib) != {"relatedIdentifierType", "relationType"}:
            continue
        if (
            rel.get("relatedIdentifierType") == "DOI"
            and rel.get("relationType") == "IsVersionOf"
            and rel.text
            and "zenodo" in rel.text
        ):
            version_ids.append(rel.text)

    if not version_ids:
        return

    version_ids = list(set(version_ids))
    if len(version_ids) > 1:
        version_ids = [min(version_ids, key=lambda doi: int(doi.split(".")[-1]))]

    version_id = version_ids[0]
    concept_num = version_id.rsplit(".", 1)[-1]
    output["versionId"] = version_id
    output["sameAs"] = [url]
    output["doi"] = version_id
    output["identifier"] = version_id.rsplit("/", 1)[-1]
    output["_id"] = "ZENODO_" + concept_num
    output["url"] = "https://zenodo.org/record/" + concept_num


def build_doc(content, record_id, types, missing_types):
    """Transform one tarball member (raw XML bytes) into an NDE doc, or None.

    Returns None when no @type can be assigned (every doc must have a type).
    """
    root = ET.fromstring(content)

    url = "https://zenodo.org/records/" + record_id
    date_modified = _date_modified(root)

    output = {
        "@context": "https://schema.org/",
        "includedInDataCatalog": {
            "@type": "DataCatalog",
            "name": "Zenodo",
            "url": "https://zenodo.org/",
            "versionDate": datetime.date.today().isoformat(),
            "archivedAt": url,
        },
        "_id": "ZENODO_" + record_id,
        "author": [],
        "identifier": "zenodo." + record_id,
        "url": url,
        "distribution": [
            {
                "contentUrl": "https://zenodo.org/api/records/" + record_id + "/files-archive",
                "dateModified": date_modified,
            }
        ],
    }
    if date_modified:
        output["dateModified"] = date_modified

    # name / title
    title = root.find(f".//{NS}titles/{NS}title")
    if title is not None and title.text:
        if len(title.text) < 30000:
            output["name"] = title.text
        else:
            logger.info("Title too long for record %s: %s", record_id, title.text[:100] + "...")


    # description
    descriptions = [
        " ".join(desc.text.split())
        for desc in root.findall(f".//{NS}descriptions/{NS}description")
        if desc.text and desc.text.strip()
    ]
    if descriptions:
        output["description"] = "\n\n".join(descriptions)

    # datePublished
    if date_published := _date_published(root):
        output["datePublished"] = date_published

    # inLanguage
    language = root.find(f".//{NS}language")
    if language is not None and language.text:
        output["inLanguage"] = {"name": language.text}

    # keywords
    keywords = [s.text for s in root.findall(f".//{NS}subjects/{NS}subject") if s.text]
    if keywords:
        output["keywords"] = keywords

    # sdPublisher: Zenodo community related identifiers
    sdps = []
    for rel in root.findall(f".//{NS}relatedIdentifier"):
        if rel.text and "https://zenodo.org/communities/" in rel.text:
            sdps.append({"name": rel.text.rsplit("/", 1)[-1], "url": rel.text})
    if sdps:
        output["sdPublisher"] = sdps

    # doi
    doi = root.find(f".//{NS}identifier[@identifierType='DOI']")
    if doi is not None and doi.text:
        output["doi"] = doi.text

    # @type
    rt = root.find(f".//{NS}resourceType[@resourceTypeGeneral]")
    if rt is not None:
        general = (rt.get("resourceTypeGeneral") or "").lower().replace(" ", "")
        specific = (rt.text or "").lower().replace(" ", "")
        if (general, specific) in types:
            output["@type"] = types[(general, specific)]
        else:
            missing_types[(general, specific)] = (general, specific)

    # authors
    for creator in root.findall(f".//{NS}creator"):
        author = {}
        name = creator.find(f"./{NS}creatorName")
        affiliation = creator.find(f"./{NS}affiliation")
        orcid = creator.find(f"./{NS}nameIdentifier[@nameIdentifierScheme='ORCID']")
        if name is not None:
            author["name"] = name.text
        # elasticsearch cannot index strings longer than 32766 chars
        if affiliation is not None and affiliation.text and len(affiliation.text) < 30000:
            author["affiliation"] = {"name": affiliation.text}
        if orcid is not None:
            author["identifier"] = orcid.text
        if author:
            output["author"].append(author)

    # license / conditionsOfAccess
    for right in root.findall(f".//{NS}rights"):
        right_uri = (right.get("rightsURI") or "").strip()
        if "openaccess" in right_uri.lower():
            output["conditionsOfAccess"] = "Open"
        elif "restrictedaccess" in right_uri.lower():
            output["conditionsOfAccess"] = "Restricted"
        elif "embargoedaccess" in right_uri.lower():
            output["conditionsOfAccess"] = "Embargoed"
        elif right_uri:
            output["license"] = right_uri

    # citedBy
    cited_by = [
        {"url": rel.text}
        for rel in root.findall(f".//{NS}relatedIdentifier")
        if rel.get("relationType") == "IsCitedBy" and rel.get("relatedIdentifierType") == "URL"
    ]
    if cited_by:
        output["citedBy"] = cited_by

    # versioning (may rewrite _id / url / identifier / doi)
    _apply_version(output, root, url)

    # funding
    fundings = []
    for contributor in root.findall(f".//{NS}contributor[@contributorType='Funder']"):
        name = contributor.find(f"./{NS}contributorName")
        if name is not None:
            fundings.append({"funder": {"name": name.text}})

    for funder in root.findall(f".//{NS}fundingReference"):
        funding = {}
        funder_identifier = funder.find(f".//{NS}funderIdentifier")
        funder_name = funder.find(f".//{NS}funderName")
        award_title = funder.find(f".//{NS}awardTitle")
        award_number = funder.find(f".//{NS}awardNumber")
        if funder_identifier is not None or funder_name is not None:
            funding = {"funder": {}}
            if funder_identifier is not None:
                funding["funder"]["identifier"] = funder_identifier.text
            if funder_name is not None:
                funding["funder"]["name"] = funder_name.text
        if award_title is not None:
            funding["name"] = award_title.text
        if award_number is not None:
            funding["identifier"] = award_number.text
        if funding:
            funding["@type"] = "MonetaryGrant"
            fundings.append(funding)
    if fundings:
        output["funding"] = fundings

    # every doc has to have a type
    return output if output.get("@type") else None


def parse(url=TARBALL_URL, limit=None):
    """Stream the export tarball and yield transformed NDE docs."""
    missing_types = {}
    types = _build_type_mapping()
    resp = requests.get(url, stream=True)
    resp.raise_for_status()
    resp.raw.decode_content = True

    with tarfile.open(fileobj=resp.raw, mode="r|gz") as tar:
        members = itertools.islice(tar, limit) if limit else tar
        for count, member in enumerate(members, start=1):
            if count % 10000 == 0:
                logger.info("Processed %s records", count)
            if not member.isfile():
                continue
            record_id = os.path.splitext(os.path.basename(member.name))[0]
            content = tar.extractfile(member).read()
            try:
                doc = build_doc(content, record_id, types, missing_types)
            except Exception as e:
                logger.exception("Error processing record %s: %s", record_id, e)
                raise
            if doc is not None:
                yield doc

    logger.info("Finished processing %s records", count)
    if missing_types:
        logger.warning("Missing type transformation: %s", list(missing_types.keys()))




