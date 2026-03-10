import datetime
import json
import logging
import re

import requests
from requests.adapters import HTTPAdapter
from sql_database import NDEDatabase
from urllib3.util.retry import Retry

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")

PROXI_BASE = "https://proteomecentral.proteomexchange.org/api/proxi/v0.1"


def _get_session():
    """Create a requests session with retry logic and connection pooling."""
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=2, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retries, pool_connections=10, pool_maxsize=10)
    session.mount("https://", adapter)
    return session


def _get_term(terms_list, term_name, field="value"):
    """Find a term by name in a list of {name, value, accession} dicts and return the specified field."""
    if not terms_list:
        return None
    for term in terms_list:
        if term.get("name") == term_name:
            val = term.get(field)
            if val and str(val).strip():
                return str(val).strip()
    return None


def _has_term(terms_list, term_name):
    """Check if a term with the given name exists in the list."""
    if not terms_list:
        return False
    return any(term.get("name") == term_name for term in terms_list)


def _parse_date(date_str):
    """Parse various date formats to YYYY-MM-DD."""
    if not date_str or date_str == "missing":
        return None
    date_str = date_str.strip()
    # Remove timezone offset like +01:00
    date_str = re.sub(r"[+-]\d{2}:\d{2}$", "", date_str)
    for fmt in ("%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.datetime.strptime(date_str, fmt).strftime("%Y-%m-%d")
        except ValueError:
            continue
    return None


def _get_latest_history(dataset):
    """Get the latest history entry from datasetHistory."""
    history = dataset.get("datasetHistory")
    if not history:
        return {}
    for entry in history:
        if entry.get("isLatestRevision") == "Y":
            return entry
    return history[-1] if history else {}


def _parse_authors_from_contacts(contacts):
    """Extract deduplicated author list from contacts array."""
    authors = []
    seen_names = set()
    for contact in contacts:
        terms = contact.get("terms", [])
        name = _get_term(terms, "contact name")
        if not name or name in seen_names:
            continue
        seen_names.add(name)
        author = {"name": name}
        affiliation = _get_term(terms, "contact affiliation")
        if affiliation:
            author["affiliation"] = {"name": affiliation}
        authors.append(author)
    return authors


def _parse_keywords_from_keyword_list(keyword_list_str):
    """Parse submitter keywords from datasetHistory keywordList string.

    Format: 'ProteomeXchange project tag: PRIME-XS Project; submitter keyword: spikes,TMT,Eriwinia;'
    """
    keywords = []
    if not keyword_list_str:
        return keywords
    match = re.search(r"submitter keyword:\s*([^;]+)", keyword_list_str)
    if match:
        kw_str = match.group(1)
        keywords = [k.strip() for k in kw_str.split(",") if k.strip()]
    return keywords


def _parse_citations_from_html(publication_html, dataset_doi=None):
    """Parse citation info from HTML publication string in datasetHistory.

    Format: '<a href="https://www.ncbi.nlm.nih.gov/pubmed/23692960" ...>text</a>; ...'
    """
    citations = []
    if not publication_html or publication_html == "no publication":
        return citations

    urls = re.findall(r'href="([^"]+)"', publication_html)
    citation = {}
    for url in urls:
        if "pubmed" in url:
            pmid_match = re.search(r"/pubmed/(\d+)", url)
            if pmid_match:
                citation["pmid"] = pmid_match.group(1)
        elif "doi.org" in url:
            doi_match = re.search(r"doi\.org/(.+?)(?:\"|$)", url)
            if doi_match:
                doi_val = doi_match.group(1).strip()
                if doi_val != dataset_doi:
                    citation["doi"] = doi_val
    if citation:
        citations.append(citation)
    return citations


def _parse_dataset(dataset):
    """Parse a single PROXI dataset record into an NDE-schema document."""
    result = {}
    result["@type"] = "Dataset"

    # ── Identifiers ──────────────────────────────────────────────────
    identifiers = dataset.get("identifiers", [])
    pxd_id = _get_term(identifiers, "ProteomeXchange accession number")
    if not pxd_id:
        return None

    result["_id"] = f"proteomexchange_{pxd_id}"
    result["identifier"] = pxd_id

    doi = _get_term(identifiers, "Digital Object Identifier (DOI)")
    if doi:
        result["doi"] = doi

    version = _get_term(identifiers, "ProteomeXchange accession number version number")

    # ── Name & Description ───────────────────────────────────────────
    name = dataset.get("title")
    description = dataset.get("description")

    # ── Authors (from contacts) ──────────────────────────────────────
    contacts = dataset.get("contacts", [])
    authors = _parse_authors_from_contacts(contacts)

    # ── Species ──────────────────────────────────────────────────────
    species_list = []
    for sp in dataset.get("species", []):
        terms = sp.get("terms", [])
        sci_name = _get_term(terms, "taxonomy: scientific name")
        tax_id = _get_term(terms, "taxonomy: NCBI TaxID")
        if sci_name:
            species_obj = {"name": sci_name}
            if tax_id:
                species_obj["identifier"] = f"taxonomy:{tax_id}"
            species_list.append(species_obj)

    # ── Instruments ──────────────────────────────────────────────────
    instrument_names = []
    for inst in dataset.get("instruments", []):
        inst_name = inst.get("name")
        if inst_name:
            instrument_names.append(inst_name)

    # ── Keywords (submitter keywords only) ───────────────────────────
    keywords = []
    for kw in dataset.get("keywords", []):
        if kw.get("name") == "submitter keyword":
            val = kw.get("value", "")
            keywords.extend([k.strip() for k in val.split(",") if k.strip()])

    # ── Citations (from publications) ────────────────────────────────
    citations = []
    for pub in dataset.get("publications", []):
        terms = pub.get("terms", [])
        # Skip "no manuscript" marker entries
        if _has_term(terms, "Dataset with no associated published manuscript"):
            continue
        citation = {}
        pmid = _get_term(terms, "PubMed identifier")
        if pmid:
            citation["pmid"] = pmid
        pub_doi = _get_term(terms, "Digital Object Identifier (DOI)")
        if pub_doi and pub_doi != doi:  # exclude the dataset's own DOI
            citation["doi"] = pub_doi
        if citation:
            citations.append(citation)

    # ── Dataset Origins (reanalysis relationships) ─────────────────────
    is_based_on = []
    is_related_to = []
    for origin in dataset.get("datasetOrigins", []):
        terms = origin.get("terms", [])
        # Skip "Original data" entries — they indicate first-hand submissions
        if _has_term(terms, "Original data"):
            continue
        # "Data derived from previous dataset" marks the dataset as a reanalysis
        if _has_term(terms, "Data derived from previous dataset"):
            is_based_on.append({"identifier": pxd_id})
        # Sibling terms carry the parent dataset identifier
        parent_pxd = _get_term(terms, "ProteomeXchange accession number")
        if parent_pxd:
            is_related_to.append({"identifier": parent_pxd})
        parent_jpost = _get_term(terms, "jPOST dataset identifier")
        if parent_jpost:
            is_related_to.append({"identifier": parent_jpost})

    # ── sdPublisher (hosting repository + PRIDE URI) ─────────────────
    summary = dataset.get("datasetSummary") or {}
    hosting_repo = summary.get("hostingRepository")

    pride_url = None
    for link in dataset.get("fullDatasetLinks", []):
        if link.get("name") == "PRIDE project URI":
            pride_url = link.get("value")
            break

    # ── Dates ────────────────────────────────────────────────────────
    announce_date = summary.get("announceDate")
    date_published = _parse_date(announce_date) if announce_date else None

    history = _get_latest_history(dataset)

    identifier_date = history.get("identifierDate")
    date_created = _parse_date(identifier_date) if identifier_date else None

    revision_date = history.get("revisionDate")
    submission_date = history.get("submissionDate")
    date_modified = _parse_date(revision_date) if revision_date else _parse_date(submission_date)

    # ── Distribution (from datasetFiles) ─────────────────────────────
    distributions = []
    for df in dataset.get("datasetFiles", []):
        dist = {}
        df_name = df.get("name")
        df_value = df.get("value")
        if df_name:
            dist["name"] = df_name
        if df_value:
            dist["contentUrl"] = df_value
        if dist:
            distributions.append(dist)

    # ══ Fallbacks from datasetHistory (only when primary is missing) ═
    if not name:
        name = history.get("title")

    if not authors:
        primary_submitter = history.get("primarySubmitter")
        if primary_submitter:
            authors = [{"name": primary_submitter}]

    if not species_list:
        hist_species = history.get("species")
        if hist_species:
            species_list = [{"name": hist_species}]

    if not instrument_names:
        hist_instrument = history.get("instrument")
        if hist_instrument:
            instrument_names = [hist_instrument]

    if not keywords:
        keywords = _parse_keywords_from_keyword_list(history.get("keywordList"))

    if not version:
        hist_version = history.get("revisionNumber")
        if hist_version:
            version = str(hist_version)

    if not citations:
        hist_pub = history.get("publication")
        if hist_pub:
            citations = _parse_citations_from_html(hist_pub, doi)

    # ══ Build final result ═══════════════════════════════════════════
    if name:
        result["name"] = name
    if description:
        result["description"] = description
    if authors:
        result["author"] = authors
    if species_list:
        result["species"] = species_list
    if instrument_names:
        result["sample"] = {"instrument": [{"name": n} for n in instrument_names]}
    if keywords:
        result["keywords"] = keywords
    if citations:
        result["citation"] = citations
    if distributions:
        result["distribution"] = distributions
    if is_based_on:
        result["isBasedOn"] = is_based_on
    if is_related_to:
        result["isRelatedTo"] = is_related_to
    # if version:
        # result["version"] = str(version)
    if date_published:
        result["datePublished"] = date_published
    if date_created:
        result["dateCreated"] = date_created
    if date_modified:
        result["dateModified"] = date_modified

    # sdPublisher
    if hosting_repo:
        sd_publisher = {"name": hosting_repo}
        if pride_url:
            sd_publisher["url"] = pride_url
        result["sdPublisher"] = sd_publisher

    # URL
    dataset_url = pride_url or f"http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={pxd_id}"
    result["url"] = dataset_url

    # includedInDataCatalog
    result["includedInDataCatalog"] = {
        "name": "ProteomeXchange",
        "archivedAt": dataset_url,
        "@type": "DataCatalog",
        "url": "https://www.proteomexchange.org/",
        "versionDate": datetime.date.today().isoformat(),
    }

    result["@type"] = "Dataset"

    # Remove any keys with None values
    return {k: v for k, v in result.items() if v is not None}


class ProteomeXchange(NDEDatabase):
    SQL_DB = "proteomexchange.db"
    EXPIRE = datetime.timedelta(days=30)

    def _get_all_pxd_ids(self, session):
        """Fetch all PXD identifiers from the listing endpoint in a single request."""
        url = f"{PROXI_BASE}/datasets?pageSize=100000&resultType=full"
        logger.info("Fetching all dataset IDs from listing endpoint")

        response = session.get(url, timeout=600)
        response.raise_for_status()
        data = response.json()

        datasets = data.get("datasets", [])
        pxd_ids = []
        for row in datasets:
            if isinstance(row, list) and row:
                pxd_ids.append(row[0])  # First column is dataset identifier

        result_set = data.get("result_set", {})
        total = result_set.get("n_available_rows", "?")
        logger.info("Collected %s PXD IDs (total available: %s)", len(pxd_ids), total)
        return pxd_ids

    def _fetch_dataset(self, session, pxd_id):
        """Fetch a single dataset from the individual PROXI endpoint."""
        url = f"{PROXI_BASE}/datasets/{pxd_id}"
        response = session.get(url, timeout=60)
        response.raise_for_status()
        return response.json()

    def load_cache(self):
        """Download all ProteomeXchange datasets and yield (_id, json_str) tuples for the cache."""
        session = _get_session()

        pxd_ids = self._get_all_pxd_ids(session)
        logger.info("Starting individual dataset fetches for %s datasets", len(pxd_ids))

        count = 0
        errors = 0
        for i, pxd_id in enumerate(pxd_ids):
            try:
                dataset = self._fetch_dataset(session, pxd_id)
                yield (pxd_id, json.dumps(dataset))
                count += 1
                if count % 500 == 0:
                    logger.info("Cached %s / %s records (errors: %s)", count, len(pxd_ids), errors)
            except requests.exceptions.RequestException as e:
                errors += 1
                logger.error("Error fetching %s: %s", pxd_id, e)
                continue
            except Exception as e:
                errors += 1
                logger.error("Error caching %s: %s", pxd_id, e)
                continue

        logger.info("Finished loading cache. Total cached: %s, Errors: %s", count, errors)

    def parse(self, records):
        """Parse cached records into NDE-schema dataset documents."""
        count = 0
        errors = 0
        for record in records:
            pxd_id = record[0]
            try:
                dataset = json.loads(record[1])
                result = _parse_dataset(dataset)
                if result:
                    yield result
                    count += 1
                    if count % 500 == 0:
                        logger.info("Parsed %s records (errors: %s)", count, errors)
            except Exception as e:
                errors += 1
                logger.error("Error parsing %s: %s", pxd_id, e)
                continue

        logger.info("Finished parsing. Total records: %s, Errors: %s", count, errors)
