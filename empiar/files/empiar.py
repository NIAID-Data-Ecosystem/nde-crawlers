import logging
import os
import re
import time
from typing import Any, Iterable

import dateutil.parser
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

API_ENTRY_URL = "https://www.ebi.ac.uk/empiar/api/entry/"
GET_ENTRY_URL_TEMPLATE = "https://www.ebi.ac.uk/empiar/api/entry/{accession}/"
CATALOG_URL = "https://www.ebi.ac.uk/empiar/"
RECORD_URL_TEMPLATE = "https://www.ebi.ac.uk/empiar/{accession}/"
DOWNLOAD_URL_TEMPLATE = "https://ftp.ebi.ac.uk/empiar/world_availability/{numeric_id}/"
STATUS_LABELS = {
    "REL": "Released",
    "OBS": "Obsolete",
}

ALWAYS_ALLOWED_KEYS = {"@context", "@type"}

CATALOG = {
    "@type": "DataCatalog",
    "name": "Electron Microscopy Public Image Archive",
    "alternateName": "EMPIAR",
    "identifier": "EMPIAR",
    "url": CATALOG_URL,
}

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "identifier": "topic_0611",
        "inDefinedTermSet": "EDAM",
        "name": "Electron microscopy",
        "url": "http://edamontology.org/topic_0611",
    },
    {
        "@type": "DefinedTerm",
        "identifier": "topic_1317",
        "inDefinedTermSet": "EDAM",
        "name": "Structural biology",
        "url": "http://edamontology.org/topic_1317",
    },
    {
        "@type": "DefinedTerm",
        "identifier": "topic_3383",
        "inDefinedTermSet": "EDAM",
        "name": "Bioimaging",
        "url": "http://edamontology.org/topic_3383",
    },
]

STATUS_MAPPINGS = {
    "REL": {"conditionsOfAccess": "Open", "creativeWorkStatus": "Active"},
    "EMDBPUB": {"conditionsOfAccess": "Embargoed", "creativeWorkStatus": "Draft"},
    "HPUB": {"conditionsOfAccess": "Embargoed", "creativeWorkStatus": "Draft"},
    "HPRE": {"conditionsOfAccess": "Embargoed", "creativeWorkStatus": "Draft"},
    "HOLD": {"conditionsOfAccess": "Embargoed", "creativeWorkStatus": "Draft"},
    "PROC": {"conditionsOfAccess": "Embargoed", "creativeWorkStatus": "Draft"},
    "WAIT": {"conditionsOfAccess": "Embargoed", "creativeWorkStatus": "Draft"},
    "OBS": {"conditionsOfAccess": "Open", "creativeWorkStatus": "Retired"},
    "WDRN": {"conditionsOfAccess": "Open", "creativeWorkStatus": "Retired"},
    "UNARCH": {"conditionsOfAccess": "Open", "creativeWorkStatus": "Retired"},
}

EXPERIMENT_MAPPINGS = {
    "ATUM-SEM": [
        {
            "@type": "DefinedTerm",
            "alternateName": "ATUM-SEM",
            "name": "Automated Tape-collecting Ultramicrotome Scanning Electron Microscopy",
        },
        {
            "@type": "DefinedTerm",
            "identifier": "D008867",
            "inDefinedTermSet": "MESH",
            "name": "Microtomy",
            "url": "http://id.nlm.nih.gov/mesh/D008867",
        },
        {
            "@type": "DefinedTerm",
            "identifier": "CHMO_0000073",
            "inDefinedTermSet": "CHMO",
            "name": "scanning electron microscopy",
            "url": "http://purl.obolibrary.org/obo/CHMO_0000073",
        },
    ],
    "CLEM": [
        {"@type": "DefinedTerm", "alternateName": "CLEM", "name": "correlative light-electron microscopy"},
        {
            "@type": "DefinedTerm",
            "identifier": "CHMO_0000068",
            "inDefinedTermSet": "CHMO",
            "name": "electron microscopy",
            "url": "http://purl.obolibrary.org/obo/CHMO_0000068",
        },
        {
            "@type": "DefinedTerm",
            "identifier": "CHMO_0000102",
            "inDefinedTermSet": "CHMO",
            "name": "light microscopy assay",
            "url": "http://purl.obolibrary.org/obo/CHMO_0000102",
        },
    ],
    "CLXM": [
        {"@type": "DefinedTerm", "alternateName": "CLXM", "name": "correlative light X-ray microscopy"},
        {
            "@type": "DefinedTerm",
            "identifier": "CHMO_0000102",
            "inDefinedTermSet": "CHMO",
            "name": "light microscopy assay",
            "url": "http://purl.obolibrary.org/obo/CHMO_0000102",
        },
        {
            "@type": "DefinedTerm",
            "identifier": "CHMO_0002312",
            "name": "X-ray microscopy",
            "url": "http://purl.obolibrary.org/obo/CHMO_0002312",
        },
    ],
    "EMDB": [
        {
            "@type": "DefinedTerm",
            "identifier": "CHMO_0000068",
            "inDefinedTermSet": "CHMO",
            "name": "electron microscopy",
            "url": "http://purl.obolibrary.org/obo/CHMO_0000068",
        }
    ],
    "FIB SEM": [
        {
            "@type": "DefinedTerm",
            "alternateName": "FIB SEM",
            "identifier": "FBbi_00050000",
            "inDefinedTermSet": "EFO",
            "name": "focussed ion beam scanning electron microscopy (FIB-SEM)",
            "url": "http://purl.obolibrary.org/obo/FBbi_00050000",
        }
    ],
    "Hard X-ray/X-ray microCT": [
        {
            "@type": "DefinedTerm",
            "alternateName": "Hard X-ray/X-ray microCT",
            "name": "Hard X-ray/X-ray micro-computed tomography",
        },
        {
            "@type": "DefinedTerm",
            "identifier": "D055114",
            "inDefinedTermSet": "MESH",
            "name": "X-Ray Microtomography",
            "url": "http://id.nlm.nih.gov/mesh/D055114",
        },
    ],
    "IHM": [
        {"@type": "DefinedTerm", "alternateName": "IHM", "name": "integrative hybrid modeling"},
        {
            "@type": "DefinedTerm",
            "identifier": "topic_2275",
            "inDefinedTermSet": "EDAM",
            "name": "Molecular Modeling",
            "url": "http://edamontology.org/topic_2275",
        },
    ],
    "MicroED": [
        {"@type": "DefinedTerm", "alternateName": "MicroED", "name": "microcrystal electron diffraction"},
        {
            "@type": "DefinedTerm",
            "identifier": "CHMO_0000142",
            "inDefinedTermSet": "CHMO",
            "name": "electron diffraction",
            "url": "http://purl.obolibrary.org/obo/CHMO_0000142",
        },
    ],
    "SBF-SEM": [
        {
            "@type": "DefinedTerm",
            "alternateName": "SBF-SEM",
            "identifier": "FBbi_00000585",
            "inDefinedTermSet": "EFO",
            "name": "serial block face SEM (SBFSEM)",
            "url": "http://purl.obolibrary.org/obo/FBbi_00000585",
        }
    ],
    "SXT": [
        {"@type": "DefinedTerm", "alternateName": "SXT", "name": "soft x-ray tomography"},
        {
            "@type": "DefinedTerm",
            "identifier": "FBbi_00001003",
            "inDefinedTermSet": "EFO",
            "name": "X-ray tomography",
            "url": "http://purl.obolibrary.org/obo/FBbi_00001003",
        },
    ],
    "ssET": [
        {"@type": "DefinedTerm", "alternateName": "ssET", "name": "serial section electron tomography"},
    ],
}

ALLOWED_OUTPUT_PATHS = {
    ("@context",),
    ("@type",),
    ("_id",),
    ("author", "affiliation", "name"),
    ("author", "email"),
    ("author", "familyName"),
    ("author", "givenName"),
    ("author", "identifier"),
    ("author", "name"),
    ("author", "url"),
    ("citation", "author", "name"),
    ("citation", "datePublished"),
    ("citation", "doi"),
    ("citation", "issueNumber"),
    ("citation", "journalName"),
    ("citation", "journalNameAbbrev"),
    ("citation", "name"),
    ("citation", "pagination"),
    ("citation", "pmid"),
    ("citation", "url"),
    ("citation", "volumeNumber"),
    ("conditionsOfAccess",),
    ("creativeWorkStatus",),
    ("dateCreated",),
    ("dateModified",),
    ("datePublished",),
    ("description",),
    ("distribution", "contentSize"),
    ("distribution", "contentUrl"),
    ("distribution", "description"),
    ("distribution", "encodingFormat"),
    ("doi",),
    ("funding", "funder", "name"),
    ("funding", "identifier"),
    ("identifier",),
    ("inLanguage",),
    ("includedInDataCatalog", "alternateName"),
    ("includedInDataCatalog", "archivedAt"),
    ("includedInDataCatalog", "identifier"),
    ("includedInDataCatalog", "name"),
    ("includedInDataCatalog", "url"),
    ("isAccessibleForFree",),
    ("isBasedOn",),
    ("isBasisFor", "identifier"),
    ("isBasisFor", "includedInDataCatalog"),
    ("isPartOf", "identifier"),
    ("isPartOf", "includedInDataCatalog"),
    ("isRelatedTo", "identifier"),
    ("isRelatedTo", "includedInDataCatalog"),
    ("keywords",),
    ("measurementTechnique",),
    ("measurementTechnique", "alternateName"),
    ("measurementTechnique", "identifier"),
    ("measurementTechnique", "inDefinedTermSet"),
    ("measurementTechnique", "name"),
    ("measurementTechnique", "url"),
    ("name",),
    ("pmids",),
    ("sameAs",),
    ("topicCategory",),
    ("url",),
}

START_ID = int(os.getenv("EMPIAR_START_ID", "10000"))
HARD_MAX_ID = int(os.getenv("EMPIAR_HARD_MAX_ID", "99999"))
BATCH_SIZE = int(os.getenv("EMPIAR_BATCH_SIZE", "500"))
EMPTY_CHUNK_LIMIT = int(os.getenv("EMPIAR_EMPTY_CHUNK_LIMIT", "3"))
SLEEP = float(os.getenv("EMPIAR_SLEEP", "0.25"))
TIMEOUT = int(os.getenv("EMPIAR_TIMEOUT", "120"))
FETCH_FULL_RECORDS = os.getenv("EMPIAR_FETCH_FULL_RECORDS", "true").casefold() not in {"0", "false", "no"}

EMPIAR_ID_RE = re.compile(r"(\d+)$")


def insert_value(d, key, value, extend=False):
    """ Insert a value into a dictionary, handling existing keys by converting to lists or extending strings as needed.
    """

    if key in d and not extend:
        if isinstance(d[key], list):
            if isinstance(value, list):
                for item in value:
                    if item not in d[key]:
                        d[key].append(item)
            elif value not in d[key]:
                d[key].append(value)
        else:
            if isinstance(value, list):
                d[key] = [d[key]] + [v for v in value if v != d[key]]
            elif d[key] != value:
                d[key] = [d[key], value]
    elif d.get(key) and extend:
        d[key] = (d.get(key) + " " + value).strip()
    else:
        d[key] = value

def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(val, ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt


def _clean_text(value):
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _as_list(value):
    if isinstance(value, list):
        return value
    if value in (None, "", []):
        return []
    return [value]


def _unique(values):
    seen = set()
    unique_values = []
    for value in values:
        if value in (None, ""):
            continue
        if isinstance(value, dict):
            key = tuple(sorted((k, str(v)) for k, v in value.items()))
        else:
            key = str(value).casefold()
        if key not in seen:
            unique_values.append(value)
            seen.add(key)
    return unique_values


def _numeric_id(accession):
    match = EMPIAR_ID_RE.search(str(accession))
    return match.group(1) if match else str(accession)


def _entry_sort_key(accession):
    try:
        return int(_numeric_id(accession))
    except ValueError:
        return 0


def _fetch_entries(session, body, retries=3):
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "nde-empiar-crawler/0.1",
    }
    for attempt in range(1, retries + 1):
        try:
            response = session.post(API_ENTRY_URL, data=str(body), headers=headers, timeout=TIMEOUT)
            if response.status_code >= 500 and attempt < retries:
                logger.warning("EMPIAR API returned %s for %s; retrying", response.status_code, body)
                time.sleep(attempt)
                continue
            response.raise_for_status()
            data = response.json()
            return data if isinstance(data, dict) else {}
        except requests.RequestException as exc:
            if attempt == retries:
                raise
            logger.warning("EMPIAR API request failed for %s: %s; retrying", body, exc)
            time.sleep(attempt)
    return {}


def _fetch_entry(session, accession, retries=3):
    headers = {
        "Accept": "application/json",
        "User-Agent": "nde-empiar-crawler/0.1",
    }
    normalized_accession = str(accession)
    if not normalized_accession.startswith("EMPIAR-"):
        normalized_accession = f"EMPIAR-{_numeric_id(normalized_accession)}"
    url = GET_ENTRY_URL_TEMPLATE.format(accession=normalized_accession)
    for attempt in range(1, retries + 1):
        try:
            response = session.get(url, headers=headers, timeout=TIMEOUT)
            if response.status_code >= 500 and attempt < retries:
                logger.warning("EMPIAR detail API returned %s for %s; retrying", response.status_code, normalized_accession)
                time.sleep(attempt)
                continue
            response.raise_for_status()
            data = response.json()
            if not isinstance(data, dict):
                return {}
            if isinstance(data.get(normalized_accession), dict):
                return data[normalized_accession]
            if isinstance(data.get(accession), dict):
                return data[accession]
            return data if data.get("title") else {}
        except requests.RequestException as exc:
            if attempt == retries:
                logger.warning("EMPIAR detail API request failed for %s: %s", normalized_accession, exc)
                return {}
            logger.warning("EMPIAR detail API request failed for %s: %s; retrying", normalized_accession, exc)
            time.sleep(attempt)
    return {}


def _with_full_record(session, accession, record):
    if not FETCH_FULL_RECORDS:
        return record
    full_record = _fetch_entry(session, accession)
    return full_record or record


def _chunked(values, chunk_size):
    for index in range(0, len(values), chunk_size):
        yield values[index:index + chunk_size]


def _has_allowed_descendant(path):
    return any(allowed_path[:len(path)] == path and len(allowed_path) > len(path) for allowed_path in ALLOWED_OUTPUT_PATHS)


def _is_allowed_or_container(path):
    return path in ALLOWED_OUTPUT_PATHS or _has_allowed_descendant(path)


def _allows_full_subtree(path):
    return path in ALLOWED_OUTPUT_PATHS and not _has_allowed_descendant(path)


def _is_empty(value):
    return value is None or value == "" or value == [] or value == {}


def _has_non_metadata_value(value):
    if not isinstance(value, dict):
        return not _is_empty(value)
    return any(key not in ALWAYS_ALLOWED_KEYS and not _is_empty(item) for key, item in value.items())


def _filter_to_mapping(value: Any, path: tuple[str, ...] = ()):
    if path and _allows_full_subtree(path):
        return value

    if isinstance(value, list):
        filtered_values = []
        for item in value:
            filtered_item = _filter_to_mapping(item, path)
            if not _is_empty(filtered_item) and _has_non_metadata_value(filtered_item):
                filtered_values.append(filtered_item)
        return filtered_values or None

    if isinstance(value, dict):
        filtered = {}
        for key, item in value.items():
            child_path = path + (key,)
            if key in ALWAYS_ALLOWED_KEYS:
                filtered[key] = item
                continue
            if not _is_allowed_or_container(child_path):
                continue
            filtered_item = _filter_to_mapping(item, child_path)
            if not _is_empty(filtered_item):
                filtered[key] = filtered_item
        if path and not _has_non_metadata_value(filtered):
            return None
        return filtered

    return value if path in ALLOWED_OUTPUT_PATHS else None


def iter_empiar_records(
    ids: Iterable[str] | None = None,
    start_id: int = START_ID,
    hard_max_id: int = HARD_MAX_ID,
    batch_size: int = BATCH_SIZE,
    empty_chunk_limit: int = EMPTY_CHUNK_LIMIT,
    sleep: float = SLEEP,
):
    """
    Yield (accession, record) pairs from the official EMPIAR entry API.

    The API exposes a documented POST endpoint that accepts comma-separated IDs
    and numeric ranges, but it does not expose a list-all endpoint. For a full
    crawl, scan numeric accession ranges and stop after several empty chunks.
    """
    session = requests.Session()

    if ids is not None:
        normalized_ids = [str(value).replace("EMPIAR-", "").strip() for value in ids if str(value).strip()]
        for chunk in _chunked(normalized_ids, batch_size):
            body = ",".join(chunk)
            entries = _fetch_entries(session, body)
            for accession in sorted(entries, key=_entry_sort_key):
                yield accession, _with_full_record(session, accession, entries[accession])
            time.sleep(sleep)
        return

    empty_chunks = 0
    current = start_id
    while current <= hard_max_id and empty_chunks < empty_chunk_limit:
        end = min(current + batch_size - 1, hard_max_id)
        body = f"{current}-{end}"
        logger.info("Fetching EMPIAR entries %s", body)
        entries = _fetch_entries(session, body)
        if entries:
            empty_chunks = 0
            logger.info("Fetched %d EMPIAR entries from %s", len(entries), body)
            for accession in sorted(entries, key=_entry_sort_key):
                yield accession, _with_full_record(session, accession, entries[accession])
        else:
            empty_chunks += 1
            logger.info("No EMPIAR entries found in %s (%d empty chunks)", body, empty_chunks)
        current = end + 1
        time.sleep(sleep)


def _person_from_author(raw_author, role=None):
    raw_author = raw_author.get("author") if isinstance(raw_author, dict) and raw_author.get("author") else raw_author
    if not isinstance(raw_author, dict):
        return None

    person = {"@type": "Person"}
    name = _clean_text(raw_author.get("name"))
    given_name = _clean_text(raw_author.get("first_name"))
    family_name = _clean_text(raw_author.get("last_name"))
    if name:
        person["name"] = name
    elif given_name or family_name:
        person["name"] = " ".join(part for part in [given_name, family_name] if part)
    if given_name:
        person["givenName"] = given_name
    if family_name:
        person["familyName"] = family_name
    if orcid := _clean_text(raw_author.get("author_orcid")):
        person["identifier"] = orcid
        person["url"] = f"https://orcid.org/{orcid}"
    if email := _clean_text(raw_author.get("email")):
        person["email"] = email
    if organization := _clean_text(raw_author.get("organization")):
        person["affiliation"] = {"@type": "Organization", "name": organization}
    if role:
        person["role"] = role
    return person if person.keys() - {"@type"} else None


def _citation_from_record(raw_citation):
    if not isinstance(raw_citation, dict):
        return None

    citation = {"@type": "ScholarlyArticle"}
    if title := _clean_text(raw_citation.get("title")):
        citation["name"] = title
    if doi := _clean_text(raw_citation.get("doi")):
        citation["doi"] = doi
        citation["identifier"] = doi
        citation["url"] = f"https://doi.org/{doi}"
    if pmid := _clean_text(raw_citation.get("pubmedid")):
        citation["pmid"] = pmid
    if year := _clean_text(raw_citation.get("year")):
        citation["datePublished"] = year if re.fullmatch(r"\d{4}", year) else _to_iso_date(year)
    if journal := _clean_text(raw_citation.get("journal")):
        citation["journalName"] = journal
    if journal_abbrev := _clean_text(raw_citation.get("journal_abbreviation")):
        citation["journalNameAbbrev"] = journal_abbrev
    if volume := _clean_text(raw_citation.get("volume")):
        citation["volumeNumber"] = volume
    if issue := _clean_text(raw_citation.get("issue")):
        citation["issueNumber"] = issue

    first_page = _clean_text(raw_citation.get("first_page"))
    last_page = _clean_text(raw_citation.get("last_page"))
    if first_page and last_page:
        citation["pagination"] = f"{first_page}-{last_page}"
    elif first_page:
        citation["pagination"] = first_page

    authors = []
    for author in _as_list(raw_citation.get("authors")):
        person = _person_from_author(author)
        if person:
            authors.append(person)
    if authors:
        citation["author"] = authors

    return citation if citation.keys() - {"@type"} else None


def _funding_from_grant(raw_grant):
    if not isinstance(raw_grant, dict):
        return None

    funding = {"@type": "MonetaryGrant"}
    if code := _clean_text(raw_grant.get("code")):
        if code.upper() not in {"N/A", "NA", "NONE"}:
            funding["identifier"] = code
    if body := _clean_text(raw_grant.get("funding_body")):
        funding["funder"] = {"@type": "Organization", "name": body}
        if country := _clean_text(raw_grant.get("country")):
            funding["funder"]["description"] = f"Country: {country}"
    return funding if funding.keys() - {"@type"} else None


def _relationship(accession, catalog_name, catalog_url, relationship, url_template=None):
    accession = accession.get("name") if isinstance(accession, dict) else accession
    accession = _clean_text(accession)
    if not accession:
        return None
    if url_template:
        url = url_template.format(accession=accession)
    else:
        url = f"{catalog_url.rstrip('/')}/{accession}/"
    return {
        "@type": "CreativeWork",
        "identifier": accession,
        "name": accession,
        "relationship": relationship,
        "url": url,
        "includedInDataCatalog": {
            "@type": "DataCatalog",
            "name": catalog_name,
            "url": catalog_url,
        },
    }


def _measurement_terms(record):
    experiment_type = _clean_text(record.get("experiment_type"))
    if not experiment_type:
        return []
    terms = EXPERIMENT_MAPPINGS.get(experiment_type)
    if terms:
        return terms
    return [{"@type": "DefinedTerm", "name": experiment_type}]


def _description(record):
    parts = []
    names = _unique(
        _clean_text(image_set.get("name"))
        for image_set in _as_list(record.get("imagesets"))
        if isinstance(image_set, dict)
    )
    parts.extend(names)
    details = _unique(
        _clean_text(image_set.get("details"))
        for image_set in _as_list(record.get("imagesets"))
        if isinstance(image_set, dict)
    )
    parts.extend(details[:3])

    return " ".join(parts) if parts else None


def _normalize_file_path(value):
    if isinstance(value, dict):
        value = value.get("path")
    path = _clean_text(value)
    if not path:
        return None
    if path.startswith("path="):
        path = path.split("=", 1)[1]
    return path.lstrip("/")


def _workflow_file(record, download_url):
    path = _normalize_file_path(record.get("workflow_file"))
    if not path:
        return None
    return {
        "@type": "CreativeWork",
        "name": "Workflow",
        "url": f"{download_url.rstrip('/')}/{path}",
    }


def _encoding_formats(record):
    formats = []
    for image_set in _as_list(record.get("imagesets")):
        if not isinstance(image_set, dict):
            continue
        for key in ("header_format", "data_format"):
            if value := _clean_text(image_set.get(key)):
                formats.append(value)
        for key in ("micrographs_file_pattern", "picked_particles_file_pattern"):
            pattern = _clean_text(image_set.get(key))
            if not pattern:
                continue
            for extension in re.findall(r"\.([A-Za-z0-9]+)(?:$|[;,\s*])", pattern):
                formats.append(extension.upper())
    return _unique(formats)


def _distribution_description(record):
    parts = []
    for image_set in _as_list(record.get("imagesets")):
        if not isinstance(image_set, dict):
            continue
        if pattern := _clean_text(image_set.get("micrographs_file_pattern")):
            parts.append(f"Micrographs file pattern: {pattern}")
        if pattern := _clean_text(image_set.get("picked_particles_file_pattern")):
            parts.append(f"Picked particles file pattern: {pattern}")
    return "; ".join(_unique(parts)) if parts else None


def _status_fields(record):
    status = _clean_text(record.get("status"))
    if not status:
        for version in reversed(_as_list(record.get("version_history"))):
            if isinstance(version, dict) and (status := _clean_text(version.get("status_code"))):
                break
    return STATUS_MAPPINGS.get(status, {"creativeWorkStatus": STATUS_LABELS.get(status, status)}) if status else {}


def _latest_modified_date(record):
    dates = [
        _to_iso_date(record.get("update_date")),
        _to_iso_date(record.get("obsolete_date")),
    ]
    for version in _as_list(record.get("version_history")):
        if isinstance(version, dict):
            dates.append(_to_iso_date(version.get("date")))
    dates = [date for date in dates if date]
    return max(dates) if dates else None


def parse(ids: Iterable[str] | None = None):
    for accession, record in iter_empiar_records(ids=ids):
        numeric_id = _numeric_id(accession)
        url = RECORD_URL_TEMPLATE.format(accession=accession)
        download_url = DOWNLOAD_URL_TEMPLATE.format(numeric_id=numeric_id)

        output = {
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": accession,
            "identifier": accession,
            "url": url,
            "includedInDataCatalog": {**CATALOG, "archivedAt": url},
            "conditionsOfAccess": "Open",
            "topicCategory": TOPIC_CATEGORY,
        }

        if name := _clean_text(record.get("title")):
            insert_value(output, "name", name)

        if description := _description(record):
            output["description"] = description

        distribution = {"@type": "DataDownload"}
        if dataset_size := _clean_text(record.get("dataset_size")):
            distribution["contentSize"] = dataset_size
        if encoding_formats := _encoding_formats(record):
            distribution["encodingFormat"] = encoding_formats[0] if len(encoding_formats) == 1 else encoding_formats
        if distribution_description := _distribution_description(record):
            distribution["description"] = distribution_description
        if distribution.keys() - {"@type"}:
            output["distribution"] = [distribution]

        if doi := _clean_text(record.get("entry_doi")):
            output["doi"] = doi
            insert_value(output, "sameAs", f"https://doi.org/{doi}")

        date_created = _to_iso_date(record.get("deposition_date"))
        date_published = _to_iso_date(record.get("release_date"))
        if date_created:
            output["dateCreated"] = date_created
        if date_published:
            output["datePublished"] = date_published
        if date_modified := _latest_modified_date(record):
            output["dateModified"] = date_modified

        for field, value in _status_fields(record).items():
            output[field] = value
        if conditions_of_access := output.get("conditionsOfAccess"):
            output["isAccessibleForFree"] = conditions_of_access == "Open"

        if measurement_terms := _measurement_terms(record):
            output["measurementTechnique"] = measurement_terms

        if workflow_file := _workflow_file(record, download_url):
            output["isBasedOn"] = workflow_file

        keywords = []
        for key in ("experiment_type", "scale"):
            if value := _clean_text(record.get(key)):
                keywords.append(value)
        for image_set in _as_list(record.get("imagesets")):
            if isinstance(image_set, dict) and (category := _clean_text(image_set.get("category"))):
                keywords.append(category)
        for citation in _as_list(record.get("citation")):
            if isinstance(citation, dict) and str(citation.get("preprint")).casefold() == "true":
                keywords.append("preprint")
        if keywords:
            output["keywords"] = _unique(keywords)

        authors = []
        for pi in _as_list(record.get("principal_investigator")):
            person = _person_from_author(pi)
            if person:
                authors.append(person)
        for corresponding_author in _as_list(record.get("corresponding_author")):
            person = _person_from_author(corresponding_author)
            if person:
                authors.append(person)
        for author in _as_list(record.get("authors")):
            person = _person_from_author(author)
            if person:
                authors.append(person)
        if authors:
            output["author"] = _unique(authors)

        citations = []
        languages = []
        for citation in _as_list(record.get("citation")):
            parsed_citation = _citation_from_record(citation)
            if parsed_citation:
                citations.append(parsed_citation)
            if isinstance(citation, dict) and (language := _clean_text(citation.get("language"))):
                languages.append(language)
        if citations:
            output["citation"] = citations
            pmids = _unique(citation.get("pmid") for citation in citations if isinstance(citation, dict))
            if pmids:
                output["pmids"] = ",".join(pmids)
        if languages:
            languages = _unique(languages)
            output["inLanguage"] = languages[0] if len(languages) == 1 else languages

        fundings = []
        for grant in _as_list(record.get("grant_references")):
            funding = _funding_from_grant(grant)
            if funding:
                fundings.append(funding)
        if fundings:
            output["funding"] = fundings

        basis_for = []
        for emdb_id in _as_list(record.get("cross_references")):
            item = _relationship(
                emdb_id,
                "Electron Microscopy Data Bank",
                "https://www.ebi.ac.uk/emdb",
                "related EMDB entry",
                url_template="https://www.ebi.ac.uk/emdb/{accession}/",
            )
            if item:
                basis_for.append(item)
        for pdb_id in _as_list(record.get("related_pdb_entries")):
            item = _relationship(
                pdb_id,
                "Protein Data Bank",
                "https://www.rcsb.org/",
                "related PDB entry",
                url_template="https://www.rcsb.org/structure/{accession}",
            )
            if item:
                basis_for.append(item)
        if basis_for:
            output["isBasisFor"] = basis_for

        part_of = []
        for biostudies_id in _as_list(record.get("biostudies_references")):
            item = _relationship(
                biostudies_id,
                "BioStudies",
                "https://www.ebi.ac.uk/biostudies/",
                "related BioStudies entry",
                url_template="https://www.ebi.ac.uk/biostudies/studies/{accession}",
            )
            if item:
                part_of.append(item)
        if part_of:
            output["isPartOf"] = part_of

        related = []
        for idr_id in _as_list(record.get("idr_references")):
            item = _relationship(
                idr_id,
                "Image Data Resource",
                "https://idr.openmicroscopy.org",
                "related IDR entry",
            )
            if item:
                related.append(item)
        for empiar_id in _as_list(record.get("empiar_references")):
            item = _relationship(
                empiar_id,
                "Electron Microscopy Public Image Archive",
                CATALOG_URL,
                "related EMPIAR entry",
                url_template="https://www.ebi.ac.uk/empiar/{accession}/",
            )
            if item:
                related.append(item)
        if related:
            output["isRelatedTo"] = related

        yield _filter_to_mapping(output)
