"""Article→dataset citations from Make Data Count (MDC).

Runs for Dataset records once the `make_data_count` dumper has loaded a release.
A record's `_id`, `identifier` and `doi`, and its URL for SRA accessions, are
looked up in that release:

  * Primary papers are added to `citation` when the record has none,
  * Secondary (reuse) papers are added to `citedBy`,
  * their number becomes an `interactionStatistic` counter.

Papers the record already lists in `citation` or `citedBy`, by DOI, PMID or
PMCID, are not added again. Citing papers are described from PubMed (PMC ID
converter, ESearch, ESummary) or else Crossref. Answers are cached next to the
release; misses are retried after 90 days.

https://github.com/digital-science/ds-make-data-count-pipeline-external
"""

import html
import os
import re
import sqlite3
import time
from collections import Counter
from datetime import date, datetime, timedelta, timezone
from itertools import batched
from urllib.parse import urlparse

import bson
import requests
from config import GEO_API_KEY, GEO_EMAIL, logger

from .cache import SqliteCache
from .common import as_list, dict_entries, sqlite

DB_DIR = "/data/nde-hub/standardizers/make_data_count"
DB_PATH = os.path.join(DB_DIR, "mdc.db")

# 202610 was built from a stale Europe PMC accession input (MDC methods, §13).
MIN_RELEASE = "202611"
# Release columns the loader reads.
_COLUMNS = frozenset({"publication_doi", "dataset_id", "dataset_id_type", "repository", "citation_type"})

# Accession families matched against NDE in nde-portal#569. Other families are not loaded.
ACCESSION_REPOSITORIES = frozenset(
    {"GEO", "PDB", "BioProject", "PRIDE", "ArrayExpress", "dbGaP", "MetaboLights", "EMPIAR", "BioStudies"}
)
# GenBank/ENA/DDBJ rows are kept for SRA studies and runs only.
SRA_REPOSITORY = "GenBank/ENA/DDBJ"
SRA_ACCESSION = re.compile(r"[SED]R[PXRS]\d+", re.IGNORECASE)

IDCONV_URL = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CROSSREF_URL = "https://api.crossref.org/works"
_TOOL = "nde-hub"

_DOI_PREFIX = re.compile(r"^(?:doi:\s*|https?://(?:dx\.)?doi\.org/)", re.IGNORECASE)
_DOI = re.compile(r"10\.\d{4,9}/\S+")
# Version suffixes, stripped on both sides so any version of a dataset matches.
_VERSIONED = (
    re.compile(r"(10\.17632/[a-z0-9]+)\.\d+"),  # Mendeley Data
    re.compile(r"(10\.\d{4,9}/\S+)\.v\d+"),  # Figshare and others
    re.compile(r"(phs\d+)\.v\d+\.p\d+"),  # dbGaP
    re.compile(r"(10\.7303/syn\d+)\.\d+"),  # Synapse
)
_PMID_IDENTIFIER = re.compile(r"pmid:\s*(\d+)", re.IGNORECASE)
_PUBMED_URL = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", re.IGNORECASE)
_PMCID = re.compile(r"\bPMC\d+\b", re.IGNORECASE)
# Crossref titles carry JATS/HTML markup and entities.
_MARKUP = re.compile(r"<[^>]+>")
_SPACES = re.compile(r"\s+")
_MONTHS = {
    name: number
    for number, name in enumerate(
        ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"), 1
    )
}

_IDCONV_SIZE = 200
_ESEARCH_SIZE = 100
_ESUMMARY_SIZE = 200
_CROSSREF_SIZE = 100
# NCBI allows 10 requests a second with an API key, 3 without. Crossref's polite
# pool (with mailto) allows 3, its public pool 1.
_NCBI_PAUSE = 0.11 if GEO_API_KEY else 0.34
_CROSSREF_PAUSE = 0.4 if GEO_EMAIL else 1.0

_MISS_DAYS = 90
_ATTEMPTS = 5
_RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
_RETRY_SECONDS = 5
_MAX_RETRY_AFTER = 120
_TIMEOUT = (10, 120)

_SQL_CHUNK_SIZE = 900
_INSERT_BATCH = 50_000
_BUSY_TIMEOUT = "PRAGMA busy_timeout = 60000"

# A record gaining more reuse papers than this is size-checked. Above the limit
# its MDC citedBy entries are dropped, so finalize does not drop the record.
_SIZE_CHECK_ENTRIES = 300
_SIZE_LIMIT = 15 * 1024 * 1024


class _Cache(SqliteCache):
    """A SqliteCache in the release database, waiting out the loader's writes."""

    def _connect(self):
        return sqlite(self.db_path, _BUSY_TIMEOUT, self._ddl)


# DOI -> {"pmid", "pmcid"}
PUBMED_IDS = _Cache(DB_PATH, "pubmed_ids", "doi", "ids", memoize=False)
# PMID -> description
PUBMED = _Cache(DB_PATH, "pubmed_descriptions", "pmid", "description", memoize=False)
# DOI -> description
CROSSREF = _Cache(DB_PATH, "crossref_descriptions", "doi", "description", memoize=False)
# DOI -> {"checked": ISO date}
PUBMED_MISSES = _Cache(DB_PATH, "pubmed_misses", "doi", "checked", memoize=False)
CROSSREF_MISSES = _Cache(DB_PATH, "crossref_misses", "doi", "checked", memoize=False)
_CACHES = (PUBMED_IDS, PUBMED, CROSSREF, PUBMED_MISSES, CROSSREF_MISSES)

# {"release", "created"} pinned for the current upload by reset_caches().
_release = None


def _connect():
    return sqlite(DB_PATH, _BUSY_TIMEOUT)


# ---------------------------------------------------------------------------
# Releases
# ---------------------------------------------------------------------------
def current_release():
    """{"release", "created"} of the newest loaded release, or None."""
    if not os.path.exists(DB_PATH):
        return None
    with _connect() as conn:
        if not conn.execute("SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'releases'").fetchone():
            return None
        row = conn.execute("SELECT release, created FROM releases ORDER BY release DESC LIMIT 1").fetchone()
    return {"release": row[0], "created": row[1]} if row else None


def lookup_file():
    """The release database, once a release has been loaded into it."""
    return DB_PATH if current_release() else None


def reset_caches():
    """Pin the newest release for this upload."""
    global _release
    _release = current_release()
    for cache in _CACHES:
        cache.reset()


def load_release(rows, release, created, row_count):
    """Load one release into `links_<release>` and make it the newest.

    `rows` are the release's CSV rows as dicts. The release is recorded only
    when `row_count` rows were read and every kept row landed. The previous
    release's table is kept for uploads still pinned to it.
    """
    if not re.fullmatch(r"\d{6}", release):
        raise ValueError(f"Unexpected Make Data Count release: {release!r}")
    os.makedirs(DB_DIR, exist_ok=True)
    table = f"links_{release}"
    staging = f"{table}_staging"
    counts = Counter()

    conn = sqlite3.connect(DB_PATH, timeout=60)
    try:
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute(
            "CREATE TABLE IF NOT EXISTS releases "
            "(release TEXT PRIMARY KEY, created TEXT, row_count INTEGER, loaded_at TEXT)"
        )
        conn.execute(f"DROP TABLE IF EXISTS {staging}")
        conn.execute(
            f"CREATE TABLE {staging} (dataset_key TEXT, publication_doi TEXT, citation_type TEXT, repository TEXT)"
        )
        conn.commit()

        for chunk in batched(_kept_links(rows, counts), _INSERT_BATCH):
            conn.executemany(f"INSERT INTO {staging} VALUES (?, ?, ?, ?)", chunk)
            conn.commit()
        if counts["read"] != row_count:
            raise RuntimeError(f"Read {counts['read']:,} rows of release {release}, expected {row_count:,}")
        if not counts["kept"]:
            raise RuntimeError(f"No usable links in release {release}")
        conn.execute(f"CREATE INDEX {table}_key_{int(time.time())} ON {staging} (dataset_key)")
        stored = conn.execute(f"SELECT COUNT(*) FROM {staging}").fetchone()[0]
        if stored != counts["kept"]:
            raise RuntimeError(f"Stored {stored:,} links of release {release}, expected {counts['kept']:,}")

        conn.execute("BEGIN IMMEDIATE")
        conn.execute(f"DROP TABLE IF EXISTS {table}")
        conn.execute(f"ALTER TABLE {staging} RENAME TO {table}")
        conn.execute(
            "INSERT OR REPLACE INTO releases VALUES (?, ?, ?, ?)",
            (release, created, row_count, datetime.now(timezone.utc).isoformat(timespec="seconds")),
        )
        conn.commit()

        # Keep this release and the one before it.
        older = [row[0] for row in conn.execute("SELECT release FROM releases ORDER BY release DESC")][2:]
        for old in older:
            conn.execute(f"DROP TABLE IF EXISTS links_{old}")
        conn.commit()
    except Exception:
        conn.rollback()
        conn.execute(f"DROP TABLE IF EXISTS {staging}")
        conn.commit()
        raise
    finally:
        conn.close()

    logger.info(
        "Loaded Make Data Count release %s (created %s): rows=%s links_kept=%s",
        release,
        created,
        counts["read"],
        counts["kept"],
    )


def _kept_links(rows, counts):
    """(dataset_key, publication_doi, citation_type, repository) of the rows the stage uses."""
    for row in rows:
        if not counts["read"] and (missing := _COLUMNS - row.keys()):
            raise ValueError(f"Make Data Count release lacks columns {sorted(missing)}; has {sorted(row)}")
        counts["read"] += 1
        citation_type = (row.get("citation_type") or "").capitalize()
        if citation_type not in ("Primary", "Secondary"):
            continue
        repository = row.get("repository") or ""
        dataset_id = (row.get("dataset_id") or "").strip()
        if (
            (row.get("dataset_id_type") or "").casefold() != "doi"
            and repository not in ACCESSION_REPOSITORIES
            and not (repository == SRA_REPOSITORY and SRA_ACCESSION.fullmatch(dataset_id))
        ):
            continue
        key = dataset_key(dataset_id)
        publication = normalize_id(row.get("publication_doi"))
        if not key or not publication or version_base(publication) == key:
            continue
        counts["kept"] += 1
        yield key, publication, citation_type, repository


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------
def normalize_id(value):
    """`value` casefolded, without a `doi:` or doi.org prefix; None if empty."""
    if isinstance(value, int) and not isinstance(value, bool):
        value = str(value)
    if not isinstance(value, str):
        return None
    return _DOI_PREFIX.sub("", value.strip()).casefold() or None


def version_base(key):
    """`key` without a known version suffix."""
    for pattern in _VERSIONED:
        if match := pattern.fullmatch(key):
            return match.group(1)
    return key


def dataset_key(value):
    """The lookup key of a dataset ID, as stored in `links_<release>`."""
    key = normalize_id(value)
    return version_base(key) if key else None


def record_keys(doc):
    """Lookup keys of a record: its `_id`, identifiers, DOIs and SRA accession URL."""
    keys = {dataset_key(value) for value in [doc.get("_id"), *as_list(doc.get("identifier")), *as_list(doc.get("doi"))]}
    for url in as_list(doc.get("url")):
        if isinstance(url, str):
            segment = urlparse(url).path.rstrip("/").rsplit("/", 1)[-1]
            if SRA_ACCESSION.fullmatch(segment):
                keys.add(segment.casefold())
    keys.discard(None)
    return keys


def _links(keys):
    """{dataset_key: [(publication_doi, citation_type, repository)]} in the pinned release."""
    table = f"links_{_release['release']}"
    links = {}
    with _connect() as conn:
        for chunk in batched(sorted(keys), _SQL_CHUNK_SIZE):
            placeholders = ",".join("?" for _ in chunk)
            for key, doi, citation_type, repository in conn.execute(
                f"SELECT dataset_key, publication_doi, citation_type, repository FROM {table} "
                f"WHERE dataset_key IN ({placeholders})",
                chunk,
            ):
                links.setdefault(key, []).append((doi, citation_type, repository))
    return links


# ---------------------------------------------------------------------------
# Papers the record already lists
# ---------------------------------------------------------------------------
def _as_doi(value):
    doi = normalize_id(value)
    return doi if doi and _DOI.fullmatch(doi) else None


def _pmcid(value):
    pmcid = str(value).strip().upper()
    return f"PMC{pmcid}" if pmcid.isdigit() else pmcid


def _known_papers(doc):
    """DOIs, PMIDs and PMCIDs of the papers in the record's `citation` and `citedBy`."""
    dois, pmids, pmcids = set(), set(), set()
    for field in ("citation", "citedBy"):
        for entry in dict_entries(doc, field):
            pmids.update(str(pmid).strip().lstrip("0") for pmid in as_list(entry.get("pmid")))
            pmcids.update(_pmcid(pmcid) for pmcid in as_list(entry.get("pmcid")))
            for value in as_list(entry.get("doi")) + as_list(entry.get("identifier")) + as_list(entry.get("url")):
                if not isinstance(value, str):
                    continue
                if doi := _as_doi(value):
                    dois.add(doi)
                if match := _PMID_IDENTIFIER.fullmatch(value.strip()) or _PUBMED_URL.search(value):
                    pmids.add(match.group(1).lstrip("0"))
                if match := _PMCID.search(value):
                    pmcids.add(match.group(0).upper())
    return dois, pmids, pmcids


def _has_citation(doc):
    """True when `citation` holds an entry with anything besides its `@type`."""
    return any(any(value for key, value in entry.items() if key != "@type") for entry in dict_entries(doc, "citation"))


# ---------------------------------------------------------------------------
# Describing citing papers
# ---------------------------------------------------------------------------
def _request_json(method, url, **kwargs):
    """JSON from `url`, retrying connection errors and 429/5xx (honoring Retry-After)."""
    for attempt in range(1, _ATTEMPTS + 1):
        try:
            response = requests.request(method, url, timeout=_TIMEOUT, **kwargs)
        except (requests.ConnectionError, requests.Timeout) as error:
            if attempt == _ATTEMPTS:
                raise
            logger.warning("%s failed (%s); retrying in %ss", url, type(error).__name__, _RETRY_SECONDS * attempt)
            time.sleep(_RETRY_SECONDS * attempt)
            continue
        if response.status_code in _RETRY_STATUSES and attempt < _ATTEMPTS:
            delay = _retry_after(response, _RETRY_SECONDS * attempt)
            logger.warning("%s returned HTTP %s; retrying in %ss", url, response.status_code, delay)
            time.sleep(delay)
            continue
        response.raise_for_status()
        return response.json()
    raise RuntimeError("Request retry loop ended unexpectedly")


def _retry_after(response, default):
    try:
        return min(int(response.headers["Retry-After"]), _MAX_RETRY_AFTER)
    except (KeyError, ValueError):
        return default


def _ncbi_params():
    params = {"tool": _TOOL}
    if GEO_EMAIL:
        params["email"] = GEO_EMAIL
    if GEO_API_KEY:
        params["api_key"] = GEO_API_KEY
    return params


def _text(value):
    """`value` without markup or HTML entities, on one line."""
    return _SPACES.sub(" ", html.unescape(_MARKUP.sub("", value))).strip()


def _partial_date(year, month=None, day=None):
    """ISO date of only the parts given, or None if they are not a real date."""
    try:
        date(year, month or 1, day or 1)
    except (TypeError, ValueError):
        return None
    if month is None:
        return f"{year:04d}"
    if day is None:
        return f"{year:04d}-{month:02d}"
    return f"{year:04d}-{month:02d}-{day:02d}"


def _pubmed_date(pubdate):
    """The parts of a PubMed date it states: '2016 Oct 14', '2019 Mar', '2019'.

    Seasons and ranges keep only the year. Nothing is filled in.
    """
    parts = pubdate.split()
    if not parts or not re.fullmatch(r"\d{4}", parts[0]):
        return None
    year = int(parts[0])
    month = _MONTHS.get(parts[1].casefold()) if len(parts) > 1 else None
    if not month:
        return _partial_date(year)
    if len(parts) > 2 and parts[2].isdigit():
        return _partial_date(year, month, int(parts[2])) or _partial_date(year, month)
    return _partial_date(year, month)


def _pubmed_description(summary):
    """pmid, doi, pmcid, name, journalName, datePublished and author of one ESummary record."""
    if not isinstance(summary, dict) or summary.get("error") or not summary.get("uid"):
        return None
    ids = {item.get("idtype"): item.get("value") for item in summary.get("articleids", []) if isinstance(item, dict)}
    description = {"pmid": str(summary["uid"])}
    if doi := ids.get("doi"):
        description["doi"] = doi.casefold()
    if pmcid := ids.get("pmc"):
        description["pmcid"] = _pmcid(pmcid)
    if title := _text(summary.get("title") or ""):
        description["name"] = title
    if journal := _text(summary.get("fulljournalname") or summary.get("source") or ""):
        description["journalName"] = journal
    if published := _pubmed_date(summary.get("pubdate") or ""):
        description["datePublished"] = published
    authors = [
        {"@type": "Organization" if author.get("authtype") == "CollectiveName" else "Person", "name": author["name"]}
        for author in summary.get("authors", [])
        if isinstance(author, dict) and author.get("name")
    ]
    if authors:
        description["author"] = authors
    return description


def _crossref_description(work):
    """name, journalName and datePublished of one Crossref work."""
    description = {}
    if title := _text((work.get("title") or [""])[0]):
        description["name"] = title
    if journal := _text((work.get("container-title") or [""])[0]):
        description["journalName"] = journal
    parts = []
    for part in ((work.get("issued") or {}).get("date-parts") or [[]])[0][:3]:
        if not isinstance(part, int):
            break
        parts.append(part)
    if parts and (published := _partial_date(*parts)):
        description["datePublished"] = published
    return description


def _pmc_ids(dois):
    """{doi: {"pmid", "pmcid"}} for the DOIs PubMed Central knows, and the DOIs whose request failed."""
    found, failed = {}, set()
    for chunk in batched(dois, _IDCONV_SIZE):
        params = {"format": "json", "tool": _TOOL, "ids": ",".join(chunk)}
        if GEO_EMAIL:
            params["email"] = GEO_EMAIL
        try:
            records = _request_json("GET", IDCONV_URL, params=params).get("records", [])
        except Exception as error:
            logger.warning("PMC ID converter failed for %s DOIs: %s", len(chunk), error)
            failed.update(chunk)
            continue
        for record in records:
            if record.get("pmid") and record.get("doi"):
                pmcid = _pmcid(record["pmcid"]) if record.get("pmcid") else None
                found[record["doi"].casefold()] = {"pmid": str(record["pmid"]), "pmcid": pmcid}
        time.sleep(_NCBI_PAUSE)
    return found, failed


def _esearch_pmids(dois):
    """PMIDs PubMed finds for `dois` searched as [DOI], and the DOIs whose request failed."""
    pmids, failed = set(), set()
    for chunk in batched([doi for doi in dois if '"' not in doi], _ESEARCH_SIZE):
        term = " OR ".join(f'"{doi}"[DOI]' for doi in chunk)
        data = {"db": "pubmed", "term": term, "retmax": 2 * len(chunk), "retmode": "json", **_ncbi_params()}
        try:
            result = _request_json("POST", EUTILS_URL + "esearch.fcgi", data=data).get("esearchresult", {})
        except Exception as error:
            logger.warning("PubMed ESearch failed for %s DOIs: %s", len(chunk), error)
            failed.update(chunk)
            continue
        pmids.update(result.get("idlist", []))
        time.sleep(_NCBI_PAUSE)
    return pmids, failed


def _pubmed_descriptions(pmids):
    """{pmid: description} from ESummary, and the PMIDs whose request failed."""
    found, failed = {}, set()
    for chunk in batched(sorted(pmids), _ESUMMARY_SIZE):
        data = {"db": "pubmed", "id": ",".join(chunk), "retmode": "json", **_ncbi_params()}
        try:
            result = _request_json("POST", EUTILS_URL + "esummary.fcgi", data=data).get("result", {})
        except Exception as error:
            logger.warning("PubMed ESummary failed for %s PMIDs: %s", len(chunk), error)
            failed.update(chunk)
            continue
        for pmid in result.get("uids", []):
            if description := _pubmed_description(result.get(pmid)):
                found[description["pmid"]] = description
        time.sleep(_NCBI_PAUSE)
    return found, failed


def _crossref_descriptions(dois):
    """{doi: description} from Crossref, the DOIs it lacks, and the DOIs whose request failed."""
    found, absent, failed = {}, set(), set()
    # Commas would split a DOI inside the filter.
    for chunk in batched([doi for doi in dois if "," not in doi], _CROSSREF_SIZE):
        params = {
            "filter": ",".join(f"doi:{doi}" for doi in chunk),
            "rows": len(chunk),
            "select": "DOI,title,container-title,issued",
        }
        if GEO_EMAIL:
            params["mailto"] = GEO_EMAIL
        try:
            message = _request_json("GET", CROSSREF_URL, params=params).get("message", {})
        except Exception as error:
            logger.warning("Crossref failed for %s DOIs: %s", len(chunk), error)
            failed.update(chunk)
            continue
        for work in message.get("items", []):
            if doi := str(work.get("DOI") or "").casefold():
                found[doi] = _crossref_description(work)
        # Only a complete page shows the rest are absent.
        if message.get("total-results", 0) <= len(chunk):
            absent.update(doi for doi in chunk if doi not in found)
        time.sleep(_CROSSREF_PAUSE)
    return found, absent, failed


def _describe(dois):
    """{doi: description} for the papers PubMed or Crossref describe, cache first.

    Also returns the DOIs left undescribed because a request failed. Those are
    not cached as misses, so the next upload asks again.
    """
    dois = sorted(set(dois))
    if not dois:
        return {}, set()
    today = date.today().isoformat()
    fresh_after = (date.today() - timedelta(days=_MISS_DAYS)).isoformat()
    failed = set()

    ids = PUBMED_IDS.get_many(dois)
    pubmed_misses = {doi for doi, miss in PUBMED_MISSES.get_many(dois).items() if miss["checked"] > fresh_after}
    unknown = [doi for doi in dois if doi not in ids and doi not in pubmed_misses]
    if unknown:
        found, uncertain = _pmc_ids(unknown)
        rest = [doi for doi in unknown if doi not in found and doi not in uncertain]
        pmids, search_failed = _esearch_pmids(rest)
        searched, summary_failed = _pubmed_descriptions(pmids)
        PUBMED.put_many(searched)
        rest_set = set(rest)
        for description in searched.values():
            if description.get("doi") in rest_set:
                found[description["doi"]] = {"pmid": description["pmid"], "pmcid": description.get("pmcid")}
        # Unmapped PMIDs can't be told apart, so a failed ESummary leaves the whole search undecided.
        uncertain |= search_failed | (rest_set if summary_failed else set())
        missing = [doi for doi in unknown if doi not in found]
        PUBMED_IDS.put_many(found)
        PUBMED_MISSES.put_many({doi: {"checked": today} for doi in missing if doi not in uncertain})
        ids.update(found)
        pubmed_misses.update(doi for doi in missing if doi not in uncertain)
        failed.update(doi for doi in missing if doi in uncertain)

    descriptions = PUBMED.get_many({entry["pmid"] for entry in ids.values()})
    if undescribed := {entry["pmid"] for entry in ids.values()} - descriptions.keys():
        fetched, pmids_failed = _pubmed_descriptions(undescribed)
        PUBMED.put_many(fetched)
        descriptions.update(fetched)
        failed.update(doi for doi, entry in ids.items() if entry["pmid"] in pmids_failed)
    described = {}
    for doi, entry in ids.items():
        if description := descriptions.get(entry["pmid"]):
            described[doi] = {**description, "pmcid": entry.get("pmcid") or description.get("pmcid")}

    crossref_dois = [doi for doi in dois if doi in pubmed_misses and doi not in ids]
    works = CROSSREF.get_many(crossref_dois)
    crossref_misses = {
        doi for doi, miss in CROSSREF_MISSES.get_many(crossref_dois).items() if miss["checked"] > fresh_after
    }
    ask = [doi for doi in crossref_dois if doi not in works and doi not in crossref_misses]
    if ask:
        fetched, absent, crossref_failed = _crossref_descriptions(ask)
        CROSSREF.put_many(fetched)
        CROSSREF_MISSES.put_many({doi: {"checked": today} for doi in absent})
        works.update(fetched)
        failed |= crossref_failed
    described.update(works)

    logger.info(
        "Make Data Count papers: dois=%s pubmed=%s crossref=%s undescribed=%s failed=%s",
        len(dois),
        len(ids),
        len(works),
        len(dois) - len(described),
        len(failed),
    )
    return described, failed


# ---------------------------------------------------------------------------
# The pipeline stage
# ---------------------------------------------------------------------------
def _paper(doi, description, with_author):
    """A citation or citedBy entry for one MDC paper; a new dict every call."""
    description = description or {}
    entry = {"@type": "ScholarlyArticle", "doi": doi}
    if pmid := description.get("pmid"):
        entry["identifier"] = f"PMID:{pmid}"
        entry["pmid"] = pmid
        entry["url"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    else:
        entry["url"] = f"https://doi.org/{doi}"
    for key in ("name", "journalName", "datePublished"):
        if value := description.get(key):
            entry[key] = value
    if with_author and description.get("author"):
        entry["author"] = [dict(author) for author in description["author"]]
    entry["fromMDC"] = True
    return entry


def _reuse_counter(count):
    """The InteractionCounter for `count` reuse papers in the pinned release."""
    counter = {
        "@type": "InteractionCounter",
        "name": "Data Reuse Citations",
        "userInteractionCount": count,
        "interactionType": {
            "@type": "Action",
            "name": "Cites as Data Source",
            "url": "http://purl.org/spar/cito/citesAsDataSource",
        },
        "owner": {"@type": "Organization", "name": "Make Data Count", "parentOrganization": "DataCite"},
        "fromMDC": True,
    }
    if _release.get("created"):
        counter["endTime"] = _release["created"]
    return counter


def _attach(doc, papers, described, failed, stats):
    """Add one record's MDC papers and reuse counter, skipping papers it already lists."""
    dois, pmids, pmcids = _known_papers(doc)
    has_citation = _has_citation(doc)
    citations, cited_by = [], []
    for doi, citation_type in papers.items():
        description = described.get(doi)
        listed = doi in dois or bool(
            description and (description.get("pmid") in pmids or description.get("pmcid") in pmcids)
        )
        if citation_type == "Primary" and has_citation:
            stats["primary_agree" if listed else "primary_disagree"] += 1
            continue
        if listed:
            continue
        if description is None:
            stats["failed_stubs" if doi in failed else "stubs"] += 1
        entries = citations if citation_type == "Primary" else cited_by
        entries.append(_paper(doi, description, with_author=citation_type == "Primary"))

    if citations:
        doc["citation"] = as_list(doc.get("citation")) + citations
        stats["citation_added"] += len(citations)
    if cited_by:
        cited_by.sort(key=lambda entry: (entry.get("datePublished", ""), entry["doi"]), reverse=True)
        doc["citedBy"] = as_list(doc.get("citedBy")) + cited_by
        stats["citedBy_added"] += len(cited_by)
    if reuse := sum(citation_type == "Secondary" for citation_type in papers.values()):
        doc["interactionStatistic"] = as_list(doc.get("interactionStatistic")) + [_reuse_counter(reuse)]
        stats["counters"] += 1

    if len(cited_by) > _SIZE_CHECK_ENTRIES and len(bson.BSON.encode(doc)) > _SIZE_LIMIT:
        added = {id(entry) for entry in cited_by}
        kept = [entry for entry in doc["citedBy"] if id(entry) not in added]
        if kept:
            doc["citedBy"] = kept
        else:
            doc.pop("citedBy")
        stats["size_trimmed"] += 1
        logger.warning(
            "%s is over %s bytes with %s Make Data Count citedBy entries; kept only its counter",
            doc.get("_id"),
            _SIZE_LIMIT,
            len(cited_by),
        )


def add_data_citations(docs):
    """Add MDC papers and reuse counters to the Dataset records in `docs`, in place.

    `docs` is one pipeline batch; every record in it is looked up at once.
    """
    docs = list(docs)
    if _release is None:
        return docs
    keys_by_doc = [record_keys(doc) for doc in docs]
    links = _links(set().union(*keys_by_doc))

    stats = Counter()
    papers_by_doc = []
    for keys in keys_by_doc:
        papers = {}
        repositories = set()
        for key in keys:
            for doi, citation_type, repository in links.get(key, ()):
                # A link from the dataset to itself under another version.
                if version_base(doi) in keys:
                    continue
                if papers.get(doi) != "Primary":
                    papers[doi] = citation_type
                repositories.add("DOI" if key.startswith("10.") else repository)
        stats.update(f"matched_{repository}" for repository in repositories)
        papers_by_doc.append(papers)

    described, failed = _describe(doi for papers in papers_by_doc for doi in papers)
    for doc, papers in zip(docs, papers_by_doc):
        if papers:
            stats["records_matched"] += 1
            _attach(doc, papers, described, failed, stats)

    logger.info(
        "Make Data Count %s: records=%s %s",
        _release["release"],
        len(docs),
        " ".join(f"{key}={value}" for key, value in sorted(stats.items())) or "no matches",
    )
    return docs
