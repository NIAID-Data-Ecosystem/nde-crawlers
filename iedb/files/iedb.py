#!/usr/bin/env python3
"""IEDB species/strain-by-tab DataCollection crawler for the NDE."""

import copy
import datetime
import json
import logging
import os
import re
import sqlite3
import tarfile
import tempfile
import time
import urllib.parse
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional

import dateutil.parser
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger("nde-logger")


def insert_value(d, key, value, extend=False):
    """Insert a value into a dictionary, promoting to a list on repeats and dropping duplicates."""

    if extend:
        d[key] = (d[key] + " " + value).strip() if d.get(key) else value
        return

    was_list = isinstance(d.get(key), list)
    merged = list(d[key]) if was_list else ([d[key]] if key in d else [])
    for item in value if isinstance(value, list) else [value]:
        if item not in merged:
            merged.append(item)
    d[key] = merged if was_list or isinstance(value, list) or len(merged) > 1 else merged[0]


def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(str(val), ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError, OverflowError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt


EXPORT_PAGE = "https://www.iedb.org/database_export_v3.php"
EXPORT_URL = "https://www.iedb.org/downloader.php?file_name=doc/iedb_export.zip"
ORGANISM_LIST_URL = "https://www.iedb.org/downloader.php?file_name=doc/OrganismList.zip"
ASSAY_TYPE_LIST_URL = "https://www.iedb.org/downloader.php?file_name=doc/AssayTypeList.zip"
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
QUERY_API_BASE = "https://query-api.iedb.org"
WEB_BASE = "https://www.iedb.org"
USER_AGENT = "nde-iedb-datacollection-crawler/0.1"

REQUEST_TIMEOUT = int(os.environ.get("IEDB_REQUEST_TIMEOUT", "180"))
REQUEST_INTERVAL = float(os.environ.get("IEDB_REQUEST_INTERVAL", "0.1"))
QUERY_PAGE_SIZE = int(os.environ.get("IEDB_QUERY_PAGE_SIZE", "1000"))
HAS_PART_LIMIT = int(os.environ.get("IEDB_HAS_PART_LIMIT", "1000"))
COMMIT_INTERVAL = int(os.environ.get("IEDB_COMMIT_INTERVAL", "50000"))

SOURCE_ORGANIZATION = {
    "@type": "Organization",
    "name": "Immune Epitope Database and Analysis Resource",
    "alternateName": "IEDB",
    "parentOrganization": "La Jolla Institute for Immunology",
    "url": "https://www.iedb.org/",
}

IEDB_OPERATOR = {
    "@type": "Organization",
    "name": "La Jolla Institute for Immunology",
    "url": "https://www.lji.org/",
}

CANONICAL_CITATION = {
    "@type": "ScholarlyArticle",
    "name": "Immune Epitope Database and Analysis Resource",
    "doi": "10.1093/nar/gkae1092",
    "pmid": "39558162",
    "url": "https://doi.org/10.1093/nar/gkae1092",
}

FUNDING = {
    "@type": "MonetaryGrant",
    "name": "IEDB support",
    "funder": {
        "@type": "Organization",
        "name": "National Institute of Allergy and Infectious Diseases",
        "alternateName": "NIAID",
        "url": "https://www.niaid.nih.gov/",
    },
}

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "name": "Immunology",
        "identifier": "topic_0804",
        "url": "http://edamontology.org/topic_0804",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Immunoinformatics",
        "identifier": "topic_3948",
        "url": "http://edamontology.org/topic_3948",
        "inDefinedTermSet": "EDAM",
    },
]

IMMUNE_EPITOPE_ASSAY = {
    "@type": "DefinedTerm",
    "name": "immune epitope assay",
    "identifier": "OBI:1110128",
    "url": "http://purl.obolibrary.org/obo/OBI_1110128",
    "inDefinedTermSet": "OBI",
}

IMMUNE_RESPONSE = {
    "@type": "DefinedTerm",
    "name": "immune response",
    "identifier": "GO:0006955",
    "url": "http://purl.obolibrary.org/obo/GO_0006955",
    "inDefinedTermSet": "GO",
}

TAB_CONFIG = {
    "epitopes": {
        "label": "epitopes",
        "unit": "epitopes",
        "about": {
            "@type": "DefinedTerm",
            "name": "Epitope",
            "url": "http://purl.obolibrary.org/obo/NCIT_C13189",
        },
        "endpoints": ("epitope_search",),
    },
    "antigens": {
        "label": "antigens",
        "unit": "antigens",
        "about": {"@type": "DefinedTerm", "name": "Antigen"},
        "endpoints": ("antigen_search",),
    },
    "assays": {
        "label": "assays",
        "unit": "assays",
        "about": {
            "@type": "DefinedTerm",
            "name": "immune epitope assay",
            "url": "http://purl.obolibrary.org/obo/OBI_1110128",
        },
        "endpoints": ("tcell_search", "bcell_search", "mhc_search"),
    },
    "receptors": {
        "label": "receptors",
        "unit": "receptors",
        "about": {"@type": "DefinedTerm", "name": "Immune receptor"},
        "endpoints": ("tcr_search", "bcr_search"),
    },
    "references": {
        "label": "references",
        "unit": "references",
        "about": {"@type": "DefinedTerm", "name": "IEDB source reference"},
        "endpoints": ("reference_search",),
    },
}

TAB_ORDER = {tab: index for index, tab in enumerate(TAB_CONFIG)}
ASSAY_TAGS = {
    "TCell": ("TCellId", "T-cell assay"),
    "BCell": ("BCellId", "B-cell assay"),
    "MhcBinding": ("MhcBindingId", "MHC binding assay"),
    "MhcLigandElution": ("MhcLigandElutionId", "MHC ligand-elution assay"),
}
MEMBERSHIP_ORGANISM_TAGS = {
    "SourceOrganismId",
    "SecondarySourceOrganismId",
    "PeptideSourceOrganismId",
    "MoleculeSourceOrganismId",
}
INVALID_XML_CONTROL_BYTES = re.compile(rb"[\x00-\x08\x0b\x0c\x0e-\x1f]")


class _SanitizedXMLStream:
    """Strip control bytes that are illegal in XML 1.0 from an archive member."""

    def __init__(self, source):
        self.source = source
        self.removed_bytes = 0

    def read(self, size: int = -1) -> bytes:
        data = self.source.read(size)
        cleaned, removed = INVALID_XML_CONTROL_BYTES.subn(b"", data)
        self.removed_bytes += removed
        return cleaned


def _is_empty(value: Any) -> bool:
    return value is None or value == "" or value == [] or value == {}


def _clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: cleaned
            for key, item in value.items()
            if not _is_empty(item) and not _is_empty(cleaned := _clean(item))
        }
    if isinstance(value, list):
        return [cleaned for item in value if not _is_empty(item) and not _is_empty(cleaned := _clean(item))]
    return value


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _element_text(element: Optional[ET.Element]) -> Optional[str]:
    if element is None or element.text is None:
        return None
    value = element.text.strip()
    return value or None


def _first_descendant_text(element: ET.Element, name: str) -> Optional[str]:
    for descendant in element.iter():
        if _local_name(descendant.tag) == name:
            return _element_text(descendant)
    return None


def _all_descendant_text(element: ET.Element, name: str) -> list[str]:
    values = []
    for descendant in element.iter():
        if _local_name(descendant.tag) == name and (value := _element_text(descendant)) and value not in values:
            values.append(value)
    return values


def _bool_env(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off"}


def _int_env(name: str) -> Optional[int]:
    value = os.environ.get(name)
    return int(value) if value else None


def _selected_organism_ids() -> Optional[set[str]]:
    value = os.environ.get("IEDB_ORGANISM_IDS")
    if not value:
        return None
    return {item for item in re.split(r"[\s,]+", value) if item}


def _get_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(
        total=6,
        connect=6,
        read=6,
        backoff_factor=2,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "HEAD"),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retries, pool_connections=4, pool_maxsize=4)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": USER_AGENT})
    session._nde_last_request = 0.0
    return session


def _wait_for_request_slot(session: requests.Session) -> None:
    elapsed = time.monotonic() - getattr(session, "_nde_last_request", 0.0)
    if REQUEST_INTERVAL > elapsed:
        time.sleep(REQUEST_INTERVAL - elapsed)


def _download(session: requests.Session, url: str, destination: Path) -> Path:
    logger.info("Downloading %s", url)
    _wait_for_request_slot(session)
    with session.get(url, stream=True, timeout=REQUEST_TIMEOUT) as response:
        session._nde_last_request = time.monotonic()
        response.raise_for_status()
        with destination.open("wb") as output:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    output.write(chunk)
    logger.info("Downloaded %s (%s bytes)", destination.name, destination.stat().st_size)
    return destination


def _resource_path(
    session: requests.Session,
    temp_dir: Path,
    env_name: str,
    url: str,
    filename: str,
) -> Path:
    if supplied := os.environ.get(env_name):
        path = Path(supplied)
        if not path.exists():
            raise FileNotFoundError(f"{env_name} does not exist: {path}")
        return path
    return _download(session, url, temp_dir / filename)


class _Store:
    def __init__(self, path: Path):
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self.connection.execute("PRAGMA temp_store=FILE")
        self._writes = 0
        self._create_schema()

    def _create_schema(self) -> None:
        self.connection.executescript(
            """
            CREATE TABLE used_organisms (
                organism_id TEXT PRIMARY KEY
            ) WITHOUT ROWID;
            CREATE TABLE organisms (
                organism_id TEXT PRIMARY KEY,
                tax_id TEXT,
                parent_tax_id TEXT,
                name TEXT,
                eligible INTEGER NOT NULL DEFAULT 0
            ) WITHOUT ROWID;
            CREATE TABLE collections (
                organism_id TEXT NOT NULL,
                tab TEXT NOT NULL,
                max_date TEXT,
                PRIMARY KEY (organism_id, tab)
            ) WITHOUT ROWID;
            CREATE TABLE members (
                organism_id TEXT NOT NULL,
                tab TEXT NOT NULL,
                member_key TEXT NOT NULL,
                member_json TEXT NOT NULL,
                PRIMARY KEY (organism_id, tab, member_key)
            ) WITHOUT ROWID;
            CREATE TABLE assay_terms (
                organism_id TEXT NOT NULL,
                assay_type_id TEXT NOT NULL,
                PRIMARY KEY (organism_id, assay_type_id)
            ) WITHOUT ROWID;
            CREATE TABLE health_conditions (
                organism_id TEXT NOT NULL,
                value TEXT NOT NULL,
                PRIMARY KEY (organism_id, value)
            ) WITHOUT ROWID;
            CREATE TABLE ncbi_taxonomy (
                tax_id TEXT PRIMARY KEY,
                parent_tax_id TEXT NOT NULL,
                rank TEXT NOT NULL
            ) WITHOUT ROWID;
            CREATE INDEX members_by_collection ON members (organism_id, tab);
            """
        )

    def _maybe_commit(self) -> None:
        self._writes += 1
        if self._writes % COMMIT_INTERVAL == 0:
            self.connection.commit()

    def add_event(self, event: dict[str, Any]) -> None:
        organism_id = str(event.get("organism_id") or "")
        tab = str(event.get("tab") or "")
        member = event.get("member") or {}
        member_key = str(event.get("member_key") or member.get("identifier") or "")
        if not organism_id or tab not in TAB_CONFIG or not member_key:
            return

        if organism := event.get("organism"):
            self.add_organism(
                organism_id,
                organism.get("tax_id"),
                organism.get("parent_tax_id"),
                organism.get("name"),
                bool(organism.get("eligible", True)),
            )

        date_modified = _to_iso_date(event.get("date_modified"))
        self.connection.execute(
            "INSERT OR IGNORE INTO used_organisms (organism_id) VALUES (?)",
            (organism_id,),
        )
        self.connection.execute(
            """
            INSERT INTO collections (organism_id, tab, max_date)
            VALUES (?, ?, ?)
            ON CONFLICT (organism_id, tab) DO UPDATE SET
                max_date = CASE
                    WHEN excluded.max_date > COALESCE(collections.max_date, '') THEN excluded.max_date
                    ELSE collections.max_date
                END
            """,
            (organism_id, tab, date_modified),
        )
        self.connection.execute(
            "INSERT OR IGNORE INTO members (organism_id, tab, member_key, member_json) VALUES (?, ?, ?, ?)",
            (organism_id, tab, member_key, json.dumps(_clean(member), sort_keys=True, separators=(",", ":"))),
        )
        if assay_type_id := event.get("assay_type_id"):
            self.connection.execute(
                "INSERT OR IGNORE INTO assay_terms (organism_id, assay_type_id) VALUES (?, ?)",
                (organism_id, str(assay_type_id)),
            )
        for condition in event.get("health_conditions") or []:
            if condition:
                self.connection.execute(
                    "INSERT OR IGNORE INTO health_conditions (organism_id, value) VALUES (?, ?)",
                    (organism_id, str(condition)),
                )
        self._maybe_commit()

    def add_organism(
        self,
        organism_id: str,
        tax_id: Optional[Any],
        parent_tax_id: Optional[Any],
        name: Optional[str],
        eligible: bool = False,
    ) -> None:
        self.connection.execute(
            """
            INSERT INTO organisms (organism_id, tax_id, parent_tax_id, name, eligible)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT (organism_id) DO UPDATE SET
                tax_id = COALESCE(excluded.tax_id, organisms.tax_id),
                parent_tax_id = COALESCE(excluded.parent_tax_id, organisms.parent_tax_id),
                name = COALESCE(excluded.name, organisms.name),
                eligible = MAX(excluded.eligible, organisms.eligible)
            """,
            (
                str(organism_id),
                str(tax_id) if tax_id not in (None, "") else None,
                str(parent_tax_id) if parent_tax_id not in (None, "") else None,
                name,
                int(eligible),
            ),
        )

    def used_organism_ids(self) -> set[str]:
        return {row[0] for row in self.connection.execute("SELECT organism_id FROM used_organisms")}

    def remove_ineligible(self) -> None:
        for table in ("collections", "members", "assay_terms", "health_conditions"):
            self.connection.execute(
                f"DELETE FROM {table} WHERE organism_id NOT IN "
                "(SELECT organism_id FROM organisms WHERE eligible = 1)"
            )
        self.connection.commit()

    def eligible_organism_maps(self) -> tuple[dict[str, list[str]], dict[str, str]]:
        by_tax_id: dict[str, list[str]] = {}
        by_iedb_id: dict[str, str] = {}
        for organism_id, tax_id in self.connection.execute(
            "SELECT organism_id, tax_id FROM organisms WHERE eligible = 1"
        ):
            by_iedb_id[organism_id] = organism_id
            if tax_id:
                by_tax_id.setdefault(tax_id, []).append(organism_id)
        return by_tax_id, by_iedb_id

    def close(self) -> None:
        self.connection.commit()
        self.connection.close()


def _source_organism_ids(structure: ET.Element) -> list[str]:
    values = []
    for descendant in structure.iter():
        if _local_name(descendant.tag) in MEMBERSHIP_ORGANISM_TAGS:
            value = _element_text(descendant)
            if value and value not in values:
                values.append(value)
    return values


def _structure_member_name(structure: ET.Element) -> Optional[str]:
    for field in ("LinearSequence", "PeptideLinearSequence", "PeptideSequence", "MoleculeName"):
        if value := _first_descendant_text(structure, field):
            return value
    return None


def _reference_member(reference: dict[str, Any]) -> Optional[dict[str, Any]]:
    reference_id = reference.get("reference_id")
    if not reference_id:
        return None
    member = {
        "@type": "ScholarlyArticle" if reference.get("reference_type") == "article" else "CreativeWork",
        "identifier": f"IEDB_REFERENCE:{reference_id}",
        "url": f"{WEB_BASE}/reference/{reference_id}",
        "name": reference.get("title") or f"IEDB reference {reference_id}",
        "pmid": reference.get("pmid"),
        "journalName": reference.get("journal"),
    }
    return _clean(member)


def _parse_article(article: ET.Element) -> dict[str, Any]:
    return {
        "reference_type": "article",
        "title": _first_descendant_text(article, "ArticleTitle"),
        "pmid": _first_descendant_text(article, "PubmedId"),
        "journal": _first_descendant_text(article, "Title"),
    }


def _parse_submission(submission: ET.Element) -> dict[str, Any]:
    return {
        "reference_type": "submission",
        "title": _first_descendant_text(submission, "SubmissionTitle"),
    }


def _assay_event(
    assay: ET.Element,
    tag: str,
    organism_id: str,
    reference: dict[str, Any],
) -> Optional[dict[str, Any]]:
    id_field, label = ASSAY_TAGS[tag]
    assay_id = _first_descendant_text(assay, id_field)
    if not assay_id:
        return None
    member = {
        "@type": "CreativeWork",
        "identifier": f"IEDB_ASSAY:{assay_id}",
        "name": f"IEDB {label} {assay_id}",
        "url": f"{WEB_BASE}/assay/{assay_id}",
        "additionalType": {
            "@type": "DefinedTerm",
            "name": label,
        },
    }
    return {
        "organism_id": organism_id,
        "tab": "assays",
        "member_key": f"IEDB_ASSAY:{assay_id}",
        "member": member,
        "date_modified": reference.get("date_modified"),
        "assay_type_id": _first_descendant_text(assay, "AssayTypeId"),
        "health_conditions": _all_descendant_text(assay, "DiseaseState"),
    }


def iter_iedb_records(
    archive_path: Path | str,
    selected_organisms: Optional[set[str]] = None,
    max_references: Optional[int] = None,
) -> Iterator[dict[str, Any]]:
    """Stream normalized member events from the reference-oriented XML archive."""
    with zipfile.ZipFile(archive_path) as archive:
        members = [info for info in archive.infolist() if info.filename.lower().endswith(".xml")]
        members.sort(key=lambda info: int(Path(info.filename).stem) if Path(info.filename).stem.isdigit() else 0)
        if max_references is not None:
            members = members[:max_references]

        for index, info in enumerate(members, 1):
            reference: dict[str, Any] = {}
            epitope: dict[str, Any] = {}
            stack: list[ET.Element] = []
            with archive.open(info) as raw_source:
                source = _SanitizedXMLStream(raw_source)
                parser = iter(ET.iterparse(source, events=("start", "end")))
                while True:
                    try:
                        action, element = next(parser)
                    except StopIteration:
                        break
                    except ET.ParseError as error:
                        logger.error(
                            "Stopped at malformed IEDB XML in %s after retaining parsed records: %s",
                            info.filename,
                            error,
                        )
                        break
                    tag = _local_name(element.tag)
                    if action == "start":
                        stack.append(element)
                        if tag == "Reference":
                            reference = {}
                        elif tag == "Epitope":
                            epitope = {"source_organism_ids": []}
                        continue

                    parent = stack[-2] if len(stack) > 1 else None
                    if tag == "ReferenceId" and not epitope:
                        reference["reference_id"] = _element_text(element)
                    elif tag == "DateLastUpdated" and not epitope:
                        reference["date_modified"] = _to_iso_date(_element_text(element))
                    elif tag == "Article":
                        reference.update(_parse_article(element))
                    elif tag == "Submission":
                        reference.update(_parse_submission(element))
                    elif tag == "EpitopeName":
                        epitope["name"] = _element_text(element)
                    elif tag == "EpitopeId":
                        epitope["epitope_id"] = _element_text(element)
                    elif tag == "EpitopeStructure":
                        epitope["source_organism_ids"] = _source_organism_ids(element)
                        epitope["structure_name"] = _structure_member_name(element)
                    elif tag in ASSAY_TAGS:
                        for organism_id in epitope.get("source_organism_ids") or []:
                            if selected_organisms is None or organism_id in selected_organisms:
                                if event := _assay_event(element, tag, organism_id, reference):
                                    yield event
                    elif tag == "Epitope":
                        epitope_id = epitope.get("epitope_id")
                        source_organism_ids = epitope.get("source_organism_ids") or []
                        if epitope_id and source_organism_ids:
                            epitope_member = _clean(
                                {
                                    "@type": "CreativeWork",
                                    "identifier": f"IEDB_EPITOPE:{epitope_id}",
                                    "name": epitope.get("name")
                                    or epitope.get("structure_name")
                                    or f"IEDB epitope {epitope_id}",
                                    "url": f"{WEB_BASE}/epitope/{epitope_id}",
                                    "additionalType": copy.deepcopy(TAB_CONFIG["epitopes"]["about"]),
                                }
                            )
                            reference_member = _reference_member(reference)
                            for organism_id in source_organism_ids:
                                if selected_organisms is not None and organism_id not in selected_organisms:
                                    continue
                                yield {
                                    "organism_id": organism_id,
                                    "tab": "epitopes",
                                    "member_key": f"IEDB_EPITOPE:{epitope_id}",
                                    "member": epitope_member,
                                    "date_modified": reference.get("date_modified"),
                                }
                                if reference_member:
                                    yield {
                                        "organism_id": organism_id,
                                        "tab": "references",
                                        "member_key": reference_member["identifier"],
                                        "member": reference_member,
                                        "date_modified": reference.get("date_modified"),
                                    }
                        epitope = {}

                    if tag in {
                        "Article",
                        "Submission",
                        "EpitopeStructure",
                        "TCell",
                        "BCell",
                        "MhcBinding",
                        "MhcLigandElution",
                        "Assays",
                        "Epitope",
                        "Epitopes",
                        "Reference",
                    }:
                        element.clear()
                        if parent is not None:
                            try:
                                parent.remove(element)
                            except ValueError:
                                pass
                    stack.pop()

                if source.removed_bytes:
                    logger.warning(
                        "Removed %s illegal XML control byte(s) from %s",
                        source.removed_bytes,
                        info.filename,
                    )

            if index % 500 == 0:
                logger.info("Parsed %s/%s IEDB reference XML files", index, len(members))


def _load_used_organisms(store: _Store, archive_path: Path) -> None:
    used = store.used_organism_ids()
    logger.info("Loading metadata for %s used IEDB organisms", len(used))
    with zipfile.ZipFile(archive_path) as archive:
        xml_members = [info for info in archive.infolist() if info.filename.lower().endswith(".xml")]
        if len(xml_members) != 1:
            raise RuntimeError(f"Expected one OrganismList XML member; found {len(xml_members)}")
        stack: list[ET.Element] = []
        with archive.open(xml_members[0]) as source:
            for action, element in ET.iterparse(source, events=("start", "end")):
                tag = _local_name(element.tag)
                if action == "start":
                    stack.append(element)
                    continue
                parent = stack[-2] if len(stack) > 1 else None
                if tag == "Organism":
                    organism_id = _first_descendant_text(element, "OrganismId")
                    if organism_id in used:
                        store.add_organism(
                            organism_id,
                            _first_descendant_text(element, "TaxId"),
                            _first_descendant_text(element, "ParentTaxId"),
                            _first_descendant_text(element, "OrganismName"),
                        )
                    element.clear()
                    if parent is not None:
                        try:
                            parent.remove(element)
                        except ValueError:
                            pass
                stack.pop()
    store.connection.commit()
    found = {row[0] for row in store.connection.execute("SELECT organism_id FROM organisms")}
    for organism_id in sorted(used - found):
        logger.warning("No OrganismList metadata for used OrganismId %s", organism_id)


def _load_ncbi_taxonomy(store: _Store, archive_path: Path) -> None:
    logger.info("Loading NCBI taxonomy ranks")
    with tarfile.open(archive_path, "r:gz") as archive:
        member = archive.getmember("nodes.dmp")
        source = archive.extractfile(member)
        if source is None:
            raise RuntimeError("nodes.dmp is missing from NCBI taxdump")
        batch = []
        for raw_line in source:
            fields = raw_line.decode("utf-8").split("\t|\t")
            if len(fields) < 3:
                continue
            batch.append((fields[0].strip(), fields[1].strip(), fields[2].strip()))
            if len(batch) >= 50000:
                store.connection.executemany(
                    "INSERT OR REPLACE INTO ncbi_taxonomy (tax_id, parent_tax_id, rank) VALUES (?, ?, ?)",
                    batch,
                )
                batch.clear()
        if batch:
            store.connection.executemany(
                "INSERT OR REPLACE INTO ncbi_taxonomy (tax_id, parent_tax_id, rank) VALUES (?, ?, ?)",
                batch,
            )
    store.connection.commit()


def _has_species_ancestor(store: _Store, tax_id: Optional[str], cache: dict[str, bool]) -> bool:
    if not tax_id:
        return False
    visited = []
    current = str(tax_id)
    result = False
    while current and current not in visited:
        if current in cache:
            result = cache[current]
            break
        visited.append(current)
        row = store.connection.execute(
            "SELECT parent_tax_id, rank FROM ncbi_taxonomy WHERE tax_id = ?",
            (current,),
        ).fetchone()
        if row is None:
            break
        parent_tax_id, rank = row
        if rank == "species":
            result = True
            break
        if parent_tax_id == current:
            break
        current = parent_tax_id
    for visited_tax_id in visited:
        cache[visited_tax_id] = result
    return result


def _mark_eligible_organisms(store: _Store) -> None:
    if _bool_env("IEDB_SKIP_RANK_FILTER", False):
        logger.warning("IEDB_SKIP_RANK_FILTER is enabled; retaining all used organism ranks")
        store.connection.execute("UPDATE organisms SET eligible = 1")
        store.connection.commit()
        return

    cache: dict[str, bool] = {}
    rows = list(store.connection.execute("SELECT organism_id, tax_id, parent_tax_id FROM organisms"))
    for organism_id, tax_id, parent_tax_id in rows:
        lineage_start = tax_id or parent_tax_id
        eligible = _has_species_ancestor(store, lineage_start, cache)
        store.connection.execute(
            "UPDATE organisms SET eligible = ? WHERE organism_id = ?",
            (int(eligible), organism_id),
        )
    store.connection.commit()
    included = store.connection.execute("SELECT COUNT(*) FROM organisms WHERE eligible = 1").fetchone()[0]
    logger.info("Retained %s/%s used organisms at species rank or below", included, len(rows))


def _load_assay_types(archive_path: Path) -> dict[str, dict[str, Optional[str]]]:
    lookup: dict[str, dict[str, Optional[str]]] = {}
    with zipfile.ZipFile(archive_path) as archive:
        xml_members = [info for info in archive.infolist() if info.filename.lower().endswith(".xml")]
        if len(xml_members) != 1:
            return lookup
        with archive.open(xml_members[0]) as source:
            for _, element in ET.iterparse(source, events=("end",)):
                if _local_name(element.tag) != "AssayType":
                    continue
                assay_type_id = _first_descendant_text(element, "AssayTypeId")
                if not assay_type_id:
                    continue
                child_name = next(
                    (_element_text(child) for child in list(element) if _local_name(child.tag) == "AssayType"),
                    None,
                )
                lookup[assay_type_id] = {
                    "name": child_name,
                    "response": _first_descendant_text(element, "Response"),
                    "category": _first_descendant_text(element, "Category"),
                }
                element.clear()
    return lookup


def _request_json(
    session: requests.Session,
    url: str,
    params: Optional[dict[str, Any]] = None,
) -> Any:
    _wait_for_request_slot(session)
    response = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
    session._nde_last_request = time.monotonic()
    response.raise_for_status()
    return response.json()


def _iter_query_api(
    session: requests.Session,
    endpoint: str,
    select: str,
    order: str,
) -> Iterator[dict[str, Any]]:
    offset = 0
    page = 0
    max_pages = _int_env("IEDB_QUERY_MAX_PAGES")
    while True:
        payload = _request_json(
            session,
            f"{QUERY_API_BASE}/{endpoint}",
            {"select": select, "order": order, "limit": QUERY_PAGE_SIZE, "offset": offset},
        )
        if not isinstance(payload, list):
            raise TypeError(f"Unexpected {endpoint} response: {type(payload).__name__}")
        yield from payload
        page += 1
        if len(payload) < QUERY_PAGE_SIZE or (max_pages is not None and page >= max_pages):
            return
        offset += QUERY_PAGE_SIZE
        if page % 50 == 0:
            logger.info("Fetched %s %s rows", offset, endpoint)


def _organism_ids_from_iris(
    iris: Any,
    by_tax_id: dict[str, list[str]],
    by_iedb_id: dict[str, str],
) -> list[str]:
    organism_ids = []
    for iri in _as_list(iris):
        value = str(iri or "")
        if value.startswith("NCBITaxon:"):
            candidates = by_tax_id.get(value.split(":", 1)[1], [])
        elif value.startswith("IEDB_ORGANISM:"):
            organism_id = value.split(":", 1)[1]
            candidates = [by_iedb_id[organism_id]] if organism_id in by_iedb_id else []
        else:
            candidates = []
        for organism_id in candidates:
            if organism_id not in organism_ids:
                organism_ids.append(organism_id)
    return organism_ids


def _supplement_antigens_and_receptors(store: _Store, session: requests.Session) -> None:
    if not _bool_env("IEDB_INCLUDE_QUERY_SUPPLEMENTS", True):
        logger.warning("Skipping Antigens and Receptors Query API supplements")
        return
    by_tax_id, by_iedb_id = store.eligible_organism_maps()

    for row in _iter_query_api(
        session,
        "antigen_search",
        "parent_source_antigen_id,parent_source_antigen_iri,parent_source_antigen_names,source_organism_iris",
        "parent_source_antigen_iri.asc",
    ):
        antigen_iri = row.get("parent_source_antigen_iri")
        if not antigen_iri:
            continue
        names = _as_list(row.get("parent_source_antigen_names"))
        member = _clean(
            {
                "@type": "CreativeWork",
                "identifier": antigen_iri,
                "name": next((str(name) for name in names if name), str(antigen_iri)),
                "url": row.get("parent_source_antigen_id"),
                "additionalType": copy.deepcopy(TAB_CONFIG["antigens"]["about"]),
            }
        )
        for organism_id in _organism_ids_from_iris(row.get("source_organism_iris"), by_tax_id, by_iedb_id):
            store.add_event(
                {
                    "organism_id": organism_id,
                    "tab": "antigens",
                    "member_key": str(antigen_iri),
                    "member": member,
                }
            )

    for endpoint, receptor_family in (("tcr_search", "TCR"), ("bcr_search", "BCR")):
        for row in _iter_query_api(
            session,
            endpoint,
            "receptor_group_id,receptor_group_iri,receptor_type,receptor_names,source_organism_iris",
            "receptor_group_id.asc",
        ):
            receptor_iri = row.get("receptor_group_iri")
            receptor_id = row.get("receptor_group_id")
            if not receptor_iri and not receptor_id:
                continue
            identifier = str(receptor_iri or f"IEDB_RECEPTOR:{receptor_id}")
            names = _as_list(row.get("receptor_names"))
            member = _clean(
                {
                    "@type": "CreativeWork",
                    "identifier": identifier,
                    "name": next((str(name) for name in names if name), f"IEDB receptor {receptor_id}"),
                    "url": f"{WEB_BASE}/receptor/{receptor_id}" if receptor_id else None,
                    "additionalType": {
                        "@type": "DefinedTerm",
                        "name": receptor_family,
                        "alternateName": row.get("receptor_type"),
                    },
                }
            )
            for organism_id in _organism_ids_from_iris(row.get("source_organism_iris"), by_tax_id, by_iedb_id):
                store.add_event(
                    {
                        "organism_id": organism_id,
                        "tab": "receptors",
                        "member_key": identifier,
                        "member": member,
                    }
                )
    store.connection.commit()


def _archive_release_date(archive_path: Path) -> str:
    with zipfile.ZipFile(archive_path) as archive:
        dates = [datetime.date(*info.date_time[:3]) for info in archive.infolist()]
    return max(dates).isoformat() if dates else datetime.date.today().isoformat()


def _query_filter(organism: dict[str, Any], endpoint: str) -> Optional[str]:
    tax_id = organism.get("tax_id")
    if not tax_id:
        return None
    field = (
        "source_organism_iri" if endpoint in {"tcell_search", "bcell_search", "mhc_search"} else "source_organism_iris"
    )
    operator = f"eq.NCBITaxon:{tax_id}" if field.endswith("_iri") else f"cs.{{NCBITaxon:{tax_id}}}"
    return f"{QUERY_API_BASE}/{endpoint}?{urllib.parse.urlencode({field: operator})}"


def _record_url(tab: str, organism: dict[str, Any]) -> str:
    endpoints = TAB_CONFIG[tab]["endpoints"]
    if len(endpoints) == 1 and (query_url := _query_filter(organism, endpoints[0])):
        return query_url
    return f"{WEB_BASE}/result_v3.php" if organism.get("tax_id") else EXPORT_PAGE


def _species_term(organism: dict[str, Any]) -> dict[str, Any]:
    tax_id = organism.get("tax_id")
    return _clean(
        {
            "@type": "DefinedTerm",
            "name": organism.get("name") or f"IEDB organism {organism['organism_id']}",
            "identifier": tax_id or f"IEDB_ORGANISM:{organism['organism_id']}",
            "url": f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={tax_id}" if tax_id else None,
            "inDefinedTermSet": "NCBI Taxonomy" if tax_id else "IEDB",
        }
    )


def _health_condition(value: str) -> dict[str, Any]:
    if match := re.fullmatch(r"DOID[:_](\d+)", value, flags=re.IGNORECASE):
        identifier = f"DOID:{match.group(1)}"
        return {
            "@type": "DefinedTerm",
            "identifier": identifier,
            "url": f"http://purl.obolibrary.org/obo/DOID_{match.group(1)}",
            "inDefinedTermSet": "DOID",
        }
    return {"@type": "DefinedTerm", "name": value}


def _assay_terms(
    store: _Store,
    organism_id: str,
    lookup: dict[str, dict[str, Optional[str]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    techniques = [copy.deepcopy(IMMUNE_EPITOPE_ASSAY)]
    variables = [copy.deepcopy(IMMUNE_RESPONSE)]
    categories = []
    for (assay_type_id,) in store.connection.execute(
        "SELECT assay_type_id FROM assay_terms WHERE organism_id = ? ORDER BY assay_type_id",
        (organism_id,),
    ):
        item = lookup.get(assay_type_id) or {}
        techniques.append(
            _clean(
                {
                    "@type": "DefinedTerm",
                    "identifier": f"IEDB_ASSAY_TYPE:{assay_type_id}",
                    "name": item.get("name") or f"IEDB assay type {assay_type_id}",
                    "inDefinedTermSet": "IEDB",
                }
            )
        )
        if response := item.get("response"):
            variables.append(
                {
                    "@type": "DefinedTerm",
                    "identifier": f"IEDB_ASSAY_RESPONSE:{assay_type_id}",
                    "name": response,
                    "inDefinedTermSet": "IEDB",
                }
            )
        if (category := item.get("category")) and category not in categories:
            categories.append(category)
    return techniques, variables, categories


def _provenance(tab: str, organism: dict[str, Any], release_date: str) -> list[dict[str, Any]]:
    query_works = []
    for endpoint in TAB_CONFIG[tab]["endpoints"]:
        if query_url := _query_filter(organism, endpoint):
            query_works.append(
                {
                    "@type": "CreativeWork",
                    "name": f"IEDB {endpoint} exact-organism query",
                    "url": query_url,
                }
            )
    return [
        {
            "@type": "Action",
            "name": "IEDB DataCollection generation process",
            "description": (
                f"Aggregated the IEDB {tab} result type by exact source OrganismId, "
                "retaining only species-rank organisms and their descendants."
            ),
            "actionProcess": {
                "@type": "HowTo",
                "step": [
                    "Stream the IEDB curation XML reference archive.",
                    "Group primary epitope-source associations by exact IEDB OrganismId.",
                    "Project the grouped result into the five IEDB result-tab entity types.",
                    "Resolve NCBI rank and retain species and strain-level collections.",
                ],
            },
        },
        {
            "@type": "ResourceCatalog",
            "name": "Immune Epitope Database and Analysis Resource",
            "url": f"{WEB_BASE}/",
            "datePublished": release_date,
        },
        {
            "@type": "CreativeWork",
            "name": "IEDB complete XML database export",
            "url": EXPORT_URL,
        },
        *query_works,
    ]


def _example_of_work(tab: str, member: dict[str, Any], query_url: str) -> dict[str, Any]:
    identifier = member.get("identifier") or member.get("url") or member.get("name")
    return _clean(
        {
            "@type": "CreativeWork",
            "about": copy.deepcopy(TAB_CONFIG[tab]["about"]),
            "encodingFormat": [
                {
                    "@type": "DefinedTerm",
                    "name": "XML",
                    "identifier": "application/xml",
                },
                {
                    "@type": "DefinedTerm",
                    "name": "JSON",
                    "identifier": "application/json",
                },
            ],
            "schemaVersion": "IEDB Curation XML v3",
            "potentialAction": {
                "@type": "SearchAction",
                "name": "Search this IEDB tab for the exact organism",
                "url": query_url,
            },
            "additionalProperty": {
                "@type": "PropertyValue",
                "name": "Representative member identifier",
                "value": identifier,
            },
        }
    )


def _build_documents(
    store: _Store,
    release_date: str,
    assay_type_lookup: Optional[dict[str, dict[str, Optional[str]]]] = None,
    has_part_limit: int = HAS_PART_LIMIT,
) -> Iterator[dict[str, Any]]:
    assay_type_lookup = assay_type_lookup or {}
    rows = list(
        store.connection.execute(
            """
            SELECT c.organism_id, c.tab, c.max_date, o.tax_id, o.parent_tax_id, o.name
            FROM collections c
            JOIN organisms o ON o.organism_id = c.organism_id
            WHERE o.eligible = 1
            """
        )
    )
    rows.sort(
        key=lambda row: (
            (0, int(row[0])) if row[0].isdigit() else (1, row[0]),
            TAB_ORDER[row[1]],
        )
    )

    for organism_id, tab, max_date, tax_id, parent_tax_id, organism_name in rows:
        organism = {
            "organism_id": organism_id,
            "tax_id": tax_id,
            "parent_tax_id": parent_tax_id,
            "name": organism_name,
        }
        config = TAB_CONFIG[tab]
        url = _record_url(tab, organism)
        member_count = store.connection.execute(
            "SELECT COUNT(*) FROM members WHERE organism_id = ? AND tab = ?",
            (organism_id, tab),
        ).fetchone()[0]
        if not member_count:
            continue
        representative_row = store.connection.execute(
            "SELECT member_json FROM members WHERE organism_id = ? AND tab = ? ORDER BY member_key LIMIT 1",
            (organism_id, tab),
        ).fetchone()
        representative = json.loads(representative_row[0])
        identifier = f"IEDB:{tab.upper()}:{organism_id}"
        name = f"IEDB {config['label']} from {organism_name or f'organism {organism_id}'}"
        output = {
            "@context": "http://schema.org/",
            "@type": "DataCollection",
            "_id": f"iedb_{tab}_{organism_id}",
            "identifier": identifier,
            "url": url,
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "Immune Epitope Database and Analysis Resource",
                "url": f"{WEB_BASE}/",
                "versionDate": release_date,
                "archivedAt": url,
            },
        }
        insert_value(output, "date", release_date)
        insert_value(output, "dateModified", max_date or release_date)
        insert_value(output, "name", name)
        insert_value(
            output,
            "description",
            (
                f"IEDB {config['label']} directly associated with {organism_name or f'IEDB organism {organism_id}'} "
                f"({f'NCBITaxon:{tax_id}' if tax_id else f'IEDB_ORGANISM:{organism_id}'}). "
                f"This collection contains {member_count:,} {config['unit']}."
            ),
            extend=True,
        )
        insert_value(
            output,
            "collectionSize",
            {"@type": "QuantitativeValue", "value": member_count, "unitText": config["unit"]},
        )
        insert_value(output, "about", copy.deepcopy(config["about"]))
        insert_value(output, "species", _species_term(organism))
        insert_value(output, "author", copy.deepcopy(IEDB_OPERATOR))
        insert_value(output, "creator", copy.deepcopy(IEDB_OPERATOR))
        insert_value(output, "sourceOrganization", copy.deepcopy(SOURCE_ORGANIZATION))
        insert_value(output, "citation", copy.deepcopy(CANONICAL_CITATION))
        insert_value(output, "conditionsOfAccess", "Open")
        insert_value(output, "isAccessibleForFree", True)
        insert_value(output, "license", "https://creativecommons.org/licenses/by/4.0/")
        insert_value(output, "funding", copy.deepcopy(FUNDING))
        insert_value(
            output,
            "keywords",
            ["IEDB", config["label"], organism_name] if organism_name else ["IEDB", config["label"]],
        )
        insert_value(output, "topicCategory", copy.deepcopy(TOPIC_CATEGORY))
        insert_value(
            output,
            "usageInfo",
            {
                "@type": "CreativeWork",
                "name": "How to cite the IEDB",
                "url": "https://discuss.iedb.org/t/citing-the-iedb/63",
            },
        )
        insert_value(output, "exampleOfWork", _example_of_work(tab, representative, url))
        insert_value(output, "isBasedOn", _provenance(tab, organism, release_date))

        if tab == "assays":
            techniques, variables, categories = _assay_terms(store, organism_id, assay_type_lookup)
            insert_value(output, "measurementTechnique", techniques)
            insert_value(output, "variableMeasured", variables)
            if categories:
                insert_value(output, "keywords", categories)
            conditions = [
                _health_condition(row[0])
                for row in store.connection.execute(
                    "SELECT value FROM health_conditions WHERE organism_id = ? ORDER BY value",
                    (organism_id,),
                )
            ]
            if conditions:
                insert_value(output, "healthCondition", conditions)

        if has_part_limit >= 0 and member_count <= has_part_limit:
            members = [
                json.loads(row[0])
                for row in store.connection.execute(
                    "SELECT member_json FROM members WHERE organism_id = ? AND tab = ? ORDER BY member_key",
                    (organism_id, tab),
                )
            ]
            insert_value(output, "hasPart", members)
        elif has_part_limit >= 0:
            logger.info(
                "Omitting hasPart for %s because %s members exceeds IEDB_HAS_PART_LIMIT=%s",
                output["_id"],
                member_count,
                has_part_limit,
            )

        yield _clean(output)


def _ingest_events(store: _Store, records: Iterable[dict[str, Any]]) -> None:
    for record in records:
        store.add_event(record)
    store.connection.commit()


def parse(
    records: Optional[Iterable[dict[str, Any]]] = None,
    *,
    release_date: Optional[str] = None,
    has_part_limit: int = HAS_PART_LIMIT,
) -> Iterator[dict[str, Any]]:
    """Yield one DataCollection for every non-empty exact-organism IEDB tab."""
    with tempfile.TemporaryDirectory(prefix="iedb-crawler-") as temp_name:
        temp_dir = Path(temp_name)
        store = _Store(temp_dir / "iedb.sqlite")
        try:
            if records is not None:
                _ingest_events(store, records)
                store.remove_ineligible()
                resolved_release_date = _to_iso_date(release_date) or datetime.date.today().isoformat()
                yield from _build_documents(
                    store,
                    resolved_release_date,
                    has_part_limit=has_part_limit,
                )
                return

            session = _get_session()
            export_path = _resource_path(
                session,
                temp_dir,
                "IEDB_EXPORT_PATH",
                EXPORT_URL,
                "iedb_export.zip",
            )
            organism_path = _resource_path(
                session,
                temp_dir,
                "IEDB_ORGANISM_LIST_PATH",
                ORGANISM_LIST_URL,
                "OrganismList.zip",
            )
            assay_type_path = _resource_path(
                session,
                temp_dir,
                "IEDB_ASSAY_TYPE_LIST_PATH",
                ASSAY_TYPE_LIST_URL,
                "AssayTypeList.zip",
            )
            taxdump_path = _resource_path(
                session,
                temp_dir,
                "IEDB_TAXDUMP_PATH",
                TAXDUMP_URL,
                "taxdump.tar.gz",
            )
            selected_organisms = _selected_organism_ids()
            _ingest_events(
                store,
                iter_iedb_records(
                    export_path,
                    selected_organisms=selected_organisms,
                    max_references=_int_env("IEDB_MAX_REFERENCES"),
                ),
            )
            _load_used_organisms(store, organism_path)
            _load_ncbi_taxonomy(store, taxdump_path)
            _mark_eligible_organisms(store)
            store.remove_ineligible()
            _supplement_antigens_and_receptors(store, session)
            resolved_release_date = _to_iso_date(release_date) or _archive_release_date(export_path)
            assay_type_lookup = _load_assay_types(assay_type_path)
            yield from _build_documents(
                store,
                resolved_release_date,
                assay_type_lookup=assay_type_lookup,
                has_part_limit=has_part_limit,
            )
        finally:
            store.close()


if __name__ == "__main__":
    for document in parse():
        print(json.dumps(document, sort_keys=True))
