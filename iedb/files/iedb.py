#!/usr/bin/env python3
"""IEDB species-by-tab DataCollection crawler for the NDE."""

import collections
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


EXPORT_URL = "https://www.iedb.org/downloader.php?file_name=doc/iedb_export.zip"
ORGANISM_LIST_URL = "https://www.iedb.org/downloader.php?file_name=doc/OrganismList.zip"
ASSAY_TYPE_LIST_URL = "https://www.iedb.org/downloader.php?file_name=doc/AssayTypeList.zip"
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
QUERY_API_BASE = "https://query-api.iedb.org"
WEB_BASE = "https://www.iedb.org"
# Host of IEDB's result-tab searches, <host>/query/<tab>?epitope_source_organism=NCBITaxon:<id>. IEDB launched
# them on https://test.iedb.org and is moving them to www.iedb.org.
SEARCH_BASE = os.environ.get("IEDB_SEARCH_BASE", WEB_BASE).rstrip("/")
USER_AGENT = "nde-iedb-datacollection-crawler/0.1"

REQUEST_TIMEOUT = int(os.environ.get("IEDB_REQUEST_TIMEOUT", "180"))
REQUEST_INTERVAL = float(os.environ.get("IEDB_REQUEST_INTERVAL", "0.1"))
QUERY_PAGE_SIZE = int(os.environ.get("IEDB_QUERY_PAGE_SIZE", "1000"))
# The Query API serves at most 10,000 rows per request
LINEAGE_PAGE_SIZE = int(os.environ.get("IEDB_LINEAGE_PAGE_SIZE", "10000"))
HAS_PART_LIMIT = int(os.environ.get("IEDB_HAS_PART_LIMIT", "1000"))
COMMIT_INTERVAL = int(os.environ.get("IEDB_COMMIT_INTERVAL", "50000"))
# How many of the organisms curated below a species to name in its records, so they are found by strain or isolate
DESCRIPTION_ORGANISM_LIMIT = 5
KEYWORD_ORGANISM_LIMIT = 10

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
            "identifier": "C13189",
            "inDefinedTermSet": "NCIT",
        },
        "endpoints": ("epitope_search",),
    },
    "antigens": {
        "label": "antigens",
        "unit": "antigens",
        "about": {
            "@type": "DefinedTerm",
            "name": "Antigen",
            "url": "http://purl.obolibrary.org/obo/NCIT_C268",
            "identifier": "C268",
            "inDefinedTermSet": "NCIT",
        },
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
# Ontology terms for the additionalType of assay and receptor members, keyed by the label they're given
ASSAY_TERMS = {
    "T-cell assay": {
        "url": "http://purl.obolibrary.org/obo/OBI_0003128",
        "identifier": "0003128",
        "inDefinedTermSet": "OBI",
    },
    "B-cell assay": {
        "url": "http://purl.obolibrary.org/obo/OBI_0003127",
        "identifier": "0003127",
        "inDefinedTermSet": "OBI",
    },
    "MHC binding assay": {
        "url": "http://purl.obolibrary.org/obo/OBI_1110126",
        "identifier": "1110126",
        "inDefinedTermSet": "OBI",
    },
    "MHC ligand-elution assay": {
        "url": "http://purl.obolibrary.org/obo/OBI_1110125",
        "identifier": "1110125",
        "inDefinedTermSet": "OBI",
    },
}
RECEPTOR_TERMS = {
    "TCR": {
        "url": "http://purl.obolibrary.org/obo/NCIT_C17065",
        "identifier": "C17065",
        "inDefinedTermSet": "NCIT",
    },
    "BCR": {
        "url": "http://purl.obolibrary.org/obo/NCIT_C107558",
        "identifier": "C107558",
        "inDefinedTermSet": "NCIT",
    },
}
MEMBERSHIP_ORGANISM_TAGS = {
    "SourceOrganismId",
    "SecondarySourceOrganismId",
    "PeptideSourceOrganismId",
    "MoleculeSourceOrganismId",
}
# Query API organism IRIs: NCBITaxon:<id> for NCBI taxa, taxon:<id> for organisms only IEDB's tree has
ORGANISM_IRI_PREFIXES = ("NCBITaxon:", "taxon:")
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


def _selected_species_ids() -> Optional[set[str]]:
    value = os.environ.get("IEDB_SPECIES_IDS")
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
        # (species_id, tab) -> how many member events each curated organism contributed, to name the most-curated
        self.organism_counts: dict[tuple[str, str], collections.Counter] = collections.defaultdict(collections.Counter)
        self._create_schema()

    def _create_schema(self) -> None:
        # taxa merges the IEDB organism tree (parents, names, IEDB-only organisms) with NCBI Taxonomy ranks; ncbi
        # marks real NCBI taxa. Every other table is keyed by the taxon a collection is grouped under.
        self.connection.executescript(
            """
            CREATE TABLE taxa (
                organism_id TEXT PRIMARY KEY,
                parent_id TEXT,
                rank TEXT,
                name TEXT,
                ncbi INTEGER NOT NULL DEFAULT 0
            ) WITHOUT ROWID;
            CREATE TABLE collections (
                species_id TEXT NOT NULL,
                tab TEXT NOT NULL,
                max_date TEXT,
                PRIMARY KEY (species_id, tab)
            ) WITHOUT ROWID;
            CREATE TABLE members (
                species_id TEXT NOT NULL,
                tab TEXT NOT NULL,
                member_key TEXT NOT NULL,
                member_json TEXT NOT NULL,
                PRIMARY KEY (species_id, tab, member_key)
            ) WITHOUT ROWID;
            CREATE TABLE assay_terms (
                species_id TEXT NOT NULL,
                assay_type_id TEXT NOT NULL,
                PRIMARY KEY (species_id, assay_type_id)
            ) WITHOUT ROWID;
            CREATE TABLE health_conditions (
                species_id TEXT NOT NULL,
                value TEXT NOT NULL,
                PRIMARY KEY (species_id, value)
            ) WITHOUT ROWID;
            """
        )

    def _maybe_commit(self) -> None:
        self._writes += 1
        if self._writes % COMMIT_INTERVAL == 0:
            self.connection.commit()

    def add_taxa(self, rows: Iterable[tuple[Any, Any, Any, Any]]) -> None:
        """Insert (organism_id, parent_id, NCBI rank, name) rows, e.g. to seed a small organism tree.

        A row without a rank is an IEDB-only organism.
        """
        self.connection.executemany(
            "INSERT OR REPLACE INTO taxa (organism_id, parent_id, rank, name, ncbi) VALUES (?, ?, ?, ?, ?)",
            [
                (*(None if value in (None, "") else str(value) for value in row), int(row[2] not in (None, "")))
                for row in rows
            ],
        )
        self.connection.commit()

    def names(self, organism_ids: Iterable[str]) -> dict[str, str]:
        found: dict[str, str] = {}
        ids = list(dict.fromkeys(str(organism_id) for organism_id in organism_ids))
        for start in range(0, len(ids), 500):
            chunk = ids[start : start + 500]
            placeholders = ",".join("?" for _ in chunk)
            for organism_id, name in self.connection.execute(
                f"SELECT organism_id, name FROM taxa WHERE organism_id IN ({placeholders}) AND name IS NOT NULL",
                chunk,
            ):
                found[organism_id] = name
        return found

    def add_event(self, species_id: str, organism_id: str, event: dict[str, Any]) -> None:
        """Add one tab member to a species collection, crediting the organism it was curated from."""
        tab = str(event.get("tab") or "")
        member = event.get("member") or {}
        member_key = str(event.get("member_key") or member.get("identifier") or "")
        if not species_id or tab not in TAB_CONFIG or not member_key:
            return

        self.organism_counts[(species_id, tab)][organism_id] += 1
        date_modified = _to_iso_date(event.get("date_modified"))
        self.connection.execute(
            """
            INSERT INTO collections (species_id, tab, max_date)
            VALUES (?, ?, ?)
            ON CONFLICT (species_id, tab) DO UPDATE SET
                max_date = CASE
                    WHEN excluded.max_date > COALESCE(collections.max_date, '') THEN excluded.max_date
                    ELSE collections.max_date
                END
            """,
            (species_id, tab, date_modified),
        )
        self.connection.execute(
            "INSERT OR IGNORE INTO members (species_id, tab, member_key, member_json) VALUES (?, ?, ?, ?)",
            (species_id, tab, member_key, json.dumps(_clean(member), sort_keys=True, separators=(",", ":"))),
        )
        if assay_type_id := event.get("assay_type_id"):
            self.connection.execute(
                "INSERT OR IGNORE INTO assay_terms (species_id, assay_type_id) VALUES (?, ?)",
                (species_id, str(assay_type_id)),
            )
        for condition in event.get("health_conditions") or []:
            if condition:
                self.connection.execute(
                    "INSERT OR IGNORE INTO health_conditions (species_id, value) VALUES (?, ?)",
                    (species_id, str(condition)),
                )
        self._maybe_commit()

    def close(self) -> None:
        self.connection.commit()
        self.connection.close()


class _GroupResolver:
    """Map a curated source organism to the species-level taxon its records are grouped under.

    The group is the highest taxon on the organism's lineage in the tree behind IEDB's own searches (``lineages``)
    that is still at or below an NCBI species. Where that tree keeps the NCBI species, as for Mycobacterium
    tuberculosis, the group is the species. Where it skips it, as for SARS-CoV-2, Influenza A virus and Dengue
    virus, which sit directly under a subgenus, genus or IEDB grouping, the group is the taxon IEDB treats as the
    species. Either way IEDB's search for the group returns the same records.

    Whether a taxon is at or below species comes from the IEDB organism list's parents and NCBI ranks (the ``taxa``
    table). An organism missing from the lineages takes the group of its nearest listed ancestor, else its NCBI
    species. Retired NCBI tax IDs, which the curation export still uses, are first mapped to their replacements.

    ``lineages`` come from an organism's own epitopes and are exact. ``shared_lineages``, for organisms only ever
    curated alongside others, can also hold those others, so only the parts that are also the organism's ancestors
    in the organism list are used.
    """

    def __init__(
        self,
        connection: sqlite3.Connection,
        lineages: Optional[dict[str, frozenset[str]]] = None,
        merged: Optional[dict[str, str]] = None,
        shared_lineages: Optional[dict[str, frozenset[str]]] = None,
    ):
        self.connection = connection
        self.exact = set(lineages or {})
        self.lineages = {**(shared_lineages or {}), **(lineages or {})}
        self.merged = merged or {}
        self._parents: dict[str, Optional[str]] = {}
        self._species: dict[str, Optional[str]] = {}
        self._depths: dict[str, int] = {}
        self._groups: dict[str, Optional[str]] = {}
        self.fallbacks = 0

    def canonical(self, organism_id: Optional[str]) -> str:
        organism_id = str(organism_id or "")
        return self.merged.get(organism_id, organism_id)

    def _parent(self, organism_id: str) -> Optional[str]:
        if organism_id not in self._parents:
            row = self.connection.execute("SELECT parent_id FROM taxa WHERE organism_id = ?", (organism_id,)).fetchone()
            parent_id = row[0] if row else None
            self._parents[organism_id] = parent_id if parent_id != organism_id else None
        return self._parents[organism_id]

    def species_id(self, organism_id: Optional[str]) -> Optional[str]:
        """The NCBI species of an organism: itself at species rank, else its nearest species ancestor."""
        if not organism_id:
            return None
        visited: list[str] = []
        current: Optional[str] = str(organism_id)
        species_id = None
        while current and current not in visited:
            if current in self._species:
                species_id = self._species[current]
                break
            row = self.connection.execute("SELECT rank FROM taxa WHERE organism_id = ?", (current,)).fetchone()
            if row is None:
                break
            visited.append(current)
            if row[0] == "species":
                species_id = current
                break
            current = self._parent(current)
        # Every organism walked through sits below the species found (or has none), so they share the answer
        for node in visited:
            self._species[node] = species_id
        return species_id

    def _depth(self, organism_id: str) -> int:
        if organism_id not in self._depths:
            depth = 0
            current = self._parent(organism_id)
            seen = {organism_id}
            while current and current not in seen:
                seen.add(current)
                depth += 1
                current = self._parent(current)
            self._depths[organism_id] = depth
        return self._depths[organism_id]

    def _path(self, organism_id: str) -> set[str]:
        """The organism and its ancestors in the organism list."""
        path: set[str] = set()
        current: Optional[str] = organism_id
        while current and current not in path:
            path.add(current)
            current = self._parent(current)
        return path

    def _above(self, upper: str, lower: str) -> bool:
        """Whether ``upper`` sits above ``lower`` in IEDB's tree, given both are on one organism's lineage."""
        if lower in self.exact:
            return upper in self.lineages[lower]
        if upper in self.exact:
            return lower not in self.lineages[upper]
        # A shared lineage holds every true ancestor plus organisms that only appear alongside it, so absence is
        # reliable and presence is not
        upper_lineage = self.lineages.get(upper)
        lower_lineage = self.lineages.get(lower)
        if lower_lineage is not None and upper not in lower_lineage:
            return False
        if upper_lineage is not None and lower not in upper_lineage:
            return True
        # No evidence, or each appears in the other's lineage: the organism list's depth stands in for IEDB's
        return (self._depth(upper), upper) < (self._depth(lower), lower)

    def _highest(self, nodes: Iterable[str]) -> Optional[str]:
        highest = None
        # Sorted first so the answer does not depend on set iteration order
        for node in sorted(nodes, key=lambda node: (self._depth(node), node)):
            if highest is None or self._above(node, highest):
                highest = node
        return highest

    def group_id(self, organism_id: Optional[str]) -> Optional[str]:
        organism_id = self.canonical(organism_id)
        if not organism_id:
            return None
        if organism_id not in self._groups:
            if (lineage := self.lineages.get(organism_id)) is not None:
                if organism_id not in self.exact:
                    lineage = lineage & self._path(organism_id)
                group = self._highest(node for node in lineage if self.species_id(node) is not None)
            elif ancestor := self._nearest_ancestor_with_lineage(organism_id):
                group = self.group_id(ancestor)
            else:
                group = self.species_id(organism_id)
                if group is not None and self.lineages:
                    self.fallbacks += 1
            self._groups[organism_id] = group
        return self._groups[organism_id]

    def _nearest_ancestor_with_lineage(self, organism_id: str) -> Optional[str]:
        """The closest ancestor IEDB's search data places, without leaving the organism's NCBI species."""
        seen = {organism_id}
        current = self._parent(organism_id)
        while current and current not in seen and self.species_id(current) is not None:
            if current in self.lineages:
                return current
            seen.add(current)
            current = self._parent(current)
        return None


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
            **ASSAY_TERMS[label],
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
    max_references: Optional[int] = None,
) -> Iterator[dict[str, Any]]:
    """Stream member events, one per curated source organism, from the reference-oriented XML archive."""
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


def _iter_organism_list(archive_path: Path) -> Iterator[tuple[str, Optional[str], Optional[str], int]]:
    """Yield (OrganismId, ParentTaxId, OrganismName, has NCBI TaxId) for every node of the IEDB organism list."""
    with zipfile.ZipFile(archive_path) as archive:
        xml_members = [info for info in archive.infolist() if info.filename.lower().endswith(".xml")]
        if len(xml_members) != 1:
            raise RuntimeError(f"Expected one OrganismList XML member; found {len(xml_members)}")
        with archive.open(xml_members[0]) as source:
            context = ET.iterparse(source, events=("start", "end"))
            _, root = next(context)
            for action, element in context:
                if action != "end" or _local_name(element.tag) != "Organism":
                    continue
                if organism_id := _first_descendant_text(element, "OrganismId"):
                    yield (
                        organism_id,
                        _first_descendant_text(element, "ParentTaxId"),
                        _first_descendant_text(element, "OrganismName"),
                        int(bool(_first_descendant_text(element, "TaxId"))),
                    )
                # The tree has millions of nodes; drop each one once read
                root.clear()


def _load_organism_list(store: _Store, archive_path: Path) -> None:
    """Load the parents and names of the IEDB organism tree, which also holds the IEDB-only organisms."""
    logger.info("Loading the IEDB organism tree")
    sql = """
        INSERT INTO taxa (organism_id, parent_id, name, ncbi) VALUES (?, ?, ?, ?)
        ON CONFLICT (organism_id) DO UPDATE SET
            parent_id = COALESCE(excluded.parent_id, taxa.parent_id),
            name = COALESCE(excluded.name, taxa.name),
            ncbi = MAX(excluded.ncbi, taxa.ncbi)
    """
    batch = []
    loaded = 0
    for row in _iter_organism_list(archive_path):
        batch.append(row)
        if len(batch) >= 50000:
            store.connection.executemany(sql, batch)
            loaded += len(batch)
            batch.clear()
    if batch:
        store.connection.executemany(sql, batch)
        loaded += len(batch)
    store.connection.commit()
    logger.info("Loaded %s IEDB organisms", f"{loaded:,}")


def _load_ncbi_taxonomy(store: _Store, archive_path: Path) -> dict[str, str]:
    """Load NCBI ranks and return retired tax IDs mapped to their replacements.

    A parent from NCBI is kept only for nodes the IEDB organism list does not place.
    """
    logger.info("Loading NCBI taxonomy ranks")
    sql = """
        INSERT INTO taxa (organism_id, parent_id, rank, ncbi) VALUES (?, ?, ?, 1)
        ON CONFLICT (organism_id) DO UPDATE SET
            parent_id = COALESCE(taxa.parent_id, excluded.parent_id),
            rank = excluded.rank,
            ncbi = 1
    """
    merged: dict[str, str] = {}
    loaded_nodes = False
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            if member.name not in {"merged.dmp", "nodes.dmp"}:
                continue
            source = archive.extractfile(member)
            if source is None:
                raise RuntimeError(f"{member.name} is missing from NCBI taxdump")
            loaded_nodes = loaded_nodes or member.name == "nodes.dmp"
            batch = []
            for raw_line in source:
                fields = raw_line.decode("utf-8").split("\t|\t")
                if member.name == "merged.dmp":
                    if len(fields) >= 2:
                        merged[fields[0].strip()] = fields[1].strip(" \t|\n")
                    continue
                if len(fields) < 3:
                    continue
                batch.append((fields[0].strip(), fields[1].strip(), fields[2].strip()))
                if len(batch) >= 50000:
                    store.connection.executemany(sql, batch)
                    batch.clear()
            if batch:
                store.connection.executemany(sql, batch)
    store.connection.commit()
    if not loaded_nodes:
        raise RuntimeError("nodes.dmp is missing from NCBI taxdump")
    if not merged:
        logger.warning("No merged.dmp in NCBI taxdump; retired tax IDs will not be mapped")
    return merged


def _organism_id_from_iri(iri: Any) -> Optional[str]:
    """An IEDB OrganismId from a Query API organism IRI: NCBITaxon:<id>, or taxon:<id> for IEDB-only organisms."""
    value = str(iri or "")
    for prefix in ORGANISM_IRI_PREFIXES:
        if value.startswith(prefix):
            organism_id = value[len(prefix) :]
            return organism_id if organism_id.isdigit() else None
    return None


def _load_iedb_lineages(session: requests.Session) -> tuple[dict[str, frozenset[str]], dict[str, frozenset[str]]]:
    """Each curated source organism's lineage in the organism tree that IEDB's searches use: (exact, shared).

    That tree is IEDB's own: it can skip NCBI ranks, e.g. SARS-CoV-2 sits directly under the subgenus Sarbecovirus
    rather than the species Betacoronavirus pandemicum. The Query API lists every epitope's source organisms with
    the union of their lineages (source_organism_iri_search). An epitope with a single source organism gives that
    organism's lineage exactly; for an organism never curated alone, intersecting the epitopes it shares leaves
    its lineage plus any organism it always appears with.
    """
    logger.info("Loading source organism lineages from the IEDB Query API")
    exact: dict[str, frozenset[str]] = {}
    shared: dict[str, frozenset[str]] = {}
    last_id = None
    rows_read = 0
    while True:
        params: dict[str, Any] = {
            "select": "structure_id,source_organism_iris,source_organism_iri_search",
            "order": "structure_id.asc",
            "limit": LINEAGE_PAGE_SIZE,
        }
        if last_id is not None:
            params["structure_id"] = f"gt.{last_id}"
        rows = _request_json(session, f"{QUERY_API_BASE}/epitope_search", params)
        if not isinstance(rows, list):
            raise TypeError(f"Unexpected epitope_search response: {type(rows).__name__}")
        for row in rows:
            ancestors = frozenset(
                organism_id
                for iri in _as_list(row.get("source_organism_iri_search"))
                if (organism_id := _organism_id_from_iri(iri))
            )
            organism_ids = _organism_ids_from_iris(row.get("source_organism_iris"))
            lineages = exact if len(organism_ids) == 1 else shared
            for organism_id in organism_ids:
                known = lineages.get(organism_id)
                lineages[organism_id] = ancestors if known is None else known & ancestors
        rows_read += len(rows)
        if len(rows) < LINEAGE_PAGE_SIZE:
            break
        last_id = rows[-1]["structure_id"]
        if rows_read % (LINEAGE_PAGE_SIZE * 50) == 0:
            logger.info("Read %s epitope_search rows", f"{rows_read:,}")
    shared = {organism_id: lineage for organism_id, lineage in shared.items() if organism_id not in exact}
    logger.info(
        "Loaded IEDB lineages for %s source organisms (%s only from shared epitopes) from %s epitopes",
        f"{len(exact) + len(shared):,}",
        f"{len(shared):,}",
        f"{rows_read:,}",
    )
    return exact, shared


def _log_skipped(store: _Store, skipped: collections.Counter, source: str) -> None:
    """Report the organisms whose members were left out because they have no species."""
    if not skipped:
        return
    top = skipped.most_common(10)
    names = store.names(organism_id for organism_id, _ in top)
    examples = "; ".join(
        f"{names.get(organism_id, 'unknown')} ({organism_id}): {count:,}" for organism_id, count in top
    )
    logger.info(
        "Skipped %s %s member events from %s organisms above species rank or missing from the IEDB organism tree. "
        "Most frequent: %s",
        f"{sum(skipped.values()):,}",
        source,
        f"{len(skipped):,}",
        examples,
    )


def _add_to_group(
    store: _Store,
    resolver: _GroupResolver,
    organism_id: Optional[str],
    event: dict[str, Any],
    skipped: collections.Counter,
    selected_groups: Optional[set[str]] = None,
) -> None:
    organism_id = resolver.canonical(organism_id)
    group_id = resolver.group_id(organism_id)
    if group_id is None:
        skipped[organism_id] += 1
    elif selected_groups is None or group_id in selected_groups:
        store.add_event(group_id, organism_id, event)


def _ingest_events(
    store: _Store,
    records: Iterable[dict[str, Any]],
    resolver: _GroupResolver,
    selected_groups: Optional[set[str]] = None,
) -> None:
    """Add each member event to the collection of the species-level group its curated source organism is in."""
    skipped: collections.Counter = collections.Counter()
    for record in records:
        _add_to_group(store, resolver, record.get("organism_id"), record, skipped, selected_groups)
    store.connection.commit()
    _log_skipped(store, skipped, "curation XML")


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


def _organism_ids_from_iris(iris: Any) -> list[str]:
    return list(dict.fromkeys(organism_id for iri in _as_list(iris) if (organism_id := _organism_id_from_iri(iri))))


def _supplement_antigens_and_receptors(
    store: _Store,
    session: requests.Session,
    resolver: _GroupResolver,
    selected_groups: Optional[set[str]] = None,
) -> None:
    if not _bool_env("IEDB_INCLUDE_QUERY_SUPPLEMENTS", True):
        logger.warning("Skipping Antigens and Receptors Query API supplements")
        return
    skipped: collections.Counter = collections.Counter()

    def add_to_groups(iris: Any, event: dict[str, Any]) -> None:
        for organism_id in _organism_ids_from_iris(iris):
            _add_to_group(store, resolver, organism_id, event, skipped, selected_groups)

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
        add_to_groups(
            row.get("source_organism_iris"),
            {"tab": "antigens", "member_key": str(antigen_iri), "member": member},
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
                        **RECEPTOR_TERMS[receptor_family],
                        "alternateName": row.get("receptor_type"),
                    },
                }
            )
            add_to_groups(
                row.get("source_organism_iris"),
                {"tab": "receptors", "member_key": identifier, "member": member},
            )
    store.connection.commit()
    _log_skipped(store, skipped, "Query API")


def _archive_release_date(archive_path: Path) -> str:
    with zipfile.ZipFile(archive_path) as archive:
        dates = [datetime.date(*info.date_time[:3]) for info in archive.infolist()]
    return max(dates).isoformat() if dates else datetime.date.today().isoformat()


def _taxon_curie(group_id: str, is_ncbi: bool) -> str:
    return f"NCBITaxon:{group_id}" if is_ncbi else f"taxon:{group_id}"


def _search_url(tab: str, group_id: str) -> str:
    """IEDB's result-tab search for an NCBI taxon and everything curated below it."""
    query = urllib.parse.urlencode({"epitope_source_organism": f"NCBITaxon:{group_id}"}, safe=":")
    return f"{SEARCH_BASE}/query/{tab}?{query}"


def _query_api_url(endpoint: str, curie: str) -> str:
    """Query API records for a taxon; source_organism_iri_search also matches everything curated below it."""
    query = urllib.parse.urlencode({"source_organism_iri_search": f"cs.{{{curie}}}"}, safe=":")
    return f"{QUERY_API_BASE}/{endpoint}?{query}"


def _record_url(tab: str, group_id: str, is_ncbi: bool) -> str:
    # IEDB's searches do not resolve IEDB-only taxa (as NCBITaxon:, taxon: or IEDBTaxon:), but the Query API does
    if is_ncbi:
        return _search_url(tab, group_id)
    return _query_api_url(TAB_CONFIG[tab]["endpoints"][0], _taxon_curie(group_id, is_ncbi))


def _species_term(group_id: str, name: str, is_ncbi: bool) -> dict[str, Any]:
    # No inDefinedTermSet: the hub's terms stage standardizes the term from its NCBI tax ID (or, for an IEDB-only
    # taxon, its name) and moves pathogens to infectiousAgent, which it skips for terms that name a term set.
    if not is_ncbi:
        return {"@type": "DefinedTerm", "name": name}
    return {
        "@type": "DefinedTerm",
        "name": name,
        "identifier": group_id,
        "url": f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={group_id}",
    }


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
    species_id: str,
    lookup: dict[str, dict[str, Optional[str]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    techniques = [copy.deepcopy(IMMUNE_EPITOPE_ASSAY)]
    variables = [copy.deepcopy(IMMUNE_RESPONSE)]
    categories = []
    for (assay_type_id,) in store.connection.execute(
        "SELECT assay_type_id FROM assay_terms WHERE species_id = ? ORDER BY assay_type_id",
        (species_id,),
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


def _provenance(tab: str, curie: str, release_date: str) -> list[dict[str, Any]]:
    query_works = [
        {
            "@type": "CreativeWork",
            "name": f"IEDB Query API {endpoint} species query",
            "url": _query_api_url(endpoint, curie),
        }
        for endpoint in TAB_CONFIG[tab]["endpoints"]
    ]
    return [
        {
            "@type": "Action",
            "name": "IEDB DataCollection generation process",
            "description": (
                f"Aggregated the IEDB {tab} result type at species level: results curated from a species, or from "
                "a strain, isolate or other taxon below it in IEDB's organism tree, count toward the species, as in "
                "IEDB's own search, and organisms curated above species rank are left out."
            ),
            "actionProcess": {
                "@type": "HowTo",
                "step": [
                    "Stream the IEDB curation XML reference archive.",
                    (
                        "Assign each primary epitope source organism to the highest taxon on its IEDB organism tree "
                        "lineage that is at or below an NCBI Taxonomy species."
                    ),
                    "Project each species' results into the five IEDB result-tab entity types.",
                    "Add Antigens and Receptors tab members from the IEDB Query API.",
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


def _example_of_work(tab: str, member: dict[str, Any], search_url: str) -> dict[str, Any]:
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
                "name": "Search this IEDB tab for the species",
                "url": search_url,
            },
            "additionalProperty": {
                "@type": "PropertyValue",
                "name": "Representative member identifier",
                "value": identifier,
            },
        }
    )


def _curated_organisms(store: _Store, group_id: str, tab: str) -> tuple[list[str], int]:
    """Names of the organisms curated below a group for one tab, most-curated first, and how many there are."""
    counts = store.organism_counts.get((group_id, tab)) or {}
    below = sorted(
        ((organism_id, count) for organism_id, count in counts.items() if organism_id != group_id),
        key=lambda item: (-item[1], int(item[0]) if item[0].isdigit() else 0, item[0]),
    )
    shown = [organism_id for organism_id, _ in below[: max(DESCRIPTION_ORGANISM_LIMIT, KEYWORD_ORGANISM_LIMIT)]]
    names = store.names(shown)
    return list(dict.fromkeys(names[organism_id] for organism_id in shown if organism_id in names)), len(below)


def _description(
    config: dict[str, Any],
    curie: str,
    name: str,
    member_count: int,
    curated_names: list[str],
    curated_total: int,
) -> str:
    # Records group by the epitope's source organism, not the assay host, so say so on every tab
    subject = "IEDB epitopes" if config["label"] == "epitopes" else f"IEDB {config['label']} for epitopes"
    text = (
        f"{subject} whose source organism is {name} ({curie}) or a strain, isolate or other taxon below it in "
        "IEDB's organism tree."
    )
    if shown := curated_names[:DESCRIPTION_ORGANISM_LIMIT]:
        listed = shown[0] if len(shown) == 1 else f"{', '.join(shown[:-1])} and {shown[-1]}"
        remaining = curated_total - len(shown)
        text += f" Organisms curated below it include {listed}"
        text += f", plus {remaining:,} more." if remaining > 0 else "."
    return f"{text} This collection contains {member_count:,} {config['unit']}."


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
            SELECT c.species_id, c.tab, c.max_date, t.name, COALESCE(t.ncbi, 1)
            FROM collections c
            LEFT JOIN taxa t ON t.organism_id = c.species_id
            """
        )
    )
    rows.sort(
        key=lambda row: (
            (0, int(row[0])) if row[0].isdigit() else (1, row[0]),
            TAB_ORDER[row[1]],
        )
    )

    for group_id, tab, max_date, group_name, is_ncbi in rows:
        is_ncbi = bool(is_ncbi)
        curie = _taxon_curie(group_id, is_ncbi)
        group_name = group_name or curie
        config = TAB_CONFIG[tab]
        url = _record_url(tab, group_id, is_ncbi)
        member_count = store.connection.execute(
            "SELECT COUNT(*) FROM members WHERE species_id = ? AND tab = ?",
            (group_id, tab),
        ).fetchone()[0]
        if not member_count:
            continue
        representative_row = store.connection.execute(
            "SELECT member_json FROM members WHERE species_id = ? AND tab = ? ORDER BY member_key LIMIT 1",
            (group_id, tab),
        ).fetchone()
        representative = json.loads(representative_row[0])
        curated_names, curated_total = _curated_organisms(store, group_id, tab)
        output = {
            "@context": "http://schema.org/",
            "@type": "DataCollection",
            "_id": f"iedb_{tab}_{group_id}",
            "identifier": f"IEDB:{tab.upper()}:{group_id}",
            "url": url,
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "Immune Epitope Database (IEDB)",
                "url": f"{WEB_BASE}/",
                "versionDate": release_date,
                "archivedAt": url,
            },
        }
        insert_value(output, "date", release_date)
        insert_value(output, "dateModified", max_date or release_date)
        insert_value(output, "name", f"IEDB {config['label']} from {group_name}")
        insert_value(
            output,
            "description",
            _description(config, curie, group_name, member_count, curated_names, curated_total),
            extend=True,
        )
        insert_value(
            output,
            "collectionSize",
            {"@type": "QuantitativeValue", "value": member_count, "unitText": config["unit"]},
        )
        insert_value(output, "about", copy.deepcopy(config["about"]))
        insert_value(output, "species", _species_term(group_id, group_name, is_ncbi))
        insert_value(output, "author", copy.deepcopy(IEDB_OPERATOR))
        insert_value(output, "creator", copy.deepcopy(IEDB_OPERATOR))
        insert_value(output, "citation", copy.deepcopy(CANONICAL_CITATION))
        insert_value(output, "conditionsOfAccess", "Open")
        insert_value(output, "isAccessibleForFree", True)
        insert_value(output, "license", "https://creativecommons.org/licenses/by/4.0/")
        insert_value(output, "funding", copy.deepcopy(FUNDING))
        insert_value(
            output,
            "keywords",
            ["IEDB", config["label"], group_name, *curated_names[:KEYWORD_ORGANISM_LIMIT]],
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
        insert_value(output, "isBasedOn", _provenance(tab, curie, release_date))

        if tab == "assays":
            techniques, variables, categories = _assay_terms(store, group_id, assay_type_lookup)
            insert_value(output, "measurementTechnique", techniques)
            insert_value(output, "variableMeasured", variables)
            if categories:
                insert_value(output, "keywords", categories)
            conditions = [
                _health_condition(row[0])
                for row in store.connection.execute(
                    "SELECT value FROM health_conditions WHERE species_id = ? ORDER BY value",
                    (group_id,),
                )
            ]
            if conditions:
                insert_value(output, "healthCondition", conditions)

        if has_part_limit >= 0 and member_count <= has_part_limit:
            members = [
                json.loads(row[0])
                for row in store.connection.execute(
                    "SELECT member_json FROM members WHERE species_id = ? AND tab = ? ORDER BY member_key",
                    (group_id, tab),
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


def parse(
    records: Optional[Iterable[dict[str, Any]]] = None,
    *,
    taxa: Iterable[tuple[Any, Any, Any, Any]] = (),
    lineages: Optional[dict[str, frozenset[str]]] = None,
    shared_lineages: Optional[dict[str, frozenset[str]]] = None,
    release_date: Optional[str] = None,
    has_part_limit: int = HAS_PART_LIMIT,
) -> Iterator[dict[str, Any]]:
    """Yield one DataCollection for every non-empty IEDB tab of every species-level taxon.

    Given ``records`` (member events), skip the downloads and group them with the organism tree in ``taxa``,
    (organism_id, parent_id, NCBI rank, name) rows, and the IEDB search lineages of curated organisms (see
    _GroupResolver).
    """
    with tempfile.TemporaryDirectory(prefix="iedb-crawler-") as temp_name:
        temp_dir = Path(temp_name)
        store = _Store(temp_dir / "iedb.sqlite")
        try:
            if records is not None:
                store.add_taxa(taxa)
                resolver = _GroupResolver(store.connection, lineages, shared_lineages=shared_lineages)
                _ingest_events(store, records, resolver)
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
            _load_organism_list(store, organism_path)
            merged = _load_ncbi_taxonomy(store, taxdump_path)
            exact_lineages, shared_lineages = _load_iedb_lineages(session)
            resolver = _GroupResolver(store.connection, exact_lineages, merged, shared_lineages)
            selected_groups = _selected_species_ids()
            _ingest_events(
                store,
                iter_iedb_records(export_path, max_references=_int_env("IEDB_MAX_REFERENCES")),
                resolver,
                selected_groups,
            )
            _supplement_antigens_and_receptors(store, session, resolver, selected_groups)
            if resolver.fallbacks:
                logger.warning(
                    "Grouped %s organisms missing from IEDB's search lineages by their NCBI species",
                    resolver.fallbacks,
                )
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
