import json
import os
import sqlite3
from typing import Iterable, List, Set

import orjson
from biothings_client import get_client

DB_PATH = "/data/nde-hub/standardizers/lineage_lookup/lineage_lookup.db"

_TAXA_CHUNK_SIZE = 1000
_BATCH_SIZE = 1000

# ---------------------------------------------------------------------------
# In-memory hot caches (populated from SQLite on first miss, then kept warm).
# ---------------------------------------------------------------------------
_taxon_lineage_cache: dict = {}   # taxid -> list of ancestor taxids
_taxon_parent_cache: dict = {}    # taxid -> parent taxid (or None)
_fetched_taxon_ids: set = set()   # taxids already in memory
_mt = None
_db_loaded = False


def _get_client():
    global _mt
    if _mt is None:
        _mt = get_client("taxon")
    return _mt


def _ensure_db():
    """Create the SQLite database and tables if they don't exist."""
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute(
            "CREATE TABLE IF NOT EXISTS taxon_lineage "
            "(taxid INTEGER PRIMARY KEY, lineage TEXT NOT NULL)"
        )
        conn.execute(
            "CREATE TABLE IF NOT EXISTS taxon_parent "
            "(taxid INTEGER PRIMARY KEY, parent_taxid INTEGER)"
        )


def _load_db_into_memory():
    """Bulk-load the entire SQLite cache into the in-memory dicts (once)."""
    global _db_loaded
    if _db_loaded:
        return
    _ensure_db()
    with sqlite3.connect(DB_PATH) as conn:
        for taxid, lineage_json in conn.execute("SELECT taxid, lineage FROM taxon_lineage"):
            _taxon_lineage_cache[taxid] = json.loads(lineage_json)
            _fetched_taxon_ids.add(taxid)
        for taxid, parent in conn.execute("SELECT taxid, parent_taxid FROM taxon_parent"):
            _taxon_parent_cache[taxid] = parent
            _fetched_taxon_ids.add(taxid)
    _db_loaded = True


def _save_to_db(lineage_rows: list, parent_rows: list):
    """Persist newly fetched taxon data to SQLite."""
    if not lineage_rows and not parent_rows:
        return
    with sqlite3.connect(DB_PATH) as conn:
        if lineage_rows:
            conn.executemany(
                "INSERT OR REPLACE INTO taxon_lineage (taxid, lineage) VALUES (?, ?)",
                lineage_rows,
            )
        if parent_rows:
            conn.executemany(
                "INSERT OR REPLACE INTO taxon_parent (taxid, parent_taxid) VALUES (?, ?)",
                parent_rows,
            )


def _chunked(iterable: Iterable[int], chunk_size: int) -> Iterable[List[int]]:
    chunk: List[int] = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _extract_taxids(record: dict) -> Set[int]:
    taxids: Set[int] = set()
    for field in ["species", "infectiousAgent"]:
        value = record.get(field)
        if value is None:
            continue
        items = value if isinstance(value, list) else [value]
        for item in items:
            if not isinstance(item, dict):
                continue
            taxid = item.get("identifier")
            if taxid is None:
                continue
            taxid_str = str(taxid)
            if taxid_str.isdigit():
                taxids.add(int(taxid_str))
    return taxids


def _fetch_taxon_info(taxon_ids: Set[int]):
    """Fetch and cache taxon lineage/parent info for any IDs not yet cached."""
    _load_db_into_memory()
    mt = _get_client()
    new_ids = taxon_ids - _fetched_taxon_ids
    if not new_ids:
        return

    lineage_rows: list = []
    parent_rows: list = []
    lineage_taxon_ids: Set[int] = set()

    for chunk in _chunked(sorted(new_ids), _TAXA_CHUNK_SIZE):
        taxon_info_list = mt.gettaxa(chunk)
        for taxon_info in taxon_info_list:
            taxid = taxon_info.get("taxid")
            lineage = taxon_info.get("lineage", [])
            parent_taxid = taxon_info.get("parent_taxid")
            if taxid is not None:
                _taxon_lineage_cache[taxid] = lineage
                parent = parent_taxid if parent_taxid else None
                _taxon_parent_cache[taxid] = parent
                lineage_rows.append((taxid, json.dumps(lineage)))
                parent_rows.append((taxid, parent))
                try:
                    lineage_taxon_ids.update(lineage)
                except TypeError:
                    pass

    _fetched_taxon_ids.update(new_ids)

    # Fetch parent info for lineage ancestors not yet in cache
    missing = lineage_taxon_ids - _fetched_taxon_ids
    if missing:
        for chunk in _chunked(sorted(missing), _TAXA_CHUNK_SIZE):
            try:
                info_list = mt.gettaxa(chunk)
                for taxon_info in info_list:
                    taxid = taxon_info.get("taxid")
                    parent_taxid = taxon_info.get("parent_taxid")
                    if taxid is not None:
                        parent = parent_taxid if parent_taxid else None
                        _taxon_parent_cache[taxid] = parent
                        parent_rows.append((taxid, parent))
            except Exception as e:
                print(f"Error fetching lineage taxon info chunk: {e}")
        _fetched_taxon_ids.update(missing)

    _save_to_db(lineage_rows, parent_rows)


def _get_lineage_entries(taxid: int):
    lineage_ids = list(_taxon_lineage_cache.get(taxid, [])) + [taxid]
    entries = []
    for taxon in lineage_ids:
        entry = {"taxon": taxon}
        parent_taxon = _taxon_parent_cache.get(taxon)
        if taxon != 1 and parent_taxon is not None:
            entry["parent_taxon"] = parent_taxon
        entries.append(entry)
    return entries


def _annotate_record(record: dict):
    """Add _meta.lineage to a single record based on its taxon IDs."""
    lineage_entries_set: set = set()
    for taxid in _extract_taxids(record):
        for entry in _get_lineage_entries(taxid):
            lineage_entries_set.add((entry["taxon"], entry.get("parent_taxon")))

    lineage_entries = []
    for taxon, parent_taxon in lineage_entries_set:
        entry = {"taxon": taxon}
        if parent_taxon is not None:
            entry["parent_taxon"] = parent_taxon
        lineage_entries.append(entry)

    record.setdefault("_meta", {})["lineage"] = lineage_entries


def _iter_docs(docs):
    """Normalise *docs* into an iterator of dicts (handles path strings too)."""
    if isinstance(docs, str):
        with open(os.path.join(docs, "data.ndjson"), "rb") as f:
            for line in f:
                yield orjson.loads(line)
    else:
        yield from docs


def process_lineage(docs):
    """Add taxonomy lineage to documents using a persistent SQLite-backed cache.

    Accepts an iterable of dicts (or a path to an ndjson directory) and yields
    each document with ``_meta.lineage`` populated.  Taxon lookups are cached in
    SQLite at ``DB_PATH`` so data persists across process restarts, and in
    in-memory dicts for fast per-record access within a run.

    Documents are collected into small internal batches so that API calls are
    amortised without ever materialising the full dataset in memory.
    """
    batch: list = []
    for doc in _iter_docs(docs):
        batch.append(doc)
        if len(batch) >= _BATCH_SIZE:
            all_ids: Set[int] = set()
            for rec in batch:
                all_ids.update(_extract_taxids(rec))
            if all_ids:
                _fetch_taxon_info(all_ids)
            for rec in batch:
                _annotate_record(rec)
                yield rec
            batch = []

    if batch:
        all_ids = set()
        for rec in batch:
            all_ids.update(_extract_taxids(rec))
        if all_ids:
            _fetch_taxon_info(all_ids)
        for rec in batch:
            _annotate_record(rec)
            yield rec
