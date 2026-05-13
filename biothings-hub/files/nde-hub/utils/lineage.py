import json
import os
import sqlite3
from typing import Iterable, List, Set

import orjson
from biothings_client import get_client

DB_PATH = "/data/nde-hub/standardizers/lineage_lookup/lineage_lookup.db"

_TAXA_CHUNK_SIZE = 1000
_BATCH_SIZE = 1000

_mt = None


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


def _load_cached_lineages(taxon_ids: Set[int], lineage_cache: dict):
    """Load SQLite-cached lineage rows into a batch-local cache."""
    missing_ids = set(taxon_ids) - set(lineage_cache)
    if not missing_ids:
        return

    _ensure_db()
    with sqlite3.connect(DB_PATH) as conn:
        for chunk in _chunked(sorted(missing_ids), _TAXA_CHUNK_SIZE):
            placeholders = ",".join("?" for _ in chunk)
            for taxid, lineage_json in conn.execute(
                f"SELECT taxid, lineage FROM taxon_lineage WHERE taxid IN ({placeholders})",
                chunk,
            ):
                lineage_cache[taxid] = json.loads(lineage_json)


def _load_cached_parents(taxon_ids: Set[int], parent_cache: dict):
    """Load SQLite-cached parent rows into a batch-local cache."""
    missing_ids = set(taxon_ids) - set(parent_cache)
    if not missing_ids:
        return

    _ensure_db()
    with sqlite3.connect(DB_PATH) as conn:
        for chunk in _chunked(sorted(missing_ids), _TAXA_CHUNK_SIZE):
            placeholders = ",".join("?" for _ in chunk)
            for taxid, parent in conn.execute(
                f"SELECT taxid, parent_taxid FROM taxon_parent WHERE taxid IN ({placeholders})",
                chunk,
            ):
                parent_cache[taxid] = parent


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


def _fetch_taxon_info(taxon_ids: Set[int], lineage_cache: dict, parent_cache: dict):
    """Fetch and persist lineage/parent info not present in the batch cache."""
    _load_cached_lineages(taxon_ids, lineage_cache)
    _load_cached_parents(taxon_ids, parent_cache)

    new_ids = set(taxon_ids) - set(lineage_cache)

    mt = None
    lineage_rows: list = []
    parent_rows: list = []

    if new_ids:
        mt = _get_client()
        for chunk in _chunked(sorted(new_ids), _TAXA_CHUNK_SIZE):
            taxon_info_list = mt.gettaxa(chunk)
            for taxon_info in taxon_info_list:
                taxid = taxon_info.get("taxid")
                lineage = taxon_info.get("lineage", [])
                parent_taxid = taxon_info.get("parent_taxid")
                if taxid is not None:
                    lineage_cache[taxid] = lineage
                    parent = parent_taxid if parent_taxid else None
                    parent_cache[taxid] = parent
                    lineage_rows.append((taxid, json.dumps(lineage)))
                    parent_rows.append((taxid, parent))

    lineage_taxon_ids: Set[int] = set(taxon_ids)
    for lineage in lineage_cache.values():
        try:
            lineage_taxon_ids.update(lineage)
        except TypeError:
            pass

    # Fetch parent info for lineage ancestors not yet in the batch cache.
    _load_cached_parents(lineage_taxon_ids, parent_cache)
    missing = lineage_taxon_ids - set(parent_cache)
    if missing:
        if mt is None:
            mt = _get_client()
        for chunk in _chunked(sorted(missing), _TAXA_CHUNK_SIZE):
            try:
                info_list = mt.gettaxa(chunk)
                for taxon_info in info_list:
                    taxid = taxon_info.get("taxid")
                    parent_taxid = taxon_info.get("parent_taxid")
                    if taxid is not None:
                        parent = parent_taxid if parent_taxid else None
                        parent_cache[taxid] = parent
                        parent_rows.append((taxid, parent))
            except Exception as e:
                print(f"Error fetching lineage taxon info chunk: {e}")

    _save_to_db(lineage_rows, parent_rows)


def _get_lineage_entries(taxid: int, lineage_cache: dict, parent_cache: dict):
    lineage_ids = list(lineage_cache.get(taxid, [])) + [taxid]
    entries = []
    for taxon in lineage_ids:
        entry = {"taxon": taxon}
        parent_taxon = parent_cache.get(taxon)
        if taxon != 1 and parent_taxon is not None:
            entry["parent_taxon"] = parent_taxon
        entries.append(entry)
    return entries


def _annotate_record(record: dict, lineage_cache: dict, parent_cache: dict):
    """Add _meta.lineage to a single record based on its taxon IDs."""
    lineage_entries_set: set = set()
    for taxid in _extract_taxids(record):
        for entry in _get_lineage_entries(taxid, lineage_cache, parent_cache):
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


def _process_batch(batch: list):
    all_ids: Set[int] = set()
    for rec in batch:
        all_ids.update(_extract_taxids(rec))

    lineage_cache: dict = {}
    parent_cache: dict = {}
    if all_ids:
        _fetch_taxon_info(all_ids, lineage_cache, parent_cache)

    for rec in batch:
        _annotate_record(rec, lineage_cache, parent_cache)
        yield rec


def process_lineage(docs):
    """Add taxonomy lineage to documents using a persistent SQLite-backed cache.

    Accepts an iterable of dicts (or a path to an ndjson directory) and yields
    each document with ``_meta.lineage`` populated.  Taxon lookups are cached in
    SQLite at ``DB_PATH`` so data persists across process restarts. In-memory
    lineage and parent dictionaries are scoped to one internal batch and are
    discarded before the next batch is processed.

    Documents are collected into small internal batches so that API calls are
    amortised without ever materialising the full dataset in memory.
    """
    batch: list = []
    for doc in _iter_docs(docs):
        batch.append(doc)
        if len(batch) >= _BATCH_SIZE:
            yield from _process_batch(batch)
            batch = []

    if batch:
        yield from _process_batch(batch)
