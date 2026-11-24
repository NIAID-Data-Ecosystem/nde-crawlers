"""Taxonomy lineage for the taxonomy browser.

Adds `_meta.lineage` to the record types the portal exposes, from the numeric
taxonomy identifiers on `species` / `infectiousAgent`. Taxon lookups are cached
(`cache.Cache`) so they persist across runs.
"""

from itertools import batched
from typing import Set

from biothings_client import get_client
from config import logger

from .cache import Cache
from .common import as_list, dict_entries

_TAXA_CHUNK_SIZE = 1000
_BATCH_SIZE = 1000
_LINEAGE_TYPES = {"Dataset", "DataCollection", "ResourceCatalog"}
_BIOTOOLS_CATALOG_NAME = "bio.tools"

_mt = None


def _get_client():
    global _mt
    if _mt is None:
        _mt = get_client("taxon")
    return _mt


# taxid -> its ancestors' taxids
LINEAGES = Cache("taxon_lineage", normalize=False)
# taxid -> its parent's taxid, or None
PARENTS = Cache("taxon_parent", normalize=False)


def _load_cached(taxon_ids: Set[int], cache: dict, table: Cache):
    """Load `table` rows for the taxon ids not already in `cache`."""
    missing_ids = set(taxon_ids) - set(cache)
    if missing_ids:
        cache.update({int(taxid): value for taxid, value in table.get_many(missing_ids).items()})


def _load_cached_lineages(taxon_ids: Set[int], lineage_cache: dict):
    _load_cached(taxon_ids, lineage_cache, LINEAGES)


def _load_cached_parents(taxon_ids: Set[int], parent_cache: dict):
    _load_cached(taxon_ids, parent_cache, PARENTS)


def _iter_string_values(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from _iter_string_values(item)
    elif isinstance(value, (list, tuple, set)):
        for item in value:
            yield from _iter_string_values(item)


def _has_biosample_additional_type(record: dict) -> bool:
    return any(value.strip() == "BioSample" for value in _iter_string_values(record.get("additionalType")))


def _iter_catalog_names(value):
    for catalog in as_list(value):
        if isinstance(catalog, dict):
            yield from _iter_string_values(catalog.get("name"))


def _is_biotools_record(record: dict) -> bool:
    return any(
        value.strip().lower() == _BIOTOOLS_CATALOG_NAME
        for value in _iter_catalog_names(record.get("includedInDataCatalog"))
    )


def _should_annotate_lineage(record: dict) -> bool:
    record_types = set(_iter_string_values(record.get("@type")))
    if record_types & _LINEAGE_TYPES:
        return True
    if "Sample" in record_types:
        return _has_biosample_additional_type(record)
    if "ComputationalTool" in record_types:
        return _is_biotools_record(record)
    return False


def _extract_taxids(record: dict) -> Set[int]:
    taxids: Set[int] = set()
    for field in ["species", "infectiousAgent"]:
        for item in dict_entries(record, field):
            taxid = str(item.get("identifier"))
            if taxid.isdigit():
                taxids.add(int(taxid))
    return taxids


def _remove_lineage(record: dict):
    meta = record.get("_meta")
    if isinstance(meta, dict):
        meta.pop("lineage", None)


def _fetch_taxon_info(taxon_ids: Set[int], lineage_cache: dict, parent_cache: dict):
    """Fetch and persist lineage/parent info not present in the batch cache."""
    _load_cached_lineages(taxon_ids, lineage_cache)
    _load_cached_parents(taxon_ids, parent_cache)

    new_ids = set(taxon_ids) - set(lineage_cache)

    mt = None
    new_lineages: dict = {}
    new_parents: dict = {}

    if new_ids:
        mt = _get_client()
        for chunk in batched(sorted(new_ids), _TAXA_CHUNK_SIZE):
            taxon_info_list = mt.gettaxa(list(chunk))
            for taxon_info in taxon_info_list:
                taxid = taxon_info.get("taxid")
                lineage = taxon_info.get("lineage", [])
                parent_taxid = taxon_info.get("parent_taxid")
                if taxid is not None:
                    lineage_cache[taxid] = lineage
                    parent = parent_taxid if parent_taxid else None
                    parent_cache[taxid] = parent
                    new_lineages[taxid] = lineage
                    new_parents[taxid] = parent

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
        for chunk in batched(sorted(missing), _TAXA_CHUNK_SIZE):
            try:
                info_list = mt.gettaxa(list(chunk))
                for taxon_info in info_list:
                    taxid = taxon_info.get("taxid")
                    parent_taxid = taxon_info.get("parent_taxid")
                    if taxid is not None:
                        parent = parent_taxid if parent_taxid else None
                        parent_cache[taxid] = parent
                        new_parents[taxid] = parent
            except Exception as e:
                logger.error("Error fetching lineage taxon info chunk: %s", e)

    LINEAGES.put_many(new_lineages)
    PARENTS.put_many(new_parents)


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

    if lineage_entries:
        record.setdefault("_meta", {})["lineage"] = lineage_entries
    else:
        _remove_lineage(record)


def _process_batch(batch: list):
    eligible_records = []
    all_ids: Set[int] = set()
    for rec in batch:
        if _should_annotate_lineage(rec):
            eligible_records.append(rec)
            all_ids.update(_extract_taxids(rec))
        else:
            _remove_lineage(rec)

    lineage_cache: dict = {}
    parent_cache: dict = {}
    if all_ids:
        _fetch_taxon_info(all_ids, lineage_cache, parent_cache)

    for rec in eligible_records:
        _annotate_record(rec, lineage_cache, parent_cache)

    yield from batch


def process_lineage(docs):
    """
        Add _meta.lineage (taxonoy lineage) to portal-visible records with numeric taxonomy IDs.
    """
    batch: list = []
    for doc in docs:
        batch.append(doc)
        if len(batch) >= _BATCH_SIZE:
            yield from _process_batch(batch)
            batch = []

    if batch:
        yield from _process_batch(batch)
