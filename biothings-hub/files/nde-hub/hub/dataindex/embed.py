"""Embedding pipeline for NDE indices.

Generates dense-vector embeddings for Dataset, ComputationalTool, and
ResourceCatalog documents and writes them back to the same ES index.
Embeddings are reused from a local SQLite cache when possible, and model
inference is delegated to a remote embed server (see embed_server.py).

Can be imported and called from the hub indexer via ``run_embeddings()``,
or executed standalone::

    python embed.py
"""

import hashlib
import logging
import os
import re
import sqlite3
import time
from array import array

import requests
from elasticsearch import Elasticsearch, helpers

logger = logging.getLogger(__name__)

# ── constants ────────────────────────────────────────────────────────────────
EMBED_BATCH_SIZE = int(os.getenv("EMBED_BATCH_SIZE", "256"))
BULK_UPDATE_SIZE = int(os.getenv("BULK_UPDATE_SIZE", "512"))
SCAN_BATCH_SIZE = int(os.getenv("ES_SCAN_BATCH_SIZE", "1000"))
DOC_TYPES = ["Dataset", "ComputationalTool", "ResourceCatalog"]

ES_COMPATIBLE_WITH = os.getenv("ES_COMPATIBLE_WITH", "9").strip().lower()
ES_REQUEST_TIMEOUT = int(os.getenv("ES_REQUEST_TIMEOUT", "120"))

EMBED_SERVER_URL = os.getenv("EMBED_SERVER_URL", "http://su10:8485").rstrip("/")
EMBED_REQUEST_TIMEOUT = int(os.getenv("EMBED_REQUEST_TIMEOUT", "300"))
EMBED_CACHE_INDEX = os.getenv("EMBED_CACHE_INDEX", "").strip()
EMBED_CACHE_DB_PATH = os.getenv(
    "EMBED_CACHE_DB_PATH",
    "/data/nde-hub/cache/embeddings/embedding_cache.db",
).strip()
# Bump these when the remote model/config or _build_text() semantics change.
EMBED_MODEL_VERSION = os.getenv(
    "EMBED_MODEL_VERSION",
    "ibm-granite/granite-embedding-125m-english;normalize=true",
).strip()
EMBED_TEXT_VERSION = os.getenv("EMBED_TEXT_VERSION", "sample-summary-v1").strip()

# ── model configs ────────────────────────────────────────────────────────────
MODEL_CONFIGS = {
    "ibm-granite": {
        "field_name": "ibmGraniteEmbedding",
        "dims": 768,
        "index_options": {"type": "int8_hnsw"},
        "model_version": EMBED_MODEL_VERSION,
        "text_version": EMBED_TEXT_VERSION,
    },
}


# ── remote embed helper ─────────────────────────────────────────────────────
def _remote_embed(texts):
    """Send texts to the remote embedding server and return vectors."""
    resp = requests.post(
        f"{EMBED_SERVER_URL}/embed",
        json={"texts": texts},
        timeout=EMBED_REQUEST_TIMEOUT,
    )
    resp.raise_for_status()
    return resp.json()["embeddings"]


# ── ES helpers ───────────────────────────────────────────────────────────────
def _build_es_headers():
    if ES_COMPATIBLE_WITH in {"", "0", "false", "off", "none"}:
        return {
            "Accept": "application/json",
            "Content-Type": "application/json",
        }
    compatible_version = ES_COMPATIBLE_WITH.lstrip("v")
    if not compatible_version.isdigit():
        raise ValueError(
            f"ES_COMPATIBLE_WITH must be an integer like '8' (or 'none'), got: {ES_COMPATIBLE_WITH!r}"
        )
    media_type = f"application/vnd.elasticsearch+json; compatible-with={compatible_version}"
    return {"Accept": media_type, "Content-Type": media_type}


def _make_es_client(es_hosts):
    """Create a synchronous Elasticsearch client."""
    if isinstance(es_hosts, str):
        es_hosts = [es_hosts]
    return Elasticsearch(
        hosts=es_hosts,
        headers=_build_es_headers(),
        request_timeout=ES_REQUEST_TIMEOUT,
    )


_field_mapping_cache: dict = {}


def _add_embedding_field(es, index, field_name, dims, index_options=None):
    mapping = _field_mapping_cache.get(index)
    if mapping is None:
        mapping = es.indices.get_mapping(index=index)
        _field_mapping_cache[index] = mapping

    properties = mapping.get(index, {}).get("mappings", {}).get("properties", {})
    if field_name in properties:
        existing = properties[field_name]
        desired_type = (index_options or {}).get("type")
        existing_type = (existing.get("index_options") or {}).get("type")
        if desired_type and existing_type and desired_type != existing_type:
            logger.warning(
                "Field '%s' already has index_options.type='%s'; cannot change to '%s' via put_mapping",
                field_name,
                existing_type,
                desired_type,
            )
        return

    vector_props = {
        "type": "dense_vector",
        "dims": dims,
        "index": True,
        "similarity": "cosine",
    }
    if index_options:
        vector_props["index_options"] = index_options

    try:
        es.indices.put_mapping(index=index, body={"properties": {field_name: vector_props}})
    finally:
        _field_mapping_cache.pop(index, None)


def _add_embedding_metadata_fields(es, index, field_names):
    mapping = _field_mapping_cache.get(index)
    if mapping is None:
        mapping = es.indices.get_mapping(index=index)
        _field_mapping_cache[index] = mapping

    properties = mapping.get(index, {}).get("mappings", {}).get("properties", {})
    missing = {
        field_name: {"type": "keyword"}
        for field_name in field_names
        if field_name not in properties
    }
    if not missing:
        return

    try:
        es.indices.put_mapping(index=index, body={"properties": missing})
    finally:
        _field_mapping_cache.pop(index, None)


def _embedding_metadata_fields(target_field):
    return {
        "hash": f"{target_field}TextHash",
        "model": f"{target_field}Model",
        "text_version": f"{target_field}TextVersion",
    }


def _embedding_text_hash(model_name, model_version, text_version, text):
    payload = "\n".join(
        [
            f"model={model_name}",
            f"model_version={model_version}",
            f"text_version={text_version}",
            text,
        ]
    )
    return "sha256:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _serialize_vector(vector):
    return array("f", (float(value) for value in vector)).tobytes()


def _deserialize_vector(payload, dims):
    vector = array("f")
    vector.frombytes(payload)
    if len(vector) != dims:
        return None
    return list(vector)


def _normalize_vector(vector, dims):
    if isinstance(vector, (str, bytes, dict)):
        return None
    try:
        normalized = [float(value) for value in vector]
    except (TypeError, ValueError):
        return None
    if len(normalized) != dims:
        return None
    return normalized


def _open_embedding_cache(db_path, log):
    db_path = (db_path or "").strip()
    if db_path.lower() in {"", "0", "false", "off", "none"}:
        return None

    try:
        db_dir = os.path.dirname(db_path)
        if db_dir:
            os.makedirs(db_dir, exist_ok=True)

        conn = sqlite3.connect(db_path, timeout=60)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.execute(
            """CREATE TABLE IF NOT EXISTS embedding_cache (
                text_hash TEXT PRIMARY KEY,
                model_name TEXT NOT NULL,
                model_version TEXT NOT NULL,
                text_version TEXT NOT NULL,
                dims INTEGER NOT NULL,
                vector BLOB NOT NULL,
                source_doc_id TEXT,
                created_at INTEGER NOT NULL,
                updated_at INTEGER NOT NULL
            )"""
        )
        conn.commit()
        log.info("Embedding SQLite cache enabled at '%s'", db_path)
        return conn
    except Exception:
        log.exception("Could not open embedding SQLite cache at '%s'; continuing without it", db_path)
        return None


def _lookup_cached_embeddings(conn, text_hashes, dims):
    if not conn or not text_hashes:
        return {}

    unique_hashes = list(dict.fromkeys(text_hashes))
    placeholders = ",".join("?" for _ in unique_hashes)
    cursor = conn.cursor()
    cursor.execute(
        f"SELECT text_hash, dims, vector FROM embedding_cache WHERE text_hash IN ({placeholders})",
        unique_hashes,
    )

    cache = {}
    for text_hash, cached_dims, payload in cursor.fetchall():
        if cached_dims != dims:
            continue
        vector = _deserialize_vector(payload, dims)
        if vector is not None:
            cache[text_hash] = vector
    return cache


def _store_cached_embeddings(conn, rows, model_name, model_version, text_version, dims, log=None):
    if not conn or not rows:
        return True

    now = int(time.time())
    deduped_rows = list({row["text_hash"]: row for row in rows}.values())
    try:
        payloads = [
            (
                row["text_hash"],
                model_name,
                model_version,
                text_version,
                dims,
                _serialize_vector(row["vector"]),
                row.get("doc_id"),
                now,
                now,
            )
            for row in deduped_rows
        ]
        conn.executemany(
            """INSERT INTO embedding_cache (
                text_hash,
                model_name,
                model_version,
                text_version,
                dims,
                vector,
                source_doc_id,
                created_at,
                updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(text_hash) DO UPDATE SET
                model_name = excluded.model_name,
                model_version = excluded.model_version,
                text_version = excluded.text_version,
                dims = excluded.dims,
                vector = excluded.vector,
                source_doc_id = excluded.source_doc_id,
                updated_at = excluded.updated_at""",
            payloads,
        )
        conn.commit()
        return True
    except Exception:
        (log or logger).exception("Embedding SQLite cache store failed; continuing without local cache")
        try:
            conn.rollback()
        except Exception:
            pass
        return False


# ── text building ────────────────────────────────────────────────────────────
_SAMPLE_TEXT_MAX_FIELD_VALUES = 6
_SAMPLE_TEXT_MAX_TOTAL_VALUES = 40
_SAMPLE_TEXT_SKIP_KEYS = {
    "@context",
    "@id",
    "@type",
    "_id",
    "additionalIdentifier",
    "identifier",
    "inDefinedTermSet",
    "sameAs",
    "termCode",
    "unitCode",
    "url",
}
_SAMPLE_TEXT_REFERENCE_KEYS = {"itemListElement"}
_SAMPLE_TEXT_QUANTITATIVE_KEYS = {"name", "value", "minValue", "maxValue", "unitText"}
_SAMPLE_TEXT_QUANTITATIVE_SIGNAL_KEYS = {"value", "minValue", "maxValue", "unitText"}
_SAMPLE_TEXT_FIELD_LABELS = {
    "additionalType": ("sample category", "Sample categories"),
    "associatedGenotype": ("associated genotype", "Associated genotypes"),
    "associatedPhenotype": ("associated phenotype", "Associated phenotypes"),
    "anatomicalStructure": ("anatomical structure", "Anatomical structures"),
    "anatomicalSystem": ("anatomical system", "Anatomical systems"),
    "cellType": ("cell type", "Cell types"),
    "developmentalStage": ("developmental stage", "Developmental stages"),
    "experimentalPurpose": ("experimental purpose", "Experimental purposes"),
    "sampleAvailability": ("sample availability", "Sample availability"),
    "sampleProcess": ("sample process", "Sample processes"),
    "sampleQuantity": ("sample quantity", "Sample quantities"),
    "sampleState": ("sample state", "Sample states"),
    "sampleStorageTemperature": ("sample storage temperature", "Sample storage temperatures"),
    "sampleType": ("sample type", "Sample types"),
    "sex": ("sex", "Sex values"),
}
_SAMPLE_TEXT_ACCESSION_RE = re.compile(
    r"^(SAMN|SRS|SRX|SRR|ERS|ERX|ERR|DRS|DRX|DRR|GSM|GSE)[A-Z0-9_.-]*\d$",
    re.IGNORECASE,
)


def _clean_text(value):
    if not isinstance(value, str):
        return None
    value = " ".join(value.split()).strip()
    return value or None


def _append_unique_text(parts, seen, value):
    value = _clean_text(value)
    if not value:
        return
    key = value.casefold()
    if key in seen:
        return
    seen.add(key)
    parts.append(value)


def _humanize_key(key):
    if not key:
        return None
    words = []
    current = []
    previous = ""
    for char in str(key).replace("_", " "):
        if char.isupper() and previous and previous not in " -" and not previous.isupper():
            words.append("".join(current))
            current = [char.lower()]
        else:
            current.append(char.lower())
        previous = char
    if current:
        words.append("".join(current))
    return " ".join(" ".join(words).split()) or None


def _format_quantitative_value(payload):
    if not isinstance(payload, dict):
        return None
    if not any(key in payload for key in _SAMPLE_TEXT_QUANTITATIVE_SIGNAL_KEYS):
        return None

    pieces = []
    for key in ("name", "value", "minValue", "maxValue", "unitText"):
        if key not in payload:
            continue
        value = payload[key]
        if isinstance(value, (str, int, float)):
            pieces.append(str(value))
    return _clean_text(" ".join(pieces))


def _scalar_to_text(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (str, int, float)):
        return _clean_text(str(value))
    return None


def _looks_like_low_semantic_identifier(value):
    value = _clean_text(value)
    if not value:
        return True
    if value.lower().startswith(("http://", "https://")):
        return True
    if _SAMPLE_TEXT_ACCESSION_RE.match(value):
        return True

    compact = value.replace("_", "").replace("-", "").replace(".", "")
    if not compact.isalnum():
        return False
    alpha_count = sum(char.isalpha() for char in compact)
    digit_count = sum(char.isdigit() for char in compact)
    has_id_separator = any(separator in value for separator in ("_", "-", "."))

    # Keep biological terms like H1N1 or COVID-19, but drop accession-like
    # strings such as 1010108_1-Pm_HC1 that otherwise dominate the embedding.
    if len(value) >= 12 and " " not in value and has_id_separator and digit_count:
        return True
    return len(value) >= 8 and has_id_separator and digit_count > alpha_count


def _append_sample_summary_value(values_by_label, value_label, value, emitted):
    if emitted["count"] >= _SAMPLE_TEXT_MAX_TOTAL_VALUES:
        return

    value = _scalar_to_text(value)
    if not value or _looks_like_low_semantic_identifier(value):
        return

    values = values_by_label.setdefault(value_label, [])
    if len(values) >= _SAMPLE_TEXT_MAX_FIELD_VALUES:
        return

    key = value.casefold()
    if key in {existing.casefold() for existing in values}:
        return

    values.append(value)
    emitted["count"] += 1


def _iter_sample_field_values(payload, field_key):
    if isinstance(payload, list):
        for item in payload:
            yield from _iter_sample_field_values(item, field_key)
        return

    if isinstance(payload, dict):
        quantitative_value = _format_quantitative_value(payload)
        if quantitative_value:
            yield quantitative_value
            return

        for preferred_key in ("name", "value"):
            value = _scalar_to_text(payload.get(preferred_key))
            if value:
                yield value

        for key, value in payload.items():
            if key in _SAMPLE_TEXT_SKIP_KEYS or key in _SAMPLE_TEXT_REFERENCE_KEYS:
                continue
            if key in _SAMPLE_TEXT_QUANTITATIVE_KEYS:
                continue
            yield from _iter_sample_field_values(value, field_key)
        return

    if isinstance(payload, bool):
        if field_key == "sampleAvailability":
            yield "sample available" if payload else "sample unavailable"
        else:
            yield "true" if payload else "false"
        return

    value = _scalar_to_text(payload)
    if value:
        yield value


def _extract_sample_count(sample):
    if isinstance(sample, list):
        return len(sample) or None
    if not isinstance(sample, dict):
        return None

    number_of_items = sample.get("numberOfItems")
    if isinstance(number_of_items, dict):
        value = number_of_items.get("value")
    else:
        value = number_of_items

    if isinstance(value, (int, float)):
        return int(value) if float(value).is_integer() else value
    if isinstance(value, str):
        value = value.strip()
        if value.isdigit():
            return int(value)

    item_list = sample.get("itemListElement")
    if isinstance(item_list, list) and item_list:
        return len(item_list)
    return None


def _collect_sample_summary(payload, field_labels, values_by_label, emitted):
    if isinstance(payload, list):
        for item in payload:
            _collect_sample_summary(item, field_labels, values_by_label, emitted)
        return

    if not isinstance(payload, dict):
        return

    for key, value in payload.items():
        if key in _SAMPLE_TEXT_SKIP_KEYS or key in _SAMPLE_TEXT_REFERENCE_KEYS:
            continue

        label_info = _SAMPLE_TEXT_FIELD_LABELS.get(key)
        if label_info:
            field_label, value_label = label_info
            field_labels.add(field_label)
            for field_value in _iter_sample_field_values(value, key):
                _append_sample_summary_value(
                    values_by_label,
                    value_label,
                    field_value,
                    emitted,
                )
        else:
            _collect_sample_summary(value, field_labels, values_by_label, emitted)


def _append_sample_text(parts, seen, sample):
    if not sample:
        return

    field_labels = set()
    values_by_label = {}
    emitted = {"count": 0}
    _collect_sample_summary(sample, field_labels, values_by_label, emitted)

    # Keep the signal compact: say that sample metadata exists, summarize the
    # useful sample dimensions, and avoid accession/reference dumps.
    _append_unique_text(parts, seen, "Dataset has sample metadata")
    _append_unique_text(parts, seen, "Dataset includes sample-level metadata")

    sample_count = _extract_sample_count(sample)
    if sample_count:
        _append_unique_text(parts, seen, f"Dataset includes {sample_count} samples")

    ordered_field_labels = [
        field_label
        for field_label, _value_label in _SAMPLE_TEXT_FIELD_LABELS.values()
        if field_label in field_labels
    ]
    if ordered_field_labels:
        _append_unique_text(
            parts,
            seen,
            "Sample metadata includes " + ", ".join(ordered_field_labels),
        )

    for _field_label, value_label in _SAMPLE_TEXT_FIELD_LABELS.values():
        values = values_by_label.get(value_label)
        if values:
            verb = "includes" if value_label == "Sample availability" else "include"
            _append_unique_text(
                parts,
                seen,
                f"{value_label} {verb} {', '.join(values)}",
            )


def _build_text(source):
    parts = []
    seen = set()

    def _append_name(payload):
        if isinstance(payload, list):
            for item in payload:
                if isinstance(item, dict) and item.get("name"):
                    _append_unique_text(parts, seen, item["name"])
                elif isinstance(item, str):
                    _append_unique_text(parts, seen, item)
        elif isinstance(payload, dict) and payload.get("name"):
            _append_unique_text(parts, seen, payload["name"])
        elif isinstance(payload, str):
            _append_unique_text(parts, seen, payload)

    _append_name(source.get("author"))
    _append_name(source.get("measurementTechnique"))
    _append_name(source.get("species"))
    _append_name(source.get("healthCondition"))
    _append_name(source.get("infectiousAgent"))

    fundings = source.get("funding")
    if fundings:
        for f in (fundings if isinstance(fundings, list) else [fundings]):
            if isinstance(f, dict):
                _append_name(f.get("funder"))

    if name := source.get("name"):
        _append_unique_text(parts, seen, name)

    _append_sample_text(parts, seen, source.get("sample"))

    if desc := source.get("description"):
        _append_unique_text(parts, seen, desc)

    return "\n\n".join(parts).strip()


# ── formatting ───────────────────────────────────────────────────────────────
def _fmt_duration(seconds):
    seconds = max(int(seconds), 0)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


def _eta(processed, total, start):
    if not total or processed <= 0:
        return "calculating..."
    elapsed = time.time() - start
    if elapsed <= 0:
        return "calculating..."
    rate = processed / elapsed
    remaining = max(total - processed, 0)
    return "00:00" if remaining == 0 else _fmt_duration(remaining / rate)


# ── main entry point ─────────────────────────────────────────────────────────
def run_embeddings(es_hosts, index_name, cache_index_name=None, log=None):
    """Generate and write embeddings for all applicable docs in *index_name*.

    Parameters
    ----------
    es_hosts : str or list[str]
        Elasticsearch host URL(s).
    index_name : str
        Target index name.
    cache_index_name : str, optional
        Previous index/cache index to reuse embeddings from when text hashes match.
    log : logging.Logger, optional
        Logger instance; falls back to the module logger.
    """
    if log is None and cache_index_name is not None and not isinstance(cache_index_name, str):
        log = cache_index_name
        cache_index_name = None

    log = log or logger
    es = _make_es_client(es_hosts)
    sqlite_cache = _open_embedding_cache(EMBED_CACHE_DB_PATH, log)
    cache_index_name = (cache_index_name or EMBED_CACHE_INDEX or "").strip() or None
    if cache_index_name == index_name:
        log.warning("Ignoring embedding cache index '%s' because it matches target index", cache_index_name)
        cache_index_name = None

    try:
        for model_name, model_info in MODEL_CONFIGS.items():
            sqlite_cache = _run_model_embeddings(
                es,
                index_name,
                model_name,
                model_info,
                sqlite_cache,
                cache_index_name,
                log,
            )
    finally:
        if sqlite_cache:
            sqlite_cache.close()


def _run_model_embeddings(es, index_name, model_name, model_info, sqlite_cache, cache_index_name, log):
    target_field = model_info["field_name"]
    metadata_fields = _embedding_metadata_fields(target_field)
    model_version = model_info.get("model_version") or model_name
    text_version = model_info.get("text_version") or "default"
    dims = model_info["dims"]
    log.info("=== Embedding model: %s  field: %s ===", model_name, target_field)

    _add_embedding_field(
        es,
        index_name,
        target_field,
        dims,
        index_options=model_info.get("index_options"),
    )
    _add_embedding_metadata_fields(es, index_name, metadata_fields.values())

    # Ensure all recently-indexed docs are visible before querying.
    es.indices.refresh(index=index_name)
    active_cache_index = None
    if cache_index_name:
        if es.indices.exists(index=cache_index_name):
            active_cache_index = cache_index_name
            log.info("Embedding ES cache enabled for '%s' from '%s'", target_field, active_cache_index)
        else:
            log.warning("Embedding ES cache index '%s' does not exist; all misses will be embedded", cache_index_name)

    update_filter = {
        "bool": {
            "filter": [{"terms": {"@type": DOC_TYPES}}],
            "should": [
                {"bool": {"must_not": [{"exists": {"field": target_field}}]}},
                {"bool": {"must_not": [{"exists": {"field": metadata_fields["hash"]}}]}},
                {"bool": {"must_not": [{"term": {metadata_fields["model"]: model_version}}]}},
                {"bool": {"must_not": [{"term": {metadata_fields["text_version"]: text_version}}]}},
            ],
            "minimum_should_match": 1,
        }
    }
    current_filter = {
        "bool": {
            "filter": [
                {"terms": {"@type": DOC_TYPES}},
                {"exists": {"field": target_field}},
                {"exists": {"field": metadata_fields["hash"]}},
                {"term": {metadata_fields["model"]: model_version}},
                {"term": {metadata_fields["text_version"]: text_version}},
            ],
        }
    }
    update_scan_query = {
        "_source": {"excludes": [target_field]},
        "query": update_filter,
    }
    current_scan_query = {
        "_source": {"excludes": [target_field]},
        "query": current_filter,
    }
    update_total = es.count(index=index_name, body={"query": update_filter}).get("count", 0)
    current_total = es.count(index=index_name, body={"query": current_filter}).get("count", 0)
    total = update_total + current_total
    if total == 0:
        log.info(
            "No documents need '%s' embedding checks (%s / %s). Skipping.",
            target_field,
            model_version,
            text_version,
        )
        return sqlite_cache
    log.info(
        "Checking %d documents for '%s' (%d current candidates, %d missing/stale metadata)",
        total,
        target_field,
        current_total,
        update_total,
    )

    embedded = 0
    reused = 0
    sqlite_reused = 0
    es_reused = 0
    cache_miss = 0
    skip_count = 0
    unchanged = 0
    cleared = 0
    checked = 0
    processed = 0
    cache_batch: list[dict] = []
    embed_batch: list[dict] = []
    bulk_actions: list[dict] = []
    seen_doc_ids = set()
    start_time = time.time()

    def _embedding_update_doc(vector, text_hash):
        return {
            target_field: vector,
            metadata_fields["hash"]: text_hash,
            metadata_fields["model"]: model_version,
            metadata_fields["text_version"]: text_version,
        }

    def _log_progress():
        log.info(
            "Checked %d/%d (%d updated: %d embedded, %d reused, %d cleared; %d unchanged; ETA %s)",
            checked,
            total,
            processed,
            embedded,
            reused,
            cleared,
            unchanged,
            _eta(checked, total, start_time),
        )

    def _flush_bulk():
        nonlocal bulk_actions
        if bulk_actions:
            helpers.bulk(es, bulk_actions, stats_only=True, request_timeout=120)
            bulk_actions = []

    def _disable_sqlite_cache():
        nonlocal sqlite_cache
        conn = sqlite_cache
        sqlite_cache = None
        if not conn:
            return
        try:
            conn.rollback()
        except Exception:
            pass
        try:
            conn.close()
        except Exception:
            pass

    def _reuse_embedding(payload, vector, source):
        nonlocal es_reused, processed, reused, sqlite_reused
        bulk_actions.append(
            {
                "_op_type": "update",
                "_index": index_name,
                "_id": payload["doc_id"],
                "doc": _embedding_update_doc(vector, payload["text_hash"]),
                "doc_as_upsert": True,
            }
        )
        reused += 1
        processed += 1
        if source == "sqlite":
            sqlite_reused += 1
        elif source == "es":
            es_reused += 1
        if processed % 1000 == 0:
            _log_progress()
        if len(bulk_actions) >= BULK_UPDATE_SIZE:
            _flush_bulk()

    def _clear_embedding(doc_id):
        nonlocal cleared, processed
        bulk_actions.append(
            {
                "_op_type": "update",
                "_index": index_name,
                "_id": doc_id,
                "script": {
                    "source": "for (field in params.fields) { ctx._source.remove(field); }",
                    "params": {"fields": [target_field, *metadata_fields.values()]},
                },
            }
        )
        cleared += 1
        processed += 1
        if processed % 1000 == 0:
            _log_progress()
        if len(bulk_actions) >= BULK_UPDATE_SIZE:
            _flush_bulk()

    def _flush_embeds():
        nonlocal embed_batch, embedded, processed, sqlite_cache
        if not embed_batch:
            return

        unique_payloads = list({item["text_hash"]: item for item in embed_batch}.values())
        texts = [item["text"] for item in unique_payloads]
        vectors = _remote_embed(texts)
        if len(vectors) != len(unique_payloads):
            raise RuntimeError(
                f"Embedding server returned {len(vectors)} vectors for {len(unique_payloads)} texts"
            )

        vectors_by_hash = {}
        for payload, vec in zip(unique_payloads, vectors):
            normalized = _normalize_vector(vec, dims)
            if normalized is None:
                raise RuntimeError(
                    f"Embedding server returned an invalid vector for document {payload['doc_id']!r}"
                )
            vectors_by_hash[payload["text_hash"]] = normalized

        cache_rows = []
        for payload in embed_batch:
            vec = vectors_by_hash[payload["text_hash"]]
            bulk_actions.append(
                {
                    "_op_type": "update",
                    "_index": index_name,
                    "_id": payload["doc_id"],
                    "doc": _embedding_update_doc(vec, payload["text_hash"]),
                    "doc_as_upsert": True,
                }
            )
            cache_rows.append(
                {
                    "doc_id": payload["doc_id"],
                    "text_hash": payload["text_hash"],
                    "vector": vec,
                }
            )
            embedded += 1
            processed += 1
            if processed % 1000 == 0:
                _log_progress()
            if len(bulk_actions) >= BULK_UPDATE_SIZE:
                _flush_bulk()
        if sqlite_cache and not _store_cached_embeddings(
            sqlite_cache,
            cache_rows,
            model_name,
            model_version,
            text_version,
            dims,
            log,
        ):
            _disable_sqlite_cache()
        embed_batch = []

    def _append_to_embed_batch(payload):
        embed_batch.append(payload)
        if len(embed_batch) >= EMBED_BATCH_SIZE:
            _flush_embeds()

    def _try_es_cache(payloads):
        nonlocal active_cache_index, cache_miss, sqlite_cache
        if not payloads:
            return

        if not active_cache_index:
            for payload in payloads:
                cache_miss += 1
                _append_to_embed_batch(payload)
            return

        unique_hashes = list(dict.fromkeys(payload["text_hash"] for payload in payloads))
        source_fields = [
            target_field,
            metadata_fields["hash"],
            metadata_fields["model"],
            metadata_fields["text_version"],
        ]
        body = {
            "_source": source_fields,
            "size": len(unique_hashes),
            "query": {
                "bool": {
                    "filter": [
                        {"terms": {metadata_fields["hash"]: unique_hashes}},
                        {"term": {metadata_fields["model"]: model_version}},
                        {"term": {metadata_fields["text_version"]: text_version}},
                        {"exists": {"field": target_field}},
                    ]
                }
            },
            "collapse": {"field": metadata_fields["hash"]},
        }
        try:
            cache_res = es.search(
                index=active_cache_index,
                body=body,
                request_timeout=ES_REQUEST_TIMEOUT,
            )
        except Exception:
            log.exception(
                "Embedding ES cache lookup failed for '%s'; disabling ES cache for this run",
                active_cache_index,
            )
            active_cache_index = None
            for payload in payloads:
                cache_miss += 1
                _append_to_embed_batch(payload)
            return

        vectors_by_hash = {}
        for cache_doc in cache_res.get("hits", {}).get("hits", []):
            cache_source = cache_doc.get("_source") or {}
            text_hash = cache_source.get(metadata_fields["hash"])
            cached_vector = _normalize_vector(cache_source.get(target_field), dims)
            if (
                text_hash in unique_hashes
                and cached_vector is not None
                and cache_source.get(metadata_fields["model"]) == model_version
                and cache_source.get(metadata_fields["text_version"]) == text_version
            ):
                vectors_by_hash[text_hash] = cached_vector

        sqlite_cache_rows = []
        for payload in payloads:
            cached_vector = vectors_by_hash.get(payload["text_hash"])
            if cached_vector is not None:
                _reuse_embedding(payload, cached_vector, "es")
                sqlite_cache_rows.append(
                    {
                        "doc_id": payload["doc_id"],
                        "text_hash": payload["text_hash"],
                        "vector": cached_vector,
                    }
                )
            else:
                cache_miss += 1
                _append_to_embed_batch(payload)

        if sqlite_cache and not _store_cached_embeddings(
            sqlite_cache,
            sqlite_cache_rows,
            model_name,
            model_version,
            text_version,
            dims,
            log,
        ):
            _disable_sqlite_cache()

    def _flush_cache_lookups():
        nonlocal cache_batch, sqlite_cache
        if not cache_batch:
            return

        remaining = cache_batch
        if sqlite_cache:
            try:
                local_hits = _lookup_cached_embeddings(
                    sqlite_cache,
                    [payload["text_hash"] for payload in cache_batch],
                    dims,
                )
            except Exception:
                log.exception("Embedding SQLite cache lookup failed; disabling SQLite cache for this run")
                _disable_sqlite_cache()
                local_hits = {}

            if local_hits:
                remaining = []
                for payload in cache_batch:
                    cached_vector = local_hits.get(payload["text_hash"])
                    if cached_vector is not None:
                        _reuse_embedding(payload, cached_vector, "sqlite")
                    else:
                        remaining.append(payload)

        _try_es_cache(remaining)
        cache_batch = []

    def _scan_for_updates(scan_query, verify_current_hash):
        nonlocal checked, skip_count, unchanged
        for doc in helpers.scan(es, query=scan_query, index=index_name, size=SCAN_BATCH_SIZE):
            doc_id = doc["_id"]
            if doc_id in seen_doc_ids:
                continue
            seen_doc_ids.add(doc_id)
            checked += 1
            source = doc.get("_source", {})

            text = _build_text(source)
            if not text:
                skip_count += 1
                _clear_embedding(doc_id)
                if checked % 10000 == 0:
                    _log_progress()
                continue

            text_hash = _embedding_text_hash(
                model_name,
                model_version,
                text_version,
                text,
            )
            if verify_current_hash and source.get(metadata_fields["hash"]) == text_hash:
                unchanged += 1
                if checked % 10000 == 0:
                    _log_progress()
                continue

            cache_batch.append(
                {
                    "doc_id": doc_id,
                    "text": text,
                    "text_hash": text_hash,
                }
            )
            if len(cache_batch) >= EMBED_BATCH_SIZE:
                _flush_cache_lookups()

            if checked % 10000 == 0:
                _log_progress()

    _scan_for_updates(current_scan_query, verify_current_hash=True)
    _scan_for_updates(update_scan_query, verify_current_hash=False)

    _flush_cache_lookups()
    _flush_embeds()
    _flush_bulk()
    elapsed = _fmt_duration(time.time() - start_time)
    log.info(
        (
            "Finished %s: %d/%d checked, %d updated (%d embedded, %d reused, %d cleared), "
            "%d unchanged, %d cache misses, %d skipped empty text (elapsed %s)"
        ),
        target_field,
        checked,
        total,
        processed,
        embedded,
        reused,
        cleared,
        unchanged,
        cache_miss,
        skip_count,
        elapsed,
    )
    log.info(
        "Cache reuse breakdown for %s: %d SQLite, %d ES",
        target_field,
        sqlite_reused,
        es_reused,
    )
    return sqlite_cache
