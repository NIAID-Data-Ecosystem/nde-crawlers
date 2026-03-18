"""Embedding pipeline for NDE indices.

Generates dense-vector embeddings for Dataset, ComputationalTool, and
ResourceCatalog documents and writes them back to the same ES index.
Model inference is delegated to a remote embed server (see embed_server.py).

Can be imported and called from the hub indexer via ``run_embeddings()``,
or executed standalone::

    python embed.py
"""

import logging
import os
import time

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

# ── model configs ────────────────────────────────────────────────────────────
MODEL_CONFIGS = {
    "ibm-granite": {
        "field_name": "ibmGraniteEmbedding",
        "dims": 768,
        "index_options": {"type": "int8_hnsw"},
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


# ── text building ────────────────────────────────────────────────────────────
def _build_text(source):
    parts = []

    def _append_name(payload):
        if isinstance(payload, list):
            for item in payload:
                if isinstance(item, dict) and item.get("name"):
                    parts.append(item["name"])
                elif isinstance(item, str):
                    parts.append(item)
        elif isinstance(payload, dict) and payload.get("name"):
            parts.append(payload["name"])
        elif isinstance(payload, str):
            parts.append(payload)

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
        parts.append(name)
    if desc := source.get("description"):
        parts.append(desc)

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
def run_embeddings(es_hosts, index_name, log=None):
    """Generate and write embeddings for all applicable docs in *index_name*.

    Parameters
    ----------
    es_hosts : str or list[str]
        Elasticsearch host URL(s).
    index_name : str
        Target index name.
    log : logging.Logger, optional
        Logger instance; falls back to the module logger.
    """
    log = log or logger
    es = _make_es_client(es_hosts)

    for model_name, model_info in MODEL_CONFIGS.items():
        target_field = model_info["field_name"]
        log.info("=== Embedding model: %s  field: %s ===", model_name, target_field)

        _add_embedding_field(
            es,
            index_name,
            target_field,
            model_info["dims"],
            index_options=model_info.get("index_options"),
        )

        query = {
            "query": {
                "bool": {
                    "filter": [{"terms": {"@type": DOC_TYPES}}],
                    "must_not": [{"exists": {"field": target_field}}],
                }
            }
        }
        total = es.count(index=index_name, body=query).get("count", 0)
        if total == 0:
            log.info("All documents already have field '%s'. Skipping.", target_field)
            continue
        log.info("Found %d documents missing '%s'", total, target_field)

        embedded = 0
        skip_count = 0
        embed_batch: list[dict] = []
        bulk_actions: list[dict] = []
        start_time = time.time()

        def _flush_bulk():
            nonlocal bulk_actions
            if bulk_actions:
                helpers.bulk(es, bulk_actions, stats_only=True, request_timeout=120)
                bulk_actions = []

        def _flush_embeds():
            nonlocal embed_batch, embedded
            if not embed_batch:
                return
            texts = [item["text"] for item in embed_batch]
            vectors = _remote_embed(texts)
            for payload, vec in zip(embed_batch, vectors):
                bulk_actions.append(
                    {
                        "_op_type": "update",
                        "_index": index_name,
                        "_id": payload["doc_id"],
                        "doc": {target_field: vec},
                        "doc_as_upsert": True,
                    }
                )
                embedded += 1
                if embedded % 1000 == 0:
                    log.info(
                        "Embedded %d/%d (ETA %s)", embedded, total, _eta(embedded, total, start_time)
                    )
                if len(bulk_actions) >= BULK_UPDATE_SIZE:
                    _flush_bulk()
            embed_batch = []

        scan_count = 0
        for doc in helpers.scan(es, query=query, index=index_name, size=SCAN_BATCH_SIZE):
            scan_count += 1
            source = doc.get("_source", {})
            if target_field in source:
                skip_count += 1
                continue

            text = _build_text(source)
            if not text:
                continue

            embed_batch.append({"doc_id": doc["_id"], "text": text})
            if len(embed_batch) >= EMBED_BATCH_SIZE:
                _flush_embeds()

            if scan_count % 10000 == 0:
                log.info("Scanned %d/%d (ETA %s)", scan_count, total, _eta(scan_count, total, start_time))

        _flush_embeds()
        _flush_bulk()
        elapsed = _fmt_duration(time.time() - start_time)
        log.info(
            "Finished %s: %d/%d embeddings added, %d skipped (elapsed %s)",
            target_field,
            embedded,
            total,
            skip_count,
            elapsed,
        )
