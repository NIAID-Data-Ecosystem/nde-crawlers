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
import re
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

        # Ensure all recently-indexed docs are visible before querying.
        es.indices.refresh(index=index_name)

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
