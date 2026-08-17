"""The single upload pipeline that every NDE source runs through.

`nde_upload_wrapper` decorates an uploader's `load_data` and runs the records
through every applicable stage. A stage is applicable when

  * the records in front of it carry the field it works on (`funding`,
    `species`, `pmids`, ...), and
  * for the stages driven by a curated file (`measurementTechnique`,
    `topicCategory`, ...), that file exists for this source.

An inapplicable stage costs one dict lookup per record and does not open its
lookup tables, a database connection or the network.

Records are processed in batches of 1000 by default.

At the end of an upload, each active stage logs repository-wide before/after
counts for the fields it manages, along with the number of records it changed.

Usage in an uploader::

    class MySourceUploader(NDESourceUploader):
        name = "my_source"

The base class already decorates `load_data`. Override it (keeping
`@nde_upload_wrapper`) when the records need custom parsing. Two optional
settings on the uploader change the pipeline:

    post_process(self, doc)      applied after every stage but before
                                 `lineage`; return the document or None to
                                 drop it
    skip_stages = ("...",)       stage names this source should not run
"""

import copy
import functools
import os
import time
from itertools import batched

import bson
from config import logger

from .common import dict_entries, supports_description_enrichment
from .corrections import apply_corrections
from .validate import add_date, add_metadata_score, check_schema, clean_description, drop_placeholder_terms

DEFAULT_BATCH_SIZE = 1000
MONGO_DOC_SIZE_LIMIT = 16 * 1024 * 1024


# ---------------------------------------------------------------------------
# Stage definition
# ---------------------------------------------------------------------------
class Stage:
    __slots__ = ("name", "tracked_fields", "_run", "_applies", "_lookup_file", "_reset")

    def __init__(self, name, run, applies, lookup_file=None, reset=None, tracked_fields=()):
        self.name = name
        self.tracked_fields = tuple(tracked_fields)
        self._run = run
        self._applies = applies
        self._lookup_file = lookup_file
        self._reset = reset

    def run(self, docs, source):
        return self._run(docs, source)

    def applies(self, batch):
        return self._applies(batch)

    def reset(self):
        """Drop any lookup data cached by a previous upload in this process."""
        if self._reset is not None:
            self._reset()

    def available(self, source):
        """False when this stage is driven by a per-source lookup file the source has none of."""
        if self._lookup_file is None:
            return True
        path = self._lookup_file(source)
        if path and os.path.exists(path):
            return True
        logger.debug("Pipeline: no %s lookup file for %s (%s), skipping stage", self.name, source, path)
        return False


# ---------------------------------------------------------------------------
# Per-stage statistics
# ---------------------------------------------------------------------------
def _new_stage_stats(stage):
    return {
        "records_seen": 0,
        "records_processed": 0,
        "records_output": 0,
        "records_changed": 0,
        "records_added": 0,
        "records_dropped": 0,
        "fields": {
            field: {
                "records_before": 0,
                "records_after": 0,
                "values_before": 0,
                "values_after": 0,
                "changed": 0,
                "added": 0,
                "modified": 0,
                "removed": 0,
            }
            for field in stage.tracked_fields
        },
    }


def _record_key(doc, occurrences):
    identifier = doc.get("_id") if isinstance(doc, dict) else None
    try:
        hash(identifier)
    except TypeError:
        identifier = repr(identifier)
    base = ("_id", identifier) if identifier is not None else ("object", id(doc))
    occurrence = occurrences.get(base, 0)
    occurrences[base] = occurrence + 1
    return base, occurrence


def _field_state(doc, field, copy_value):
    value = doc
    for part in field.split("."):
        if not isinstance(value, dict) or part not in value:
            return False, None
        value = value[part]
    return True, copy.deepcopy(value) if copy_value else value


def _snapshot_docs(docs, fields, *, copy_values):
    occurrences = {}
    return {
        _record_key(doc, occurrences): {field: _field_state(doc, field, copy_values) for field in fields}
        for doc in docs
    }


def _value_count(state):
    present, value = state
    if not present or not value:
        return 0
    return len(value) if isinstance(value, (list, tuple, set)) else 1


def _update_stage_stats(stats, before, after, *, applied):
    stats["records_seen"] += len(before)
    if applied:
        stats["records_processed"] += len(before)
        stats["records_output"] += len(after)

    before_keys = set(before)
    after_keys = set(after)
    stats["records_added"] += len(after_keys - before_keys)
    stats["records_dropped"] += len(before_keys - after_keys)

    changed_records = before_keys ^ after_keys
    for key in before_keys & after_keys:
        if any(before[key][field] != after[key][field] for field in stats["fields"]):
            changed_records.add(key)
    stats["records_changed"] += len(changed_records)

    absent = (False, None)
    for field, field_stats in stats["fields"].items():
        for key in before_keys | after_keys:
            before_state = before.get(key, {}).get(field, absent)
            after_state = after.get(key, {}).get(field, absent)
            before_count = _value_count(before_state)
            after_count = _value_count(after_state)

            field_stats["records_before"] += bool(before_count)
            field_stats["records_after"] += bool(after_count)
            field_stats["values_before"] += before_count
            field_stats["values_after"] += after_count

            if before_state == after_state:
                continue
            field_stats["changed"] += 1
            if not before_count and after_count:
                field_stats["added"] += 1
            elif before_count and not after_count:
                field_stats["removed"] += 1
            else:
                field_stats["modified"] += 1


def _log_stage_stats(source, stages, stage_stats):
    repository = source or "unknown"
    for stage in stages:
        stats = stage_stats[stage.name]
        logger.info(
            "Pipeline util stats: repository=%s util=%s records_seen=%s records_processed=%s "
            "records_output=%s records_changed=%s records_added=%s records_dropped=%s",
            repository,
            stage.name,
            stats["records_seen"],
            stats["records_processed"],
            stats["records_output"],
            stats["records_changed"],
            stats["records_added"],
            stats["records_dropped"],
        )
        for field, field_stats in stats["fields"].items():
            logger.info(
                "Pipeline util stats: repository=%s util=%s field=%s records=%s->%s values=%s->%s "
                "changed=%s added=%s modified=%s removed=%s",
                repository,
                stage.name,
                field,
                field_stats["records_before"],
                field_stats["records_after"],
                field_stats["values_before"],
                field_stats["values_after"],
                field_stats["changed"],
                field_stats["added"],
                field_stats["modified"],
                field_stats["removed"],
            )


def _any_doc(predicate):
    def applies(batch):
        return any(predicate(doc) for doc in batch)

    return applies


def _always(batch):
    return True


# ---------------------------------------------------------------------------
# When each stage applies
# ---------------------------------------------------------------------------
def _needs_citations(doc):
    # Key presence, not truthiness: the stage is also what strips `pmids` and
    # `pmcs` off the record, and neither belongs in the index.
    return "pmids" in doc or "pmcs" in doc or any(entry.get("doi") for entry in dict_entries(doc, "citation"))


def _needs_funding(doc):
    return bool(doc.get("funding"))


def _needs_terms(doc):
    return bool(doc.get("species") or doc.get("infectiousAgent") or doc.get("healthCondition"))


def _needs_descriptions(doc):
    # EXTRACT only mines a record for the entities it is missing, so a record
    # that already carries both taxonomy and health conditions has nothing to gain.
    if not supports_description_enrichment(doc) or not doc.get("description"):
        return False
    return ("species" not in doc and "infectiousAgent" not in doc) or "healthCondition" not in doc


# ---------------------------------------------------------------------------
# Stage runners. Each imports its module on call, so a hub process only loads
# text2term, pandas, rdflib and Bio for the stages a source actually runs.
# ---------------------------------------------------------------------------
def _run_citations(docs, source):
    from .citations import add_citations

    return add_citations(docs)


def _run_funding(docs, source):
    from .funding import standardize_funding

    return standardize_funding(docs)


def _run_terms(docs, source):
    from .terms import standardize_terms

    return standardize_terms(docs)


def _run_descriptions(docs, source):
    from .descriptions import augment_from_descriptions

    # Stage.applies only gates the batch. Filter again here so an eligible
    # Dataset in a mixed batch does not send neighbouring ineligible records
    # (for example ordinary Samples) to EXTRACT.
    doc_list = list(docs)
    eligible = [doc for doc in doc_list if _needs_descriptions(doc)]
    if eligible:
        augment_from_descriptions(eligible)
    return doc_list


def _run_measurement_technique(docs, source):
    from .measurement_technique import process_measurement_technique

    return process_measurement_technique(docs, source)


def _run_nctid(docs, source):
    from .nctid import add_nct_measurement_techniques

    return add_nct_measurement_techniques(docs)


def _run_topic_category(docs, source):
    from .topic_category import add_topic_category

    return add_topic_category(docs, source)


def _run_disambiguating_description(docs, source):
    from .disambiguating_description import add_disambiguating_description

    return add_disambiguating_description(docs, source)


def _run_lineage(docs, source):
    from .lineage import process_lineage

    return process_lineage(docs)


def _reset_terms():
    from .terms import reset_caches

    reset_caches()


def _reset_descriptions():
    from .descriptions import reset_caches

    reset_caches()


def _measurement_technique_file(source):
    from .measurement_technique import lookup_file

    return lookup_file(source)


def _nctid_file(source):
    from .nctid import lookup_file

    return lookup_file()


def _topic_category_file(source):
    from .topic_category import lookup_file

    return lookup_file(source)


def _disambiguating_description_file(source):
    from .disambiguating_description import lookup_file

    return lookup_file(source)


# The pipeline, in order. Citations run first because they add funding, species
# and health conditions the later stages then standardize.
STAGES = (
    Stage(
        "citations",
        _run_citations,
        _any_doc(_needs_citations),
        tracked_fields=("citation", "funding", "species", "infectiousAgent", "healthCondition"),
    ),
    Stage("funding", _run_funding, _any_doc(_needs_funding), tracked_fields=("funding",)),
    Stage(
        "terms",
        _run_terms,
        _any_doc(_needs_terms),
        reset=_reset_terms,
        tracked_fields=("species", "infectiousAgent", "healthCondition"),
    ),
    Stage(
        "descriptions",
        _run_descriptions,
        _any_doc(_needs_descriptions),
        reset=_reset_descriptions,
        tracked_fields=("species", "infectiousAgent", "healthCondition"),
    ),
    Stage(
        "measurement_technique",
        _run_measurement_technique,
        _any_doc(lambda doc: bool(doc.get("measurementTechnique"))),
        _measurement_technique_file,
        tracked_fields=("measurementTechnique", "keywords"),
    ),
    Stage(
        "nctid",
        _run_nctid,
        _any_doc(lambda doc: bool(doc.get("nctid"))),
        _nctid_file,
        tracked_fields=("measurementTechnique",),
    ),
    Stage(
        "topic_category",
        _run_topic_category,
        _always,
        _topic_category_file,
        tracked_fields=("topicCategory",),
    ),
    Stage(
        "disambiguating_description",
        _run_disambiguating_description,
        _always,
        _disambiguating_description_file,
        tracked_fields=("disambiguatingDescription",),
    ),
    Stage("lineage", _run_lineage, _always, tracked_fields=("_meta.lineage",)),
)

STAGE_NAMES = tuple(stage.name for stage in STAGES)


def _post_process_stage(post_process):
    """Wrap an uploader's `post_process` method as a stage.

    It runs after every augmentation but before `lineage`, so a source can still
    add, correct or drop taxonomy and have the lineage reflect it.
    """

    def run(docs, source):
        for doc in docs:
            doc = post_process(doc)
            if doc is not None:
                yield doc

    return Stage("post_process", run, _always)


# ---------------------------------------------------------------------------
# The pipeline
# ---------------------------------------------------------------------------
def run_pipeline(docs, source=None, skip=(), batch_size=None, post_process=None):
    """Run `docs` through every applicable stage, then finish each document.

    Yields documents ready for MongoDB. Oversized documents are dropped with a
    warning; everything else is validated by `check_schema` and will raise if
    the source emitted something invalid.
    """
    unknown = set(skip) - set(STAGE_NAMES)
    if unknown:
        raise ValueError(
            "Unknown pipeline stage(s) in skip_stages: %s. Known stages: %s" % (sorted(unknown), list(STAGE_NAMES))
        )

    stages = []
    skipped = {}
    for stage in STAGES:
        if stage.name in skip:
            skipped[stage.name] = "skip_stages"
        elif not stage.available(source):
            skipped[stage.name] = "no lookup file"
        else:
            stages.append(stage)

    for stage in stages:
        stage.reset()

    if post_process is not None:
        lineage_index = next((i for i, stage in enumerate(stages) if stage.name == "lineage"), len(stages))
        stages.insert(lineage_index, _post_process_stage(post_process))

    logger.info("Pipeline for %s: %s", source, ", ".join(stage.name for stage in stages) or "no stages")

    batch_size = batch_size or DEFAULT_BATCH_SIZE
    started = time.monotonic()
    ran = {}
    stage_stats = {stage.name: _new_stage_stats(stage) for stage in stages}
    total = 0
    yielded = 0

    try:
        for batch in batched(docs, batch_size):
            total += len(batch)
            for stage in stages:
                applied = stage.applies(batch)
                before = _snapshot_docs(batch, stage.tracked_fields, copy_values=applied)
                if applied:
                    batch = list(stage.run(batch, source))
                    ran[stage.name] = ran.get(stage.name, 0) + 1
                    after = _snapshot_docs(batch, stage.tracked_fields, copy_values=False)
                else:
                    after = before
                _update_stage_stats(stage_stats[stage.name], before, after, applied=applied)

            for doc in batch:
                doc = finalize(doc)
                if doc is not None:
                    yielded += 1
                    yield doc

            logger.info("Pipeline: %s documents processed (%s emitted)", total, yielded)
    finally:
        # Also report partial statistics when validation or upload consumption
        # stops the generator before the repository has finished.
        _log_stage_stats(source, stages, stage_stats)

    # A stage that was active but never applied had no record that needed it.
    skipped.update({stage.name: "no matching records" for stage in stages if stage.name not in ran})

    logger.info(
        "Pipeline for %s finished: %s documents, %s emitted, %.1fs",
        source,
        total,
        yielded,
        time.monotonic() - started,
    )
    logger.info("Pipeline for %s ran: %s", source, ran or "no stages")
    logger.info("Pipeline for %s skipped: %s", source, skipped or "nothing")


def finalize(doc):
    """Apply the per-document finishing touches. Returns None for oversized documents."""
    # Apply sourceOrganization corrections on the fly, using the index cached
    # once per process. Matches by _id (records.txt) AND funding.identifier.
    apply_corrections(doc)
    add_date(doc)
    add_metadata_score(doc)
    clean_description(doc)
    drop_placeholder_terms(doc)

    # final validation and size check
    check_schema(doc)
    doc["_id"] = doc["_id"].casefold()
    bson_size = len(bson.BSON.encode(doc))
    if bson_size >= MONGO_DOC_SIZE_LIMIT:
        logger.warning("Document %s exceeds MongoDB's size limit: %s bytes", doc["_id"], bson_size)
        return None
    return doc


def nde_upload_wrapper(func=None, *, skip=(), batch_size=None):
    """Run an uploader's `load_data` output through the upload pipeline.

    Usable bare (`@nde_upload_wrapper`) or configured
    (`@nde_upload_wrapper(skip=("descriptions",))`). Configuration also comes
    from the uploader itself, so sources that don't override `load_data` can
    still tune the pipeline:

        skip_stages       -- stage names this source should not run
        post_process(doc) -- last word on each document; return None to drop it
    """

    def decorate(load_data):
        @functools.wraps(load_data)
        def wrapper(*args, **kwargs):
            uploader = args[0] if args else None
            source = getattr(uploader, "main_source", None) or getattr(uploader, "name", None)
            skipped = set(skip) | set(_source_attr(uploader, "skip_stages") or ())
            yield from run_pipeline(
                load_data(*args, **kwargs),
                source=source,
                skip=skipped,
                batch_size=batch_size,
                post_process=_source_attr(uploader, "post_process"),
            )

        return wrapper

    return decorate(func) if callable(func) else decorate


def _source_attr(uploader, name):
    """Read a pipeline setting the source itself declares.

    Only looks at classes above biothings' own uploader base classes, so a
    setting is never picked up from whatever biothings happens to define.
    """
    for klass in type(uploader).__mro__:
        if klass.__module__.startswith("biothings"):
            return None
        if name in klass.__dict__:
            return getattr(uploader, name)
    return None
