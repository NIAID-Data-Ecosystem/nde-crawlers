"""The single upload pipeline that every NDE source runs through.

`nde_upload_wrapper` decorates an uploader's `load_data` and runs the records
through every applicable stage. A stage is applicable when

  * the records in front of it carry the field it works on (`funding`,
    `species`, `pmids`, ...), and
  * for the stages driven by a curated file (`measurementTechnique`,
    `topicCategory`, ...), that file exists for this source.

An inapplicable stage costs one dict lookup per record and does not open its
lookup tables, a database connection or the network.

Records flow through in batches, so the batch-oriented stages have enough
records to make their bulk lookups worthwhile while memory stays bounded
regardless of source size.

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
    """One augmentation step: what it is called, how to run it, and when it applies."""

    __slots__ = ("name", "_run", "_applies", "_lookup_file", "_reset")

    def __init__(self, name, run, applies, lookup_file=None, reset=None):
        self.name = name
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
    Stage("citations", _run_citations, _any_doc(_needs_citations)),
    Stage("funding", _run_funding, _any_doc(_needs_funding)),
    Stage("terms", _run_terms, _any_doc(_needs_terms), reset=_reset_terms),
    Stage("descriptions", _run_descriptions, _any_doc(_needs_descriptions), reset=_reset_descriptions),
    Stage(
        "measurement_technique",
        _run_measurement_technique,
        _any_doc(lambda doc: bool(doc.get("measurementTechnique"))),
        _measurement_technique_file,
    ),
    Stage("nctid", _run_nctid, _any_doc(lambda doc: bool(doc.get("nctid"))), _nctid_file),
    Stage("topic_category", _run_topic_category, _always, _topic_category_file),
    Stage("disambiguating_description", _run_disambiguating_description, _always, _disambiguating_description_file),
    Stage("lineage", _run_lineage, _always),
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
    total = 0
    yielded = 0

    for batch in batched(docs, batch_size):
        total += len(batch)
        for stage in stages:
            if not stage.applies(batch):
                continue
            batch = list(stage.run(batch, source))
            ran[stage.name] = ran.get(stage.name, 0) + 1

        for doc in batch:
            doc = finalize(doc)
            if doc is not None:
                yielded += 1
                yield doc

        logger.info("Pipeline: %s documents processed (%s emitted)", total, yielded)

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
