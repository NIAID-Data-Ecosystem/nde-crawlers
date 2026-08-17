"""Shared utilities for the NDE hub.

Uploaders need exactly one thing from here: `nde_upload_wrapper`, which runs an
uploader's records through the whole standardization pipeline. See
`utils/pipeline.py` for what the pipeline does and how a source can tune it.

The individual stages live in their own modules and are picked by the
pipeline for each record.
"""

from .common import as_list, iter_ndjson, retry
from .pipeline import STAGE_NAMES, finalize, nde_upload_wrapper, run_pipeline

__all__ = [
    "STAGE_NAMES",
    "as_list",
    "finalize",
    "iter_ndjson",
    "nde_upload_wrapper",
    "retry",
    "run_pipeline",
]
