"""Prewarm description-extraction caches from an NDJSON file.

This runs the same description augmentation used by the uploader but discards
the transformed documents. The subsequent upload follows its normal code path
and finds the EXTRACT, species, disease, and negative lookup results cached.
"""

import argparse
import os
import time
from itertools import batched, islice

from config import logger

from .common import iter_ndjson
from .descriptions import augment_from_descriptions, reset_caches


def _augment_for_prewarm(docs):
    """Run description enrichment without the uploader's record-type filter."""
    return augment_from_descriptions(docs, filter_supported_types=False)


def prewarm_descriptions(path, batch_size=1000, limit=None, augment=_augment_for_prewarm):
    """Populate description-related caches from `path`; return documents read."""
    path = os.path.abspath(os.fspath(path))
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    if limit is not None and limit < 1:
        raise ValueError("limit must be at least 1")

    docs = iter_ndjson(os.path.dirname(path), os.path.basename(path))
    if limit is not None:
        docs = islice(docs, limit)

    reset_caches()
    started = time.monotonic()
    total = 0
    for batch in batched(docs, batch_size):
        augment(batch)
        total += len(batch)
        logger.info("Description cache prewarm: %s documents processed", total)

    logger.info(
        "Description cache prewarm finished: path=%s documents=%s elapsed=%.1fs",
        path,
        total,
        time.monotonic() - started,
    )
    return total


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Populate the existing EXTRACT and term-standardization caches from an NDJSON file. "
            "The input file is never modified."
        )
    )
    parser.add_argument("ndjson_path", help="Path to the NDJSON records that will enter the upload pipeline")
    parser.add_argument("--batch-size", type=int, default=1000, help="Documents per committed batch (default: 1000)")
    parser.add_argument("--limit", type=int, help="Only process the first N documents (useful for a trial run)")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    prewarm_descriptions(args.ndjson_path, batch_size=args.batch_size, limit=args.limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
