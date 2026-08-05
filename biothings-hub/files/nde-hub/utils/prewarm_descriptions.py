"""Prewarm description-extraction caches from an NDJSON file.

This runs the same description augmentation used by the uploader but discards
the transformed documents. The subsequent upload follows its normal code path
and finds the EXTRACT, species, disease, and negative lookup results cached.
"""

import argparse
import os
import time
from itertools import batched, islice

import orjson
from config import logger

from .descriptions import augment_from_descriptions, reset_caches


def iter_ndjson(path):
    """Yield documents from the NDJSON file at `path`."""
    with open(path, "rb") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                yield orjson.loads(line)
            except orjson.JSONDecodeError as e:
                raise ValueError(f"Invalid JSON on line {line_number} of {path}: {e}") from e


def prewarm_descriptions(path, batch_size=1000, limit=None, augment=augment_from_descriptions):
    """Populate description-related caches from `path`; return documents read."""
    path = os.path.abspath(os.fspath(path))
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    if limit is not None and limit < 1:
        raise ValueError("limit must be at least 1")

    docs = iter_ndjson(path)
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
