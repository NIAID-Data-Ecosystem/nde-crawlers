"""Make Data Count links for the make_data_count pipeline stage.

One document per kept link, `_id` "<dataset_key>|<publication_doi>|<citation_type>",
indexed by `dataset_key`.
"""

import csv
import gzip
import json
import os
import re
from collections import Counter

import biothings
import config

biothings.config_for_app(config)
from biothings.hub.dataload.uploader import BaseSourceUploader
from biothings.utils.mongo import get_src_db
from biothings.utils.storage import BasicStorage
from config import logger
from utils.make_data_count import COLLECTION, kept_links

from .dumper import RELEASE_FILE


class MakeDataCountUploader(BaseSourceUploader):
    name = COLLECTION
    storage_class = BasicStorage
    keep_archive = 1

    __metadata__ = {
        "src_meta": {
            "url": "https://makedatacount.org/",
            "description": "Make Data Count article-to-dataset citation links",
        }
    }

    def load_data(self, data_folder):
        get_src_db()[self.temp_collection_name].create_index("dataset_key")
        with open(os.path.join(data_folder, RELEASE_FILE)) as handle:
            release = json.load(handle)
        if not re.fullmatch(r"\d{6}", release["release"]):
            raise ValueError(f"Unexpected Make Data Count release: {release['release']!r}")

        counts = Counter()
        seen = set()
        path = os.path.join(data_folder, f"{release['table']}.csv.gz")
        with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
            for key, doi, citation_type, repository in kept_links(csv.DictReader(handle), counts):
                _id = f"{key}|{doi}|{citation_type}"
                if _id in seen:
                    counts["repeated"] += 1
                    continue
                seen.add(_id)
                yield {
                    "_id": _id,
                    "dataset_key": key,
                    "publication_doi": doi,
                    "citation_type": citation_type,
                    "repository": repository,
                    "release": release["release"],
                    "created": release["created"],
                }

        if counts["read"] != release["row_count"]:
            raise RuntimeError(
                f"Read {counts['read']:,} rows of release {release['release']}, expected {release['row_count']:,}"
            )
        if not counts["kept"]:
            raise RuntimeError(f"No usable links in release {release['release']}")
        logger.info(
            "Make Data Count release %s (created %s): rows=%s links_kept=%s repeated=%s",
            release["release"],
            release["created"],
            counts["read"],
            counts["kept"],
            counts["repeated"],
        )
