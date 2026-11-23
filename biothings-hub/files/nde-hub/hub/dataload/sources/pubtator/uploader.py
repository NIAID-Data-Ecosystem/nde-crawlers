"""PubTator3 annotations by PMID, read by the citations stage.

One document per PMID:
{"_id": <int pmid>, "annotations": [{"id": <concept id>, "mentions": [...]}]}.
"""

import csv
import gzip
import os
import re
from collections import Counter

import biothings
import config

biothings.config_for_app(config)
from biothings.hub.dataload.uploader import BaseSourceUploader
from biothings.utils.storage import MergerStorage
from config import logger

# Tab-separated columns: PMID, Type, Concept ID, Mentions ("|"-joined), Resource.
# Fields are not quoted; a mention can start with a literal '"'.
_PMID, _CONCEPT_ID, _MENTIONS = 0, 2, 3
_DIGITS = re.compile(r"[0-9]+")


def load_annotations(path):
    """Yield one document per run of consecutive rows with the same PMID.

    A concept repeated within a run keeps its last row.
    """
    counts = Counter()
    pmid, annotations = None, {}
    with gzip.open(path, "rt") as handle:
        for row in csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE):
            counts["rows"] += 1
            if len(row) <= _MENTIONS:
                counts["short_rows"] += 1
                continue
            if not _DIGITS.fullmatch(row[_PMID]):
                counts["bad_pmids"] += 1
                continue
            if row[_PMID] != pmid:
                if annotations:
                    counts["documents"] += 1
                    yield {"_id": int(pmid), "annotations": list(annotations.values())}
                pmid, annotations = row[_PMID], {}
            annotations[row[_CONCEPT_ID]] = {"id": row[_CONCEPT_ID], "mentions": row[_MENTIONS].split("|")}
    if annotations:
        counts["documents"] += 1
        yield {"_id": int(pmid), "annotations": list(annotations.values())}
    logger.info("PubTator %s: %s", os.path.basename(path), " ".join(f"{k}={v}" for k, v in sorted(counts.items())))


class PubTatorStorage(MergerStorage):
    """Merges a PMID that reappears in the file, keeping the later annotation per concept ID."""

    merges = 0

    @staticmethod
    def merge_func(new, existing, **kwargs):
        PubTatorStorage.merges += 1
        annotations = {annotation["id"]: annotation for annotation in existing.get("annotations", [])}
        annotations.update((annotation["id"], annotation) for annotation in new.get("annotations", []))
        existing["annotations"] = list(annotations.values())
        return existing

    def process(self, *args, **kwargs):
        PubTatorStorage.merges = 0
        total = super().process(*args, **kwargs)
        self.logger.info("Merged %s reappearing PMIDs", PubTatorStorage.merges)
        return total


class PubTatorUploader(BaseSourceUploader):
    main_source = "pubtator"
    storage_class = PubTatorStorage
    keep_archive = 0  # a reload can come from the dump folder
    filename = None

    __metadata__ = {
        "src_meta": {
            "url": "https://www.ncbi.nlm.nih.gov/research/pubtator3/",
            "description": "PubTator3 species and disease annotations of PubMed articles",
        }
    }

    def load_data(self, data_folder):
        yield from load_annotations(os.path.join(data_folder, self.filename))


class PubTatorSpeciesUploader(PubTatorUploader):
    name = "pubtator_species"
    filename = "species2pubtator3.gz"


class PubTatorDiseaseUploader(PubTatorUploader):
    name = "pubtator_disease"
    filename = "disease2pubtator3.gz"
