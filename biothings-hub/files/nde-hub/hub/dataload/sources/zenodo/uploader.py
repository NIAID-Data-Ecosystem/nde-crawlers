import os
import sqlite3

import orjson
from config import logger
from hub.dataload.nde import NDESourceUploader
from utils import as_list, iter_ndjson, nde_upload_wrapper


def _merge_versions(doc, other):
    """Keep the later-published version and append the other's sameAs and distribution to it."""
    # datePublished is an ISO date, so string order is date order. A missing date sorts oldest.
    if (doc.get("datePublished") or "") > (other.get("datePublished") or ""):
        newer, older = doc, other
    else:
        newer, older = other, doc
    for field in ("sameAs", "distribution"):
        if values := as_list(newer.get(field)) + as_list(older.get(field)):
            newer[field] = values
    return newer


class ZenodoUploader(NDESourceUploader):
    name = "zenodo"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        """Yield Zenodo Datasets, with the versions of each record merged into one document.

        Versions share a versionId (the concept DOI) and are merged in a scratch
        SQLite database. Records of any other @type are dropped before the pipeline.
        """
        db_path = os.path.join(data_folder, "zenodo.db")
        if os.path.exists(db_path):
            os.remove(db_path)
        con = sqlite3.connect(db_path)
        try:
            con.execute("CREATE TABLE zenodo (versionId TEXT NOT NULL PRIMARY KEY, doc BLOB NOT NULL)")

            count_total = count_unversioned = count_versions = count_merged = count_datasets = 0
            for count_total, doc in enumerate(iter_ndjson(data_folder), start=1):
                if count_total % 10000 == 0:
                    logger.info("Looping through ndjson: %s records", count_total)

                version_id = doc.pop("versionId", None)
                if not version_id:
                    count_unversioned += 1
                    if doc.get("@type") == "Dataset":
                        count_datasets += 1
                        yield doc
                    continue

                count_versions += 1
                row = con.execute("SELECT doc FROM zenodo WHERE versionId = ?", (version_id,)).fetchone()
                if row:
                    doc = _merge_versions(doc, orjson.loads(row[0]))
                else:
                    count_merged += 1
                con.execute(
                    "INSERT INTO zenodo VALUES (?, ?) ON CONFLICT(versionId) DO UPDATE SET doc = excluded.doc",
                    (version_id, orjson.dumps(doc)),
                )
            con.commit()

            logger.info(
                "Total records: %s. Records without a versionId: %s. Versions: %s, merged into %s records.",
                count_total,
                count_unversioned,
                count_versions,
                count_merged,
            )
            logger.info("Retrieving merged records from zenodo database...")

            # The merged record takes the @type of its latest version, so filter only after merging.
            for count_db, (blob,) in enumerate(con.execute("SELECT doc FROM zenodo"), start=1):
                if count_db % 10000 == 0:
                    logger.info("Retrieving records from zenodo database: %s records", count_db)
                doc = orjson.loads(blob)
                if doc.get("@type") == "Dataset":
                    count_datasets += 1
                    yield doc

            logger.info(
                "Finished. Uploaded %s Datasets. Dropped %s non-Dataset records.",
                count_datasets,
                count_unversioned + count_merged - count_datasets,
            )
        finally:
            con.close()
            os.remove(db_path)
