import json
import os
import sqlite3
from datetime import datetime

import orjson
from config import logger
from hub.dataload.nde import NDESourceUploader
from utils.utils import nde_upload_wrapper

# from utils.utils import zenodo_upload_wrapper


class ZenodoUploader(NDESourceUploader):
    name = "zenodo"

    # TODO RERUN ZENODO CRAWLER atm some do not have a type
    # @nde_upload_wrapper
    # @zenodo_upload_wrapper
    # def load_data(self, data_folder):
    #     with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
    #         for line in f:
    #             doc = orjson.loads(line)
    #             yield doc

    @nde_upload_wrapper
    def load_data(self, data_folder):
        with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
            # connect to database
            con = sqlite3.connect(data_folder + "/zenodo.db")
            c = con.cursor()
            c.execute("DROP TABLE IF EXISTS zenodo")
            c.execute(
                """CREATE TABLE zenodo (
                        versionId text NOT NULL PRIMARY KEY,
                        doc text NOT NULL
                        )"""
            )
            con.commit()

            count_uploaded = 0
            for count_total, line in enumerate(f, start=1):
                new_doc = orjson.loads(line)
                if version_id := new_doc.pop("versionId", None):
                    new_doc_str = json.dumps(new_doc)
                    is_newer = True
                    # check if there is a need to compare date_published.
                    c.execute("SELECT doc from zenodo WHERE versionId=(?)", (version_id,))
                    current_doc = c.fetchone()
                    if current_doc:
                        current_doc = json.loads(current_doc[0])

                        # Append sameAs to either new doc or current doc before upserting depending on datePublished
                        if datetime.fromisoformat(new_doc["datePublished"]) > datetime.fromisoformat(
                            current_doc["datePublished"]
                        ):
                            new_doc["sameAs"] += current_doc.get("sameAs")
                            new_doc_str = json.dumps(new_doc)
                        else:
                            current_doc["sameAs"] += new_doc.get("sameAs")
                            current_doc_str = json.dumps(current_doc)
                            is_newer = False

                    # insert into the database
                    if is_newer:
                        c.execute(
                            """INSERT INTO zenodo VALUES(?, ?)
                                        ON CONFLICT(versionId) DO UPDATE SET doc=excluded.doc
                                    """,
                            (version_id, new_doc_str),
                        )
                        con.commit()
                    else:
                        c.execute(
                            """INSERT INTO zenodo VALUES(?, ?)
                                        ON CONFLICT(versionId) DO UPDATE SET doc=excluded.doc
                                    """,
                            (version_id, current_doc_str),
                        )
                        con.commit()

                else:
                    count_uploaded += 1
                    # does not have a versionId just yield
                    yield new_doc

                if count_total % 10000 == 0:
                    logger.info("Looping through ndjson: %s records", count_total)
            count_inserted = count_total - count_uploaded
            logger.info(
                "Total records: %s. Uploaded records %s. Records inserted into zenodo database %s.",
                count_total,
                count_uploaded,
                count_inserted,
            )

            logger.info("Retrieving dumped records from zenodo database...")

            # loop through the database and upload remaining records
            c.execute("SELECT * from zenodo")
            for count_db, record in enumerate(c, start=1):
                if count_db % 10000 == 0:
                    logger.info("Retrieving records from zenodo sqlitedb. %s records", count_db)
                yield json.loads(record[1])
            logger.info(
                "Finished. Records retrieved and uploaded from database: %s. Total duplicate records %s",
                count_db,
                count_inserted - count_db,
            )
