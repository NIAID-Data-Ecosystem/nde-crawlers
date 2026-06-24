import datetime
import logging
import os
import sqlite3

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")


class NDEDatabase:
    EXPIRE = datetime.timedelta(days=30)
    SQL_DB = None
    NO_CACHE = False

    def __init__(self, sql_db=None):
        self.SQL_DB = self.SQL_DB or sql_db
        self.DIR_NAME = self.SQL_DB.split(".")[0]
        self.path = os.path.join("/cache/", self.DIR_NAME)
        os.makedirs(self.path, exist_ok=True)

    @property
    def db_path(self):
        return os.path.join(self.path, self.SQL_DB)

    def _metadata_value(self, cursor, name):
        cursor.execute("SELECT date FROM metadata WHERE name=?", (name,))
        row = cursor.fetchone()
        return row[0] if row else None

    def metadata_value(self, name):
        with sqlite3.connect(self.db_path) as con:
            return self._metadata_value(con.cursor(), name)

    def _upsert_metadata(self, name, value):
        with sqlite3.connect(self.db_path) as con:
            con.execute(
                """INSERT INTO metadata VALUES(?, ?)
                   ON CONFLICT(name) DO UPDATE SET date=excluded.date""",
                (name, value),
            )

    def is_cache_expired(self):
        with sqlite3.connect(self.db_path) as con:
            c = con.cursor()
            c.execute("""SELECT name FROM sqlite_master WHERE type='table' AND name='metadata' COLLATE NOCASE""")
            if not c.fetchone():
                return True

            if self._metadata_value(c, "load_complete") != "true":
                return True

            date_created = self._metadata_value(c, "date_created")
            if not date_created:
                return True

            today = datetime.date.today()
            return today - datetime.date.fromisoformat(date_created) >= self.EXPIRE

    def new_cache(self):
        logger.info("Cache does not exist or is expired. Creating new cache.")
        today = datetime.date.today().isoformat()

        with sqlite3.connect(self.db_path) as con:
            c = con.cursor()
            c.execute(
                """CREATE TABLE IF NOT EXISTS metadata (
                    name text NOT NULL PRIMARY KEY,
                    date text NOT NULL
                )"""
            )
            c.execute(
                """INSERT INTO metadata VALUES(?, ?)
                   ON CONFLICT(name) DO UPDATE SET date=excluded.date""",
                ("date_created", today),
            )
            c.execute(
                """INSERT INTO metadata VALUES(?, ?)
                   ON CONFLICT(name) DO UPDATE SET date=excluded.date""",
                ("date_updated", today),
            )
            c.execute(
                """INSERT INTO metadata VALUES(?, ?)
                   ON CONFLICT(name) DO UPDATE SET date=excluded.date""",
                ("load_complete", "false"),
            )
            c.execute("DROP TABLE IF EXISTS cache")
            c.execute(
                """CREATE TABLE cache (
                    _id text NOT NULL PRIMARY KEY,
                    data text NOT NULL
                )"""
            )

    def mark_cache_complete(self):
        self._upsert_metadata("load_complete", "true")

    def load_cache(self):
        raise NotImplementedError("Define in subclass")

    def update_cache(self):
        return None

    def retreive_last_updated(self):
        with sqlite3.connect(self.db_path) as con:
            c = con.cursor()
            last_updated = self._metadata_value(c, "date_updated")
            if not last_updated:
                raise RuntimeError("No date_updated value found in cache metadata.")
            return last_updated

    def insert_last_updated(self):
        today = datetime.date.today().isoformat()
        self._upsert_metadata("date_updated", today)
        logger.info("Cache last_updated on: %s", today)

    def dump(self, records):
        logger.info("Dumping records...")
        count = 0
        with sqlite3.connect(self.db_path) as con:
            c = con.cursor()
            for _id, data in records:
                if not (isinstance(_id, str) and isinstance(data, str)):
                    raise TypeError("_id and data must be a string")
                c.execute(
                    """INSERT INTO cache VALUES(?, ?)
                       ON CONFLICT(_id) DO UPDATE SET data=excluded.data""",
                    (_id, data),
                )
                count += 1
                if count % 1000 == 0:
                    con.commit()
            con.commit()
        logger.info("Finished dumping %s records.", count)

    def retreive_cache(self):
        con = sqlite3.connect(self.db_path)
        c = con.cursor()
        logger.info("Retrieving dumped records from cache...")
        c.execute("SELECT * FROM cache")
        count = 0
        try:
            for record in c:
                count += 1
                if count % 10000 == 0:
                    logger.info("Retrieving and parsing cache. Parsed %s records", count)
                yield record
        finally:
            logger.info("Finished Parsing. Total Records: %s", count)
            con.close()

    def parse(self, records):
        raise NotImplementedError("Define in subclass")

    def upload(self):
        if self.is_cache_expired() or self.NO_CACHE:
            self.new_cache()
            self.dump(self.load_cache())
            self.mark_cache_complete()
        else:
            records = self.update_cache()
            if records:
                logger.info("Updating cache with new records.")
                self.dump(records)
            else:
                logger.info("No new records to update cache.")

        return self.parse(self.retreive_cache())
