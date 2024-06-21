import datetime
import functools
import logging
import os
import sqlite3
import time
import traceback

import dateutil
import requests
from dateutil.relativedelta import relativedelta
from sickle import Sickle

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


def handle_error(err, retry_sleep_sec):
    """Handle error for retrying on 502 and 504 errors"""
    if err.response.status_code in [502, 504]:
        logger.error("Could not connect to server. Retrying in %s seconds.", retry_sleep_sec)
        logger.error(err)
        logger.error(traceback.format_exc())
        time.sleep(retry_sleep_sec)
    else:
        raise err


def retry_oai(retry_num, retry_sleep_sec):
    """
    retry help decorator for OAI-PMH requests
    We only want to retry on 502 and 504 errors and the interval is days=1
    :param retry_num: the retry num; retry sleep sec
    :return: decorator
    """

    def decorator(func):
        """decorator"""

        # preserve information about the original function, or the func name will be "wrapper" not "func"
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """wrapper"""
            for attempt in range(retry_num):
                try:
                    yield from func(*args, **kwargs)  # need to yeild from the generator
                    break
                except requests.exceptions.HTTPError as err:
                    interval = kwargs.get("interval") or args[3]
                    granularity = kwargs.get("granularity") or args[4]
                    if granularity != "YYYY-MM-DDThh:mm:ssZ":
                        if interval and interval.get("days") == 1:
                            handle_error(err, retry_sleep_sec)
                        else:
                            raise err
                    else:
                        if interval and interval.get("hours") == 1:
                            handle_error(err, retry_sleep_sec)
                        else:
                            raise err
                if attempt + 1 == retry_num:
                    logger.error("func %s retry failed", func)
                    raise Exception("Exceed max retry num: {} failed".format(retry_num))

        return wrapper

    return decorator


class NDEDatabase:
    # override in subclass
    # how many days before cache expires
    EXPIRE = datetime.timedelta(days=30)
    # name of database
    SQL_DB = None
    # Never uses cache, defaults as false
    NO_CACHE = False

    def __init__(self, sql_db=None):
        # database name example: zenodo.db
        self.SQL_DB = self.SQL_DB or sql_db
        # parses directory name from database name
        self.DIR_NAME = self.SQL_DB.split(".")[0]
        # puts directory into cache folder
        self.path = os.path.join("/cache/", self.DIR_NAME)
        # make the directory
        os.makedirs(self.path, exist_ok=True)

    def is_cache_expired(self):
        """Uses the SQL_DB to connect to a sqlite db to check if the cache is expired in metadata table
        Returns:
            True: if cache does not exist or is expired, False otherwise
        """

        # connect to database
        con = sqlite3.connect(self.path + "/" + self.SQL_DB)
        c = con.cursor()

        # SQLite table names are case insensitive, but comparison is case sensitive by default.
        # To make this work properly in all cases you need to add COLLATE NOCASE.
        c.execute("""SELECT name FROM sqlite_master WHERE type='table' AND name='metadata' COLLATE NOCASE""")

        # if the table exists check date_created if not return True
        if c.fetchone():
            c.execute("""SELECT date FROM metadata WHERE name='date_created'""")
            date_created = c.fetchall()
            assert len(date_created) <= 1, "There is more than one date_created."
            date_created = datetime.date.fromisoformat(date_created[0][0])
            today = datetime.date.today()
            # compare if date is expired
            if today - date_created < self.EXPIRE:
                con.close()
                return False
        con.close()
        return True

    def new_cache(self):
        """Creates a new tables: metadata and cache. Upserts two entries: date_created, date_updated in metadata table"""

        # Read comments in base class for the first part
        con = sqlite3.connect(self.path + "/" + self.SQL_DB)
        c = con.cursor()

        logger.info("Cache does not exist or is expired. Creating new cache.")

        c.execute(
            """CREATE TABLE IF NOT EXISTS metadata (
                name text NOT NULL PRIMARY KEY,
                date text NOT NULL
                )"""
        )

        today = datetime.date.today().isoformat()

        # used for testing
        # today = datetime.date(2022, 6, 7).isoformat()

        # upserting in sqlite https://www.sqlite.org/lang_UPSERT.html
        # https://stackoverflow.com/questions/62274285/sqlite3-programmingerror-incorrect-number-of-bindings-supplied-the-current-stat
        c.execute(
            """INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """,
            ("date_created", today),
        )

        # add date_updated
        c.execute(
            """INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """,
            ("date_updated", today),
        )
        con.commit()

        c.execute("DROP TABLE IF EXISTS cache")
        c.execute(
            """CREATE TABLE cache (
                      _id text NOT NULL PRIMARY KEY,
                      data text NOT NULL
                     )"""
        )
        con.commit()

        con.close()

    def load_cache(self):
        """Download the request information"""
        raise NotImplementedError("Define in subclass")

    def update_cache(self):
        """Update cache with new data"""
        pass

    def retreive_last_updated(self):
        """Helper method for update_cache. Gets the date_updated value from metadata table."""

        # connect to database
        con = sqlite3.connect(self.path + "/" + self.SQL_DB)
        c = con.cursor()

        # checks date_updated
        c.execute("SELECT date from metadata WHERE name='date_updated'")
        last_updated = c.fetchall()
        assert len(last_updated) <= 1, "There is more than one last_updated."
        last_updated = last_updated[0][0]

        con.close()

        return last_updated

    def insert_last_updated(self):
        """Helper method for update_cache. Changes the date_updated in the metadata table to today"""
        # connect to database
        con = sqlite3.connect(self.path + "/" + self.SQL_DB)
        c = con.cursor()

        today = datetime.date.today().isoformat()
        # upsert date_updated to today
        c.execute(
            """INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """,
            ("date_updated", today),
        )

        logger.info("Cache last_updated on: %s", today)
        con.commit()

        con.close()

    def dump(self, records):
        """Stores raw data from a list or generator into the cache table. Each value of the list contains (_id, data)"""

        # connect to database
        con = sqlite3.connect(self.path + "/" + self.SQL_DB)
        c = con.cursor()

        logger.info("Dumping records...")
        for record in records:
            _id = record[0]
            data = record[1]
            if not (isinstance(_id, str) and isinstance(data, str)):
                raise TypeError("_id and data must be a string")
            # insert doc into cache table _id first column second column data as a json string
            # we upsert here in case there is repeat ids from existing cache
            c.execute(
                """INSERT INTO cache VALUES(?, ?)
                            ON CONFLICT(_id) DO UPDATE SET data=excluded.data
                        """,
                (_id, data),
            )
            con.commit()

        logger.info("Finished dumping records.")

        con.close()

    def retreive_cache(self):
        con = sqlite3.connect(self.path + "/" + self.SQL_DB)
        c = con.cursor()
        logger.info("Retrieving dumped records from cache...")
        c.execute("SELECT * from cache")
        count = 0
        for record in c:
            count += 1
            if count % 10000 == 0:
                logger.info("Retrieving and parsing cache. Parsed %s records", count)
            yield record
        logger.info("Finished Parsing. Total Records: %s", count)
        con.close()

    def parse(self, records):
        """Parse the request information"""
        raise NotImplementedError("Define in subclass")

    def upload(self):
        """Checks if cache is expired or NO_CACHE is True if so make a new cache and dump else update existing cache"""
        if self.is_cache_expired() or self.NO_CACHE:
            # make new cache
            self.new_cache()
            # api request to get all the data in a generator. Each entry (_id, data) -> (str, str)
            records = self.load_cache()
            # dumps all records into the cache table of sql database
            self.dump(records)
            # retreives all the records from the cache table
            records = self.retreive_cache()
            # pipeline to transform data to put into the ndjson file
            return self.parse(records)
        else:
            # gets new records
            records = self.update_cache()
            # dumps all new records into cache table that has existing data
            if records:
                logger.info("Updating cache with new records.")
                self.dump(records)
            else:
                logger.info("No new records to update cache.")
            # retreives all the records from the cache table
            records = self.retreive_cache()
            # pipeline to transform data to put into the ndjson file
            return self.parse(records)


class OAIDatabase(NDEDatabase):
    """Subclass of NDEDatabase that uses OAI-PMH to get data. Uses sickle.ListRecords to get data.
    START: datetime.date object of the start date
    END: datetime.date object of the end date.
    METADATA_PREFIX: str of the metadata prefix
    INTERVAL: dictionary of the interval to query by can be years or days
    HOST: str of the host url
    SLEEP_COUNT: int of the number of records to load before sleeping
    SLEEP_LENGTH: int of the number of seconds to sleep
    """

    START = None
    # end date needs to be tomorrow since it is not inclusive
    END = datetime.datetime.now().replace(microsecond=0) + datetime.timedelta(days=1)
    METADATA_PREFIX = None
    INTERVAL = None
    HOST = None
    SLEEP_COUNT = 100
    SLEEP_LENGTH = 1

    def __init__(self, sql_db=None):
        super().__init__(sql_db)
        self.sickle = Sickle(self.HOST, max_retries=4, default_retry_after=10)
        self.identify = self.sickle.Identify()
        self.granularity = self.identify.granularity
        self.start = self.START or self.identify.earliestDatestamp
        self.end = self.END
        assert (
            "hours" or "days" or "years" in self.INTERVAL.keys()
        ), "INTERVAL format only permits hours, days, or years"
        if self.granularity != "YYYY-MM-DDThh:mm:ssZ":
            self.start = self.start.date()
            self.end = self.end.date()
            if "hours" in self.INTERVAL.keys():
                raise AssertionError("OAI format does not permit hourly intervals")

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, value):
        if isinstance(value, str):
            self._start = dateutil.parser.parse(value, ignoretz=True)
        else:
            self._start = value

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, value):
        if value:
            if isinstance(value, str):
                self._end = dateutil.parser.parse(value, ignoretz=True)
            else:
                self._end = value
        else:
            self._end = datetime.datetime.now().replace(microsecond=0) + datetime.timedelta(days=1)

    def record_data(self, record):
        """Choose what information in the request record to yield into the cache table
        Args:   records: generator of records
        """
        raise NotImplementedError("Define in Subclass")

    @retry_oai(retry_num=5, retry_sleep_sec=10)
    def request_data(self, start, until, interval, granularity):
        """Request data from the sickle object
        Args:
            sickle: sickle object
            start: datetime object of the start date
            until: datetime object of the end date
        """
        records = self.sickle.ListRecords(
            **{
                "metadataPrefix": self.METADATA_PREFIX,
                "from": start.isoformat(),
                "until": until.isoformat(),
                "ignore_deleted": True,
            }
        )
        for count, record in enumerate(records, start=1):
            if count % self.SLEEP_COUNT == 0:
                time.sleep(self.SLEEP_LENGTH)
                logger.info("Loading cache... Loaded %s records", count)
            yield self.record_data(record)
        logger.info("Total Records for this interval: %s", count)

    def load_cache(self, start=None):
        if start:
            self.start = start
        # since dictionaries are mutable we need to make a copy
        interval = self.INTERVAL.copy()

        while self.start < self.end:
            try:
                until = min(self.start + relativedelta(**interval), self.end)
                logger.info("Start: %s, Until: %s, Interval: %s", self.start, until, interval)
                records = self.request_data(self.start, until, interval, self.granularity)
                for record in records:
                    yield record
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 422 and "UNPROCESSABLE ENTITY" in e.response.reason:
                    logger.info("No records found for this url %s. Continue to next interval.", e.response.url)
                    self.start = until
                    interval = self.INTERVAL.copy()
                else:
                    # break the while loop failed to many times
                    if self.granularity != "YYYY-MM-DDThh:mm:ssZ":
                        if interval.get("days") and interval.get("days") // 2 == 0:
                            logger.error(traceback.format_exc())
                            raise e
                    elif interval.get("hours") and interval.get("hours") // 2 == 0:
                        logger.error(traceback.format_exc())
                        raise e
                    # when querying by year fails change to days
                    if "years" in interval.keys():
                        if interval["years"] // 2 == 0:
                            interval["days"] = 365
                            del interval["years"]
                    if "days" in interval.keys():
                        if interval["days"] // 2 == 0:
                            interval["hours"] = 24
                            del interval["days"]
                    # cut interval in half every time it fails
                    interval = {k: max(1, v // 2) for k, v in interval.items()}
            else:
                self.start = until
                interval = self.INTERVAL.copy()

    def update_cache(self):
        last_updated = self.retreive_last_updated()

        # get all the records since last_updated to add into current cache
        logger.info("Updating cache from %s", last_updated)

        records = self.load_cache(start=last_updated)
        for record in records:
            yield record

        self.insert_last_updated()
