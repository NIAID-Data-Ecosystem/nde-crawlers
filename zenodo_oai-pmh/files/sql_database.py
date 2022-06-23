import sqlite3
import datetime
import logging
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


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
        self.DIR_NAME = self.SQL_DB.split('.')[0]
        # puts directory into cache folder
        self.path = os.path.join('/cache/', self.DIR_NAME)
        # make the directory
        os.makedirs(self.path, exist_ok = True)

    def is_cache_expired(self):
        """ Uses the SQL_DB to connect to a sqlite db to check if the cache is expired in metadata table
            Returns:
                True: if cache does not exist or is expired, False otherwise
        """

        # connect to database
        con = sqlite3.connect(self.path + '/' + self.SQL_DB)
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
        con = sqlite3.connect(self.path + '/' + self.SQL_DB)
        c = con.cursor()

        logger.info("Cache does not exist or is expired. Creating new cache.")

        c.execute("""CREATE TABLE IF NOT EXISTS metadata (
                name text NOT NULL PRIMARY KEY,
                date text NOT NULL
                )""")

        today = datetime.date.today().isoformat()

        # used for testing
        # today = datetime.date(2022, 6, 7).isoformat()

        # upserting in sqlite https://www.sqlite.org/lang_UPSERT.html
        # https://stackoverflow.com/questions/62274285/sqlite3-programmingerror-incorrect-number-of-bindings-supplied-the-current-stat
        c.execute("""INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """, ('date_created', today))
        
        # add date_updated
        c.execute("""INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """, ('date_updated', today))
        con.commit()

        c.execute("DROP TABLE IF EXISTS cache")
        c.execute("""CREATE TABLE cache (
                      _id text NOT NULL PRIMARY KEY,
                      data text NOT NULL
                     )""")
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
        con = sqlite3.connect(self.path + '/' + self.SQL_DB)
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
        con = sqlite3.connect(self.path + '/' + self.SQL_DB)
        c = con.cursor()

        today = datetime.date.today().isoformat()
        # upsert date_updated to today
        c.execute("""INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """, ('date_updated', today))

        logger.info("Cache last_updated on: %s", today)
        con.commit()

        con.close()

    def dump(self, records):
        """Stores raw data from a list or generator into the cache table. Each value of the list contains (_id, data)"""
        
        # connect to database
        con = sqlite3.connect(self.path + '/' + self.SQL_DB)
        c = con.cursor()

        logger.info("Dumping records...")
        for record in records:
            _id = record[0]
            data = record[1]
            if not (isinstance(_id, str) and isinstance(data, str)):
                raise TypeError("_id and data must be a string")
            # insert doc into cache table _id first column second column data as a json string
            # we upsert here in case there is repeat ids from existing cache
            c.execute("""INSERT INTO cache VALUES(?, ?)
                            ON CONFLICT(_id) DO UPDATE SET data=excluded.data
                        """, (_id, data))
            con.commit()

        logger.info("Finished dumping records.")

        con.close()

    def retreive_cache(self):
        con = sqlite3.connect(self.path + '/' + self.SQL_DB)
        c = con.cursor()
        c.execute("SELECT * from cache")
        records = c.fetchall()
        con.close()
        return records

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
            self.dump(records)
            # retreives all the records from the cache table
            records = self.retreive_cache()
            # pipeline to transform data to put into the ndjson file
            return self.parse(records)
        