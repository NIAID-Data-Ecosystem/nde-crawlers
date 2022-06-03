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
    DBM_NAME = None
    # Never uses cache, defaults as false
    NO_CACHE = False

    def __init__(self, dbm_name=None):
        # database name example: zenodo.db
        self.DBM_NAME = self.DBM_NAME or dbm_name
        # parses directory name from database name
        self.DIR_NAME = self.DBM_NAME.split('.')[0]
        # puts directory into cache folder
        self.path = os.path.join('/cache/', self.DIR_NAME)
        # make the directory
        os.makedirs(self.path, exist_ok = True)

    def is_cache_expired(self):
        """ Uses the DBM_NAME to connect to a sqlite db to check if the cache is expired in metadata table
            Returns:
                True: if cache does not exist or is expired, False otherwise
        """

        # connect to database
        con = sqlite3.connect(self.path + '/' + self.DBM_NAME)
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
        """Creates a new table metadata in the sqlite db using DBM_NAME. Upserts an entry: date_created"""

        # connect to database
        con = sqlite3.connect(self.path + '/' + self.DBM_NAME)
        c = con.cursor()

        # Creates a new table metadata if it does not exist
        c.execute("""CREATE TABLE IF NOT EXISTS metadata (
                name text NOT NULL PRIMARY KEY,
                date text NOT NULL
                )""")

        today = datetime.date.today().isoformat()

        # used for testing
        # today= datetime.date(2022, 4, 10).isoformat()

        # upserting in sqlite https://www.sqlite.org/lang_UPSERT.html
        # https://stackoverflow.com/questions/62274285/sqlite3-programmingerror-incorrect-number-of-bindings-supplied-the-current-stat
        c.execute("""INSERT INTO metadata VALUES(?, ?)
                        ON CONFLICT(name) DO UPDATE SET date=excluded.date
                  """, ('date_created', today))

        con.commit()
        
        # delete if cache table exists and add a new cache table
        c.execute("DROP TABLE IF EXISTS cache")
        c.execute("""CREATE TABLE cache (
                      _id text NOT NULL PRIMARY KEY,
                      data text NOT NULL
                     )""")
        con.commit()
        
        con.close()

    def update_cache(self):
        """Update cache with new data"""
        pass

    def dump(self):
        """Download the request information"""
        raise NotImplementedError("Define in subclass")

    def parse(self):
        """Parse the request information"""
        raise NotImplementedError("Define in subclass")

    def upload(self):
        """Checks if cache is expired or NO_CACHE is True if so make a new cache and dump else update existing cache"""
        if self.is_cache_expired() or self.NO_CACHE:
            self.new_cache()
            self.dump()
            return self.parse()
        else:
            self.update_cache()
            return self.parse()





