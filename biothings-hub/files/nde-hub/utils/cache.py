"""API answers cached across uploads, in one SQLite file.

Each table is (key, value, updated). `value` is JSON, or NULL for a miss: the
question was asked and got no answer. `updated` is the write time in epoch
seconds. Reads go to the file every time, so a row written by one batch or
worker is seen by the next.

A connection is kept per process and thread: worker processes are forked, and
an SQLite connection must not be used on both sides of a fork or across threads.
"""

import json
import os
import sqlite3
import threading
import time
from itertools import batched

from config import logger

DB_PATH = "/data/nde-hub/standardizers/cache.db"

# SQLite allows 999 bound variables by default; keep headroom.
_SQL_CHUNK_SIZE = 900
_BUSY_TIMEOUT_SECONDS = 60

_local = threading.local()
# Connections inherited from a parent process. Closing one could disturb the
# parent's, so they are kept referenced and never used.
_inherited = []


def _connection():
    """(connection, tables created on it) for this process and thread."""
    state = getattr(_local, "state", None)
    if state is not None:
        pid, path, conn, tables = state
        if pid == os.getpid() and path == DB_PATH:
            return conn, tables
        if pid == os.getpid():
            conn.close()
        else:
            _inherited.append(conn)
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    conn = sqlite3.connect(DB_PATH, timeout=_BUSY_TIMEOUT_SECONDS, isolation_level=None)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    _local.state = (os.getpid(), DB_PATH, conn, set())
    return conn, _local.state[3]


class Cache:
    """One table of key -> JSON answer, with misses and write times.

    `normalize` lower-cases and strips keys, on write and on read.
    `max_age` (seconds) limits a read to rows written that recently.
    """

    def __init__(self, table, normalize=True):
        self.table = table
        self.normalize = normalize

    def key(self, key):
        """`key` as stored."""
        return str(key).lower().strip() if self.normalize else str(key)

    def _conn(self):
        conn, tables = _connection()
        if self.table not in tables:
            conn.execute(
                f"CREATE TABLE IF NOT EXISTS {self.table} (key TEXT PRIMARY KEY, value TEXT, updated REAL NOT NULL)"
            )
            tables.add(self.table)
        return conn

    def _select(self, keys, condition, max_age):
        """(key, value) rows for `keys` matching `condition`, an SQL expression."""
        keys = {self.key(key) for key in keys}
        if not keys:
            return []
        cutoff = () if max_age is None else (time.time() - max_age,)
        if cutoff:
            condition += " AND updated > ?"
        rows = []
        try:
            conn = self._conn()
            for chunk in batched(sorted(keys), _SQL_CHUNK_SIZE):
                placeholders = ",".join("?" * len(chunk))
                rows += conn.execute(
                    f"SELECT key, value FROM {self.table} WHERE key IN ({placeholders}) AND {condition}",
                    (*chunk, *cutoff),
                ).fetchall()
        except (sqlite3.Error, OSError) as e:
            logger.error("Error reading %s: %s", self.table, e)
        return rows

    def _write(self, sql, rows):
        if not rows:
            return
        try:
            conn = self._conn()
            conn.execute("BEGIN IMMEDIATE")
            try:
                conn.executemany(sql, rows)
                conn.execute("COMMIT")
            except BaseException:
                if conn.in_transaction:
                    conn.execute("ROLLBACK")
                raise
        except (sqlite3.Error, OSError) as e:
            logger.error("Error writing to %s: %s", self.table, e)

    def get(self, key, max_age=None):
        """The answer for `key`, or None."""
        return self.get_many([key], max_age).get(self.key(key))

    def get_many(self, keys, max_age=None):
        """{key: answer} for the `keys` answered, keys as stored."""
        return {key: json.loads(value) for key, value in self._select(keys, "value IS NOT NULL", max_age)}

    def misses(self, keys, max_age=None):
        """The `keys` recorded as misses, keys as stored."""
        return {key for key, _ in self._select(keys, "value IS NULL", max_age)}

    def put(self, key, value):
        self.put_many({key: value})

    def put_many(self, values):
        """Store answers. An answer replaces a miss."""
        now = time.time()
        rows = [(self.key(key), json.dumps(value), now) for key, value in values.items()]
        self._write(f"INSERT OR REPLACE INTO {self.table} (key, value, updated) VALUES (?, ?, ?)", rows)

    def miss_many(self, keys):
        """Record misses. A miss never replaces an answer; a repeated miss gets a new time."""
        now = time.time()
        rows = [(key, now) for key in {self.key(key) for key in keys if key is not None} if key]
        self._write(
            f"INSERT INTO {self.table} (key, value, updated) VALUES (?, NULL, ?) "
            "ON CONFLICT(key) DO UPDATE SET updated = excluded.updated WHERE value IS NULL",
            rows,
        )

    def items(self):
        """[(key, answer)] for every answer, oldest row first."""
        try:
            rows = self._conn().execute(f"SELECT key, value FROM {self.table} WHERE value IS NOT NULL ORDER BY rowid")
            return [(key, json.loads(value)) for key, value in rows.fetchall()]
        except (sqlite3.Error, OSError) as e:
            logger.error("Error reading %s: %s", self.table, e)
            return []
