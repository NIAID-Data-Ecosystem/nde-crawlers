"""SQLite-backed caches shared by the standardization stages.

Every stage that resolves a term against an external service keeps the answer in
a SQLite table so later runs don't ask again, with an in-memory cache in front of
it for the current run. These two classes are that pattern:

    SqliteCache   -- key -> JSON value (a standardized term, a grant, a citation)
    SqliteKeySet  -- keys known to resolve to nothing, so we stop retrying them

Caching strategy per instance:

    preload=True                 the whole table is read on first use. For tables
                                 small enough to hold and looked up constantly.
    preload=False, memoize=True  rows are read by key and kept, so a repeated
                                 lookup is free.
    preload=False, memoize=False rows are read by key and not kept. For tables
                                 too large to hold, where a batch asks for its
                                 own keys and moves on.

Values are held as JSON text and parsed on read, so callers always get their own
object and can never mutate another record's copy -- or the cache itself.

A connection is opened per operation: the species resolvers write from a thread
pool, and SQLite connections cannot be shared across threads.
"""

import json
from itertools import batched

from config import logger

from .common import sqlite

# SQLite allows 999 bound variables by default; keep headroom.
_SQL_CHUNK_SIZE = 900


class _SqliteBacked:
    """Shared plumbing: the table, its DDL, key normalization and chunked reads."""

    def __init__(self, db_path, table, key_column, columns, preload=False, normalize=True):
        self.db_path = db_path
        self.table = table
        self.key_column = key_column
        self._ddl = f"CREATE TABLE IF NOT EXISTS {table} ({columns})"
        self.preload = preload
        self.normalize = normalize
        self._loaded = False

    def _connect(self):
        return sqlite(self.db_path, self._ddl)

    def _key(self, key):
        return str(key).lower().strip() if self.normalize else str(key)

    def _select_in(self, conn, columns, keys):
        """Yield rows for `keys`, a chunk of bound variables at a time."""
        for chunk in batched(sorted(keys), _SQL_CHUNK_SIZE):
            placeholders = ",".join("?" for _ in chunk)
            yield from conn.execute(
                f"SELECT {columns} FROM {self.table} WHERE {self.key_column} IN ({placeholders})",
                chunk,
            )

    def reset(self):
        """Forget everything held in memory, so the next read sees current data."""
        self._loaded = False


class SqliteCache(_SqliteBacked):
    """A SQLite table of key -> JSON value, fronted by an in-memory cache."""

    def __init__(
        self,
        db_path,
        table,
        key_column="original_name",
        value_column="standard_dict",
        preload=False,
        memoize=True,
        normalize=True,
    ):
        super().__init__(
            db_path,
            table,
            key_column,
            f"{key_column} TEXT PRIMARY KEY, {value_column} TEXT",
            preload=preload,
            normalize=normalize,
        )
        self.value_column = value_column
        self.memoize = memoize or preload
        self._memo = {}

    def reset(self):
        super().reset()
        self._memo.clear()

    def _load_all(self):
        if self._loaded:
            return
        self._loaded = True
        try:
            with self._connect() as conn:
                rows = conn.execute(f"SELECT {self.key_column}, {self.value_column} FROM {self.table}").fetchall()
            self._memo.update({self._key(key): value for key, value in rows if value})
            logger.info("Loaded %s rows from %s", len(self._memo), self.table)
        except Exception as e:
            logger.error("Error loading %s: %s", self.table, e)

    def get(self, key):
        """The cached value for `key`, or None."""
        if self.preload:
            self._load_all()
        key = self._key(key)
        if key in self._memo:
            return json.loads(self._memo[key])
        if self.preload:
            return None

        try:
            with self._connect() as conn:
                row = conn.execute(
                    f"SELECT {self.value_column} FROM {self.table} WHERE {self.key_column} = ?", (key,)
                ).fetchone()
        except Exception as e:
            logger.error("Error reading %s: %s", self.table, e)
            return None

        if not row or not row[0]:
            return None
        if self.memoize:
            self._memo[key] = row[0]
        return json.loads(row[0])

    def get_many(self, keys):
        """The cached values for `keys`, as {normalized key: value}, skipping misses."""
        keys = {self._key(key) for key in keys}
        if not keys:
            return {}
        if self.preload:
            self._load_all()

        found = {key: json.loads(self._memo[key]) for key in keys if key in self._memo}
        missing = keys - found.keys()
        if self.preload or not missing:
            return found

        try:
            with self._connect() as conn:
                for key, value in self._select_in(conn, f"{self.key_column}, {self.value_column}", missing):
                    if not value:
                        continue
                    if self.memoize:
                        self._memo[self._key(key)] = value
                    found[self._key(key)] = json.loads(value)
        except Exception as e:
            logger.error("Error reading %s: %s", self.table, e)
        return found

    def put(self, key, value):
        """Store `value` under `key`, in the table and in memory."""
        self.put_many({key: value})

    def put_many(self, values):
        """Store several values in one SQLite transaction."""
        rows = [(self._key(key), json.dumps(value)) for key, value in values.items()]
        if not rows:
            return
        if self.memoize:
            self._memo.update(rows)
        try:
            with self._connect() as conn:
                conn.executemany(
                    f"INSERT OR REPLACE INTO {self.table} ({self.key_column}, {self.value_column}) VALUES (?, ?)",
                    rows,
                )
        except Exception as e:
            logger.error("Error writing to %s: %s", self.table, e)


class SqliteKeySet(_SqliteBacked):
    """A SQLite table of keys, fronted by an in-memory set.

    Used as a negative cache: a key that is present resolved to nothing before,
    so there is no point asking again.
    """

    def __init__(self, db_path, table, key_column="original_name", preload=False, normalize=True):
        super().__init__(db_path, table, key_column, f"{key_column} TEXT PRIMARY KEY", preload=preload, normalize=normalize)
        self._keys = set()

    def reset(self):
        super().reset()
        self._keys.clear()

    def _load_all(self):
        if self._loaded:
            return
        self._loaded = True
        try:
            with self._connect() as conn:
                rows = conn.execute(f"SELECT {self.key_column} FROM {self.table}").fetchall()
            self._keys.update(self._key(row[0]) for row in rows if row and row[0])
            logger.info("Loaded %s keys from %s", len(self._keys), self.table)
        except Exception as e:
            logger.error("Error loading %s: %s", self.table, e)

    def known(self, keys):
        """The subset of `keys` already recorded, as normalized keys."""
        if self.preload:
            self._load_all()
        keys = {self._key(key) for key in keys}
        hits = keys & self._keys
        remaining = keys - hits
        if self.preload or not remaining:
            return hits

        try:
            with self._connect() as conn:
                for (key,) in self._select_in(conn, self.key_column, remaining):
                    key = self._key(key)
                    self._keys.add(key)
                    hits.add(key)
        except Exception as e:
            logger.error("Error reading %s: %s", self.table, e)
        return hits

    def __contains__(self, key):
        return bool(self.known([key]))

    def add(self, key):
        """Record `key` as resolving to nothing."""
        self.add_many([key])

    def add_many(self, keys):
        """Record several keys in one SQLite transaction."""
        keys = {self._key(key) for key in keys if key is not None}
        keys.discard("")
        if not keys:
            return
        self._keys.update(keys)
        try:
            with self._connect() as conn:
                conn.executemany(
                    f"INSERT OR IGNORE INTO {self.table} ({self.key_column}) VALUES (?)",
                    ((key,) for key in keys),
                )
        except Exception as e:
            logger.error("Error writing to %s: %s", self.table, e)
