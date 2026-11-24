"""Copy the per-module SQLite caches into cache.db, the `Cache` layout.

Run once, with the hub stopped, as the hub user:

    python -m utils.migrate_caches [--dry-run] [--force] [--from DIR] [--to PATH]

The old files are opened read-only and left as they are. cache.db is built
under a temporary name and renamed when complete. An existing cache.db stops
the run unless --force is given.

Row rules:
  * A key keeps the position of its first row and the value of its last
    answer, so an answer replaces a miss and a miss never replaces an answer.
  * Rows with no known write time get one spread over the 90 days before the
    run, fixed per key.
"""

import argparse
import json
import os
import sqlite3
import sys
import time
import zlib
from collections import Counter
from datetime import date

from .cache import DB_PATH

STANDARDIZERS = "/data/nde-hub/standardizers"
OLD_FILES = {
    "pmid_lookup": f"{STANDARDIZERS}/pmid_lookup/pmid_lookup.db",
    "pubtator_lookup": f"{STANDARDIZERS}/pubtator_lookup/pubtator_lookup.db",
    "extract_lookup": f"{STANDARDIZERS}/extract_lookup/extract_lookup.db",
    "funding_lookup": f"{STANDARDIZERS}/funding_lookup/funding_lookup.db",
    "lineage_lookup": f"{STANDARDIZERS}/lineage_lookup/lineage_lookup.db",
    "nctid_lookup": "/nvme/nde-hub/standardizers/nctid_lookup/nctid_lookup.db",
    "mdc": f"{STANDARDIZERS}/make_data_count/mdc.db",
}
# Old tables that are not caches: the bulk dumps now uploaded to Mongo.
NOT_CACHES = {"pmid_lookup": ("species_data", "disease_data"), "mdc": ("releases", "links_")}

_SPREAD_SECONDS = 90 * 86400
_EXTRACT_FAILURE_TTL_SECONDS = 3600
_EXTRACT_TABLES = {"-2": "extract_species", "-26": "extract_disease"}
_BATCH = 50_000


def _midnight(iso_date):
    return time.mktime(date.fromisoformat(iso_date).timetuple())


# A rule turns one old row into (new table, key, JSON value or None for a miss,
# updated or None for a spread time), or None to skip the row.
def _answer(table, encode=None):
    """An answer. Stored JSON is kept as is; `encode` makes JSON of a plain value."""

    def rule(key, value):
        if encode:
            return table, key, encode(value), None
        return (table, key, value, None) if value else None

    return rule


def _miss(table):
    return lambda key, *_: (table, key, None, None)


def _doi_pmid(key, pmid):
    return "doi_pmid", key, None if pmid is None else json.dumps(pmid), None


def _extract_failure(key, retry_after, entity_type):
    return _EXTRACT_TABLES[entity_type], key, None, retry_after - _EXTRACT_FAILURE_TTL_SECONDS


def _trial_design(key, entry):
    entry = json.loads(entry)
    return "trial_designs", key, json.dumps(entry["design"]), _midnight(entry["fetched"])


def _dated_miss(table):
    return lambda key, checked: (table, key, None, _midnight(json.loads(checked)["checked"]))


# (old file, old table, selected columns, rule). Misses come before the answers they merge with.
MAPPING = [
    ("pmid_lookup", "eutils_cache", "pmid, data", _answer("pubmed_efetch")),
    ("pmid_lookup", "doi_cache", "doi, pmid", _doi_pmid),
    ("pubtator_lookup", "species", "original_name, standard_dict", _answer("pubtator_species_terms")),
    ("pubtator_lookup", "health_conditions", "original_name, standard_dict", _answer("health_condition_terms")),
    ("pubtator_lookup", "health_conditions_negative", "original_name", _miss("health_condition_rejects")),
    ("extract_lookup", "extract_failures", "ndeid, retry_after, entity_type", _extract_failure),
    ("extract_lookup", "species", "ndeid, text_response", _answer("extract_species", json.dumps)),
    ("extract_lookup", "disease", "ndeid, text_response", _answer("extract_disease", json.dumps)),
    ("extract_lookup", "species_negative", "original_name", _miss("species_details")),
    ("extract_lookup", "species_details", "original_name, standard_dict", _answer("species_details")),
    ("extract_lookup", "extract_taxon_failures", "identifier", _miss("taxon_details")),
    ("extract_lookup", "extract_taxon_details", "identifier, standard_dict", _answer("taxon_details")),
    ("extract_lookup", "description_species_negative", "lookup_key", _miss("description_species_rejects")),
    ("funding_lookup", "funding_lookup", "funding_id, funding", _answer("nih_grants")),
    ("funding_lookup", "funder_cache", "funder_name, funder_data", _answer("crossref_funders")),
    ("lineage_lookup", "taxon_lineage", "taxid, lineage", _answer("taxon_lineage")),
    ("lineage_lookup", "taxon_parent", "taxid, parent_taxid", _answer("taxon_parent", json.dumps)),
    ("nctid_lookup", "trial_designs", "nctid, design", _trial_design),
    ("nctid_lookup", "ncit_names", "iri, name", _answer("ncit_names")),
    ("mdc", "pubmed_misses", "doi, checked", _dated_miss("doi_pubmed_ids")),
    ("mdc", "pubmed_ids", "doi, ids", _answer("doi_pubmed_ids")),
    ("mdc", "pubmed_descriptions", "pmid, description", _answer("pubmed_summaries")),
    ("mdc", "crossref_misses", "doi, checked", _dated_miss("crossref_works")),
    ("mdc", "crossref_descriptions", "doi, description", _answer("crossref_works")),
]

_UPSERT = (
    "INSERT INTO {table} (key, value, updated) VALUES (?, ?, ?) "
    "ON CONFLICT(key) DO UPDATE SET value = excluded.value, updated = excluded.updated "
    "WHERE excluded.value IS NOT NULL"
)


def _tables(conn):
    return {name for (name,) in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}


def _uncovered(old_files):
    covered = {(source, table) for source, table, _, _ in MAPPING}
    return sorted(
        f"{source}.{table}"
        for source, conn in old_files.items()
        for table in _tables(conn)
        if (source, table) not in covered and not table.startswith(("sqlite_", *NOT_CACHES.get(source, ())))
    )


class _Writer:
    """Upserts converted rows into the new file, a batch at a time."""

    def __init__(self, path):
        self.conn = sqlite3.connect(path, isolation_level=None)
        self.conn.execute("PRAGMA journal_mode = OFF")
        self.conn.execute("PRAGMA synchronous = OFF")
        self.conn.execute("PRAGMA cache_size = -1000000")
        self.conn.execute(
            "CREATE TABLE migrated_from (old TEXT, new TEXT, rows INTEGER, misses INTEGER, skipped INTEGER)"
        )
        self.created = set()
        self.batch = []

    def add(self, table, key, value, updated):
        self.batch.append((table, key, value, updated))
        if len(self.batch) >= _BATCH:
            self.flush()

    def flush(self):
        for table in sorted({row[0] for row in self.batch}):
            if table not in self.created:
                self.conn.execute(f"CREATE TABLE {table} (key TEXT PRIMARY KEY, value TEXT, updated REAL NOT NULL)")
                self.created.add(table)
            self.conn.executemany(_UPSERT.format(table=table), [row[1:] for row in self.batch if row[0] == table])
        self.batch.clear()


def migrate(old_paths, target, dry_run=False):
    started = time.time()
    old_files = {}
    for source, path in old_paths.items():
        if os.path.exists(path):
            old_files[source] = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        else:
            print(f"skipped {source}: {path} not found")
    if uncovered := _uncovered(old_files):
        sys.exit(f"Tables with no migration rule: {', '.join(uncovered)}. Nothing was written.")

    building = target + ".migrating"
    writer = None
    if not dry_run:
        for suffix in ("", "-journal", "-wal", "-shm"):
            if os.path.exists(building + suffix):
                os.remove(building + suffix)
        writer = _Writer(building)

    converted = Counter()
    for source, table, columns, rule in MAPPING:
        if source not in old_files or table not in _tables(old_files[source]):
            continue
        counts = Counter()
        new_tables = set()
        if writer:
            writer.conn.execute("BEGIN")
        for row in old_files[source].execute(f"SELECT {columns} FROM {table} ORDER BY rowid"):
            counts["rows"] += 1
            out = rule(str(row[0]), *row[1:]) if row[0] is not None else None
            if out is None:
                counts["skipped"] += 1
                continue
            new_table, key, value, updated = out
            counts["misses"] += value is None
            converted[new_table] += 1
            new_tables.add(new_table)
            if writer:
                writer.add(
                    new_table,
                    key,
                    value,
                    started - zlib.crc32(key.encode()) % _SPREAD_SECONDS if updated is None else updated,
                )
        if writer:
            writer.flush()
            writer.conn.execute(
                "INSERT INTO migrated_from VALUES (?, ?, ?, ?, ?)",
                (
                    f"{source}.{table}",
                    ",".join(sorted(new_tables)),
                    counts["rows"],
                    counts["misses"],
                    counts["skipped"],
                ),
            )
            writer.conn.execute("COMMIT")
        print(
            f"{source}.{table} -> {','.join(sorted(new_tables)) or '-'}: rows={counts['rows']:,} "
            f"misses={counts['misses']:,} skipped={counts['skipped']:,}",
            flush=True,
        )
    for conn in old_files.values():
        conn.close()

    if not writer:
        print(f"dry run, nothing written: {time.time() - started:.0f}s")
        return
    for table in sorted(writer.created):
        hits, misses = writer.conn.execute(f"SELECT COUNT(value), COUNT(*) - COUNT(value) FROM {table}").fetchone()
        merged = converted[table] - hits - misses
        print(f"  {table}: {hits:,} answers, {misses:,} misses, {merged:,} rows merged into an earlier key")
    writer.conn.execute("PRAGMA journal_mode = WAL")
    writer.conn.close()
    for suffix in ("-wal", "-shm"):
        if os.path.exists(target + suffix):
            os.remove(target + suffix)
    os.replace(building, target)
    print(f"wrote {target}: {os.path.getsize(target) / 2**30:.2f} GB in {time.time() - started:.0f}s")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dry-run", action="store_true", help="read and convert everything, write nothing")
    parser.add_argument("--force", action="store_true", help="replace an existing cache.db")
    parser.add_argument("--from", dest="from_dir", help="read the old files from this directory instead")
    parser.add_argument("--to", default=DB_PATH, help=f"the file to write (default {DB_PATH})")
    args = parser.parse_args(argv)

    old_paths = OLD_FILES
    if args.from_dir:
        old_paths = {source: os.path.join(args.from_dir, os.path.basename(path)) for source, path in OLD_FILES.items()}
    if os.path.exists(args.to) and not args.dry_run and not args.force:
        sys.exit(f"{args.to} exists: it was migrated already or the new code has run. --force replaces it.")
    migrate(old_paths, args.to, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
