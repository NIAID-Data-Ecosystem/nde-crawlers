"""Make Data Count data-citation releases for the make_data_count pipeline stage.

Dumper only: `post_dump` loads each release into the stage's SQLite database.
There is no uploader, and the source is never part of a build. Releases are
BigQuery tables, read as in MDC's own downloader:
https://github.com/digital-science/ds-make-data-count-pipeline-external/blob/main/examples/download_releases.py
"""

import csv
import gzip
import json
import os

import biothings
import config

biothings.config_for_app(config)
from biothings.hub.dataload.dumper import BaseDumper
from config import DATA_ARCHIVE_ROOT, MDC_GOOGLE_CREDENTIALS
from utils.make_data_count import MIN_RELEASE, current_release, load_release

DATA_PROJECT = "ds-make-data-count"
DATA_DATASET = "public"
CHANGELOG_TABLE = "release_changelog"
PAGE_SIZE = 50_000
REQUEST_TIMEOUT = 120
RELEASE_FILE = "release.json"


def _bigquery_client():
    from google.cloud import bigquery
    from google.oauth2 import service_account

    credentials = service_account.Credentials.from_service_account_file(
        MDC_GOOGLE_CREDENTIALS, scopes=["https://www.googleapis.com/auth/bigquery.readonly"]
    )
    # The project only labels the client; reading table rows runs no query jobs.
    return bigquery.Client(project=DATA_PROJECT, credentials=credentials)


class MakeDataCountDumper(BaseDumper):
    SRC_NAME = "make_data_count"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SCHEDULE = "0 12 3 * *"  # 3rd of each month at 12:00; MDC publishes on the 2nd at 06:00 UTC
    ARCHIVE = False  # releases are kept in the stage's database
    AUTO_UPLOAD = False  # no uploader; post_dump loads the release
    DISABLED = not MDC_GOOGLE_CREDENTIALS

    @property
    def client(self):
        # BigQuery clients are made where they are used; download runs in another process.
        return None

    def prepare_client(self):
        pass

    def release_client(self):
        pass

    def need_prepare(self):
        return False

    def create_todump_list(self, force=False, **kwargs):
        """Queue the newest release in the changelog when it is newer than the loaded one."""
        client = _bigquery_client()
        changelog = client.get_table(f"{DATA_PROJECT}.{DATA_DATASET}.{CHANGELOG_TABLE}")
        snapshots = {}
        for row in client.list_rows(changelog, timeout=REQUEST_TIMEOUT):
            row = dict(row.items())
            if row.get("release") and row.get("snapshot_table") and str(row["release"]) >= MIN_RELEASE:
                # snapshot_table is recorded as "public.data_citations_YYYYMM"
                snapshots[str(row["release"])] = row["snapshot_table"].split(".")[-1]

        loaded = current_release()
        if not snapshots or (not force and loaded and loaded["release"] >= max(snapshots)):
            self.to_dump = []
            return
        self.release = max(snapshots)
        snapshot = snapshots[self.release]
        self.to_dump = [{"remote": snapshot, "local": os.path.join(self.new_data_folder, f"{snapshot}.csv.gz")}]

    def download(self, remotefile, localfile):
        """Write snapshot table `remotefile` to `localfile` as gzipped CSV, then record its release."""
        self.prepare_local_folders(localfile)
        client = _bigquery_client()
        table = client.get_table(f"{DATA_PROJECT}.{DATA_DATASET}.{remotefile}")
        partial = f"{localfile}.part"
        rows = 0
        try:
            with gzip.open(partial, "wt", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow([field.name for field in table.schema])
                for row in client.list_rows(table, page_size=PAGE_SIZE, timeout=REQUEST_TIMEOUT):
                    writer.writerow(["" if value is None else value for value in row.values()])
                    rows += 1
            if rows != table.num_rows:
                raise RuntimeError(f"{remotefile}: downloaded {rows:,} rows, expected {table.num_rows:,}")
            os.replace(partial, localfile)
        finally:
            if os.path.exists(partial):
                os.remove(partial)

        release = {
            "release": remotefile.rsplit("_", 1)[-1],
            "table": remotefile,
            "created": table.created.date().isoformat() if table.created else None,
            "row_count": rows,
        }
        with open(os.path.join(os.path.dirname(localfile), RELEASE_FILE), "w") as handle:
            json.dump(release, handle)

    def post_dump(self, *args, **kwargs):
        folder = self.new_data_folder
        with open(os.path.join(folder, RELEASE_FILE)) as handle:
            release = json.load(handle)
        path = os.path.join(folder, f"{release['table']}.csv.gz")
        with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
            load_release(csv.DictReader(handle), release["release"], release["created"], release["row_count"])
