import os
import tempfile
from datetime import datetime, timezone

import biothings
import biothings.hub.dataload.dumper as dumper
import config
import requests
from Bio import Entrez
from config import GEO_API_KEY, GEO_EMAIL
from utils import retry

biothings.config_for_app(config)


class NCBI_Geo_Dumper(dumper.BaseDumper):

    Entrez.email = GEO_EMAIL
    Entrez.api_key = GEO_API_KEY
    SRC_NAME = "ncbi_geo"
    SRC_ROOT_FOLDER = os.path.join(config.DATA_ARCHIVE_ROOT, SRC_NAME)
    SERIES_TERM = "GSE[ETYP]"
    # gds UIDs encode the accession number: GSE1234 -> 200001234
    GSE_UID_OFFSET = 200000000
    SCHEDULE = "0 17 1 1,4,7,10 *"  # 1st day of Jan, Apr, Jul, Oct at 17:00 UTC
    MAX_PARALLEL_DUMP = 5

    def prepare_client(self):
        # one session per dump process, so acc.cgi connections get reused
        self.client = requests.Session()

    def release_client(self):
        if self._state["client"]:
            self.client.close()
            self.client = None

    def remote_is_better(self, remotefile, localfile):
        return True

    def set_release(self):
        self.release = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    @retry(5, 5)
    def query_acc(self, retstart, retmax):
        """
        Fetch one page of series Accession numbers.
        The accession is derived from the gds UID, so no ESummary call is needed.
        """
        handle = Entrez.esearch(db="gds", term=self.SERIES_TERM, retstart=retstart, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        return [f"GSE{int(uid) - self.GSE_UID_OFFSET}" for uid in record["IdList"]]

    @retry(5, 5)
    def fetch_family(self, acc):
        """
        Fetch the full SOFT family for a series: the SERIES record, its PLATFORM, and every SAMPLE.
        view=brief omits the data tables.
        """
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=all&form=text&view=brief"
        resp = self.client.get(url, timeout=300)
        resp.raise_for_status()
        return resp.text

    @retry(5, 5)
    def count_all_acc(self):
        """Total number of "GSE[ETYP]" records in NCBI GEO."""
        handle = Entrez.esearch(db="gds", term=self.SERIES_TERM, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])

    def create_todump_list(self, force=False, **kwargs):
        self.set_release()  # so we can generate new_data_folder
        # only series are enumerated; their samples come from the family files
        total_count = self.count_all_acc()
        self.logger.info(f"Total {self.SERIES_TERM} records to download: {total_count}")
        retmax = 9999
        for retstart in range(0, total_count, retmax):
            self.logger.info(f"Preparing to dump families to {self.new_data_folder}, start at {retstart}")
            self.to_dump.append({"remote": [retstart, retmax], "local": self.new_data_folder})

    def create_subdir(self, localfile, acc):
        # Extract prefix (GSE/GSM) and numeric part
        prefix = acc[:3]
        num = acc[3:]
        # Pad numeric part to at least 3 digits for nnn, or 6 for full
        padded = num.zfill(6)
        subdir = prefix + padded[:3] + "nnn"
        subdir_path = os.path.join(localfile, subdir)
        os.makedirs(subdir_path, exist_ok=True)
        return subdir_path

    def split_family(self, text):
        """
        Split a SOFT family stream into its records, yielding (accession, record_text).
        Each record runs from its "^TYPE = ACC" header to the next one.
        """
        acc, lines = None, []
        for line in text.splitlines():
            if line.startswith("^"):
                if acc:
                    yield acc, "\n".join(lines) + "\n"
                _, _, value = line.partition("=")
                acc = value.strip()
                lines = [line]
            elif acc:
                lines.append(line)
        if acc:
            yield acc, "\n".join(lines) + "\n"

    def write_record(self, path, text):
        """
        Write via a uniquely named temp file in the destination dir, then rename.
        A sample shared by several series is written by several dump processes at once;
        os.replace is atomic, so readers only ever see a complete file.
        The temp file is not named .txt, so one orphaned by a killed process is not parsed.
        """
        fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), prefix=".tmp-", suffix=".part")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                f.write(text)
            os.replace(tmp, path)
        except BaseException:
            if os.path.exists(tmp):
                os.unlink(tmp)
            raise

    def download(self, remoteurl, localfile):
        retstart, retmax = remoteurl
        dirs = {"GSE": os.path.join(localfile, "gse"), "GSM": os.path.join(localfile, "gsm")}

        for acc in self.query_acc(retstart, retmax):
            self.logger.info(f"Downloading SOFT family for: {acc}")
            try:
                family = self.fetch_family(acc)
            except Exception as e:
                self.logger.error(f"Error downloading family for {acc}: {e}")
                continue
            for rec_acc, rec_text in self.split_family(family):
                # skip the PLATFORM and DATABASE records in the family
                base = dirs.get(rec_acc[:3])
                if not base:
                    continue
                subdir = self.create_subdir(base, rec_acc)
                try:
                    self.write_record(os.path.join(subdir, f"{rec_acc}.txt"), rec_text)
                except Exception as e:
                    self.logger.error(f"Error writing {rec_acc} from {acc}: {e}")
