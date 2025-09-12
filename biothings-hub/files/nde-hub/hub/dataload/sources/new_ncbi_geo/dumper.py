import os
from datetime import datetime, timezone

import biothings
import biothings.hub.dataload.dumper as dumper
import config
import requests
from Bio import Entrez
from config import GEO_API_KEY, GEO_EMAIL
from utils.utils import retry

biothings.config_for_app(config)


class New_NCBI_Geo_Dumper(dumper.BaseDumper):

    Entrez.email = GEO_EMAIL
    Entrez.api_key = GEO_API_KEY
    MAX_PARALLEL_DUMP = 5
    SRC_NAME = "new_ncbi_geo"
    SRC_ROOT_FOLDER = os.path.join(config.DATA_ARCHIVE_ROOT, SRC_NAME)

    def prepare_client(self):
        pass

    def release_client(self):
        pass

    def remote_is_better(self, remotefile, localfile):
        return True

    def set_release(self):
        self.release = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    def query_acc(self, term, retmax, retstart):
        handle = Entrez.esearch(db="gds", term=term, usehistory="y")
        record = Entrez.read(handle)
        handle.close()
        webenv = record["WebEnv"]
        query_key = record["QueryKey"]
        handle = Entrez.esummary(
            db="gds",
            query_key=query_key,
            webenv=webenv,
            retmode="xml",
            retmax=retmax,
            retstart=retstart,
        )

        records = Entrez.read(handle)
        handle.close()
        accs = [rec["Accession"] for rec in records]
        return accs

    def fetch_all_acc(self, term):
        """
        Fetch all Accession numbers from NCBI GEO using Biopython Entrez ESearch with pagination.
        In our case terms would be "GSE[ETYP]" or "GSM[ETYP]".
        """

        # retmax = 9999  # NCBI max per request
        retmax = 10  # for testing
        # First, get total count
        handle = Entrez.esearch(db="gds", term=term, usehistory="y")
        record = Entrez.read(handle)
        handle.close()
        total = int(record["Count"])
        self.logger.info(f"Total {term} records to download: {total}")

        for restart in range(0, total, retmax):
            accs = self.query_acc(term, retmax, restart)
            yield accs

    def create_todump_list(self, force=False, **kwargs):
        self.set_release()  # so we can generate new_data_folder
        for term in ["GSE[ETYP]", "GSM[ETYP]"]:
            for accs in self.fetch_all_acc(term):
                for acc in accs:
                    if term.startswith("GSE"):
                        new_localfile = os.path.join(self.new_data_folder, f"gse/{acc}.txt")
                    else:
                        new_localfile = os.path.join(self.new_data_folder, f"gsm/{acc}.txt")
                    remoteurl = (
                        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief"
                    )
                    self.to_dump.append({"remote": remoteurl, "local": new_localfile})
                break  # for testing, remove this to fetch all

    @retry(3, 5)
    def download(self, remoteurl, localfile):
        self.prepare_local_folders(localfile)

        self.logger.info(f"Downloading SOFT file from: {remoteurl}")
        data = requests.get(remoteurl)
        with open(localfile, "w") as f:
            f.write(data.text)
