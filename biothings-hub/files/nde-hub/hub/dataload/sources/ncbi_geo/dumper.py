import os
import subprocess
from datetime import datetime, timezone

import biothings
import biothings.hub.dataload.dumper as dumper
import config
import requests
from Bio import Entrez
from config import GEO_API_KEY, GEO_EMAIL
from utils.utils import retry

biothings.config_for_app(config)


class NCBI_Geo_Dumper(dumper.BaseDumper):

    Entrez.email = GEO_EMAIL
    Entrez.api_key = GEO_API_KEY
    SRC_NAME = "ncbi_geo"
    SRC_ROOT_FOLDER = os.path.join(config.DATA_ARCHIVE_ROOT, SRC_NAME)
    MAX_PARALLEL_DUMP = 3

    def prepare_client(self):
        pass

    def release_client(self):
        pass

    def remote_is_better(self, remotefile, localfile):
        return True

    def set_release(self):
        self.release = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    def query_acc(self, term, retstart, retmax):
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

    # def fetch_all_acc(self, term):
    #     """
    #     Fetch all Accession numbers from NCBI GEO using Biopython Entrez ESearch with pagination.
    #     In our case terms would be "GSE[ETYP]" or "GSM[ETYP]".
    #     """

    #     # retmax = 9999  # NCBI max per request
    #     retmax = 10  # for testing
    #     # First, get total count
    #     handle = Entrez.esearch(db="gds", term=term, usehistory="y")
    #     record = Entrez.read(handle)
    #     handle.close()
    #     total = int(record["Count"])
    #     self.logger.info(f"Total {term} records to download: {total}")

    #     for retstart in range(0, total, retmax):
    #         accs = self.query_acc(term, retstart, retmax)
    #         yield accs

    # def create_todump_list(self, force=False, **kwargs):
    #     self.set_release()  # so we can generate new_data_folder
    #     for term in ["GSE[ETYP]", "GSM[ETYP]"]:
    #         if term.startswith("GSE"):
    #             new_localfile = os.path.join(self.new_data_folder, "gse/")
    #         else:
    #             new_localfile = os.path.join(self.new_data_folder, "gsm/")

    #         self.logger.info(f"Preparing to dump {term} files to {new_localfile}")
    #         self.to_dump.append({"remote": term, "local": new_localfile})  # placeholder


    @retry(3, 5)
    def wget_download(self, url, output_path):
        subprocess.run(["wget", url, "-O", output_path], check=True)

    # def download(self, remoteurl, localfile):
    #     self.prepare_local_folders(localfile)

    #     for accs in self.fetch_all_acc(remoteurl):
    #         for acc in accs:
    #             if remoteurl.startswith("GSE"):
    #                 new_localfile = os.path.join(localfile, f"{acc}.txt")
    #             else:
    #                 new_localfile = os.path.join(localfile, f"{acc}.txt")

    #             remoteurl = (
    #                 f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief"
    #             )
    #             self.logger.info(f"Downloading SOFT file from: {remoteurl}")
    #             self.wget_download(remoteurl, new_localfile)
    #         break  # for testing, remove this to fetch all


    def count_all_acc(self, term):
        """
        Fetch all Accession numbers from NCBI GEO using Biopython Entrez ESearch with pagination.
        In our case terms would be "GSE[ETYP]" or "GSM[ETYP]".
        """


        # get total count
        handle = Entrez.esearch(db="gds", term=term, usehistory="y")
        record = Entrez.read(handle)
        handle.close()
        total = int(record["Count"])


        return total


    def create_todump_list(self, force=False, **kwargs):
        self.set_release()  # so we can generate new_data_folder
        for term in ["GSE[ETYP]", "GSM[ETYP]"]:
            if term.startswith("GSE"):
                new_localfile = os.path.join(self.new_data_folder, "gse/")
            else:
                new_localfile = os.path.join(self.new_data_folder, "gsm/")

            total_count = self.count_all_acc(term)
            self.logger.info(f"Total {term} records to download: {total_count}")
            retmax = 9999
            count = 0   # for testing
            for retstart in range(0, total_count, retmax):
                self.logger.info(f"Preparing to dump {term} files to {new_localfile}, start at {retstart}")
                remote = [term, retstart, retmax]
                self.to_dump.append({"remote": remote, "local": new_localfile})


    def create_subdir(self, localfile, acc):
        # Extract prefix (GSE/GSM) and numeric part
        prefix = acc[:3]
        num = acc[3:]
        # Pad numeric part to at least 3 digits for nnn, or 6 for full
        padded = num.zfill(6)
        subdir = prefix + padded[:3] + "nnn"
        subdir_path = os.path.join(localfile, subdir)
        if not os.path.exists(subdir_path):
            os.makedirs(subdir_path)
        return subdir_path



    def download(self, remoteurl, localfile):
        self.prepare_local_folders(localfile)

        accs = self.query_acc(remoteurl[0], remoteurl[1], remoteurl[2])
        for acc in accs:
            if remoteurl[0].startswith("GSE"):
                subdir = self.create_subdir(localfile, acc)
                new_localfile = os.path.join(subdir, f"{acc}.txt")
            else:
                subdir = self.create_subdir(localfile, acc)
                new_localfile = os.path.join(subdir, f"{acc}.txt")

            remoteurl = (
                f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief"
            )
            self.logger.info(f"Downloading SOFT file from: {remoteurl}")
            try:
                self.wget_download(remoteurl, new_localfile)
            except Exception as e:
                print(f"Error downloading {remoteurl}: {e}")
