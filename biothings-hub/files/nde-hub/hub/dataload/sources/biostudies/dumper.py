import os
from datetime import datetime, timezone

import biothings
import biothings.hub.dataload.dumper as dumper
import config
import requests
from utils.utils import retry

biothings.config_for_app(config)


class Biostudies_Dumper(dumper.BaseDumper):

    MAX_PARALLEL_DUMP = 3
    SRC_NAME = "biostudies"
    SRC_ROOT_FOLDER = os.path.join(config.DATA_ARCHIVE_ROOT, SRC_NAME)

    def prepare_client(self):
        pass

    def release_client(self):
        pass

    def remote_is_better(self, remotefile, localfile):
        return True

    def set_release(self):
        self.release = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    def create_todump_list(self, force=False, **kwargs):
        self.set_release()  # so we can generate new_data_folder
        # Initial request to get the facets
        request = requests.get("https://www.ebi.ac.uk/biostudies/api/v1/public/facets/facet.collection/").json()

        # Extract facets, excluding "europepmc"
        facets = [child.get("value") for child in request.get("children") if child.get("value") != "europepmc"]
        # Base URL
        base_url = "https://www.ebi.ac.uk/biostudies/api/v1/search?"
        # Construct the URL with all facets combined
        facet_string = "&".join([f"facet.collection={facet}" for facet in facets])
        url = f"{base_url}{facet_string}"
        # Make the request
        hits = requests.get(url).json().get("totalHits")
        # Calculate the number of pages
        pages = (hits // 100) + 1 if hits % 100 != 0 else hits // 100
        # form urls to dump TODO change to pages + 1 after testing
        for page in range(1, pages + 1):
            remoteurl = f"{url}&pageSize=100&page={page}"
            new_localfile = os.path.join(self.new_data_folder, f"biostudies_{page}.txt")
            self.to_dump.append({"remote": remoteurl, "local": new_localfile})

    @retry(3, 5)
    def download(self, remoteurl, localfile):
        self.prepare_local_folders(localfile)
        self.logger.info("Downloading accno from %s to %s", remoteurl, localfile)
        data = requests.get(remoteurl).json()
        with open(localfile, "w") as f:
            for hit in data.get("hits"):
                accession = hit.get("accession")
                if accession:
                    f.write(accession + "\n")
