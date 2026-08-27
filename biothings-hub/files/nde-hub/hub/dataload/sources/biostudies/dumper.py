import os
from datetime import datetime, timezone
from urllib.parse import urlencode

import biothings
import biothings.hub.dataload.dumper as dumper
import config
import requests
from utils import retry

biothings.config_for_app(config)


class Biostudies_Dumper(dumper.BaseDumper):

    MAX_PARALLEL_DUMP = 3
    SRC_NAME = "biostudies"
    SRC_ROOT_FOLDER = os.path.join(config.DATA_ARCHIVE_ROOT, SRC_NAME)

    FACETS_URL = "https://www.ebi.ac.uk/biostudies/api/v1/public/facets/facet.collection/"
    SEARCH_URL = "https://www.ebi.ac.uk/biostudies/api/v1/search"
    # hits per search request, and accessions per output file. The uploader makes one
    # job per file and the parser issues one study request per accession under a 1800s
    # timeout; 1000 studies takes roughly 750s.
    PAGE_SIZE = 1000

    def __getstate__(self):
        # do_dump pickles self to send download() to a worker process. _state holds a
        # pymongo client, which is unpicklable, and BaseDumper's lazy properties
        # re-attach it on any read after unprepare() -- including do_dump's own error
        # handler, which then breaks every job submitted after the first failure.
        # The worker rebuilds what it needs on first use.
        state = self.__dict__.copy()
        state["_state"] = dict.fromkeys(self._state)
        return state

    def prepare_client(self):
        pass

    def release_client(self):
        pass

    def remote_is_better(self, remotefile, localfile):
        return True

    def set_release(self):
        self.release = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    @retry(3, 5)
    def collections(self):
        """Collections available in the BioStudies search index exclude europepmc."""
        resp = requests.get(self.FACETS_URL, timeout=300)
        resp.raise_for_status()
        return [child["value"] for child in resp.json()["children"] if child["value"] != "europepmc"]

    @retry(3, 5)
    def total_hits(self, collection):
        """Studies the search index will return for a collection."""
        resp = requests.get(
            self.SEARCH_URL,
            params={"facet.collection": collection, "pageSize": 1, "pagination": "cursor"},
            timeout=300,
        )
        resp.raise_for_status()
        return resp.json().get("totalHits")

    def search_url(self, collection):
        params = {
            "facet.collection": collection,
            "pageSize": self.PAGE_SIZE,
            "pagination": "cursor",
        }
        return f"{self.SEARCH_URL}?{urlencode(params)}"

    @retry(3, 5)
    def fetch(self, url, cursor=None):
        resp = requests.get(url, params={"cursor": cursor} if cursor else None, timeout=300)
        resp.raise_for_status()
        return resp.json()

    def pages(self, url):
        """
        Walk cursor pagination, yielding each page's accessions.
        """
        cursor = None
        seen_cursors = set()
        while True:
            data = self.fetch(url, cursor)
            hits = data.get("hits") or []
            accessions = [accession for hit in hits if (accession := hit.get("accession"))]
            if accessions:
                yield accessions
            cursor = data.get("nextCursor")
            if not cursor:
                return
            if not hits:
                self.logger.warning("empty page with cursor still set, stopping walk of %s", url)
                return
            if cursor in seen_cursors:
                self.logger.warning("cursor %s repeated, stopping walk of %s", cursor, url)
                return
            seen_cursors.add(cursor)

    def create_todump_list(self, force=False, **kwargs):
        self.set_release()  # so we can generate new_data_folder
        # cursor pagination is sequential and cursors expire, so only the cursor-less
        # first URL can be queued here and a collection cannot be split across jobs;
        # one job per collection still runs MAX_PARALLEL_DUMP collections at once
        for collection in self.collections():
            remoteurl = self.search_url(collection)
            new_localfile = os.path.join(self.new_data_folder, f"biostudies_{collection}.txt")
            self.to_dump.append({"remote": remoteurl, "local": new_localfile})
            self.logger.info("%s: %s studies to dump", collection, self.total_hits(collection))

    def download(self, remoteurl, localfile):
        self.prepare_local_folders(localfile)
        base = os.path.splitext(localfile)[0]
        self.logger.info("Downloading accnos from %s to %s_*.txt", remoteurl, base)
        total = 0
        page = 0
        for page, accessions in enumerate(self.pages(remoteurl), start=1):
            total += len(accessions)
            with open(f"{base}_{page}.txt", "w") as f:
                f.write("".join(accession + "\n" for accession in accessions))
        self.logger.info("Wrote %s accnos to %s files from %s", total, page, remoteurl)
