import os
from datetime import date, datetime, timedelta, timezone
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
    PAGE_SIZE = 100
    # the search API returns HTTP 500 for any request where pageSize * page > 20000,
    # so no single query can be paged past this many results
    MAX_OFFSET = 20000
    # lower bound for release_date bisection; the earliest study is from 2001
    EPOCH = date(1990, 1, 1)

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

    def search_url(self, params, **extra):
        return f"{self.SEARCH_URL}?{urlencode({**params, **extra})}"

    @retry(3, 5)
    def count(self, params):
        """totalHits for a query, without paging through it. Counts are exact."""
        resp = requests.get(self.search_url(params, pageSize=1, page=1), timeout=300)
        resp.raise_for_status()
        return resp.json()["totalHits"]

    @retry(3, 5)
    def collections(self):
        """(name, hits) per collection. europepmc is citations, not datasets."""
        resp = requests.get(self.FACETS_URL, timeout=300)
        resp.raise_for_status()
        return [
            (child["value"], child["hits"])
            for child in resp.json()["children"]
            if child["value"] != "europepmc"
        ]

    def windows(self, collection, hits):
        """Yield (params, hits) for queries that each stay within MAX_OFFSET."""
        params = {"facet.collection": collection}
        if hits <= self.MAX_OFFSET:
            yield params, hits
            return
        yield from self.bisect(params, self.EPOCH, date.today() + timedelta(days=1))

    def bisect(self, params, start, end):
        """
        Split a query on release_date until each half fits under MAX_OFFSET.
        The bounds are inclusive and the halves are adjacent days, so every
        accession falls in exactly one window.
        """
        window = {**params, "release_date": f"[{start} TO {end}]"}
        hits = self.count(window)
        if hits <= self.MAX_OFFSET:
            if hits:
                yield window, hits
            return
        if start == end:
            self.logger.error("%s studies released on %s exceed the offset cap; truncating", hits, start)
            yield window, self.MAX_OFFSET
            return
        mid = start + (end - start) // 2
        yield from self.bisect(params, start, mid)
        yield from self.bisect(params, mid + timedelta(days=1), end)

    def create_todump_list(self, force=False, **kwargs):
        self.set_release()  # so we can generate new_data_folder
        for collection, hits in self.collections():
            windows = list(self.windows(collection, hits))
            covered = sum(window_hits for _, window_hits in windows)
            if covered != hits:
                self.logger.warning("%s: windows cover %s of %s studies", collection, covered, hits)
            n = 0
            for params, window_hits in windows:
                pages = -(-window_hits // self.PAGE_SIZE)
                for page in range(1, pages + 1):
                    n += 1
                    remoteurl = self.search_url(params, pageSize=self.PAGE_SIZE, page=page)
                    new_localfile = os.path.join(self.new_data_folder, f"biostudies_{collection}_{n}.txt")
                    self.to_dump.append({"remote": remoteurl, "local": new_localfile})
            self.logger.info("%s: %s studies, %s windows, %s pages to dump", collection, hits, len(windows), n)

    @retry(3, 5)
    def download(self, remoteurl, localfile):
        self.prepare_local_folders(localfile)
        self.logger.info("Downloading accno from %s to %s", remoteurl, localfile)
        resp = requests.get(remoteurl, timeout=300)
        resp.raise_for_status()
        with open(localfile, "w") as f:
            for hit in resp.json().get("hits") or []:
                accession = hit.get("accession")
                if accession:
                    f.write(accession + "\n")
