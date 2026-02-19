import datetime
import json
import logging
import os
import time
from html.parser import HTMLParser

import psutil
import requests
from api_secret import TOKEN as api_key
from sql_database import NDEDatabase

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")

MAX_RETRIES = 5
BASE_DELAY = 2  # 2 seconds
MAX_DELAY = 60  # 1 minute
REQUEST_DELAY = float(os.getenv("DATAVERSE_REQUEST_DELAY", "1.0"))  # seconds between API calls
COOLDOWN_THRESHOLD = 10  # consecutive export failures before triggering cooldown
COOLDOWN_SECONDS = float(os.getenv("DATAVERSE_COOLDOWN_SECONDS", "60"))  # cooldown sleep duration


class Dataverse(NDEDatabase):
    SQL_DB = "dataverse.db"
    EXPIRE = datetime.timedelta(days=90)

    DATAVERSE_SERVER = "https://dataverse.harvard.edu/api"
    EXPORT_URL = f"{DATAVERSE_SERVER}/datasets/export?exporter=schema.org"
    DATAVERSE_JSON_EXPORT_URL = f"{DATAVERSE_SERVER}/datasets/export?exporter=dataverse_json"

    def __init__(self, sql_db=None):
        super().__init__(sql_db)
        self._session = requests.Session()
        self._headers = {"User-Agent": "nde-crawlers/dataverse"}
        self._headers["X-Dataverse-key"] = api_key

    def _sleep_for_rate_limit(self, response, backoff_time=None):
        """Sleep when rate-limited. Handles both 429 and 403 responses."""
        status = getattr(response, "status_code", None) if response else None
        if status not in (429, 403):
            return
        retry_after = response.headers.get("Retry-After") if response else None
        try:
            seconds = int(retry_after) if retry_after else (backoff_time or MAX_DELAY)
        except ValueError:
            seconds = backoff_time or MAX_DELAY
        seconds = max(1, min(MAX_DELAY, seconds))
        logger.warning(f"Rate limited ({status}). Sleeping {seconds}s before retrying...")
        time.sleep(seconds)

    @staticmethod
    def _first_str(value):
        if isinstance(value, str):
            return value
        if isinstance(value, list):
            for item in value:
                if isinstance(item, str):
                    return item
                if isinstance(item, dict):
                    for key in ("@id", "identifier", "url", "value"):
                        if isinstance(item.get(key), str):
                            return item.get(key)
        if isinstance(value, dict):
            for key in ("@id", "identifier", "url", "value"):
                if isinstance(value.get(key), str):
                    return value.get(key)
        return None

    @staticmethod
    def _select_dataset_from_jsonld(parsed):
        """Pick the most likely Dataset object from parsed JSON-LD.

        Dataverse landing pages may expose JSON-LD as:
        - a single dict
        - a list of dicts
        - a dict with an '@graph' list
        """

        def iter_nodes(obj):
            if isinstance(obj, dict):
                yield obj
                graph = obj.get("@graph")
                if isinstance(graph, list):
                    for node in graph:
                        if isinstance(node, dict):
                            yield node
            elif isinstance(obj, list):
                for item in obj:
                    yield from iter_nodes(item)

        def is_dataset(node: dict) -> bool:
            t = node.get("@type")
            if isinstance(t, str):
                return t.lower() == "dataset"
            if isinstance(t, list):
                return any(isinstance(x, str) and x.lower() == "dataset" for x in t)
            # Some pages only have schema context and an identifier/url; treat as a candidate.
            return bool(node.get("identifier") or node.get("@id") or node.get("url"))

        candidates = [node for node in iter_nodes(parsed) if isinstance(node, dict) and is_dataset(node)]
        if not candidates:
            return parsed if isinstance(parsed, dict) else None

        # Prefer a true Dataset type if present.
        for node in candidates:
            t = node.get("@type")
            if (isinstance(t, str) and t.lower() == "dataset") or (
                isinstance(t, list) and any(isinstance(x, str) and x.lower() == "dataset" for x in t)
            ):
                return node
        return candidates[0]

    @staticmethod
    def _funder_has_identifier(funder) -> bool:
        if not isinstance(funder, dict):
            return False
        for key in ("identifier", "@id", "sameAs"):
            value = funder.get(key)
            if isinstance(value, str) and value.strip():
                return True
            if isinstance(value, list) and any(isinstance(v, str) and v.strip() for v in value):
                return True
            if isinstance(value, dict) and any(
                isinstance(v, str) and v.strip() for v in value.values() if v is not None
            ):
                return True
        return False

    @classmethod
    def _has_funder_identifiers(cls, document: dict) -> bool:
        if not isinstance(document, dict):
            return False
        funders = []
        funding = document.get("funding")
        if isinstance(funding, dict) and isinstance(funding.get("funder"), list):
            funders.extend(funding.get("funder"))
        if isinstance(document.get("funder"), list):
            funders.extend(document.get("funder"))
        if isinstance(document.get("funder"), dict):
            funders.append(document.get("funder"))
        return any(cls._funder_has_identifier(f) for f in funders)

    @classmethod
    def _extract_funders(cls, document: dict) -> list[dict]:
        """Extract funder objects from a JSON-LD document.

        Returns a list of dict funders (normalizing strings to {'name': str}).
        """
        if not isinstance(document, dict):
            return []

        funders = []
        funding = document.get("funding")
        if isinstance(funding, dict) and isinstance(funding.get("funder"), list):
            funders.extend(funding.get("funder"))

        direct = document.get("funder")
        if isinstance(direct, list):
            funders.extend(direct)
        elif isinstance(direct, dict) or isinstance(direct, str):
            funders.append(direct)

        normalized: list[dict] = []
        for funder in funders:
            if isinstance(funder, dict):
                normalized.append(funder)
            elif isinstance(funder, str) and funder.strip():
                normalized.append({"name": funder.strip()})
        return normalized

    @staticmethod
    def _norm_name(value) -> str:
        if not value:
            return ""
        if isinstance(value, str):
            return " ".join(value.split()).casefold()
        return ""

    @classmethod
    def _merge_funder_identifiers(cls, schema_doc: dict, jsonld_doc: dict) -> bool:
        """Merge funder identifiers from landing-page JSON-LD into schema.org export.

        Returns True if schema_doc was modified.
        """
        if not isinstance(schema_doc, dict) or not isinstance(jsonld_doc, dict):
            return False

        jsonld_funders = [f for f in cls._extract_funders(jsonld_doc) if cls._funder_has_identifier(f)]
        if not jsonld_funders:
            return False

        schema_funders_raw = schema_doc.get("funder")
        if not isinstance(schema_funders_raw, list):
            schema_funders: list = []
        else:
            schema_funders = schema_funders_raw

        changed = False
        schema_by_name = {}
        for funder in schema_funders:
            if isinstance(funder, dict):
                name = cls._norm_name(funder.get("name"))
                if name:
                    schema_by_name[name] = funder

        for jsonld_funder in jsonld_funders:
            jsonld_name = cls._norm_name(jsonld_funder.get("name"))
            if jsonld_name and jsonld_name in schema_by_name:
                target = schema_by_name[jsonld_name]
                for key in ("identifier", "@id", "sameAs"):
                    if key in jsonld_funder and not target.get(key):
                        target[key] = jsonld_funder.get(key)
                        changed = True
            else:
                schema_funders.append(jsonld_funder)
                changed = True

        if changed:
            schema_doc["funder"] = schema_funders
        return changed

    @staticmethod
    def _extract_dates_from_dataverse_json(dv_json):
        """Extract dateCreated and dateModified from a dataverse_json export.

        Returns a dict with 'dateCreated' and/or 'dateModified' in YYYY-MM-DD format.
        """
        dates = {}
        version = dv_json.get("datasetVersion", {}) if isinstance(dv_json, dict) else {}
        if not isinstance(version, dict):
            return dates

        create_time = version.get("createTime")
        if isinstance(create_time, str) and create_time:
            try:
                dates["dateCreated"] = (
                    datetime.datetime.fromisoformat(create_time.replace("Z", "+00:00")).strftime("%Y-%m-%d")
                )
            except (ValueError, TypeError):
                pass

        update_time = version.get("lastUpdateTime")
        if isinstance(update_time, str) and update_time:
            try:
                dates["dateModified"] = (
                    datetime.datetime.fromisoformat(update_time.replace("Z", "+00:00")).strftime("%Y-%m-%d")
                )
            except (ValueError, TypeError):
                pass

        return dates

    @staticmethod
    def _extract_funding_from_dataverse_json(dv_json):
        """Extract funding info from the grantNumber fields of a dataverse_json export.

        Returns a list of funding dicts (Schema.org Grant-style), or an empty list.
        """
        if not isinstance(dv_json, dict):
            return []

        version = dv_json.get("datasetVersion", {})
        if not isinstance(version, dict):
            return []

        metadata_blocks = version.get("metadataBlocks", {})
        citation = metadata_blocks.get("citation", {}) if isinstance(metadata_blocks, dict) else {}
        fields = citation.get("fields", []) if isinstance(citation, dict) else []

        grant_field = None
        for field in fields:
            if isinstance(field, dict) and field.get("typeName") == "grantNumber":
                grant_field = field
                break

        if not grant_field or not isinstance(grant_field.get("value"), list):
            return []

        funding = []
        for grant in grant_field["value"]:
            if not isinstance(grant, dict):
                continue

            entry = {}
            funder = {}

            agency = grant.get("grantNumberAgency", {})
            if isinstance(agency, dict):
                # Use the expanded value for the human-readable funder name
                expanded = agency.get("expandedvalue", {})
                if isinstance(expanded, dict) and expanded.get("termName"):
                    funder["name"] = expanded["termName"]
                # The agency value field often contains a ROR URL
                agency_value = agency.get("value", "")
                if isinstance(agency_value, str) and agency_value.strip():
                    if agency_value.startswith("http"):
                        funder["@id"] = agency_value.strip()
                    elif not funder.get("name"):
                        # Fall back to raw value as funder name
                        funder["name"] = agency_value.strip()

            grant_number = grant.get("grantNumberValue", {})
            if isinstance(grant_number, dict) and grant_number.get("value"):
                entry["identifier"] = grant_number["value"]

            if funder:
                entry["funder"] = funder

            if entry:
                funding.append(entry)

        return funding

    def extract_schema_json(self, url):
        """Fallback: scrape schema.org JSON-LD from a dataset landing page.

        Looks for a <script type="application/ld+json"> tag and returns a parsed dict.
        """

        class SchemaScraper(HTMLParser):
            def __init__(self):
                super().__init__()
                self._reading_schema = False
                self.schemas = []
                self._buffer = []

            def handle_starttag(self, tag, attrs):
                attrs = dict(attrs)
                if tag == "script" and attrs.get("type") == "application/ld+json":
                    self._reading_schema = True
                    self._buffer = []

            def handle_data(self, data):
                if self._reading_schema:
                    self._buffer.append(data)

            def handle_endtag(self, tag):
                if tag == "script" and self._reading_schema:
                    raw = "".join(self._buffer).strip()
                    if raw:
                        self.schemas.append(raw)
                    self._reading_schema = False
                    self._buffer = []

        retries = 0
        backoff_time = BASE_DELAY
        logger.info(f"Scraping {url} for schema representation")

        while retries <= MAX_RETRIES:
            try:
                if REQUEST_DELAY > 0:
                    time.sleep(REQUEST_DELAY)
                req = self._session.get(url, headers=self._headers, timeout=30)
                if req.status_code == 429:
                    self._sleep_for_rate_limit(req, backoff_time)
                    retries += 1
                    backoff_time = min(MAX_DELAY, backoff_time * 2)
                    continue
                if not req.ok:
                    logger.info(f"unable to retrieve {url}, status code: {req.status_code}")
                    return False
                parser = SchemaScraper()
                parser.feed(req.text)
                if not parser.schemas:
                    return False

                for schema_text in parser.schemas:
                    raw = schema_text.strip().replace("\r", "").replace("\t", "")
                    try:
                        parsed = json.loads(raw)
                    except json.decoder.JSONDecodeError:
                        # Some pages include multiple JSON objects or formatting quirks; skip invalid blocks.
                        continue
                    selected = self._select_dataset_from_jsonld(parsed)
                    if isinstance(selected, dict):
                        return selected
                logger.info(f"scraped schema blocks were not valid JSON-LD for {url}")
                return False
            except requests.RequestException as request_exception:
                retries += 1
                if retries > MAX_RETRIES:
                    logger.info(
                        f"Failed to scrape {url} after {MAX_RETRIES} attempts due to {request_exception}"
                    )
                    return False
                logger.info(f"Scraping failed due to {request_exception}, retrying in {backoff_time} seconds...")
                time.sleep(backoff_time)
                backoff_time = min(MAX_DELAY, backoff_time * 2)

    def run_dataverse_schema_export(self, gid, url, verbose=False):
        """
        Run the dataverse schema.org export on a harvard dataset
        Returns the schema.org json object
        """
        dv_schema_export = f"{self.EXPORT_URL}&persistentId={gid}"
        logger.info(f"Running export using dataverse api - {dv_schema_export}")

        retries = 0
        backoff_time = BASE_DELAY

        while retries <= MAX_RETRIES:
            try:
                if REQUEST_DELAY > 0:
                    time.sleep(REQUEST_DELAY)
                req = self._session.get(dv_schema_export, headers=self._headers, timeout=30)
                if req.status_code == 429:
                    self._sleep_for_rate_limit(req, backoff_time)
                    retries += 1
                    backoff_time = min(MAX_DELAY, backoff_time * 2)
                    continue

                # Dataverse occasionally returns HTML/empty bodies (e.g., transient upstream errors)
                # which will raise JSONDecodeError. Log minimal diagnostics to aid debugging.
                if not req.ok:
                    # Retry transient server errors; do not spin forever on deterministic client errors.
                    if 500 <= req.status_code < 600:
                        req.raise_for_status()
                    logger.info(
                        "schema.org export returned non-OK status=%s content-type=%s url=%s",
                        req.status_code,
                        req.headers.get("Content-Type"),
                        dv_schema_export,
                    )
                    document = self.extract_schema_json(url)
                    if document and isinstance(document, dict):
                        logger.info(f"successfully scraped schema JSON-LD, {url}")
                        return document
                    return False

                content_type = (req.headers.get("Content-Type") or "").lower()
                if "json" not in content_type:
                    logger.info(
                        "schema.org export returned non-JSON content-type=%s status=%s url=%s",
                        req.headers.get("Content-Type"),
                        req.status_code,
                        dv_schema_export,
                    )
                    document = self.extract_schema_json(url)
                    if document and isinstance(document, dict):
                        logger.info(f"successfully scraped schema JSON-LD, {url}")
                        return document

                try:
                    res = req.json()
                except json.decoder.JSONDecodeError as e:
                    snippet = (req.text or "").strip().replace("\n", " ")[:200]
                    logger.info(
                        "schema.org export JSON decode failed status=%s content-type=%s url=%s snippet=%r",
                        req.status_code,
                        req.headers.get("Content-Type"),
                        dv_schema_export,
                        snippet,
                    )
                    raise e
                # odd case of a response with a status of ERROR
                if res.get("status") and res.get("status") == "ERROR":
                    document = self.extract_schema_json(url)
                    if document and isinstance(document, dict):
                        logger.info(f"successfully scraped schema JSON-LD, {url}")
                        return document
                    logger.info(
                        f"schema.org export failed on {dv_schema_export}, {url}, {req.status_code}, {res.get('status')}, {res.get('message')}"
                    )
                    return False
                else:
                    if isinstance(res, dict):
                        # Dataverse schema.org exporter may omit funder identifiers.
                        # Always require/augment funder IDs by scraping landing-page JSON-LD and merging.
                        if not self._has_funder_identifiers(res):
                            document = self.extract_schema_json(url)
                            if document and isinstance(document, dict):
                                if self._merge_funder_identifiers(res, document):
                                    logger.info(f"Augmented schema.org export with funder IDs from landing JSON-LD, {url}")
                        return res
                    return False
            except (requests.RequestException, json.decoder.JSONDecodeError) as e:
                retries += 1
                if retries > MAX_RETRIES:
                    logger.info(f"Failed to get {dv_schema_export} after {MAX_RETRIES} attempts due to {e}")
                    return False
                logger.info(f"Request failed due to {e}, retrying in {backoff_time} seconds...")
                time.sleep(backoff_time)
                # double the wait time, but cap at MAX_DELAY
                backoff_time = min(MAX_DELAY, backoff_time * 2)

    def fetch_dataverse_json(self, gid):
        """Fetch the dataverse_json export for a dataset.

        Returns the parsed JSON dict, or None if unavailable (e.g. harvested records).
        """
        dv_json_url = f"{self.DATAVERSE_JSON_EXPORT_URL}&persistentId={gid}"
        logger.info(f"Fetching dataverse_json export - {dv_json_url}")

        retries = 0
        backoff_time = BASE_DELAY

        while retries <= MAX_RETRIES:
            try:
                if REQUEST_DELAY > 0:
                    time.sleep(REQUEST_DELAY)
                req = self._session.get(dv_json_url, headers=self._headers, timeout=30)
                if req.status_code == 429:
                    self._sleep_for_rate_limit(req, backoff_time)
                    retries += 1
                    backoff_time = min(MAX_DELAY, backoff_time * 2)
                    continue
                if not req.ok:
                    logger.info(f"dataverse_json export returned status={req.status_code} for {gid}")
                    return None
                res = req.json()
                if isinstance(res, dict) and res.get("status") == "ERROR":
                    logger.info(f"dataverse_json export returned ERROR for {gid}: {res.get('message')}")
                    return None
                return res
            except (requests.RequestException, json.decoder.JSONDecodeError) as e:
                retries += 1
                if retries > MAX_RETRIES:
                    logger.info(f"Failed to fetch dataverse_json for {gid} after {MAX_RETRIES} attempts: {e}")
                    return None
                logger.info(f"dataverse_json fetch failed: {e}, retrying in {backoff_time}s")
                time.sleep(backoff_time)
                backoff_time = min(MAX_DELAY, backoff_time * 2)
        return None

    def compile_paginated_data(self, query_endpoint, per_page=400, start=0):
        logger.info(f"Compiling paginated data from {query_endpoint}")
        continue_paging = True

        while continue_paging:
            pager = per_page
            response = None
            retries = 0
            backoff_time = BASE_DELAY
            while retries <= MAX_RETRIES and not response:
                try:
                    url = f"{query_endpoint}&per_page={pager}&start={start}"
                    logger.info(f"Requesting {url}")
                    if REQUEST_DELAY > 0:
                        time.sleep(REQUEST_DELAY)
                    req = self._session.get(url, headers=self._headers, timeout=30)
                    if req.status_code == 429:
                        self._sleep_for_rate_limit(req, backoff_time)
                        retries += 1
                        backoff_time = min(MAX_DELAY, backoff_time * 2)
                        continue
                    req.raise_for_status()  # Raise an exception for HTTP errors (including 403)
                    response = req.json()
                    retries = 0  # Reset the retry counter after a successful request
                except (requests.RequestException, ValueError) as e:
                    pager = max(1, pager // 2)
                    logger.error(f"Error accessing {url}: {str(e)}")
                    retries += 1
                    if retries > MAX_RETRIES:
                        logger.error(f"Max retries reached for {url}. Skipping.")
                        return
                    logger.info(f"Retrying in {backoff_time} seconds...")
                    time.sleep(backoff_time)
                    # double the wait time, but cap at MAX_DELAY
                    backoff_time = min(MAX_DELAY, backoff_time * 2)
            if response:
                try:
                    total = response.get("data").get("total_count")
                    page_data = response.get("data").get("items")
                    logger.info(f"Retrieved {len(page_data)} items")
                    start += pager
                    for page in page_data:
                        # dataset pages use "global_id" key
                        if "global_id" in page.keys():
                            yield (page["global_id"], page["url"], page)
                        # dataverse pages have "identifier" key
                        elif "identifier" in page.keys():
                            # EXTRACT ID HERE
                            yield (page["identifier"], page["url"], page)
                    # data_pages.extend(page_data)
                    continue_paging = total and start < total
                except Exception as exception:
                    logger.info(f"passing datapage because of exception: {exception}")
            else:
                logger.error(f"Failed to get response for start={start} after {MAX_RETRIES} retries. Stopping pagination.")
                return

    def log_memory_usage(self):
        process = psutil.Process()
        memory_use = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        logging.info(f"Current memory usage: {memory_use:.2f} MB")

    def load_cache(self):
        """
        # For Dev
        # Testing ~ set initial_start_page to view different sources
        # initial_start_page=0 -- for harvard sources
        #   - https://dataverse.harvard.edu/api/search?q=*&type=dataset&per_page=100&start=0
        # initial_start_page=94000 -- for external sources (this may increase as more data is added to dataverse)
        #   - https://dataverse.harvard.edu/api/search?q=*&type=dataset&per_page=5&start=94000
        # initial_start_page=155000 -- for handle urls
        #   - https://dataverse.harvard.edu/api/search?q=*&type=dataset&per_page=5&start=155000
        """
        logger.info("Runing load_cache process....")

        # Initial memory usage log
        self.log_memory_usage()

        # Dev helpers: allow starting later and limiting record count for local runs.
        # Example: DATAVERSE_START=0 DATAVERSE_MAX_RECORDS=200
        initial_page_start = int(os.getenv("DATAVERSE_START", "0"))
        max_records = os.getenv("DATAVERSE_MAX_RECORDS")
        max_records = int(max_records) if max_records else None
        records_processed = 0
        handle_url_ct = 0
        schemas_gathered_ct = 0
        sleep_every = int(os.getenv("DATAVERSE_SLEEP_EVERY", "1000"))
        sleep_time = float(os.getenv("DATAVERSE_SLEEP_SECONDS", "5"))
        consecutive_failures = 0

        counters = {
            "harvard_attempted_export": 0,
            "harvard_export_ok": 0,
            "harvard_export_missing_at_id": 0,
            "harvard_export_failed": 0,
            "harvard_fallback_to_search": 0,
            "harvard_dv_json_ok": 0,
            "harvard_dv_json_failed": 0,
            "harvard_dv_json_funding_replaced": 0,
            "handle_records": 0,
            "external_search_records": 0,
        }

        query_endpoint = "https://dataverse.harvard.edu/api/search?q=*&type=dataset"

        # Iterate through the paginated data starting from the adjusted initial position
        for global_id, data_url, data_page in self.compile_paginated_data(
            query_endpoint, per_page=400, start=initial_page_start
        ):
            records_processed += 1

            if max_records and records_processed > max_records:
                logger.info(f"Reached DATAVERSE_MAX_RECORDS={max_records}; stopping early.")
                break

            if "https://hdl.handle.net/" in data_url:
                if "https://hdl.handle.net/1902.4" in data_url:
                    # https://hdl.handle.net/1902.4 : empty/dead-links (404 error) -- maybe handle 404 instead?
                    pass
                elif "https://hdl.handle.net/" in data_url:
                    handle_url_ct += 1
                    if data_page:
                        yield (data_url, json.dumps(data_page))
                        schemas_gathered_ct += 1
                        counters["handle_records"] += 1

            elif global_id.startswith("doi:10.7910"):
                # case: harvard data - doi:10.7910
                # run the built-in dataverse schema export on harvard dataverse sources
                counters["harvard_attempted_export"] += 1

                # Cooldown: if many consecutive exports fail, the server may be
                # throttling us.  Pause before making more requests.
                if consecutive_failures >= COOLDOWN_THRESHOLD:
                    logger.warning(
                        f"{consecutive_failures} consecutive export failures — "
                        f"cooling down for {COOLDOWN_SECONDS}s to avoid rate limiting"
                    )
                    time.sleep(COOLDOWN_SECONDS)
                    consecutive_failures = 0

                schema_record = self.run_dataverse_schema_export(global_id, data_url)
                if schema_record and isinstance(schema_record, dict):
                    consecutive_failures = 0  # reset on success
                    # Enrich with dataverse_json export (dates + funding)
                    dv_json = self.fetch_dataverse_json(global_id)
                    if dv_json:
                        counters["harvard_dv_json_ok"] += 1
                        # Merge dateCreated and dateModified
                        dates = self._extract_dates_from_dataverse_json(dv_json)
                        if "dateCreated" in dates:
                            schema_record["dateCreated"] = dates["dateCreated"]
                        if "dateModified" in dates:
                            schema_record["dateModified"] = dates["dateModified"]
                        # Merge funding: overwrite if dataverse_json has it
                        dv_funding = self._extract_funding_from_dataverse_json(dv_json)
                        if dv_funding:
                            schema_record.pop("funder", None)
                            schema_record["funding"] = dv_funding
                            counters["harvard_dv_json_funding_replaced"] += 1
                            logger.info(f"Replaced funding with dataverse_json data for {global_id}")
                    else:
                        counters["harvard_dv_json_failed"] += 1
                    if schema_record.get("@id"):
                        counters["harvard_export_ok"] += 1
                        logger.info(f"schema export passed on {data_url}")
                        yield (schema_record["@id"], json.dumps(schema_record))
                        schemas_gathered_ct += 1
                    else:
                        counters["harvard_export_missing_at_id"] += 1
                        # Critical fallback: do not drop Harvard records.
                        if data_page:
                            counters["harvard_fallback_to_search"] += 1
                            yield (data_url, json.dumps(data_page))
                            schemas_gathered_ct += 1
                        else:
                            counters["harvard_export_failed"] += 1
                else:
                    counters["harvard_export_failed"] += 1
                    consecutive_failures += 1
                    # Critical fallback: do not drop Harvard records.
                    if data_page:
                        counters["harvard_fallback_to_search"] += 1
                        yield (data_url, json.dumps(data_page))
                        schemas_gathered_ct += 1

            else:
                # case: outside data (non-harvard registered--therefore not standardized)
                if data_page:
                    # yield data available through dataverse api, not extracting from url
                    yield (data_url, json.dumps(data_page))
                    schemas_gathered_ct += 1
                    counters["external_search_records"] += 1

            if sleep_every > 0 and records_processed % sleep_every == 0:
                logger.info(
                    f"Processed {records_processed} datasets, going to sleep for {sleep_time} seconds to manage load..."
                )
                time.sleep(sleep_time)

            # Optional - log memory usage periodically, e.g., every 100 records
            if records_processed % 100 == 0:
                self.log_memory_usage()

        # Final memory usage log
        self.log_memory_usage()

        harvard_attempted = counters.get("harvard_attempted_export", 0) or 0
        harvard_export_ok = counters.get("harvard_export_ok", 0) or 0
        harvard_fallback = counters.get("harvard_fallback_to_search", 0) or 0
        harvard_missing_at_id = counters.get("harvard_export_missing_at_id", 0) or 0
        harvard_failed = counters.get("harvard_export_failed", 0) or 0

        if harvard_attempted > 0:
            logger.info(
                "Harvard (doi:10.7910) schema.org export summary: "
                f"attempted={harvard_attempted}, full_schema_ok={harvard_export_ok} "
                f"({(harvard_export_ok / harvard_attempted) * 100:.1f}%), "
                f"fell_back_to_search_record={harvard_fallback} "
                f"({(harvard_fallback / harvard_attempted) * 100:.1f}%), "
                f"export_missing_@id={harvard_missing_at_id}, export_failed={harvard_failed}"
            )

        overall_yielded = (
            (counters.get("handle_records", 0) or 0)
            + (counters.get("external_search_records", 0) or 0)
            + harvard_export_ok
            + harvard_fallback
        )
        if overall_yielded > 0:
            logger.info(
                "Overall yielded record composition: "
                f"full_schema_ok={harvard_export_ok} ({(harvard_export_ok / overall_yielded) * 100:.1f}%), "
                f"search_record={overall_yielded - harvard_export_ok} "
                f"({((overall_yielded - harvard_export_ok) / overall_yielded) * 100:.1f}%)"
            )
        logger.info(
            f"Processed {records_processed} datasets and yielded {schemas_gathered_ct} cached records. "
            f"Counters: {counters}"
        )

    def parse(self, records):
        start_time = time.process_time()
        parse_ct = 0
        counters = {
            "parsed_schema": 0,
            "parsed_search": 0,
            "skipped_schema": 0,
            "skipped_search": 0,
        }
        logger.info("Starting metadata parser...")
        # rec = ('doi:10.18738/T8/YJMLKO', '{"name": "ChIP-seq peak calls for epigenetic marks in GBM tumors", "type": "dataset", "url": "https://doi.org/10.18738/T8/YJMLKO", "global_id": "doi:10.18738/T8/YJMLKO", "description": "MACS2 narrowPeak files from ChIP-seq experiments for 11 primary GBM tumors, each targeting CTCF transcription factor marks and H3K27Ac, H3K27Me3, H3K4Me1, H3K4Me3, H3K9Ac, and H3K9Me3 histone modifications. See Methods section of doi:10.1158/0008-5472.CAN-17-1724 for more information.", "published_at": "2018-11-05T05:17:42Z", "publisher": "Texas Data Repository Harvested Dataverse", "citationHtml": "Battenhouse, Anna; Hall, Amelia Weber, 2018, \\"ChIP-seq peak calls for epigenetic marks in GBM tumors\\", <a href=\\"https://doi.org/10.18738/T8/YJMLKO\\" target=\\"_blank\\">https://doi.org/10.18738/T8/YJMLKO</a>, Texas Data Repository Dataverse", "identifier_of_dataverse": "tdr_harvested", "name_of_dataverse": "Texas Data Repository Harvested Dataverse", "citation": "Battenhouse, Anna; Hall, Amelia Weber, 2018, \\"ChIP-seq peak calls for epigenetic marks in GBM tumors\\", https://doi.org/10.18738/T8/YJMLKO, Texas Data Repository Dataverse", "storageIdentifier": "s3://10.18738/T8/YJMLKO", "keywords": ["Medicine, Health and Life Sciences", "glioblastoma", "bivalent", "enhancer", "epigenetic", "histone modification"], "subjects": [], "fileCount": 84, "versionId": 146549, "versionState": "RELEASED", "createdAt": "2018-11-05T05:17:42Z", "updatedAt": "2018-11-05T05:17:42Z", "contacts": [{"name": "", "affiliation": ""}], "authors": ["Battenhouse, Anna", "Hall, Amelia Weber"]}')
        # records = [rec]
        # logger.info(f"Starting parsing of {len(records)} records...")
        for record in records:
            dataset = None
            is_schema = False
            try:
                dataset = json.loads(record[1])
                is_schema = isinstance(dataset, dict) and ("@context" in dataset or "@type" in dataset)
            except Exception as error:
                counters["skipped_search"] += 1
                logger.info(f"skipping (invalid JSON) - {error}, record id - {record[0]}")
                continue

            try:
                # parse the schema.org export document
                if is_schema:
                    dataset_url = (
                        self._first_str(dataset.get("identifier"))
                        or self._first_str(dataset.get("@id"))
                        or self._first_str(dataset.get("url"))
                        or record[0]
                    )
                    dataset["url"] = dataset_url

                    doi_url = self._first_str(dataset.get("identifier")) or self._first_str(dataset.get("@id"))
                    if doi_url:
                        dataset["doi"] = doi_url
                        dataset.pop("identifier", None)
                    else:
                        dataset["doi"] = dataset_url
                    dataset["identifier"] = (
                        dataset.get("doi", "")
                        .replace("https://doi.org/", "")
                        .replace("http://doi.org/", "")
                        .replace("doi:", "")
                    )
                    dataset["_id"] = (
                        dataset.get("doi", "")
                        .replace("https://doi.org", "Dataverse")
                        .replace("http://doi.org", "Dataverse")
                        .replace("doi:", "Dataverse_")
                        .replace("/", "_")
                    )
                    dataset["dateModified"] = (
                        datetime.datetime.strptime(dataset["dateModified"], "%Y-%m-%d").date().isoformat()
                    )

                    if "datePublished" in dataset:
                        dataset["datePublished"] = (
                            datetime.datetime.strptime(dataset["datePublished"], "%Y-%m-%d").date().isoformat()
                        )

                    if "dateCreated" in dataset:
                        dataset["dateCreated"] = (
                            datetime.datetime.strptime(dataset["dateCreated"], "%Y-%m-%d").date().isoformat()
                        )

                    if dataset.get("author") is None:
                        dataset.pop("author", None)
                        if dataset.get("creator"):
                            dataset["author"] = dataset.pop("creator")
                            for data_dict in dataset["author"]:
                                if "affiliation" in data_dict.keys() and isinstance(data_dict["affiliation"], str):
                                    data_dict["affiliation"] = {"name": data_dict.pop("affiliation")}
                        else:
                            dataset.pop("creator", None)
                    else:
                        for data_dict in dataset["author"]:
                            if "affiliation" in data_dict.keys() and isinstance(data_dict["affiliation"], str):
                                data_dict["affiliation"] = {"name": data_dict.pop("affiliation")}
                        dataset.pop("creator", None)

                    if dataset.get("publisher") is None:
                        dataset.pop("publisher", None)
                        if dataset.get("provider"):
                            dataset["sdPublisher"] = [dataset.pop("provider")]
                        else:
                            dataset["sdPublisher"] = [dataset.pop("provider", None)]
                    else:
                        dataset["sdPublisher"] = [dataset.pop("publisher")]
                        dataset.pop("provider", None)

                    if "funder" in dataset and "funding" not in dataset:
                        dataset["funding"] = {"funder": dataset.pop("funder")}
                        for data_dict in dataset["funding"]["funder"]:
                            data_dict.pop("@type", None)
                    elif "funder" in dataset and "funding" in dataset:
                        # funding already set (e.g. from dataverse_json); drop raw funder
                        dataset.pop("funder", None)

                    if "description" in dataset:
                        description = dataset.pop("description")
                        if isinstance(description, list):
                            dataset["description"] = " ".join(description)
                        else:
                            dataset["description"] = description

                    if "keywords" in dataset and dataset.get("keywords"):
                        dataset["keywords"] = dataset.pop("keywords")[0]

                    if dataset.get("license"):
                        if type(dataset["license"]) is str:
                            pass
                        elif type(dataset["license"]) is dict:
                            if "url" in dataset["license"].keys():
                                dataset["license"] = dataset["license"]["url"]
                            elif "text" in dataset["license"].keys():
                                dataset["license"] = dataset["license"]["text"]
                            else:
                                dataset.pop("license", None)

                    if "temporalCoverage" in dataset and dataset["temporalCoverage"]:
                        dataset["temporalCoverage"] = [
                            {"@type": "TemporalInterval", "duration": dataset["temporalCoverage"][0]}
                        ]

                    if "spatialCoverage" in dataset:
                        dataset["spatialCoverage"] = [{"name": dataset.pop("spatialCoverage")[0]}]

                    if "distribution" in dataset:
                        for data_dict in dataset["distribution"]:
                            if "@id" in data_dict:
                                data_dict["@id"] = data_dict["@id"].strip("https://doi.org/")
                            if "fileFormat" in data_dict:
                                data_dict["encodingFormat"] = data_dict.pop("fileFormat")
                            if "identifier" in data_dict:
                                data_dict["contentUrl"] = data_dict.pop("identifier")
                            if "contentSize" in data_dict:
                                data_dict.pop("contentSize")
                            if "description" in data_dict and "https://" in data_dict["description"]:
                                data_dict.pop("description")

                    if "citation" in dataset:
                        cit_list = []
                        for data_dict in dataset["citation"]:
                            if "text" in data_dict:
                                cit_list.append({"citation": data_dict["text"]})
                        if cit_list:
                            dataset["citation"] = cit_list
                        else:
                            dataset.pop("citation", None)

                    dataset["includedInDataCatalog"] = {
                        "@type": "DataCatalog",
                        "name": "Harvard Dataverse",
                        "url": "https://dataverse.harvard.edu/",
                        "versionDate": datetime.datetime.today().strftime("%Y-%m-%d"),
                        "archivedAt": dataset_url,
                    }

                    dataset.pop("@id", None)
                    dataset.pop("version", None)

                    yield dataset
                    parse_ct += 1
                    counters["parsed_schema"] += 1
                else:
                    # parse abridged search metadata
                    dataset["@context"] = "http://schema.org"
                    if "type" in dataset:
                        dataset["@type"] = dataset.pop("type")
                    else:
                        dataset["@type"] = "Dataset"
                    dataset_url = dataset.get("url") or dataset.get("doi") or record[0]
                    if "global_id" in dataset:
                        dataset["doi"] = dataset.pop("global_id")
                    dataset["identifier"] = dataset.get("doi", "").replace("doi:", "")
                    dataset["_id"] = dataset.get("doi", "").replace("doi:", "Dataverse_").replace("/", "_")
                    dataset["includedInDataCatalog"] = {
                        "@type": "DataCatalog",
                        "name": "Harvard Dataverse",
                        "url": "https://dataverse.harvard.edu/",
                        "versionDate": datetime.datetime.today().strftime("%Y-%m-%d"),
                        "archivedAt": dataset_url
                    }

                    # have to strip Z from time to get into correct format
                    if "updatedAt" in dataset:
                        dataset["dateModified"] = dataset["updatedAt"].strip("Z")
                        dataset.pop("updatedAt")
                        dataset["dateModified"] = datetime.datetime.strptime(
                            dataset["dateModified"], "%Y-%m-%dT%H:%M:%S"
                        ).strftime("%Y-%m-%d")

                    if "published_at" in dataset:
                        dataset["datePublished"] = dataset["published_at"].strip("Z")
                        dataset.pop("published_at")
                        dataset["datePublished"] = datetime.datetime.strptime(
                            dataset["datePublished"], "%Y-%m-%dT%H:%M:%S"
                        ).strftime("%Y-%m-%d")

                    if "identifier_of_dataverse" in dataset:
                        dataset["sdPublisher"] = {
                            "identifier": dataset.pop("identifier_of_dataverse"),
                            "name": dataset.pop("name_of_dataverse"),
                        }

                    if "subjects" in dataset:
                        if not dataset["subjects"]:
                            pass
                        elif len(dataset["subjects"]) > 1:
                            dataset["topicCategory"] = []
                            for subject in dataset["subjects"]:
                                dataset["topicCategory"].append({"description": subject})
                        else:
                            dataset["topicCategory"] = {"description": dataset["subjects"][0]}
                        dataset.pop("subjects")

                    if "keywords" in dataset:
                        dataset["keywords"] = ",".join(dataset["keywords"])

                    if "createdAt" in dataset:
                        dataset["dateCreated"] = dataset["createdAt"].strip("Z")
                        dataset.pop("createdAt")
                        dataset["dateCreated"] = datetime.datetime.strptime(
                            dataset["dateCreated"], "%Y-%m-%dT%H:%M:%S"
                        ).strftime("%Y-%m-%d")

                    if "authors" in dataset:
                        if not dataset["authors"]:
                            dataset.pop("authors")
                        elif len(dataset["authors"]) > 1:
                            dataset["author"] = []
                            for author in dataset["authors"]:
                                author_dict = {"name": author}
                                dataset["author"].append(author_dict)
                            dataset.pop("authors")
                        else:
                            dataset["author"] = {"name": dataset.pop("authors")[0]}
                    if "citation" in dataset:
                        dataset["citation"] = {"citation": dataset.pop("citation")}

                    keys_to_remove = [
                        "dataSources",
                        "geographicCoverage",
                        "majorVersion",
                        "producers",
                        "publications",
                        "relatedMaterial",
                        "topicCategory",
                        "publisher",
                        "citationHtml",
                        "storageIdentifier",
                        "fileCount",
                        "versionId",
                        "versionState",
                        "contacts",
                    ]

                    for key in keys_to_remove:
                        dataset.pop(key, None)

                    yield dataset
                    parse_ct += 1
                    counters["parsed_search"] += 1
            except Exception as error:
                if is_schema:
                    counters["skipped_schema"] += 1
                    logger.info(f"skipping (schema parse error) - {error}, record id - {record[0]}")
                else:
                    counters["skipped_search"] += 1
                    logger.info(f"skipping (search parse error) - {error}, record id - {record[0]}")

        process_time = time.process_time() - start_time
        logger.info(
            f"Completed parsing individual metadata, {parse_ct} records parsed in {process_time:.2f} seconds. "
            f"Counters: {counters}"
        )
        logger.info("Document parsing complete")
