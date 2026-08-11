import itertools
import json
import logging
import os
import sqlite3
import sys
import tempfile
import types
import unittest
from contextlib import contextmanager
from itertools import islice
from pathlib import Path
from unittest.mock import Mock, call, patch

# The checked-in development virtualenv is Python 3.10, while production uses
# Python 3.12. Supply itertools.batched so these tests can run in either.
if not hasattr(itertools, "batched"):

    def batched(iterable, n):
        iterator = iter(iterable)
        while batch := tuple(islice(iterator, n)):
            yield batch

    itertools.batched = batched

HUB_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HUB_ROOT))

# Unit tests do not need the hub's filesystem-creating configuration module.
config = types.ModuleType("config")
config.logger = logging.getLogger("test_descriptions")
config.token = ""
sys.modules["config"] = config

from utils import descriptions
from utils import prewarm_descriptions
from utils import terms
from utils.cache import SqliteCache, SqliteKeySet


SPECIES_RESPONSE = "Plasmodium falciparum\t-2\t5833\nhumans\t-2\t9606\nhumans\t-2\t9605"
DISEASE_RESPONSE = "malaria\t-26\tDOID:12365"
COMBINED_RESPONSE = f"{SPECIES_RESPONSE}\n{DISEASE_RESPONSE}"


class FakeResponse:
    def __init__(self, status_code, text="", headers=None):
        self.status_code = status_code
        self.text = text
        self.headers = headers or {}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise descriptions.requests.HTTPError(f"HTTP {self.status_code}", response=self)


class DescriptionExtractionTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = str(Path(self.temp_dir.name) / "extract.db")
        self.db_path_patch = patch.object(descriptions, "DB_PATH", self.db_path)
        self.db_path_patch.start()
        descriptions._reset_extract_circuit()

    def tearDown(self):
        descriptions._reset_extract_circuit()
        self.db_path_patch.stop()
        self.temp_dir.cleanup()

    @contextmanager
    def connect_cache(self):
        connection = sqlite3.connect(self.db_path)
        try:
            with connection:
                for statement in descriptions._RESPONSE_CACHE_DDL:
                    connection.execute(statement)
                yield connection
        finally:
            connection.close()

    def test_combines_only_types_the_document_needs_and_splits_cache_rows(self):
        docs = [
            {"_id": "both", "description": "Plasmodium falciparum causes malaria in humans."},
            {
                "_id": "disease-only",
                "description": "Patients with COVID-19.",
                "infectiousAgent": [{"name": "SARS-CoV-2", "identifier": "2697049"}],
            },
            {
                "_id": "species-only",
                "description": "SARS-CoV-2 infection.",
                "healthCondition": [{"name": "COVID-19", "identifier": "0080600"}],
            },
        ]

        def query(_description, entity_types):
            key = tuple(entity_types)
            return {
                (descriptions._SPECIES_TYPE, descriptions._DISEASE_TYPE): COMBINED_RESPONSE,
                (descriptions._DISEASE_TYPE,): "COVID-19\t-26\tDOID:0080600",
                (descriptions._SPECIES_TYPE,): "SARS-CoV-2\t-2\t2697049",
            }[key]

        with patch.object(descriptions, "query_extract_api", side_effect=query) as query_mock:
            descriptions._extract_entities(docs)

        self.assertEqual(
            query_mock.call_args_list,
            [
                call(docs[0]["description"], [descriptions._SPECIES_TYPE, descriptions._DISEASE_TYPE]),
                call(docs[1]["description"], [descriptions._DISEASE_TYPE]),
                call(docs[2]["description"], [descriptions._SPECIES_TYPE]),
            ],
        )
        self.assertEqual(docs[0]["species"][0]["name"], "Plasmodium falciparum")
        self.assertEqual(
            docs[0]["healthCondition"],
            [{"name": "malaria", "identifier": "DOID:12365", "fromEXTRACT": True}],
        )
        self.assertEqual(
            docs[1]["healthCondition"],
            [{"name": "COVID-19", "identifier": "DOID:0080600", "fromEXTRACT": True}],
        )
        self.assertEqual(
            docs[2]["infectiousAgent"] if "infectiousAgent" in docs[2] else docs[2]["species"],
            [{"name": "SARS-CoV-2", "identifier": "2697049", "fromEXTRACT": True}],
        )

        with self.connect_cache() as connection:
            both_species = connection.execute("SELECT text_response FROM species WHERE ndeid = 'both'").fetchone()[0]
            both_disease = connection.execute("SELECT text_response FROM disease WHERE ndeid = 'both'").fetchone()[0]
            self.assertEqual(both_species, SPECIES_RESPONSE)
            self.assertEqual(both_disease, DISEASE_RESPONSE)
            self.assertIsNone(connection.execute("SELECT 1 FROM species WHERE ndeid = 'disease-only'").fetchone())
            self.assertIsNone(connection.execute("SELECT 1 FROM disease WHERE ndeid = 'species-only'").fetchone())

    def test_existing_cache_rows_keep_the_old_shape_and_avoid_requests(self):
        with self.connect_cache() as connection:
            connection.execute("INSERT INTO species VALUES (?, ?)", ("cached", SPECIES_RESPONSE))
            connection.execute("INSERT INTO disease VALUES (?, ?)", ("cached", DISEASE_RESPONSE))

        doc = {"_id": "CACHED", "description": "Plasmodium falciparum causes malaria."}
        with patch.object(descriptions, "query_extract_api") as query_mock:
            descriptions._extract_entities([doc])

        query_mock.assert_not_called()
        self.assertEqual(
            doc["healthCondition"],
            [{"name": "malaria", "identifier": "DOID:12365", "fromEXTRACT": True}],
        )
        self.assertEqual(doc["species"][0]["identifier"], "5833")

    def test_successful_empty_response_is_cached(self):
        doc = {"_id": "empty", "description": "No recognized disease.", "species": [{"name": "existing"}]}
        with patch.object(descriptions, "query_extract_api", return_value=""):
            descriptions._extract_entities([doc])

        with self.connect_cache() as connection:
            row = connection.execute("SELECT text_response FROM disease WHERE ndeid = 'empty'").fetchone()
        self.assertEqual(row, ("",))

    def test_failed_request_is_temporarily_cached(self):
        doc = {"_id": "failure", "description": "Transient failure.", "species": [{"name": "existing"}]}
        with patch.object(descriptions, "query_extract_api", side_effect=descriptions.requests.Timeout("timed out")):
            descriptions._extract_entities([doc])

        with self.connect_cache() as connection:
            response_row = connection.execute("SELECT text_response FROM disease WHERE ndeid = 'failure'").fetchone()
            failure_row = connection.execute(
                "SELECT entity_type, retry_after, error FROM extract_failures WHERE ndeid = 'failure'"
            ).fetchone()
        self.assertIsNone(response_row)
        self.assertEqual(failure_row[0], descriptions._DISEASE_TYPE)
        self.assertGreater(failure_row[1], descriptions.time.time())
        self.assertIn("Timeout", failure_row[2])

        with patch.object(descriptions, "query_extract_api") as query_mock:
            descriptions._extract_entities([doc])
        query_mock.assert_not_called()

    def test_query_retries_transient_status_and_reuses_supplied_session(self):
        session = Mock()
        session.post.side_effect = [
            FakeResponse(503, headers={"Retry-After": "0"}),
            FakeResponse(200, COMBINED_RESPONSE),
        ]

        with patch.object(descriptions.time, "sleep") as sleep_mock:
            result = descriptions.query_extract_api("description", ["-2", "-26"], session=session)

        self.assertEqual(result, COMBINED_RESPONSE)
        self.assertEqual(session.post.call_count, 2)
        self.assertEqual(session.post.call_args.kwargs["data"]["entity_types"], "-2 -26")
        self.assertEqual(session.post.call_args.kwargs["timeout"], descriptions._EXTRACT_TIMEOUT)
        sleep_mock.assert_called_once_with(0.0)

    def test_query_does_not_retry_non_transient_client_error(self):
        session = Mock()
        session.post.return_value = FakeResponse(400)
        with self.assertRaises(descriptions.requests.HTTPError):
            descriptions.query_extract_api("description", "-26", session=session)
        self.assertEqual(session.post.call_count, 1)

    def test_query_uses_post_when_the_encoded_url_would_be_too_long(self):
        session = Mock()
        session.post.return_value = FakeResponse(200, DISEASE_RESPONSE)

        result = descriptions.query_extract_api("x" * 8000, "-26", session=session)

        self.assertEqual(result, DISEASE_RESPONSE)
        self.assertEqual(session.post.call_count, 1)
        self.assertEqual(session.post.call_args.kwargs["data"]["entity_types"], "-26")
        self.assertEqual(session.post.call_args.kwargs["timeout"], descriptions._EXTRACT_TIMEOUT)

    def test_forbidden_response_opens_circuit_without_repeating_requests(self):
        session = Mock()
        session.post.return_value = FakeResponse(403)

        with self.assertRaises(descriptions.requests.HTTPError):
            descriptions.query_extract_api("description", "-26", session=session)
        with self.assertRaises(descriptions.ExtractCircuitOpen):
            descriptions.query_extract_api("another description", "-26", session=session)

        self.assertEqual(session.post.call_count, 1)

    @unittest.skipUnless(os.environ.get("RUN_EXTRACT_LIVE_TESTS") == "1", "set RUN_EXTRACT_LIVE_TESTS=1")
    def test_live_combined_responses_equal_separate_responses(self):
        sample_descriptions = [
            "Plasmodium falciparum causes malaria in humans.",
            "Mus musculus models are used to study Alzheimer disease and Parkinson disease.",
            "Human lung epithelial cells infected with SARS-CoV-2 from patients with COVID-19.",
        ]
        for description in sample_descriptions:
            with self.subTest(description=description):
                species = descriptions.query_extract_api(description, descriptions._SPECIES_TYPE)
                disease = descriptions.query_extract_api(description, descriptions._DISEASE_TYPE)
                combined = descriptions.query_extract_api(
                    description,
                    [descriptions._SPECIES_TYPE, descriptions._DISEASE_TYPE],
                )
                self.assertEqual(species, descriptions._response_for_entity_type(combined, descriptions._SPECIES_TYPE))
                self.assertEqual(disease, descriptions._response_for_entity_type(combined, descriptions._DISEASE_TYPE))


class DiseaseIdentifierLookupTests(unittest.TestCase):
    def test_extract_ontology_id_resolves_directly_without_text_search(self):
        response = {
            "hits": [
                {
                    "_id": "DOID:9191",
                    "label": "diabetic macular edema",
                    "synonym": {},
                }
            ]
        }

        with patch.object(terms, "_retry_request", return_value=response) as request:
            result = terms.query_condition("DME", ontology_id="DOID:9191")

        request.assert_called_once_with("https://biothings.transltr.io/doid/query?q=_id%3A%22DOID%3A9191%22&limit=1")
        self.assertEqual(result["identifier"], "9191")
        self.assertEqual(result["inDefinedTermSet"], "DOID")

    def test_extract_ontology_id_does_not_fall_back_to_ambiguous_name_search(self):
        empty_response = {"hits": []}

        with patch.object(terms, "_retry_request", return_value=empty_response) as request:
            result = terms.query_condition("DME", ontology_id="DOID:9191")

        self.assertIsNone(result)
        request.assert_called_once()


class DescriptionPrewarmTests(unittest.TestCase):
    def test_prewarm_processes_inference_records_without_type_filter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "records.ndjson"
            path.write_text(
                json.dumps(
                    {
                        "_id": "gxa_record",
                        "@type": "Inference",
                        "description": "Malaria differential gene expression.",
                        "species": [{"name": "Homo sapiens"}],
                    }
                )
                + "\n"
            )

            with (
                patch.object(prewarm_descriptions, "reset_caches"),
                patch.object(prewarm_descriptions, "augment_from_descriptions") as augment,
            ):
                total = prewarm_descriptions.prewarm_descriptions(path)

            self.assertEqual(total, 1)
            batch = augment.call_args.args[0]
            self.assertEqual(batch[0]["@type"], "Inference")
            self.assertEqual(augment.call_args.kwargs, {"filter_supported_types": False})

    def test_prewarm_streams_batches_without_modifying_input(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "records.ndjson"
            original = b"".join(
                json.dumps({"_id": str(index), "description": f"description {index}"}).encode() + b"\n"
                for index in range(5)
            )
            path.write_bytes(original)
            batches = []

            def augment(batch):
                batches.append([doc["_id"] for doc in batch])

            with patch.object(prewarm_descriptions, "reset_caches"):
                total = prewarm_descriptions.prewarm_descriptions(path, batch_size=2, augment=augment)

            self.assertEqual(total, 5)
            self.assertEqual(batches, [["0", "1"], ["2", "3"], ["4"]])
            self.assertEqual(path.read_bytes(), original)

    def test_prewarm_limit(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "records.ndjson"
            path.write_text("".join(json.dumps({"_id": str(index)}) + "\n" for index in range(5)))
            batches = []
            with patch.object(prewarm_descriptions, "reset_caches"):
                total = prewarm_descriptions.prewarm_descriptions(
                    path,
                    batch_size=2,
                    limit=3,
                    augment=lambda batch: batches.append(len(batch)),
                )
            self.assertEqual(total, 3)
            self.assertEqual(batches, [2, 1])


class SqliteBatchCacheTests(unittest.TestCase):
    def test_cache_put_many_persists_all_values(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = str(Path(temp_dir) / "cache.db")
            cache = SqliteCache(path, "details")
            cache.put_many({"Alpha": {"identifier": "1"}, "Beta": {"identifier": "2"}})

            reloaded = SqliteCache(path, "details", memoize=False)
            self.assertEqual(
                reloaded.get_many(["alpha", "beta"]),
                {
                    "alpha": {"identifier": "1"},
                    "beta": {"identifier": "2"},
                },
            )

    def test_key_set_add_many_persists_all_keys(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = str(Path(temp_dir) / "cache.db")
            keys = SqliteKeySet(path, "negative")
            keys.add_many(["Alpha", "Beta", "Alpha"])

            reloaded = SqliteKeySet(path, "negative")
            self.assertEqual(reloaded.known(["alpha", "beta", "gamma"]), {"alpha", "beta"})


if __name__ == "__main__":
    unittest.main()
