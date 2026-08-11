import gzip
import importlib.util
import itertools
import json
import logging
import sqlite3
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import Mock, patch


if not hasattr(itertools, "batched"):

    def batched(iterable, n):
        iterator = iter(iterable)
        while batch := tuple(itertools.islice(iterator, n)):
            yield batch

    itertools.batched = batched


HUB_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = HUB_ROOT / "utils" / "citations.py"


def _install_optional_dependency_stubs():
    try:
        import orjson  # noqa: F401
    except ImportError:
        orjson = types.ModuleType("orjson")
        orjson.loads = json.loads
        orjson.dumps = lambda value: json.dumps(value).encode()
        sys.modules["orjson"] = orjson

    try:
        import requests  # noqa: F401
    except ImportError:
        requests = types.ModuleType("requests")

        class RequestException(Exception):
            pass

        class HTTPError(RequestException):
            pass

        class Timeout(RequestException):
            pass

        class ConnectionError(RequestException):
            pass

        requests.RequestException = RequestException
        requests.HTTPError = HTTPError
        requests.Timeout = Timeout
        requests.ConnectionError = ConnectionError
        requests.exceptions = types.SimpleNamespace(RequestException=RequestException)
        requests.Session = type("Session", (), {})
        sys.modules["requests"] = requests

    try:
        import Bio  # noqa: F401
    except ImportError:
        bio = types.ModuleType("Bio")
        bio.Entrez = types.SimpleNamespace()
        bio.Medline = types.SimpleNamespace()
        sys.modules["Bio"] = bio


def _load_citations_module():
    _install_optional_dependency_stubs()

    config = types.ModuleType("config")
    config.GEO_API_KEY = ""
    config.GEO_EMAIL = "test@example.org"
    config.logger = logging.getLogger("test_citations")
    sys.modules["config"] = config

    package_name = "citation_test_utils"
    package = types.ModuleType(package_name)
    package.__path__ = [str(HUB_ROOT / "utils")]
    sys.modules[package_name] = package

    common = types.ModuleType(f"{package_name}.common")
    common.as_list = lambda value: [] if value is None else value if isinstance(value, list) else [value]
    common.dict_entries = lambda doc, field: (
        entry for entry in common.as_list(doc.get(field)) if isinstance(entry, dict)
    )
    common.retry = lambda _attempts, _delay: lambda function: function
    sys.modules[common.__name__] = common

    funding = types.ModuleType(f"{package_name}.funding")
    funding.standardize_funder = lambda name: {"name": name}
    sys.modules[funding.__name__] = funding

    terms = types.ModuleType(f"{package_name}.terms")
    terms.DB_PATH = "/tmp/test-pubtator.db"
    terms.get_species_details = lambda _name, _identifier: None
    terms.query_condition = lambda _name, _identifier=None: None
    sys.modules[terms.__name__] = terms

    spec = importlib.util.spec_from_file_location(f"{package_name}.citations", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


citations = _load_citations_module()


class MentionMatchingTests(unittest.TestCase):
    def test_matches_complete_mentions_not_substrings(self):
        self.assertTrue(citations._mentioned_in("SAM", ("the sam cohort",)))
        self.assertFalse(citations._mentioned_in("SAM", ("biological samples",)))
        self.assertFalse(citations._mentioned_in("MS", ("clinical symptoms",)))
        self.assertFalse(citations._mentioned_in("PP", ("plasmodium spp. infections",)))

    def test_rejects_ambiguous_short_mentions_but_keeps_safe_terms(self):
        for mention in ("ACTT-1", "CMV", "IGM", "IPA", "LNS", "MS", "NON", "PP", "SAM"):
            with self.subTest(mention=mention):
                self.assertTrue(citations._is_ambiguous_short_mention(mention))

        for mention in ("HIV", "COVID-19", "malaria"):
            with self.subTest(mention=mention):
                self.assertFalse(citations._is_ambiguous_short_mention(mention))

    def test_requires_resolved_disease_label_to_match_mention(self):
        self.assertFalse(
            citations._term_matches_mention(
                {
                    "name": "Renal cell carcinoma 1",
                    "alternateName": ["Hereditary clear cell renal carcinoma"],
                    "originalName": "ACTT-1",
                },
                "ACTT-1",
            )
        )
        self.assertFalse(
            citations._term_matches_mention(
                {"name": "Doxorubicin/Semustine/Streptozocin Regimen"},
                "SAM",
            )
        )
        self.assertTrue(citations._term_matches_mention({"name": "malaria"}, "Malaria"))
        self.assertTrue(citations._term_matches_mention({"name": "HIV infectious disease"}, "HIV"))
        self.assertTrue(
            citations._term_matches_mention(
                {
                    "name": "Severe acute respiratory syndrome coronavirus 2",
                    "alternateName": ["SARS-CoV2"],
                },
                "SARS-CoV-2",
            )
        )

    def test_matches_species_common_and_display_names(self):
        human = {
            "name": "Homo sapiens",
            "commonName": "Human",
            "displayName": "Human | Homo sapiens",
        }
        self.assertTrue(citations._term_matches_mention(human, "human"))
        self.assertTrue(citations._term_matches_mention(human, "Homo sapiens"))
        self.assertFalse(citations._term_matches_mention(human, "water"))

    def test_matches_simple_plural_species_mentions(self):
        cases = [
            ({"name": "Capra hircus", "commonName": "Goat"}, "goats"),
            ({"name": "Homo sapiens", "commonName": "Human"}, "humans"),
            ({"name": "Mus musculus", "commonName": "Mouse"}, "mice"),
            ({"name": "Plasmodium vivax", "commonName": "malaria parasite P. vivax"}, "P. vivax"),
            (
                {"name": "Plasmodium falciparum", "commonName": "malaria parasite P. falciparum"},
                "P.falciparum",
            ),
        ]
        for term, mention in cases:
            with self.subTest(mention=mention):
                self.assertTrue(citations._species_term_matches_mention(term, mention))

        self.assertFalse(
            citations._species_term_matches_mention(
                {"name": "Wallaconchis ater", "alternateName": ["Onchidium ater"]},
                "water",
            )
        )
        self.assertFalse(citations._species_term_matches_mention({"name": "Danaea vivax"}, "P. vivax"))
        self.assertFalse(citations._species_term_matches_mention({"name": "Phyllostachys vivax"}, "P. vivax"))
        self.assertFalse(citations._species_term_matches_mention({"name": "Phyllostachys falciparum"}, "P.falciparum"))

    def test_rejects_non_condition_ontology_concepts(self):
        for name in (
            "Death Domain",
            "Dead",
            "Death",
            "Feces",
            "How Often Experienced Abdominal Discomfort",
            "Mortality Rate",
            "Newborn",
            "Nutrition",
            "obsolete spontaneous abortion",
            "Protozoal",
        ):
            with self.subTest(name=name):
                self.assertTrue(citations._is_non_condition_term({"name": name}))

        self.assertFalse(citations._is_non_condition_term({"name": "malaria"}))


class RecordAugmentationTests(unittest.TestCase):
    def setUp(self):
        citations._incompatible_disease_terms.clear()

    def test_disease_augmentation_rejects_observed_false_positive_mentions(self):
        cases = [
            ("Biological samples were collected.", "SAM"),
            ("Adaptive COVID-19 Treatment Trial (ACTT-1)", "ACTT-1"),
            ("A daily supplement (LNS) was provided.", "LNS"),
            ("This was a clinical trial.", "clinical"),
        ]
        for description, mention in cases:
            with self.subTest(mention=mention):
                record = {"_id": mention, "description": description}
                with patch.object(citations, "get_disease_details") as get_details:
                    citations.update_record_disease(record, {"MESH:C000000": [mention]})
                get_details.assert_not_called()
                self.assertNotIn("healthCondition", record)

    def test_disease_augmentation_keeps_a_specific_matching_condition(self):
        record = {"_id": "malaria", "description": "Participants were treated for malaria."}
        term = {"name": "malaria", "identifier": "0005136", "fromPMID": True}
        with patch.object(citations, "get_disease_details", return_value=term):
            citations.update_record_disease(record, {"MESH:D008288": ["malaria"]})
        self.assertEqual(record["healthCondition"], [term])

    def test_species_augmentation_rejects_ambiguous_acronyms(self):
        record = {"_id": "cmv", "description": "Frequent cytomegalovirus (CMV) infection was observed."}
        with patch.object(citations, "get_species_details") as get_details:
            citations.update_record_species(record, {"12305": ["CMV"]})
        get_details.assert_not_called()
        self.assertNotIn("species", record)
        self.assertNotIn("infectiousAgent", record)

    def test_species_augmentation_rejects_incompatible_cached_mapping(self):
        record = {
            "_id": "wash-bangladesh",
            "description": "The intervention improved drinking water quality.",
        }
        wallaconchis = {
            "name": "Wallaconchis ater",
            "identifier": "2231505",
            "alternateName": ["Onchidium ater"],
            "classification": "host",
        }
        with (
            patch.object(citations, "pubtator_lookup", return_value=wallaconchis),
            patch.object(citations, "get_species_details") as get_details,
        ):
            citations.update_record_species(record, {"2231505": ["water"]})
        get_details.assert_not_called()
        self.assertNotIn("species", record)
        self.assertNotIn("infectiousAgent", record)

    def test_species_augmentation_keeps_matching_common_name(self):
        record = {"_id": "human-study", "description": "Samples were collected from human participants."}
        human = {
            "name": "Homo sapiens",
            "commonName": "Human",
            "identifier": "9606",
            "classification": "host",
        }
        with patch.object(citations, "pubtator_lookup", return_value=human):
            citations.update_record_species(record, {"9606": ["human"]})
        self.assertEqual(record["species"][0]["name"], "Homo sapiens")
        self.assertTrue(record["species"][0]["fromPMID"])

    def test_disease_augmentation_rejects_generic_non_conditions(self):
        for mention in (
            "death",
            "dead",
            "die",
            "died",
            "deaths",
            "dying",
            "feces",
            "food insecurity",
            "infected",
            "infections",
            "infectious",
            "inflammatory",
            "mortality",
            "newborn",
            "nutrition",
            "protozoal",
            "weight gain",
        ):
            with self.subTest(mention=mention):
                record = {"_id": mention, "description": f"The study measured {mention}."}
                with patch.object(citations, "get_disease_details") as get_details:
                    citations.update_record_disease(record, {"MESH:C000000": [mention]})
                get_details.assert_not_called()
                self.assertNotIn("healthCondition", record)

    def test_incompatible_cached_disease_mapping_is_ignored(self):
        bad_term = {
            "name": "Renal cell carcinoma 1",
            "identifier": "C538557",
            "originalName": "ACTT-1",
        }
        response = Mock()
        response.json.return_value = {"label": {"@value": "Renal cell carcinoma 1"}}
        with (
            patch.object(citations, "pubtator_lookup", return_value=bad_term),
            patch.object(citations, "query_condition", return_value=bad_term),
            patch.object(citations.requests, "get", return_value=response, create=True),
        ):
            self.assertIsNone(citations.get_disease_details("MESH:C538557", "ACTT-1"))
        self.assertIn(("C538557", "actt-1"), citations._incompatible_disease_terms)

    def test_multilingual_mesh_label_is_parsed_and_safely_rejected(self):
        response = Mock(status_code=200)
        response.json.return_value = {
            "identifier": "D002493",
            "label": [
                {"@language": "es", "@value": "Enfermedades del Sistema Nervioso Central"},
                {"@language": "en", "@value": "Central Nervous System Diseases"},
            ],
        }

        with (
            patch.object(citations, "pubtator_lookup", return_value=None),
            patch.object(citations, "query_condition", return_value=None),
            patch.object(citations, "pubtator_add") as add,
            patch.object(citations.requests, "get", return_value=response) as get,
        ):
            result = citations.get_disease_details("MESH:D002493", "neurological complications")

        self.assertIsNone(result)
        self.assertIn(("D002493", "neurological complications"), citations._incompatible_disease_terms)
        get.assert_called_once_with(
            "https://id.nlm.nih.gov/mesh/D002493.json",
            timeout=citations._MESH_TIMEOUT,
        )
        add.assert_not_called()

    def test_transient_mesh_failure_is_retried(self):
        response = Mock(status_code=200)
        response.json.return_value = {"label": {"@language": "en", "@value": "Malaria"}}

        with (
            patch.object(citations, "pubtator_lookup", return_value=None),
            patch.object(citations, "query_condition", return_value=None),
            patch.object(citations, "pubtator_add"),
            patch.object(
                citations.requests,
                "get",
                side_effect=[citations.requests.Timeout("timed out"), response],
            ) as get,
            patch.object(citations.time, "sleep") as sleep,
        ):
            result = citations.get_disease_details("MESH:D008288", "malaria")

        self.assertEqual(result["name"], "Malaria")
        self.assertTrue(result["fromPMID"])
        self.assertEqual(get.call_count, 2)
        sleep.assert_called_once_with(citations._MESH_REQUEST_RETRY_SECONDS)

    def test_invalid_mesh_json_is_not_retried(self):
        response = Mock(status_code=200)
        response.json.side_effect = ValueError("invalid JSON")

        with (
            patch.object(citations, "pubtator_lookup", return_value=None),
            patch.object(citations, "query_condition", return_value=None),
            patch.object(citations.requests, "get", return_value=response) as get,
            patch.object(citations.time, "sleep") as sleep,
        ):
            with self.assertRaisesRegex(ValueError, "invalid JSON"):
                citations.get_disease_details("MESH:D002493", "neurological complications")

        get.assert_called_once()
        sleep.assert_not_called()

    def test_non_transient_mesh_http_error_is_not_retried(self):
        response = Mock(status_code=404)
        response.raise_for_status.side_effect = citations.requests.HTTPError(response=response)

        with (
            patch.object(citations, "pubtator_lookup", return_value=None),
            patch.object(citations, "query_condition", return_value=None),
            patch.object(citations.requests, "get", return_value=response) as get,
            patch.object(citations.time, "sleep") as sleep,
        ):
            with self.assertRaises(citations.requests.HTTPError):
                citations.get_disease_details("MESH:D002493", "neurological complications")

        get.assert_called_once()
        sleep.assert_not_called()


class PubTatorCacheReplacementTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.connection = sqlite3.connect(Path(self.temp_dir.name) / "pmid.db")
        self.connection.execute(
            "CREATE TABLE disease_data " "(pmid TEXT, entity_id TEXT, names TEXT, PRIMARY KEY (pmid, entity_id))"
        )
        self.connection.execute("INSERT INTO disease_data VALUES ('old-pmid', 'old-id', 'old-name')")
        self.connection.commit()
        self.original_connection = citations._pmid_conn
        citations._pmid_conn = self.connection

    def tearDown(self):
        citations._pmid_conn = self.original_connection
        self.connection.close()
        self.temp_dir.cleanup()

    def test_complete_dump_atomically_replaces_old_rows(self):
        dump_path = Path(self.temp_dir.name) / "disease2pubtator3.gz"
        with gzip.open(dump_path, "wt") as handle:
            handle.write("new-pmid\tDisease\tnew-id\tfirst-name\tPubTator3\n")
            handle.write("new-pmid\tDisease\tnew-id\tnew-name\tPubTator3\n")

        with patch.object(citations, "PUBTATOR_DIR", self.temp_dir.name):
            citations._stream_and_store(dump_path.name, "disease")

        rows = self.connection.execute("SELECT pmid, entity_id, names FROM disease_data").fetchall()
        self.assertEqual(rows, [("new-pmid", "new-id", "new-name")])

    def test_failed_dump_keeps_active_table(self):
        with patch.object(citations, "_dump_rows", side_effect=RuntimeError("broken dump")):
            with self.assertRaisesRegex(RuntimeError, "broken dump"):
                citations._stream_and_store("unused.gz", "disease")

        rows = self.connection.execute("SELECT pmid, entity_id, names FROM disease_data").fetchall()
        self.assertEqual(rows, [("old-pmid", "old-id", "old-name")])
        staging = self.connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='disease_data_staging'"
        ).fetchone()
        self.assertIsNone(staging)


if __name__ == "__main__":
    unittest.main()
