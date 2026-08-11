import importlib.util
import itertools
import logging
import sys
import types
import unittest
from itertools import islice
from pathlib import Path
from unittest.mock import Mock, patch

# The checked-in development virtualenv is Python 3.10, while production uses
# Python 3.12. Supply itertools.batched so these tests can run in either.
if not hasattr(itertools, "batched"):

    def batched(iterable, n):
        iterator = iter(iterable)
        while batch := tuple(islice(iterator, n)):
            yield batch

    itertools.batched = batched


HUB_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = HUB_ROOT / "utils" / "descriptions.py"


class DummyCache:
    def __init__(self, *_args, **_kwargs):
        self.values = {}

    def get(self, key):
        return self.values.get(key)

    def get_many(self, keys):
        return {key: self.values[key] for key in keys if key in self.values}

    def put(self, key, value):
        self.values[key] = value

    def put_many(self, values):
        self.values.update(values)

    def reset(self):
        self.values.clear()


class DummyKeySet(DummyCache):
    def known(self, keys):
        return set(keys) & self.values.keys()

    def __contains__(self, key):
        return key in self.values

    def add(self, key):
        self.values[key] = True

    def add_many(self, keys):
        self.values.update({key: True for key in keys})


def _load_descriptions_module():
    config = types.ModuleType("config")
    config.logger = logging.getLogger("test_description_filters")
    sys.modules["config"] = config

    package_name = "description_filter_test_utils"
    package = types.ModuleType(package_name)
    package.__path__ = [str(HUB_ROOT / "utils")]
    sys.modules[package_name] = package

    cache = types.ModuleType(f"{package_name}.cache")
    cache.SqliteCache = DummyCache
    cache.SqliteKeySet = DummyKeySet
    sys.modules[cache.__name__] = cache

    common = types.ModuleType(f"{package_name}.common")
    common.as_list = lambda value: [] if value is None else value if isinstance(value, list) else [value]
    common.sqlite = Mock()
    common.supports_description_enrichment = Mock(return_value=True)
    sys.modules[common.__name__] = common

    taxonomy = types.ModuleType(f"{package_name}.taxonomy")
    taxonomy.classify_from_lineage = lambda _lineage: None
    sys.modules[taxonomy.__name__] = taxonomy

    terms = types.ModuleType(f"{package_name}.terms")
    terms.DB_PATH = "/tmp/test-pubtator.db"
    terms.SPECIES_CACHE_DB_PATH = "/tmp/test-species.db"
    terms.fetch_taxon = Mock()
    terms.normalize_taxon_id = lambda identifier: (
        str(identifier).split("*")[-1].split(":")[-1].strip()
        if identifier and str(identifier).split("*")[-1].split(":")[-1].strip().isdigit()
        else None
    )
    terms.query_condition = Mock()
    sys.modules[terms.__name__] = terms

    spec = importlib.util.spec_from_file_location(f"{package_name}.descriptions", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


descriptions = _load_descriptions_module()


class DescriptionFilterTests(unittest.TestCase):
    def setUp(self):
        descriptions.SPECIES_DETAILS = DummyCache()
        descriptions.TAXON_DETAILS = DummyCache()
        descriptions.FAILED_TAXA = DummyKeySet()
        descriptions.REJECTED_SPECIES = DummyKeySet()
        descriptions.HEALTH_CONDITIONS = DummyCache()
        descriptions.NEGATIVE_DISEASES = DummyKeySet()

    def test_complete_mention_check_rejects_substrings(self):
        self.assertFalse(descriptions.mentioned_in("arge", ("A Large-scale malaria study",)))
        self.assertTrue(descriptions.mentioned_in("P. vivax", ("Severe P. vivax malaria",)))

    def test_generic_vector_taxon_is_dropped(self):
        response = "\n".join(
            [
                "vector\t-2\t2971083",
                "arge\t-2\t95269",
                "metagenome\t-2\t256318",
                "age strata\t-2\t1208515",
                "Matara\t-2\t2612617",
                "Kerala\t-2\t2249684",
                "transformation\t-2\t2839062",
                "scleroderma\t-2\t68787",
                "syncope\t-2\t1271638",
                "Latina\t-2\t1325907",
                "Human Microbiome\t-2\t646099",
                "Venus\t-2\t55714",
                "Napo\t-2\t706958",
                "human\t-2\t9606",
            ]
        )
        self.assertEqual(
            descriptions._tagged_entities(response, descriptions._SPECIES_TYPE),
            [("human", "9606")],
        )

    def test_extract_candidate_ids_are_preserved_for_later_disambiguation(self):
        response = "\n".join(
            [
                "SCV\t-2\t135656",
                "HCV\t-2\t3052230",
                "FDH\t-26\tDOID:0050475",
                "HNPP\t-26\tDOID:0050448",
                "HIV\t-26\tDOID:526",
            ]
        )

        self.assertEqual(
            descriptions._tagged_entities(response, descriptions._SPECIES_TYPE),
            [("SCV", "135656"), ("HCV", "3052230")],
        )
        self.assertEqual(
            descriptions._tagged_entities(response, descriptions._DISEASE_TYPE),
            [
                ("FDH", "DOID:0050475"),
                ("HNPP", "DOID:0050448"),
                ("HIV", "DOID:526"),
            ],
        )

    def test_drop_rules_only_remove_extract_terms(self):
        docs = [
            {
                "_id": "extracted",
                "species": [
                    {
                        "name": "Agestrata",
                        "identifier": "1208515",
                        "fromEXTRACT": True,
                    }
                ],
            },
            {
                "_id": "source-provided",
                "species": [
                    {
                        "name": "Agestrata",
                        "identifier": "1208515",
                        "fromEXTRACT": True,
                    },
                    {
                        "name": "Agestrata",
                        "identifier": "1208515",
                        "isCurated": True,
                    },
                ],
            },
        ]

        descriptions._dedupe_species(docs)

        self.assertNotIn("species", docs[0])
        self.assertEqual(docs[1]["species"][0]["identifier"], "1208515")

    def test_ambiguous_taxon_acronym_requires_authoritative_casing(self):
        self.assertFalse(
            descriptions._species_candidate_matches_mention(
                {"name": "Uga", "identifier": "2710078"},
                "UGA",
            )
        )
        self.assertFalse(
            descriptions._species_candidate_matches_mention(
                {
                    "name": "Strawberry crinkle virus",
                    "alternateName": ["SCV"],
                    "identifier": "135656",
                },
                "scV",
            )
        )
        self.assertTrue(
            descriptions._species_candidate_matches_mention(
                {
                    "name": "Strawberry crinkle virus",
                    "alternateName": ["SCV"],
                    "identifier": "135656",
                },
                "SCV",
            )
        )
        self.assertTrue(
            descriptions._species_candidate_matches_mention(
                {
                    "name": "Hepacivirus hominis",
                    "alternateName": ["HCV"],
                    "identifier": "3052230",
                },
                "HCV",
            )
        )

    def test_known_good_taxa_are_not_filtered(self):
        docs = [
            {
                "_id": "known-good",
                "species": [
                    {
                        "name": "Bulinus truncatus",
                        "identifier": "55810",
                        "classification": "host",
                        "fromEXTRACT": True,
                    }
                ],
                "infectiousAgent": [
                    {
                        "name": "Plasmodium vivax",
                        "identifier": "5855",
                        "classification": "infectiousAgent",
                        "fromEXTRACT": True,
                    }
                ],
            }
        ]

        descriptions._dedupe_species(docs)

        self.assertEqual(docs[0]["species"][0]["identifier"], "55810")
        self.assertEqual(docs[0]["infectiousAgent"][0]["identifier"], "5855")

    def test_incompatible_cached_mapping_is_removed_if_no_safe_resolution(self):
        descriptions.SPECIES_DETAILS.put(
            "P. vivax",
            {"name": "Danaea vivax", "identifier": "491846", "classification": "host"},
        )
        docs = [{"_id": "malaria", "species": [{"name": "P. vivax", "identifier": "756229", "fromEXTRACT": True}]}]

        with patch.object(descriptions, "_resolve_missing_species", return_value=[]) as resolve:
            descriptions._standardize_extracted_species(docs)

        self.assertNotIn("species", docs[0])
        self.assertEqual(resolve.call_args.args[0], ["P. vivax"])
        self.assertEqual(resolve.call_args.args[1], {"P. vivax": ["756229"]})

    def test_ambiguous_species_stub_is_removed_without_taxonomy_calls(self):
        descriptions.SPECIES_DETAILS.put(
            "SCV",
            {
                "name": "Strawberry crinkle virus",
                "alternateName": ["SCV"],
                "identifier": "135656",
                "classification": "infectiousAgent",
            },
        )
        docs = [
            {
                "_id": "small-colony-variant",
                "species": [{"name": "SCV", "identifier": "135656", "fromEXTRACT": True}],
            }
        ]

        with patch.object(descriptions, "_resolve_missing_species") as resolve:
            descriptions._standardize_extracted_species(docs)

        resolve.assert_not_called()
        self.assertNotIn("species", docs[0])

    def test_ambiguous_species_acronym_is_kept_when_context_spells_it_out(self):
        descriptions.SPECIES_DETAILS.put(
            "SCV",
            {
                "name": "Strawberry crinkle virus",
                "alternateName": ["SCV"],
                "identifier": "135656",
                "classification": "infectiousAgent",
            },
        )
        docs = [
            {
                "_id": "plant-virus",
                "description": "Strawberry crinkle virus (SCV) was detected.",
                "infectiousAgent": [{"name": "SCV", "identifier": "135656", "fromEXTRACT": True}],
            }
        ]

        descriptions._standardize_extracted_species(docs)

        self.assertEqual(docs[0]["infectiousAgent"][0]["identifier"], "135656")

    def test_extract_candidates_select_compatible_abbreviated_taxon(self):
        candidates = {
            "756229": {"name": "Phyllostachys vivax", "identifier": "756229"},
            "5855": {
                "name": "Plasmodium vivax",
                "commonName": "malaria parasite P. vivax",
                "identifier": "5855",
                "classification": "infectiousAgent",
                "originalName": "P. vivax",
            },
        }

        with patch.object(
            descriptions,
            "get_species_details",
            side_effect=lambda original_name, identifier: dict(candidates[identifier], originalName=original_name),
        ):
            resolved = descriptions._resolve_missing_species(
                ["P. vivax"],
                {"P. vivax": ["756229", "5855"]},
            )

        self.assertEqual([term["identifier"] for term in resolved], ["5855"])
        self.assertEqual(descriptions.SPECIES_DETAILS.get("P. vivax")["name"], "Plasmodium vivax")

    def test_extract_candidates_accept_abbreviation_without_space_after_period(self):
        candidate = {
            "name": "Plasmodium falciparum",
            "commonName": "malaria parasite P. falciparum",
            "identifier": "5833",
            "classification": "infectiousAgent",
        }

        with patch.object(
            descriptions,
            "get_species_details",
            side_effect=lambda original_name, _identifier: dict(candidate, originalName=original_name),
        ):
            resolved = descriptions._resolve_missing_species(
                ["P.falciparum"],
                {"P.falciparum": ["5833"]},
            )

        self.assertEqual([term["identifier"] for term in resolved], ["5833"])
        self.assertEqual(descriptions.SPECIES_DETAILS.get("P.falciparum")["name"], "Plasmodium falciparum")

    def test_rejected_extract_candidates_are_not_retried(self):
        candidate = {
            "name": "Uga",
            "identifier": "2710078",
            "classification": "host",
        }

        with patch.object(
            descriptions,
            "get_species_details",
            return_value=candidate,
        ) as get_details:
            first = descriptions._resolve_missing_species(["UGA"], {"UGA": ["2710078"]})
            second = descriptions._resolve_missing_species(["UGA"], {"UGA": ["2710078"]})

        self.assertEqual(first, [])
        self.assertEqual(second, [])
        get_details.assert_called_once()

    def test_new_candidate_set_retries_a_previously_rejected_mention(self):
        candidates = {
            "756229": {
                "name": "Phyllostachys vivax",
                "identifier": "756229",
                "classification": "host",
            },
            "5855": {
                "name": "Plasmodium vivax",
                "commonName": "malaria parasite P. vivax",
                "identifier": "5855",
                "classification": "infectiousAgent",
            },
        }

        with patch.object(
            descriptions,
            "get_species_details",
            side_effect=lambda _original_name, identifier: candidates[identifier],
        ) as get_details:
            first = descriptions._resolve_missing_species(
                ["P. vivax"],
                {"P. vivax": ["756229"]},
            )
            second = descriptions._resolve_missing_species(
                ["P. vivax"],
                {"P. vivax": ["756229", "5855"]},
            )

        self.assertEqual(first, [])
        self.assertEqual([term["identifier"] for term in second], ["5855"])
        self.assertEqual(get_details.call_count, 2)

    def test_transient_taxonomy_failure_is_retryable(self):
        candidate = {
            "name": "Plasmodium falciparum",
            "identifier": "5833",
            "classification": "infectiousAgent",
        }

        with patch.object(
            descriptions,
            "get_species_details",
            side_effect=[descriptions.requests.Timeout("timed out"), candidate],
        ) as get_details:
            first = descriptions._resolve_missing_species(
                ["Plasmodium falciparum"],
                {"Plasmodium falciparum": ["5833"]},
            )
            second = descriptions._resolve_missing_species(
                ["Plasmodium falciparum"],
                {"Plasmodium falciparum": ["5833"]},
            )

        self.assertEqual(first, [])
        self.assertEqual([term["identifier"] for term in second], ["5833"])
        self.assertEqual(get_details.call_count, 2)

    def test_incompatible_extracted_disease_mapping_is_removed(self):
        descriptions.HEALTH_CONDITIONS.put(
            "hypersensitivity",
            {"name": "drug allergy", "identifier": "0000775", "alternateName": ["allergy to drug"]},
        )
        docs = [{"_id": "skin", "healthCondition": [{"name": "hypersensitivity"}]}]

        with patch.object(
            descriptions,
            "query_condition",
            return_value={"name": "drug allergy", "identifier": "0000775"},
        ):
            descriptions._standardize_extracted_diseases(docs)

        self.assertNotIn("healthCondition", docs[0])

    def test_disease_candidate_ids_and_context_disambiguate_acronyms(self):
        docs = [
            {
                "_id": "chemistry",
                "description": "FIDDLEHEAD (FDH) regulates Arabidopsis cuticle biogenesis.",
                "healthCondition": [
                    {
                        "name": "FDH",
                        "identifier": "DOID:2120",
                        "fromEXTRACT": True,
                    }
                ],
            },
            {
                "_id": "neuropathy",
                "description": "Hereditary neuropathy with liability to pressure palsies (HNPP).",
                "healthCondition": [
                    {
                        "name": "HNPP",
                        "identifier": "DOID:0060843",
                        "fromEXTRACT": True,
                    }
                ],
            },
        ]
        resolved = {
            "DOID:2120": {
                "name": "focal dermal hypoplasia",
                "alternateName": ["FDH"],
                "identifier": "2120",
                "inDefinedTermSet": "DOID",
            },
            "DOID:0060843": {
                "name": "hereditary neuropathy with liability to pressure palsies",
                "alternateName": ["HNPP"],
                "identifier": "0060843",
                "inDefinedTermSet": "DOID",
            },
        }

        with patch.object(
            descriptions,
            "query_condition",
            side_effect=lambda _name, *, ontology_id: resolved[ontology_id],
        ) as query_condition:
            descriptions._standardize_extracted_diseases(docs)

        self.assertEqual(query_condition.call_count, 2)
        self.assertNotIn("healthCondition", docs[0])
        self.assertEqual(docs[1]["healthCondition"][0]["identifier"], "0060843")

    def test_extract_disease_id_bypasses_incompatible_name_only_cache_entry(self):
        descriptions.HEALTH_CONDITIONS.put(
            "dme",
            {
                "name": "Dimethyl Ether",
                "alternateName": ["DME"],
                "identifier": "C1645",
                "inDefinedTermSet": "NCIT",
            },
        )
        docs = [
            {
                "_id": "retina",
                "description": "Diabetic macular edema (DME) was assessed.",
                "healthCondition": [
                    {
                        "name": "DME",
                        "identifier": "DOID:9191",
                        "fromEXTRACT": True,
                    }
                ],
            }
        ]
        resolved = {
            "name": "diabetic macular edema",
            "alternateName": ["DME"],
            "identifier": "9191",
            "inDefinedTermSet": "DOID",
        }

        with patch.object(descriptions, "query_condition", return_value=resolved) as query_condition:
            descriptions._standardize_extracted_diseases(docs)

        query_condition.assert_called_once_with("DME", ontology_id="DOID:9191")
        self.assertEqual(docs[0]["healthCondition"][0]["identifier"], "9191")

    def test_full_and_acronym_mentions_share_the_candidate_id_cache(self):
        docs = [
            {
                "_id": "retina",
                "description": "Diabetic macular edema (DME) was assessed.",
                "healthCondition": [
                    {
                        "name": "diabetic macular edema",
                        "identifier": "DOID:9191",
                        "fromEXTRACT": True,
                    },
                    {
                        "name": "DME",
                        "identifier": "DOID:9191",
                        "fromEXTRACT": True,
                    },
                ],
            }
        ]
        resolved = {
            "name": "diabetic macular edema",
            "identifier": "9191",
            "inDefinedTermSet": "DOID",
        }

        with patch.object(descriptions, "query_condition", return_value=resolved) as query_condition:
            descriptions._standardize_extracted_diseases(docs)

        query_condition.assert_called_once()
        self.assertEqual(
            [condition["identifier"] for condition in docs[0]["healthCondition"]],
            ["9191", "9191"],
        )

    def test_negative_disease_result_removes_uncurated_stub(self):
        docs = [{"_id": "negative", "healthCondition": [{"name": "not a condition"}]}]

        with patch.object(descriptions, "query_condition", return_value=None) as query_condition:
            descriptions._standardize_extracted_diseases(docs)

        query_condition.assert_called_once_with("not a condition")
        self.assertIn("not a condition", descriptions.NEGATIVE_DISEASES)
        self.assertNotIn("healthCondition", docs[0])

    def test_negative_disease_cache_hit_removes_uncurated_stub_without_lookup(self):
        descriptions.NEGATIVE_DISEASES.add("not a condition")
        docs = [{"_id": "negative", "healthCondition": [{"name": "not a condition"}]}]

        with patch.object(descriptions, "query_condition") as query_condition:
            descriptions._standardize_extracted_diseases(docs)

        query_condition.assert_not_called()
        self.assertNotIn("healthCondition", docs[0])

    def test_safe_disease_acronym_is_still_standardized(self):
        docs = [{"_id": "virology", "healthCondition": [{"name": "HIV"}]}]
        resolved = {
            "name": "human immunodeficiency virus infectious disease",
            "alternateName": ["HIV"],
            "identifier": "DOID:526",
        }

        with patch.object(descriptions, "query_condition", return_value=resolved) as query_condition:
            descriptions._standardize_extracted_diseases(docs)

        query_condition.assert_called_once_with("HIV")
        self.assertEqual(docs[0]["healthCondition"][0]["identifier"], "DOID:526")


if __name__ == "__main__":
    unittest.main()
