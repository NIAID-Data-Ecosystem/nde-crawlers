import pathlib
import sys
import unittest

FILES_DIR = pathlib.Path(__file__).resolve().parents[1] / "files"
sys.path.insert(0, str(FILES_DIR))

import pgs_catalog  # noqa: E402


RAW_RECORD = {
    "info": {
        "latest_release": {"date": "2026-07-29"},
        "rest_api": {"version": "1.8.6"},
        "terms_of_use": "https://www.ebi.ac.uk/about/terms-of-use",
        "ensembl_version": 105,
        "citation": {
            "title": "Enhancing the Polygenic Score Catalog",
            "doi": "10.1038/s41588-024-01937-x",
            "PMID": 39327485,
            "journal": "Nature Genetics",
        },
    },
    "trait": {
        "id": "MONDO_0005041",
        "label": "glaucoma",
        "description": "Source trait description.",
        "url": "http://purl.obolibrary.org/obo/MONDO_0005041",
        "trait_categories": ["Other trait"],
        "trait_synonyms": ["glaucoma", "glaucoma (disease)"],
        "trait_mapped_terms": ["DOID:1686", "NCIT:C26782"],
        "associated_pgs_ids": ["PGS000137"],
        "child_associated_pgs_ids": ["PGS000350"],
    },
    "scores": [
        {
            "id": "PGS000137",
            "name": "MTAG_glaucoma",
            "ftp_scoring_file": "https://example.org/PGS000137.txt.gz",
            "date_release": "2020-03-27",
            "trait_reported": "Glaucoma",
            "method_name": "Clumping and Thresholding (C+T)",
            "method_params": "r2 = 0.1; p < 0.001",
            "variants_number": 2673,
            "license": "Score-specific terms",
            "samples_training": [
                {
                    "sample_number": 4672,
                    "sample_cases": 1734,
                    "sample_controls": 2938,
                    "phenotyping_free": "Advanced glaucoma",
                    "ancestry_country": "Australia, New Zealand",
                }
            ],
        }
    ],
    "performance": [
        {
            "id": "PPM000421",
            "associated_pgs_id": "PGS000137",
            "phenotyping_reported": "Primary open-angle glaucoma (POAG)",
            "sampleset": {
                "id": "PSS000245",
                "samples": [
                    {
                        "sample_number": 1795,
                        "sample_cases": 74,
                        "sample_controls": 1721,
                        "sample_percent_male": 43.0,
                        "sample_age": {
                            "estimate_type": "mean",
                            "estimate": 64.02,
                            "variability_type": "sd",
                            "variability": 8.24,
                            "unit": "years",
                        },
                        "ancestry_broad": "European",
                        "ancestry_country": "Austalia",
                    }
                ],
            },
        }
    ],
}


class PGSCatalogMappingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = next(pgs_catalog.parse([RAW_RECORD]))

    def test_reviewed_mapping_decisions(self):
        self.assertEqual(self.doc["_id"], "pgs_catalog_MONDO_0005041")
        self.assertEqual(self.doc["identifier"], "MONDO_0005041")
        self.assertNotIn("sameAs", self.doc["healthCondition"])
        self.assertNotIn("analyticalMethod", self.doc)
        self.assertNotIn("license", self.doc)
        self.assertNotIn("version", self.doc)
        self.assertIn("Other trait", self.doc["keywords"])
        self.assertEqual(
            [term["identifier"] for term in self.doc["topicCategory"]],
            ["topic_3053", "topic_0622", "topic_0625"],
        )

    def test_direct_and_child_scores_are_members(self):
        self.assertEqual(self.doc["collectionSize"]["value"], 2)
        self.assertEqual(
            [member["identifier"] for member in self.doc["hasPart"]],
            ["PGS000137", "PGS000350"],
        )

    def test_evaluation_sample_preserves_source_values(self):
        sample_set = next(item for item in self.doc["sample"] if item["identifier"] == "PSS000245")
        properties = sample_set["aggregateElement"]["additionalProperty"]
        country = next(item for item in properties if item["name"] == "Country of recruitment")
        self.assertEqual(country["value"], "Austalia")

    def test_every_output_object_is_typed(self):
        missing = []

        def walk(value, path="doc"):
            if isinstance(value, list):
                for index, item in enumerate(value):
                    walk(item, f"{path}[{index}]")
            elif isinstance(value, dict):
                if "@type" not in value:
                    missing.append(path)
                for key, item in value.items():
                    walk(item, f"{path}.{key}")

        walk(self.doc)
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
