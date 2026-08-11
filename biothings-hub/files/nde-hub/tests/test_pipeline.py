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
sys.path.insert(0, str(HUB_ROOT))

# Unit tests do not need the hub's filesystem-creating configuration module.
config = types.ModuleType("config")
config.logger = logging.getLogger("test_pipeline")
config.token = ""
sys.modules.setdefault("config", config)

from utils import pipeline


def candidate(**updates):
    doc = {"_id": "record", "@type": "Dataset", "description": "A useful description"}
    doc.update(updates)
    return doc


class DescriptionEligibilityTests(unittest.TestCase):
    def test_resource_types_are_eligible(self):
        for record_type in ("Dataset", "DataCollection", "ResourceCatalog"):
            with self.subTest(record_type=record_type):
                self.assertTrue(pipeline._needs_descriptions(candidate(**{"@type": record_type})))

    def test_only_biosample_samples_are_eligible(self):
        self.assertTrue(
            pipeline._needs_descriptions(candidate(**{"@type": "Sample", "additionalType": "BioSample"}))
        )
        self.assertTrue(
            pipeline._needs_descriptions(
                candidate(**{"@type": ["Thing", "Sample"], "additionalType": ["Specimen", "BioSample"]})
            )
        )

        for additional_type in (None, "Specimen", "biosample"):
            with self.subTest(additional_type=additional_type):
                self.assertFalse(
                    pipeline._needs_descriptions(
                        candidate(**{"@type": "Sample", "additionalType": additional_type})
                    )
                )

    def test_other_types_and_complete_records_are_ineligible(self):
        self.assertFalse(pipeline._needs_descriptions(candidate(**{"@type": "ComputationalTool"})))
        self.assertFalse(pipeline._needs_descriptions(candidate(description="")))
        self.assertFalse(
            pipeline._needs_descriptions(
                candidate(
                    species=[{"name": "Homo sapiens"}],
                    healthCondition=[{"name": "malaria"}],
                )
            )
        )

    def test_runner_filters_mixed_batches_and_preserves_order(self):
        dataset = candidate(_id="dataset")
        biosample = candidate(_id="biosample", **{"@type": "Sample", "additionalType": "BioSample"})
        ordinary_sample = candidate(_id="sample", **{"@type": "Sample"})
        tool = candidate(_id="tool", **{"@type": "ComputationalTool"})
        docs = [ordinary_sample, dataset, tool, biosample]

        descriptions = types.ModuleType("utils.descriptions")
        descriptions.augment_from_descriptions = Mock(side_effect=lambda records: list(records))
        with patch.dict(sys.modules, {"utils.descriptions": descriptions}):
            result = pipeline._run_descriptions(docs, "test_source")

        self.assertEqual(result, docs)
        self.assertIs(result[0], ordinary_sample)
        eligible = descriptions.augment_from_descriptions.call_args.args[0]
        self.assertEqual([doc["_id"] for doc in eligible], ["dataset", "biosample"])


if __name__ == "__main__":
    unittest.main()
