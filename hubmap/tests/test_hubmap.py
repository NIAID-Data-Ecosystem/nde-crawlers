import copy
import pathlib
import sys
import unittest
from unittest import mock

FILES_DIR = pathlib.Path(__file__).resolve().parents[1] / "files"
sys.path.insert(0, str(FILES_DIR))

import hubmap  # noqa: E402

DONOR = {
    "mapped_metadata": {
        "sex": ["Female"],
        "age_value": [30],
        "age_unit": ["years"],
        "race": ["Unknown"],
        "abo_blood_group_system": ["Blood Type O"],
    }
}
ORGANS = ["Kidney (Left)"]


def source_sample(uuid, hubmap_id, category):
    return {"uuid": uuid, "hubmap_id": hubmap_id, "entity_type": "Sample", "mapped_sample_category": category}


SECTION = source_sample("s3", "HBM333.CCCC.333", "Section")
SECTION_2 = source_sample("s4", "HBM444.DDDD.444", "Section")

# Organ -> Block -> Section -> raw Dataset -> processed Dataset, plus a dataset assayed on two sections.
DATASETS = [
    {
        "uuid": "d1",
        "hubmap_id": "HBM111.AAAA.111",
        "source_samples": [SECTION],
        "origin_samples_unique_mapped_organs": ORGANS,
        "donor": DONOR,
    },
    {
        "uuid": "d2",
        "hubmap_id": "HBM222.BBBB.222",
        "immediate_ancestors": [{"uuid": "d1", "entity_type": "Dataset"}],
        "source_samples": [SECTION],
        "origin_samples_unique_mapped_organs": ORGANS,
        "donor": DONOR,
    },
    {
        "uuid": "d3",
        "hubmap_id": "HBM555.EEEE.555",
        "source_samples": [SECTION, SECTION_2, SECTION_2],
        "origin_samples_unique_mapped_organs": ORGANS,
        "donor": DONOR,
    },
]

# Nearest first, as HuBMAP lists them.
UPSTREAM = [
    {"uuid": "s2", "entity_type": "Sample", "sample_category": "block"},
    {"uuid": "s1", "entity_type": "Sample", "sample_category": "organ"},
    {"uuid": "donor-1", "entity_type": "Donor"},
]

SAMPLES = {
    "s1": {
        "uuid": "s1",
        "hubmap_id": "HBM000.ORGN.000",
        "mapped_sample_category": "Organ",
        "description": "Left kidney",
    },
    "s2": {
        "uuid": "s2",
        "hubmap_id": "HBM999.BLCK.999",
        "mapped_sample_category": "Block",
        "description": "Block embedded in OCT.",
        "metadata": {"notes": "Stored at -80C"},
    },
    "s3": {
        "uuid": "s3",
        "hubmap_id": "HBM333.CCCC.333",
        "mapped_sample_category": "Section",
        "description": "Section for imaging",
        "ancestors": UPSTREAM,
    },
    "s4": {
        "uuid": "s4",
        "hubmap_id": "HBM444.DDDD.444",
        "mapped_sample_category": "Section",
        # Same text as its block, so the block's description isn't repeated.
        "description": "Block embedded in OCT",
        "ancestors": UPSTREAM,
    },
}


class HubmapSampleFlatteningTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with mock.patch.multiple(
            hubmap,
            get_ids=lambda: DATASETS,
            get_sample_ids=lambda: list(SAMPLES),
            _fetch_browse_json=lambda _entity_type, uuid, timeout=20: SAMPLES[uuid],
        ):
            docs = list(hubmap.parse())
        cls.datasets = {d["name"]: d for d in docs if d["@type"] == "Dataset"}
        cls.samples = {d["name"]: d for d in docs if d["@type"] == "Sample"}

    def test_dataset_lists_and_counts_source_samples(self):
        sample = self.datasets["HBM111.AAAA.111"]["sample"]
        self.assertEqual(sample["@type"], "SampleCollection")
        self.assertEqual(sample["numberOfItems"], {"@type": "QuantitativeValue", "value": 1, "unitText": "sample"})
        self.assertEqual(
            sample["itemListElement"],
            [
                {
                    "@type": "Sample",
                    "identifier": "HBM333.CCCC.333",
                    "url": "https://portal.hubmapconsortium.org/browse/sample/s3",
                }
            ],
        )
        self.assertEqual(self.datasets["HBM555.EEEE.555"]["sample"]["numberOfItems"]["value"], 2)

    def test_processed_dataset_gets_the_raw_datasets_samples(self):
        items = self.datasets["HBM222.BBBB.222"]["sample"]["itemListElement"]
        self.assertEqual([item["identifier"] for item in items], ["HBM333.CCCC.333"])

    def test_aggregate_element_summarizes_donor_and_organ(self):
        aggregate = self.datasets["HBM111.AAAA.111"]["sample"]["aggregateElement"]
        self.assertEqual(aggregate["additionalType"], "ExperimentalRunSample")
        self.assertEqual(aggregate["sex"], "Female")
        self.assertEqual(aggregate["sampleType"], [{"@type": "DefinedTerm", "name": "Section"}])
        self.assertEqual(aggregate["anatomicalStructure"], [{"@type": "DefinedTerm", "name": "Kidney (Left)"}])
        # "Unknown" race is dropped; the blood type remains.
        self.assertEqual(aggregate["associatedPhenotype"], [{"@type": "DefinedTerm", "name": "Blood Type O"}])

    def test_only_assayed_samples_are_records(self):
        self.assertEqual(sorted(self.samples), ["HBM333.CCCC.333", "HBM444.DDDD.444"])
        self.assertEqual({s["additionalType"] for s in self.samples.values()}, {"ExperimentalRunSample"})

    def test_upstream_sample_text_is_carried_on_the_assayed_sample(self):
        self.assertEqual(
            self.samples["HBM333.CCCC.333"]["description"],
            "Section for imaging. Block: Block embedded in OCT. Block notes: Stored at -80C. Organ: Left kidney",
        )
        self.assertEqual(
            self.samples["HBM444.DDDD.444"]["description"],
            "Block embedded in OCT. Block notes: Stored at -80C. Organ: Left kidney",
        )

    def test_assayed_sample_is_basis_for_every_dataset_it_feeds(self):
        section = self.samples["HBM333.CCCC.333"]
        self.assertEqual(
            [(entry["@type"], entry["identifier"]) for entry in section["isBasisFor"]],
            [("Dataset", "HBM111.AAAA.111"), ("Dataset", "HBM222.BBBB.222"), ("Dataset", "HBM555.EEEE.555")],
        )

    def test_portal_lookup_matches_item_list(self):
        # The portal lists a Dataset's samples with isBasisFor.identifier:"<Dataset identifier>".
        for dataset in self.datasets.values():
            linked = {
                sample["identifier"]
                for sample in self.samples.values()
                if any(entry["identifier"] == dataset["identifier"] for entry in sample.get("isBasisFor", []))
            }
            listed = {item["identifier"] for item in dataset["sample"]["itemListElement"]}
            self.assertEqual(linked, listed, dataset["identifier"])


class HubmapSearchPagingTest(unittest.TestCase):
    def test_pages_past_the_first_response_with_search_after(self):
        pages = [["a", "b"], ["c", "d"], ["e"]]
        payloads = []

        def fake_post(url, json, timeout):
            payloads.append(copy.deepcopy(json))
            hits = [{"_id": uuid, "sort": [uuid]} for uuid in pages[len(payloads) - 1]]
            return mock.Mock(**{"json.return_value": {"hits": {"hits": hits}}})

        with mock.patch.object(hubmap.requests, "post", side_effect=fake_post):
            ids = hubmap._search_entity_ids("Sample", page_size=2)

        self.assertEqual(ids, ["a", "b", "c", "d", "e"])
        self.assertEqual([payload.get("search_after") for payload in payloads], [None, ["b"], ["d"]])


if __name__ == "__main__":
    unittest.main()
