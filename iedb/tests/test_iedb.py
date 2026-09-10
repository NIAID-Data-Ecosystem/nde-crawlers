import pathlib
import sys
import tempfile
import unittest
import zipfile
from unittest import mock

FILES_DIR = pathlib.Path(__file__).resolve().parents[1] / "files"
sys.path.insert(0, str(FILES_DIR))

import iedb  # noqa: E402


ORGANISM = {
    "organism_id": "227984",
    "tax_id": "227984",
    "parent_tax_id": "2901879",
    "name": "SARS coronavirus Tor2",
    "eligible": True,
}


def event(tab, identifier, name, *, assay_type_id=None, conditions=None):
    return {
        "organism_id": ORGANISM["organism_id"],
        "organism": ORGANISM,
        "tab": tab,
        "member_key": identifier,
        "member": {
            "@type": "ScholarlyArticle" if tab == "references" else "CreativeWork",
            "identifier": identifier,
            "name": name,
            "url": f"https://example.org/{identifier}",
        },
        "date_modified": "2026-08-30-07:00",
        "assay_type_id": assay_type_id,
        "health_conditions": conditions or [],
    }


RAW_EVENTS = [
    event("epitopes", "IEDB_EPITOPE:68", "Epitope 68"),
    event("antigens", "taxon_protein:10002316-other", "Spike glycoprotein"),
    event("assays", "IEDB_ASSAY:7127", "T-cell assay 7127", assay_type_id="34", conditions=["DOID:2945"]),
    event("assays", "IEDB_ASSAY:972", "B-cell assay 972", assay_type_id="34", conditions=["DOID:2945"]),
    event("receptors", "IEDB_RECEPTOR:183965", "BCR 183965"),
    event("references", "IEDB_REFERENCE:406", "SARS coronavirus epitope study"),
]


SYNTHETIC_XML = b"""<?xml version="1.0" encoding="UTF-8"?>
<References xmlns="http://www.iedb.org/schema/CurationSchema">
  <Reference>
    <ReferenceId>406</ReferenceId>
    <DateLastUpdated>2026-08-30-07:00</DateLastUpdated>
    <Article>
      <PubmedId>15356154</PubmedId>
      <ArticleTitle>SARS coronavirus epitope study</ArticleTitle>
      <Journal><Title>Journal of Immunology</Title></Journal>
    </Article>
    <Epitopes>
      <Epitope>
        <EpitopeName>PYRVVVLSF</EpitopeName>
        <EpitopeStructure>
          <FragmentOfANaturalSequenceMolecule>
            <LinearSequence>PYRVVVLSF</LinearSequence>
            <SourceOrganismId>227984</SourceOrganismId>
          </FragmentOfANaturalSequenceMolecule>
        </EpitopeStructure>
        <EpitopeId>68</EpitopeId>
        <Assays>
          <TCell>
            <TCellId>7127</TCellId>
            <Immunization>
              <HostOrganism><OrganismId>9606</OrganismId></HostOrganism>
              <FirstInVivoProcess><DiseaseState>DOID:2945</DiseaseState></FirstInVivoProcess>
              <ImmunizationComments>invalid\x01control byte from the real export</ImmunizationComments>
            </Immunization>
            <AssayInformation><AssayTypeId>34</AssayTypeId></AssayInformation>
          </TCell>
        </Assays>
      </Epitope>
    </Epitopes>
  </Reference>
</References>
"""


class IEDBXMLIteratorTest(unittest.TestCase):
    def test_streaming_iterator_uses_source_organism_not_host(self):
        with tempfile.TemporaryDirectory() as temp_name:
            archive_path = pathlib.Path(temp_name) / "iedb_export.zip"
            with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as archive:
                archive.writestr("406.xml", SYNTHETIC_XML)

            records = list(iedb.iter_iedb_records(archive_path))

        self.assertEqual({record["organism_id"] for record in records}, {"227984"})
        self.assertNotIn("9606", {record["organism_id"] for record in records})
        self.assertEqual([record["tab"] for record in records], ["assays", "epitopes", "references"])
        assay = records[0]
        self.assertEqual(assay["member_key"], "IEDB_ASSAY:7127")
        self.assertEqual(assay["assay_type_id"], "34")
        self.assertEqual(assay["health_conditions"], ["DOID:2945"])


class IEDBGroupingTest(unittest.TestCase):
    def test_rank_filter_retains_species_and_descendants_only(self):
        with tempfile.TemporaryDirectory() as temp_name:
            store = iedb._Store(pathlib.Path(temp_name) / "iedb.sqlite")
            try:
                store.add_organism("genus", "100", "1", "Example genus")
                store.add_organism("species", "101", "100", "Example species")
                store.add_organism("strain", "102", "101", "Example strain")
                store.connection.executemany(
                    "INSERT INTO ncbi_taxonomy (tax_id, parent_tax_id, rank) VALUES (?, ?, ?)",
                    [
                        ("1", "1", "no rank"),
                        ("100", "1", "genus"),
                        ("101", "100", "species"),
                        ("102", "101", "strain"),
                    ],
                )
                iedb._mark_eligible_organisms(store)
                eligibility = dict(store.connection.execute("SELECT organism_id, eligible FROM organisms"))
            finally:
                store.close()

        self.assertEqual(eligibility, {"genus": 0, "species": 1, "strain": 1})

    def test_query_supplements_use_source_organism_iris(self):
        def query_rows(_session, endpoint, _select, _order):
            if endpoint == "antigen_search":
                yield {
                    "parent_source_antigen_id": "https://ontology.iedb.org/taxon-protein/10002316-other",
                    "parent_source_antigen_iri": "taxon_protein:10002316-other",
                    "parent_source_antigen_names": ["Spike glycoprotein"],
                    "source_organism_iris": ["NCBITaxon:227984"],
                }
            elif endpoint == "bcr_search":
                yield {
                    "receptor_group_id": 183965,
                    "receptor_group_iri": "IEDB_RECEPTOR:183965",
                    "receptor_type": "heavylight",
                    "receptor_names": ["Example BCR"],
                    "source_organism_iris": ["NCBITaxon:227984"],
                }

        with tempfile.TemporaryDirectory() as temp_name:
            store = iedb._Store(pathlib.Path(temp_name) / "iedb.sqlite")
            try:
                store.add_organism("227984", "227984", "2901879", "SARS coronavirus Tor2", True)
                with mock.patch.object(iedb, "_iter_query_api", new=query_rows):
                    iedb._supplement_antigens_and_receptors(store, object())
                counts = dict(store.connection.execute("SELECT tab, COUNT(*) FROM members GROUP BY tab ORDER BY tab"))
            finally:
                store.close()

        self.assertEqual(counts, {"antigens": 1, "receptors": 1})


class IEDBMappingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.docs = list(iedb.parse(RAW_EVENTS, release_date="2026-08-30", has_part_limit=10))
        cls.by_tab = {doc["_id"].split("_")[1]: doc for doc in cls.docs}

    def test_one_record_per_non_empty_tab(self):
        self.assertEqual(
            [doc["_id"] for doc in self.docs],
            [
                "iedb_epitopes_227984",
                "iedb_antigens_227984",
                "iedb_assays_227984",
                "iedb_receptors_227984",
                "iedb_references_227984",
            ],
        )

    def test_collection_size_units_match_ui_tabs(self):
        expected = {
            "epitopes": (1, "epitopes"),
            "antigens": (1, "antigens"),
            "assays": (2, "assays"),
            "receptors": (1, "receptors"),
            "references": (1, "references"),
        }
        for tab, (value, unit) in expected.items():
            self.assertEqual(self.by_tab[tab]["collectionSize"]["value"], value)
            self.assertEqual(self.by_tab[tab]["collectionSize"]["unitText"], unit)

    def test_assay_metadata_is_tab_specific(self):
        assay = self.by_tab["assays"]
        self.assertEqual(assay["healthCondition"][0]["identifier"], "DOID:2945")
        self.assertEqual(assay["measurementTechnique"][0]["identifier"], "OBI:1110128")
        self.assertNotIn("measurementTechnique", self.by_tab["references"])
        self.assertNotIn("healthCondition", self.by_tab["references"])

    def test_first_pass_omissions(self):
        for doc in self.docs:
            self.assertNotIn("infectiousAgent", doc)
            self.assertNotIn("sample", doc)
            self.assertNotIn("distribution", doc)
            self.assertNotIn("version", doc)
            self.assertEqual(doc["conditionsOfAccess"], "Open")

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

        for doc in self.docs:
            walk(doc)
        self.assertEqual(missing, [])

    def test_has_part_is_omitted_not_truncated_above_limit(self):
        records = [
            event("epitopes", "IEDB_EPITOPE:1", "Epitope 1"),
            event("epitopes", "IEDB_EPITOPE:2", "Epitope 2"),
        ]
        doc = next(iedb.parse(records, release_date="2026-08-30", has_part_limit=1))
        self.assertEqual(doc["collectionSize"]["value"], 2)
        self.assertNotIn("hasPart", doc)


if __name__ == "__main__":
    unittest.main()
