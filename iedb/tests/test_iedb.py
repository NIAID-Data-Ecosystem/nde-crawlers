import io
import pathlib
import sys
import tarfile
import tempfile
import unittest
import zipfile
from unittest import mock

FILES_DIR = pathlib.Path(__file__).resolve().parents[1] / "files"
sys.path.insert(0, str(FILES_DIR))

import iedb  # noqa: E402

# (organism_id, parent_id, NCBI rank, name): IEDB tree parents with NCBI ranks; IEDB-only organisms have no rank
TAXA = [
    ("1", "1", "no rank", "root"),
    ("1763", "1", "genus", "Mycobacterium"),
    ("77643", "1763", "species group", "Mycobacterium tuberculosis complex"),
    ("1773", "77643", "species", "Mycobacterium tuberculosis"),
    ("83332", "1773", "strain", "Mycobacterium tuberculosis H37Rv"),
    ("1765", "1773", "biotype", "Mycobacterium tuberculosis variant bovis"),
    ("10000323", "1765", None, "Mycobacterium bovis T/91/1378"),
    ("10000324", "10000323", None, "Mycobacterium bovis T/91/1378 sub-isolate"),
    ("1769", "1763", "species", "Mycobacterium leprae"),
    # NCBI puts SARS-CoV-1 and SARS-CoV-2 under the species Betacoronavirus pandemicum; IEDB's search tree does not
    ("10239", "1", "superkingdom", "Viruses"),
    ("694002", "10239", "genus", "Betacoronavirus"),
    ("2509511", "694002", "subgenus", "Sarbecovirus"),
    ("3418604", "2509511", "species", "Betacoronavirus pandemicum"),
    ("2697049", "3418604", "no rank", "Severe acute respiratory syndrome coronavirus 2"),
    ("9999999", "2697049", "no rank", "SARS-CoV-2 variant newer than IEDB's search data"),
    ("694009", "3418604", "no rank", "Severe acute respiratory syndrome-related coronavirus"),
    ("10002316", "694009", None, "SARS-CoV1"),
    ("2901879", "694009", "no rank", "Severe acute respiratory syndrome coronavirus"),
    ("227984", "2901879", "no rank", "SARS coronavirus Tor2"),
    # Shallower than SARS-CoV1 in the organism list, but under it in IEDB's search tree
    ("2708335", "3418604", "no rank", "Pangolin coronavirus"),
    ("10002040", "10239", None, "RNA virus"),
    # A species whose strains are only ever curated together, so their lineages contain each other
    ("583", "1", "genus", "Proteus"),
    ("102862", "583", "species", "Proteus penneri"),
    ("10001836", "102862", None, "Proteus penneri 14"),
    ("10001837", "102862", None, "Proteus penneri 15 (O52)"),
    # A species only ever curated alongside another one
    ("8638", "1", "genus", "Naja"),
    ("8649", "8638", "species", "Naja kaouthia"),
    ("35670", "8638", "species", "Naja naja"),
]

# Lineages in IEDB's search tree, as source_organism_iri_search lists them (sets, so in no particular order)
LINEAGES = {
    "83332": frozenset({"1", "1763", "77643", "1773", "83332"}),
    "1763": frozenset({"1", "1763"}),
    "2697049": frozenset({"1", "10239", "10002040", "694002", "2509511", "2697049"}),
    "10002316": frozenset({"1", "10239", "10002040", "694002", "2509511", "10002316"}),
    "227984": frozenset({"1", "10239", "10002040", "694002", "2509511", "10002316", "227984"}),
    "2708335": frozenset({"1", "10239", "10002040", "694002", "2509511", "10002316", "2708335"}),
    "8649": frozenset({"1", "8638", "8649"}),
}
# Lineages of organisms never curated alone: intersections that also hold the organisms they appear with
SHARED_LINEAGES = {
    "102862": frozenset({"1", "583", "102862", "10001836", "10001837"}),
    "10001836": frozenset({"1", "583", "102862", "10001836", "10001837"}),
    "10001837": frozenset({"1", "583", "102862", "10001836", "10001837"}),
    "35670": frozenset({"1", "8638", "8649", "35670"}),
}


def event(organism_id, tab, identifier, name, *, assay_type_id=None, conditions=None):
    return {
        "organism_id": organism_id,
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
    # The same epitope curated from the species and from a strain belongs to the species collection once
    event("1773", "epitopes", "IEDB_EPITOPE:68", "Epitope 68"),
    event("83332", "epitopes", "IEDB_EPITOPE:68", "Epitope 68"),
    event("83332", "epitopes", "IEDB_EPITOPE:69", "Epitope 69"),
    event("10000324", "epitopes", "IEDB_EPITOPE:70", "Epitope 70"),
    # Curated above species rank, so no collection
    event("1763", "epitopes", "IEDB_EPITOPE:71", "Epitope 71"),
    event("83332", "antigens", "taxon_protein:10002316-other", "Antigen"),
    event("83332", "assays", "IEDB_ASSAY:7127", "T-cell assay 7127", assay_type_id="34", conditions=["DOID:399"]),
    event("1765", "assays", "IEDB_ASSAY:972", "B-cell assay 972", assay_type_id="34", conditions=["DOID:399"]),
    event("1773", "receptors", "IEDB_RECEPTOR:183965", "BCR 183965"),
    event("83332", "references", "IEDB_REFERENCE:406", "Tuberculosis epitope study"),
    event("1769", "epitopes", "IEDB_EPITOPE:72", "Epitope 72"),
    event("2697049", "epitopes", "IEDB_EPITOPE:80", "Epitope 80"),
    event("9999999", "epitopes", "IEDB_EPITOPE:81", "Epitope 81"),
    event("227984", "epitopes", "IEDB_EPITOPE:82", "Epitope 82"),
    event("2708335", "epitopes", "IEDB_EPITOPE:83", "Epitope 83"),
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


ORGANISM_LIST_XML = b"""<?xml version="1.0" encoding="UTF-8"?>
<OrganismList xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <Organism>
        <OrganismId>1773</OrganismId>
        <TaxId>1773</TaxId>
        <ParentTaxId>77643</ParentTaxId>
        <OrganismName>Mycobacterium tuberculosis</OrganismName>
    </Organism>
    <Organism>
        <OrganismId>10000323</OrganismId>
        <ParentTaxId>1765</ParentTaxId>
        <OrganismName>Mycobacterium bovis T/91/1378</OrganismName>
    </Organism>
</OrganismList>
"""


def seeded_store(temp_name):
    store = iedb._Store(pathlib.Path(temp_name) / "iedb.sqlite")
    store.add_taxa(TAXA)
    return store


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


class IEDBSpeciesGroupingTest(unittest.TestCase):
    def test_organisms_roll_up_to_their_ncbi_species(self):
        with tempfile.TemporaryDirectory() as temp_name:
            store = seeded_store(temp_name)
            try:
                resolver = iedb._GroupResolver(store.connection)
                species = {
                    organism_id: resolver.species_id(organism_id)
                    for organism_id in ("1773", "83332", "1765", "10000323", "10000324", "1769", "1763", "77643")
                }
                species["missing"] = resolver.species_id("424242")
                # Without IEDB search lineages, the group is the NCBI species
                groups = {organism_id: resolver.group_id(organism_id) for organism_id in ("83332", "227984", "1763")}
            finally:
                store.close()

        self.assertEqual(
            species,
            {
                "1773": "1773",
                "83332": "1773",
                "1765": "1773",
                # IEDB-only organisms reach the species through their IEDB-only parents
                "10000323": "1773",
                "10000324": "1773",
                "1769": "1769",
                "1763": None,
                "77643": None,
                "missing": None,
            },
        )
        self.assertEqual(groups, {"83332": "1773", "227984": "3418604", "1763": None})

    def test_groups_follow_iedb_search_tree(self):
        with tempfile.TemporaryDirectory() as temp_name:
            store = seeded_store(temp_name)
            try:
                resolver = iedb._GroupResolver(store.connection, LINEAGES, {"11623": "2697049"}, SHARED_LINEAGES)
                groups = {
                    organism_id: resolver.group_id(organism_id)
                    for organism_id in (
                        "83332",
                        "1763",
                        "2697049",
                        "10002316",
                        "227984",
                        "2708335",
                        "9999999",
                        "11623",
                        "10001837",
                        "35670",
                    )
                }
            finally:
                store.close()

        self.assertEqual(
            groups,
            {
                # IEDB's tree keeps the NCBI species
                "83332": "1773",
                "1763": None,
                # IEDB's tree skips Betacoronavirus pandemicum, so SARS-CoV-2 and IEDB's SARS-CoV1 are the groups
                "2697049": "2697049",
                "10002316": "10002316",
                "227984": "10002316",
                # IEDB's tree, not the organism list's depth, decides which taxon is higher
                "2708335": "10002316",
                # Not in IEDB's search data: grouped with its nearest ancestor that is
                "9999999": "2697049",
                # A retired NCBI tax ID is read as its replacement
                "11623": "2697049",
                # Organisms that only share epitopes are not each other's parents
                "10001837": "102862",
                "35670": "35670",
            },
        )
        self.assertEqual(resolver.fallbacks, 0)

    def test_lineages_are_intersected_across_epitopes_with_keyset_paging(self):
        pages = [
            [
                {
                    "structure_id": 1,
                    "source_organism_iris": ["NCBITaxon:2697049", "taxon:10002316"],
                    "source_organism_iri_search": [
                        "NCBITaxon:1",
                        "NCBITaxon:2509511",
                        "NCBITaxon:2697049",
                        "OBI:0100026",
                        "taxon:10002316",
                    ],
                },
                {
                    "structure_id": 2,
                    "source_organism_iris": ["NCBITaxon:2697049"],
                    "source_organism_iri_search": ["NCBITaxon:1", "NCBITaxon:2509511", "NCBITaxon:2697049"],
                },
            ],
            [
                {
                    "structure_id": 3,
                    "source_organism_iris": ["taxon:10002316"],
                    "source_organism_iri_search": ["NCBITaxon:1", "NCBITaxon:2509511", "taxon:10002316"],
                },
                {
                    "structure_id": 4,
                    "source_organism_iris": ["taxon:10001836", "taxon:10001837"],
                    "source_organism_iri_search": [
                        "NCBITaxon:1",
                        "NCBITaxon:102862",
                        "taxon:10001836",
                        "taxon:10001837",
                    ],
                },
            ],
        ]
        calls = []

        def request_json(_session, url, params=None):
            calls.append(dict(params))
            return pages[len(calls) - 1] if len(calls) <= len(pages) else []

        with mock.patch.object(iedb, "LINEAGE_PAGE_SIZE", 2), mock.patch.object(iedb, "_request_json", request_json):
            exact, shared = iedb._load_iedb_lineages(object())

        # Their own epitopes win over the one they share
        self.assertEqual(
            exact,
            {
                "2697049": frozenset({"1", "2509511", "2697049"}),
                "10002316": frozenset({"1", "2509511", "10002316"}),
            },
        )
        # Never curated alone, so each keeps the strain it always appears with
        self.assertEqual(
            shared,
            {
                "10001836": frozenset({"1", "102862", "10001836", "10001837"}),
                "10001837": frozenset({"1", "102862", "10001836", "10001837"}),
            },
        )
        self.assertNotIn("structure_id", calls[0])
        self.assertEqual([call["structure_id"] for call in calls[1:]], ["gt.2", "gt.4"])

    def test_organism_tree_parents_win_over_ncbi_and_ranks_come_from_ncbi(self):
        with tempfile.TemporaryDirectory() as temp_name:
            organism_list = pathlib.Path(temp_name) / "OrganismList.zip"
            with zipfile.ZipFile(organism_list, "w", zipfile.ZIP_DEFLATED) as archive:
                archive.writestr("OrganismList.xml", ORGANISM_LIST_XML)
            store = iedb._Store(pathlib.Path(temp_name) / "iedb.sqlite")
            try:
                iedb._load_organism_list(store, organism_list)
                # Rows as _load_ncbi_taxonomy upserts them, after the organism tree
                store.connection.executemany(
                    """
                    INSERT INTO taxa (organism_id, parent_id, rank) VALUES (?, ?, ?)
                    ON CONFLICT (organism_id) DO UPDATE SET
                        parent_id = COALESCE(taxa.parent_id, excluded.parent_id),
                        rank = excluded.rank
                    """,
                    [("1773", "1763", "species"), ("1765", "1773", "biotype")],
                )
                rows = {
                    row[0]: row[1:]
                    for row in store.connection.execute("SELECT organism_id, parent_id, rank, name, ncbi FROM taxa")
                }
                species = iedb._GroupResolver(store.connection).species_id("10000323")
            finally:
                store.close()

        self.assertEqual(rows["1773"], ("77643", "species", "Mycobacterium tuberculosis", 1))
        self.assertEqual(rows["10000323"], ("1765", None, "Mycobacterium bovis T/91/1378", 0))
        self.assertEqual(rows["1765"], ("1773", "biotype", None, 0))
        self.assertEqual(species, "1773")

    def test_ncbi_taxdump_loads_ranks_and_retired_ids(self):
        with tempfile.TemporaryDirectory() as temp_name:
            taxdump = pathlib.Path(temp_name) / "taxdump.tar.gz"
            files = {
                "merged.dmp": b"11623\t|\t3052303\t|\n",
                "names.dmp": b"1\t|\troot\t|\t\t|\tscientific name\t|\n",
                "nodes.dmp": b"1\t|\t1\t|\tno rank\t|\t\t|\n3052303\t|\t1\t|\tspecies\t|\t\t|\n",
            }
            with tarfile.open(taxdump, "w:gz") as archive:
                for name, data in files.items():
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    archive.addfile(info, io.BytesIO(data))
            store = iedb._Store(pathlib.Path(temp_name) / "iedb.sqlite")
            try:
                merged = iedb._load_ncbi_taxonomy(store, taxdump)
                ranks = dict(store.connection.execute("SELECT organism_id, rank FROM taxa"))
                ncbi = dict(store.connection.execute("SELECT organism_id, ncbi FROM taxa"))
            finally:
                store.close()

        self.assertEqual(merged, {"11623": "3052303"})
        self.assertEqual(ranks, {"1": "no rank", "3052303": "species"})
        self.assertEqual(ncbi, {"1": 1, "3052303": 1})

    def test_query_supplements_read_ncbi_and_iedb_only_iris(self):
        def query_rows(_session, endpoint, _select, _order):
            if endpoint == "antigen_search":
                yield {
                    "parent_source_antigen_id": "https://ontology.iedb.org/taxon-protein/10002316-other",
                    "parent_source_antigen_iri": "taxon_protein:10002316-other",
                    "parent_source_antigen_names": ["Antigen 85B"],
                    "source_organism_iris": ["NCBITaxon:83332", "taxon:10000323", "NCBITaxon:1763"],
                }
            elif endpoint == "bcr_search":
                yield {
                    "receptor_group_id": 183965,
                    "receptor_group_iri": "IEDB_RECEPTOR:183965",
                    "receptor_type": "heavylight",
                    "receptor_names": ["Example BCR"],
                    "source_organism_iris": ["taxon:10000324", "NCBITaxon:1769"],
                }

        with tempfile.TemporaryDirectory() as temp_name:
            store = seeded_store(temp_name)
            try:
                resolver = iedb._GroupResolver(store.connection)
                with mock.patch.object(iedb, "_iter_query_api", new=query_rows):
                    iedb._supplement_antigens_and_receptors(store, object(), resolver)
                members = sorted(store.connection.execute("SELECT species_id, tab, member_key FROM members"))
                antigen_organisms = dict(store.organism_counts[("1773", "antigens")])
            finally:
                store.close()

        self.assertEqual(
            members,
            [
                ("1769", "receptors", "IEDB_RECEPTOR:183965"),
                ("1773", "antigens", "taxon_protein:10002316-other"),
                ("1773", "receptors", "IEDB_RECEPTOR:183965"),
            ],
        )
        self.assertEqual(antigen_organisms, {"83332": 1, "10000323": 1})

    def test_species_selection_limits_ingestion(self):
        with tempfile.TemporaryDirectory() as temp_name:
            store = seeded_store(temp_name)
            try:
                iedb._ingest_events(store, RAW_EVENTS, iedb._GroupResolver(store.connection), {"1769"})
                species = {row[0] for row in store.connection.execute("SELECT species_id FROM members")}
            finally:
                store.close()

        self.assertEqual(species, {"1769"})


class IEDBMappingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.docs = list(
            iedb.parse(
                RAW_EVENTS,
                taxa=TAXA,
                lineages=LINEAGES,
                shared_lineages=SHARED_LINEAGES,
                release_date="2026-08-30",
                has_part_limit=10,
            )
        )
        cls.by_id = {doc["_id"]: doc for doc in cls.docs}
        cls.by_tab = {doc["_id"].split("_")[1]: doc for doc in cls.docs if doc["_id"].endswith("_1773")}

    def test_one_record_per_non_empty_tab_of_each_species(self):
        self.assertEqual(
            [doc["_id"] for doc in self.docs],
            [
                "iedb_epitopes_1769",
                "iedb_epitopes_1773",
                "iedb_antigens_1773",
                "iedb_assays_1773",
                "iedb_receptors_1773",
                "iedb_references_1773",
                "iedb_epitopes_2697049",
                "iedb_epitopes_10002316",
            ],
        )
        self.assertEqual(self.by_tab["epitopes"]["identifier"], "IEDB:EPITOPES:1773")

    def test_groups_where_iedb_skips_the_ncbi_species(self):
        sars2 = self.by_id["iedb_epitopes_2697049"]
        self.assertEqual(sars2["name"], "IEDB epitopes from Severe acute respiratory syndrome coronavirus 2")
        self.assertEqual(sars2["collectionSize"]["value"], 2)
        self.assertEqual(sars2["url"], f"{iedb.SEARCH_BASE}/query/epitopes?epitope_source_organism=NCBITaxon:2697049")
        self.assertEqual(sars2["species"]["identifier"], "2697049")
        self.assertIn("include SARS-CoV-2 variant newer than IEDB's search data.", sars2["description"])

        # IEDB's searches do not resolve IEDB-only taxa, so the record links to the Query API instead
        sars1 = self.by_id["iedb_epitopes_10002316"]
        self.assertEqual(sars1["collectionSize"]["value"], 2)
        query = f"{iedb.QUERY_API_BASE}/epitope_search?source_organism_iri_search=cs.%7Btaxon:10002316%7D"
        self.assertEqual(sars1["url"], query)
        self.assertEqual(sars1["includedInDataCatalog"]["archivedAt"], query)
        self.assertEqual(sars1["species"], {"@type": "DefinedTerm", "name": "SARS-CoV1"})
        self.assertIn("IEDB epitopes whose source organism is SARS-CoV1 (taxon:10002316)", sars1["description"])
        self.assertIn(query, [work.get("url") for work in sars1["isBasedOn"]])

    def test_members_are_deduplicated_across_the_species_subtree(self):
        epitopes = self.by_tab["epitopes"]
        self.assertEqual(epitopes["collectionSize"]["value"], 3)
        self.assertEqual(
            [member["identifier"] for member in epitopes["hasPart"]],
            ["IEDB_EPITOPE:68", "IEDB_EPITOPE:69", "IEDB_EPITOPE:70"],
        )

    def test_collection_size_units_match_ui_tabs(self):
        expected = {
            "epitopes": (3, "epitopes"),
            "antigens": (1, "antigens"),
            "assays": (2, "assays"),
            "receptors": (1, "receptors"),
            "references": (1, "references"),
        }
        for tab, (value, unit) in expected.items():
            self.assertEqual(self.by_tab[tab]["collectionSize"]["value"], value)
            self.assertEqual(self.by_tab[tab]["collectionSize"]["unitText"], unit)

    def test_urls_use_iedb_species_searches(self):
        for tab, doc in self.by_tab.items():
            expected = f"{iedb.SEARCH_BASE}/query/{tab}?epitope_source_organism=NCBITaxon:1773"
            self.assertEqual(doc["url"], expected)
            self.assertEqual(doc["includedInDataCatalog"]["archivedAt"], expected)
            self.assertEqual(doc["exampleOfWork"]["potentialAction"]["url"], expected)

        assay_queries = [
            work["url"] for work in self.by_tab["assays"]["isBasedOn"] if "Query API" in work.get("name", "")
        ]
        self.assertEqual(
            assay_queries,
            [
                f"{iedb.QUERY_API_BASE}/{endpoint}?source_organism_iri_search=cs.%7BNCBITaxon:1773%7D"
                for endpoint in ("tcell_search", "bcell_search", "mhc_search")
            ],
        )

    def test_species_term_is_left_for_the_hub_to_standardize(self):
        self.assertEqual(
            self.by_tab["epitopes"]["species"],
            {
                "@type": "DefinedTerm",
                "name": "Mycobacterium tuberculosis",
                "identifier": "1773",
                "url": "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1773",
            },
        )

    def test_organisms_curated_below_the_species_are_named(self):
        epitopes = self.by_tab["epitopes"]
        self.assertEqual(epitopes["name"], "IEDB epitopes from Mycobacterium tuberculosis")
        self.assertIn("(NCBITaxon:1773)", epitopes["description"])
        self.assertIn(
            "Organisms curated below it include Mycobacterium tuberculosis H37Rv and "
            "Mycobacterium bovis T/91/1378 sub-isolate.",
            epitopes["description"],
        )
        self.assertEqual(
            epitopes["keywords"],
            [
                "IEDB",
                "epitopes",
                "Mycobacterium tuberculosis",
                "Mycobacterium tuberculosis H37Rv",
                "Mycobacterium bovis T/91/1378 sub-isolate",
            ],
        )
        # Only the species itself was curated for receptors, so there is nothing to name
        self.assertNotIn("curated below", self.by_tab["receptors"]["description"])
        self.assertEqual(self.by_tab["receptors"]["keywords"], ["IEDB", "receptors", "Mycobacterium tuberculosis"])

    def test_assay_metadata_is_tab_specific(self):
        assay = self.by_tab["assays"]
        self.assertEqual(assay["healthCondition"][0]["identifier"], "DOID:399")
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
            event("1773", "epitopes", "IEDB_EPITOPE:1", "Epitope 1"),
            event("83332", "epitopes", "IEDB_EPITOPE:2", "Epitope 2"),
        ]
        doc = next(iedb.parse(records, taxa=TAXA, release_date="2026-08-30", has_part_limit=1))
        self.assertEqual(doc["collectionSize"]["value"], 2)
        self.assertNotIn("hasPart", doc)


if __name__ == "__main__":
    unittest.main()
