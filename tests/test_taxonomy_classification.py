import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load_taxonomy():
    path = ROOT / "biothings-hub/files/nde-hub/utils/taxonomy.py"
    spec = importlib.util.spec_from_file_location("taxonomy_for_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _lineage(*names):
    return [{"scientificName": name} for name in names]


@pytest.mark.parametrize(
    ("scientific_name", "lineage", "expected"),
    [
        ("Archaea", [], "host"),
        ("Mollusca", [], "host"),
        ("Test mite", _lineage("Arthropoda", "Acari"), "infectiousAgent"),
        ("Test tick", _lineage("Arthropoda", "Acari", "Ixodida"), "host"),
        ("Ixodes example", _lineage("Arthropoda", "Acari", "Ixodes"), "host"),
        ("Test virus", _lineage("Viruses"), "infectiousAgent"),
        ("Test plant", _lineage("Embryophyta"), "host"),
        ("Test parasite", _lineage("Embryophyta", "Cuscuta"), "infectiousAgent"),
        ("Test vertebrate", _lineage("Deuterostomia"), "host"),
        ("Unknown taxon", _lineage("Cellular organisms"), "infectiousAgent"),
    ],
)
def test_classify_from_lineage(scientific_name, lineage, expected):
    taxonomy = _load_taxonomy()

    assert taxonomy.classify_from_lineage(scientific_name, lineage) == expected
