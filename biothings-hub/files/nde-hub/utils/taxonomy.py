"""Shared taxonomy classification helpers."""

_HOST_GROUPS = frozenset({"Deuterostomia", "Embryophyta", "Arthropoda", "Archaea", "Mollusca"})
_PARASITIC_PLANTS = frozenset({"Arceuthobium", "Cuscuta", "Orobanche", "Striga", "Phoradendron"})


def classify_from_lineage(scientific_name, lineage):
    """Classify a taxon as ``host`` or ``infectiousAgent`` from its UniProt lineage."""
    if scientific_name in _HOST_GROUPS:
        return "host"

    scientific_names = {item["scientificName"] for item in lineage}
    if "Viruses" in scientific_names:
        return "infectiousAgent"
    if scientific_names & {"Archaea", "Mollusca", "Deuterostomia"}:
        return "host"
    if "Embryophyta" in scientific_names and not scientific_names & _PARASITIC_PLANTS:
        return "host"
    if "Arthropoda" in scientific_names:
        if "Acari" in scientific_names and not scientific_names & {"Ixodida", "Ixodes"}:
            return "infectiousAgent"
        return "host"
    return "infectiousAgent"
