#!/usr/bin/env python3
"""
DBAASP DataCollection crawler for the NIAID Data Ecosystem.

Pulls every peptide from the DBAASP API, groups peptide records by their
target species (from peptide-level targetActivities), and produces one
DataCollection record per target species.

Data source
-----------
https://dbaasp.org/
https://dbaasp.org/v3/api-docs
"""

import copy
import datetime as dt
import logging
import re
import time
import urllib.parse
from typing import Any, Iterator, Optional

import requests

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")

API_BASE = "https://dbaasp.org"
USER_AGENT = "nde-dbaasp-crawler/0.1"
LIST_PAGE_SIZE = 500
REQUEST_TIMEOUT = 60
MAX_RETRIES = 4
RETRY_BACKOFF = 5

SOURCE_METADATA = {
    "name": "DBAASP (as DataCollection records)",
    "url": "https://dbaasp.org/",
    "description": (
        "Database of Antimicrobial Activity and Structure of Peptides "
        "(DBAASP) records grouped by target organism as DataCollections."
    ),
    "license": "https://creativecommons.org/licenses/by-nc/4.0/",
}

ABOUT_DEFINED_TERM = {
    "@type": "DefinedTerm",
    "description": "A subclass of MolecularEntity with the same properties as Protein",
    "displayName": "Peptide",
    "name": "Peptide",
    "url": "http://purl.obolibrary.org/obo/NCIT_C735",
}

INCLUDED_IN_DATA_CATALOG = {
    "@type": "DataCatalog",
    "name": "Database of Antimicrobial Activity and Structure of Peptides",
    "alternateName": "DBAASP",
    "url": "https://dbaasp.org/",
    "archivedAt": "https://data.niaid.nih.gov/resources?id=dde_cc13744ba5c15cca",
}

SPECIES = {"@type": "DefinedTerm", "name": "Homo sapiens"}

# Curated from DDE record dde_cc13744ba5c15cca.
DDE_CURATED_BY = {
    "name": "Data Discovery Engine",
    "url": "https://discovery.biothings.io/",
    "dateModified": "2026-02-28",
}

AUTHORS = [
    {
        "@type": "Person",
        "name": "Malak Pirtskhalava",
        "givenName": "Malak",
        "familyName": "Pirtskhalava",
        "affiliation": {
            "@type": "Organization",
            "name": "Ivane Beritashvili Center of Experimental Biomedicine",
        },
    },
    {
        "@type": "Person",
        "name": "Anthony A Amstrong",
        "givenName": "Anthony",
        "familyName": "Amstrong",
        "affiliation": {
            "@type": "Organization",
            "name": "Office of Cyber Infrastructure and Computational Biology",
            "parentOrganization": "National Institute of Allergy and Infectious Diseases",
        },
    },
    {
        "@type": "Person",
        "name": "Maia Grigolava",
        "givenName": "Maia",
        "familyName": "Grigolava",
        "affiliation": {
            "@type": "Organization",
            "name": "Ivane Beritashvili Center of Experimental Biomedicine",
        },
    },
    {
        "@type": "Person",
        "name": "Mindia Chubinidze",
        "givenName": "Mindia",
        "familyName": "Chubinidze",
        "affiliation": {
            "@type": "Organization",
            "name": "Ivane Beritashvili Center of Experimental Biomedicine",
        },
    },
    {
        "@type": "Person",
        "name": "Evgenia Alimbarashvili",
        "givenName": "Evgenia",
        "familyName": "Alimbarashvili",
        "affiliation": {
            "@type": "Organization",
            "name": "Ivane Beritashvili Center of Experimental Biomedicine",
        },
    },
    {
        "@type": "Person",
        "name": "Boris Vishnepolsky",
        "givenName": "Boris",
        "familyName": "Vishnepolsky",
        "affiliation": {
            "@type": "Organization",
            "name": "Ivane Beritashvili Center of Experimental Biomedicine",
        },
    },
    {
        "@type": "Person",
        "name": "Andrei Gabrielian",
        "givenName": "Andrei",
        "familyName": "Gabrielian",
        "affiliation": {
            "@type": "Organization",
            "name": "Office of Cyber Infrastructure and Computational Biology",
            "parentOrganization": "National Institute of Allergy and Infectious Diseases",
        },
    },
    {
        "@type": "Person",
        "name": "Alex Rosenthal",
        "givenName": "Alex",
        "familyName": "Rosenthal",
        "affiliation": {
            "@type": "Organization",
            "name": "Office of Cyber Infrastructure and Computational Biology",
            "parentOrganization": "National Institute of Allergy and Infectious Diseases",
        },
    },
    {
        "@type": "Person",
        "name": "Darrell E Hurt",
        "givenName": "Darrell",
        "familyName": "Hurt",
        "affiliation": {
            "@type": "Organization",
            "name": "Office of Cyber Infrastructure and Computational Biology",
            "parentOrganization": "National Institute of Allergy and Infectious Diseases",
        },
    },
    {
        "@type": "Person",
        "name": "Michael Tartakovsky",
        "givenName": "Michael",
        "familyName": "Tartakovsky",
        "affiliation": {
            "@type": "Organization",
            "name": "Office of Cyber Infrastructure and Computational Biology",
            "parentOrganization": "National Institute of Allergy and Infectious Diseases",
        },
    },
]

CREATOR = [
    {
        "@type": "Organization",
        "name": "Ivane Beritashvili Center of Experimental Biomedicine",
    },
    {
        "@type": "Organization",
        "name": "Office of Cyber Infrastructure and Computational Biology",
        "parentOrganization": "National Institute of Allergy and Infectious Diseases",
    },
]

CITATION = [
    {
        "@type": "ScholarlyArticle",
        "name": (
            "DBAASP v3: database of antimicrobial/cytotoxic activity and "
            "structure of peptides as a resource for development of new "
            "therapeutics."
        ),
        "author": [
            {"@type": "Person", "name": "Pirtskhalava M"},
            {"@type": "Person", "name": "Amstrong AA"},
            {"@type": "Person", "name": "Grigolava M"},
            {"@type": "Person", "name": "Chubinidze M"},
            {"@type": "Person", "name": "Alimbarashvili E"},
            {"@type": "Person", "name": "Vishnepolsky B"},
            {"@type": "Person", "name": "Gabrielian A"},
            {"@type": "Person", "name": "Rosenthal A"},
            {"@type": "Person", "name": "Hurt DE"},
            {"@type": "Person", "name": "Tartakovsky M"},
        ],
        "datePublished": "2021-01-08",
        "journalName": "Nucleic acids research",
        "identifier": "PMID:33151284",
        "pmid": "33151284",
        "url": "https://pubmed.ncbi.nlm.nih.gov/33151284/",
    },
]

CREDIT_TEXT = (
    "Please cite Malak Pirtskhalava, Anthony A Amstrong, Maia Grigolava, "
    "Mindia Chubinidze, Evgenia Alimbarashvili, Boris Vishnepolsky, "
    "Andrei Gabrielian, Alex Rosenthal, Darrell E Hurt, Michael Tartakovsky, "
    "DBAASP v3: database of antimicrobial/cytotoxic activity and structure "
    "of peptides as a resource for development of new therapeutics, "
    "Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, "
    "Pages D288-D297, https://doi.org/10.1093/nar/gkaa991"
)

FUNDING = [
    {
        "@type": "MonetaryGrant",
        "identifier": "HHSN316201300006W",
        "funder": {
            "@type": "Organization",
            "name": "National Institute of Allergy and Infectious Diseases",
            "alternateName": "NIAID",
            "parentOrganization": "NIH",
        },
    },
    {
        "@type": "MonetaryGrant",
        "identifier": "HHSN27200002",
        "funder": {
            "@type": "Organization",
            "name": "National Institute of Allergy and Infectious Diseases",
            "alternateName": "NIAID",
            "parentOrganization": "NIH",
        },
    },
]

MEASUREMENT_TECHNIQUE = [
    {
        "@type": "DefinedTerm",
        "name": "Minimum Inhibitory Concentration Test",
        "url": "http://purl.obolibrary.org/obo/NCIT_C128985",
        "inDefinedTermSet": "GENEPIO",
        "alternateName": [
            "MIC",
            "Minimum Inhibitory Concentration",
            "Minimum Inhibitory Concentration Test",
        ],
        "isCurated": True,
        "curatedBy": DDE_CURATED_BY,
    },
    {
        "@type": "DefinedTerm",
        "name": "cytotoxicity assay",
        "url": "http://www.bioassayontology.org/bao#BAO_0002993",
        "inDefinedTermSet": "ENM",
        "alternateName": [
            "cytotox assay",
            "cytotoxic activity",
            "in vitro cytotoxic activity",
        ],
        "isCurated": True,
        "curatedBy": DDE_CURATED_BY,
    },
    {
        "@type": "DefinedTerm",
        "name": "3D structure prediction",
        "url": "http://www.bioassayontology.org/bao#BAO_0002216",
        "inDefinedTermSet": "BAO",
        "isCurated": True,
        "curatedBy": DDE_CURATED_BY,
    },
]

VARIABLE_MEASURED = [
    {
        "@type": "DefinedTerm",
        "name": "MIC",
        "url": "http://www.bioassayontology.org/bao#BAO_0002146",
        "inDefinedTermSet": "ENM",
        "alternateName": ["minimum inhibitory concentration"],
        "isCurated": True,
        "curatedBy": DDE_CURATED_BY,
    },
    {
        "@type": "DefinedTerm",
        "name": "Hemolysis",
        "url": "http://purl.obolibrary.org/obo/NCIT_C37965",
        "inDefinedTermSet": "NCIT",
        "alternateName": ["Hemolysis", "Hemolytic"],
        "isCurated": True,
        "curatedBy": DDE_CURATED_BY,
    },
    {
        "@type": "DefinedTerm",
        "name": "cytotoxicity measurement",
        "url": "http://www.ebi.ac.uk/efo/EFO_0006952",
        "inDefinedTermSet": "EFO",
        "isCurated": True,
        "curatedBy": DDE_CURATED_BY,
    },
    {
        "@type": "DefinedTerm",
        "name": "Peptide Conformation",
        "url": "http://purl.obolibrary.org/obo/NCIT_C13407",
        "inDefinedTermSet": "NCIT",
        "alternateName": ["Peptide Conformation"],
        "isCurated": True,
        "curatedBy": DDE_CURATED_BY,
    },
]

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "name": "Chemoinformatics",
        "identifier": "topic_2814",
        "url": "http://edamontology.org/topic_2814",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Medicinal chemistry",
        "identifier": "topic_0209",
        "url": "http://edamontology.org/topic_0209",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Microbiology",
        "identifier": "topic_3301",
        "url": "http://edamontology.org/topic_3301",
        "inDefinedTermSet": "EDAM",
    },
]

USAGE_INFO = {
    "@type": "CreativeWork",
    "name": "Conditions of use",
    "url": "https://dbaasp.org/DBAASP_Terms_And_Conditions.pdf",
}

EXAMPLE_OF_WORK = {
    "@type": "CreativeWork",
    "about": [
        {
            "@type": "DefinedTerm",
            "inDefinedTermSet": "NCIT",
            "name": "Antimicrobial Peptide",
            "displayName": "Antimicrobial Peptide",
            "url": "http://purl.obolibrary.org/obo/NCIT_C123795",
            "termCode": "NCIT:C123795",
        }
    ],
    "encodingFormat": [
        {
            "@type": "DefinedTerm",
            "name": "CSV",
            "url": "http://edamontology.org/format_3752",
            "inDefinedTermSet": "EDAM",
        },
        {
            "@type": "DefinedTerm",
            "name": "FASTA",
            "url": "http://edamontology.org/format_1332",
            "inDefinedTermSet": "EDAM",
        },
    ],
    "schemaVersion": "https://dbaasp.org/v3/api-docs",
}

HOWTO_STEPS = [
    "Step 1. Batch pull all peptides in DBAASP via the API",
    (
        "Step 2. Group the peptide records by the targetSpecies.value to "
        "obtain the available organism names, and the count of peptides per "
        "organism"
    ),
    (
        "Step 3. Use each organism to generate a DataCollection url to search "
        "DBAASP by organism and get the count of the records for "
        "'collectionSize'. Generate the values for the 'species', "
        "'infectiousAgent', 'url', 'name', and 'description' properties based "
        "on the information retrieved for the organism and any templated "
        "text."
    ),
    (
        "Step 4. Fill in manually curated fields from the DBAASP resource "
        "catalog record such as 'measurementTechnique', 'variableMeasured', "
        "'topicCategory', 'conditionsOfAccess', 'usageInfo'."
    ),
]

DESCRIPTION_TEMPLATE = (
    "Amino acid sequences, chemical modifications, 3D structures, "
    "bioactivities and toxicities of peptides that possess antimicrobial "
    "properties targeting {infectiousAgent.name}. For more details, visit: "
    "{url}"
)

NAME_TEMPLATE = "Antimicrobial peptide records targeting {infectiousAgent.name} at DBAASP"

_SUBSPECIES_MARKERS = {
    "subsp", "subsp.", "subspecies", "var", "var.", "serovar",
    "pathovar", "pv", "pv.", "biovar", "bv", "bv.", "sv", "sv.",
}
_LOWER_EPITHET = re.compile(r"^[a-z]+$")


def _http_get_json(url: str, params: Optional[dict] = None) -> Any:
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    last_exc: Optional[Exception] = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(
                url,
                params=params,
                headers=headers,
                timeout=REQUEST_TIMEOUT,
            )
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_exc = e
            if attempt == MAX_RETRIES:
                break
            sleep_for = RETRY_BACKOFF * attempt
            logger.warning(
                "GET %s failed (attempt %d/%d): %s -- retrying in %ds",
                url, attempt, MAX_RETRIES, e, sleep_for,
            )
            time.sleep(sleep_for)
    raise RuntimeError(f"GET {url} failed after {MAX_RETRIES} attempts: {last_exc}")


def normalize_species(name: str) -> str:
    """Collapse strain/ATCC-style suffixes so entries group by genus+species.

    Rules:
    - "Genus species [strain/ATCC ...]" -> "Genus species"
    - Preserve "subsp./var./serovar <epithet>" annotations after the species
    - Virus-like first tokens (contain hyphen/digit) -> just the first token
    - Anything else is kept as-is
    """
    if not name:
        return ""
    cleaned = re.sub(r"\s+", " ", name.strip(" ,;"))
    if not cleaned:
        return ""
    tokens = cleaned.split()
    if len(tokens) < 2:
        return cleaned

    first, second = tokens[0], tokens[1]
    if first[:1].isupper() and _LOWER_EPITHET.fullmatch(second):
        base = f"{first} {second}"
        i = 2
        while i + 1 < len(tokens):
            marker = tokens[i].lower().rstrip(".")
            if marker in {m.rstrip(".") for m in _SUBSPECIES_MARKERS}:
                if _LOWER_EPITHET.fullmatch(tokens[i + 1]):
                    base += f" {tokens[i]} {tokens[i + 1]}"
                    i += 2
                    continue
            break
        return base

    if re.search(r"[-\d]", first):
        return first

    return cleaned


def _iter_peptide_stubs(limit: Optional[int] = None) -> Iterator[dict]:
    """Page through /peptides. Each yielded dict has id, dbaaspId, name."""
    offset = 0
    seen = 0
    while True:
        payload = _http_get_json(
            f"{API_BASE}/peptides",
            params={"limit": LIST_PAGE_SIZE, "offset": offset},
        )
        data = payload.get("data") or []
        total = payload.get("totalCount") or 0
        if not data:
            return
        for row in data:
            yield row
            seen += 1
            if limit is not None and seen >= limit:
                return
        offset += len(data)
        if offset >= total:
            return


def _fetch_peptide_detail(peptide_id: int) -> Optional[dict]:
    try:
        return _http_get_json(f"{API_BASE}/peptides/{peptide_id}")
    except Exception as e:
        logger.warning("Could not fetch peptide %s: %s", peptide_id, e)
        return None


def _extract_target_species_names(detail: dict) -> set[str]:
    names: set[str] = set()
    if not detail:
        return names
    queues: list[dict] = [detail]
    for monomer in detail.get("monomers") or []:
        if isinstance(monomer, dict):
            queues.append(monomer)
    for node in queues:
        activities = node.get("targetActivities") or []
        for a in activities:
            ts = a.get("targetSpecies") if isinstance(a, dict) else None
            if isinstance(ts, dict):
                raw = ts.get("name") or ts.get("value")
                if raw:
                    names.add(raw.strip())
            elif isinstance(ts, str) and ts.strip():
                names.add(ts.strip())
    return names


def _make_record_id(species: str) -> str:
    raw = species.lower()
    raw = re.sub(r"[^a-z0-9]+", "_", raw).strip("_")
    raw = re.sub(r"_+", "_", raw)
    return f"dbaasp_{raw}" if raw else "dbaasp_unknown"


def _species_url(species: str) -> str:
    encoded = urllib.parse.quote(species)
    return (
        "https://dbaasp.org/search?id.value=&name.value=&sequence.value="
        "&sequence.option=full&sequenceLength.value=&complexity.value="
        "&synthesisType.value=&nTerminus.value=&cTerminus.value="
        "&unusualAminoAcid.value=&intraChainBond.value="
        "&intraChainBond.chainParticipating=&intraChainBond.cycleType="
        "&interChainBond.value=&interChainBond.chainParticipating="
        "&interChainBond.participatingInCycle=&coordinationBond.value="
        "&uniprot.value=&pubchem.value=&pubchemId.value="
        "&threeDStructure.value=&kingdom.value=&source.value="
        f"&targetSpecies.value={encoded}"
        "&synergy.value=&articleAuthor.value=&articleJournal.value="
        "&articleYear.value=&articleVolume.value=&articlePages.value="
        "&articleTitle.value="
    )


def _example_action(species: str) -> dict:
    encoded = urllib.parse.quote(species)
    return {
        "@type": "Action",
        "name": "Use API call for an example record",
        "target": (
            f"https://dbaasp.org/peptides?targetSpecies.value={encoded}"
            "&limit=500&offset=0"
        ),
    }


def _build_example_of_work(species: str) -> dict:
    eow = copy.deepcopy(EXAMPLE_OF_WORK)
    eow["potentialAction"] = _example_action(species)
    return eow


def _build_is_based_on(species: str) -> list[dict]:
    return [
        {
            "@type": "Action",
            "name": (
                "Process for creating the DBAASP Peptide Data Collection "
                "in the NIAID Data Ecosystem Discovery Portal"
            ),
            "description": (
                "How this DBAASP Antimicrobial Peptide DataCollection Record "
                "was generated for the NIAID Data Ecosystem."
            ),
            "actionProcess": {
                "@type": "HowTo",
                "step": list(HOWTO_STEPS),
            },
        },
        {
            "@type": "ResourceCatalog",
            "name": "Database of Antimicrobial Activity and Structure of Peptides",
            "url": "https://data.niaid.nih.gov/resources?id=dde_cc13744ba5c15cca",
        },
    ]


def _dbaasp_version() -> Optional[str]:
    try:
        info = _http_get_json(f"{API_BASE}/v3/api-docs")
        return (info.get("info") or {}).get("version")
    except Exception as e:
        logger.warning("Could not fetch DBAASP API version: %s", e)
        return None


def _build_record(
    species: str,
    peptide_count: int,
    version: Optional[str],
    today: str,
) -> dict[str, Any]:
    url = _species_url(species)
    infectious_agent = {
        "@type": "DefinedTerm",
        "name": species,
    }
    description = DESCRIPTION_TEMPLATE.replace(
        "{infectiousAgent.name}", species
    ).replace("{url}", url)
    name = NAME_TEMPLATE.replace("{infectiousAgent.name}", species)

    record: dict[str, Any] = {
        "_id": _make_record_id(species),
        "@type": "DataCollection",
        "about": copy.deepcopy(ABOUT_DEFINED_TERM),
        "includedInDataCatalog": [
            {
                **copy.deepcopy(INCLUDED_IN_DATA_CATALOG),
                "versionDate": today,
            }
        ],
        "name": name,
        "url": url,
        "description": description,
        "collectionSize": {
            "minValue": peptide_count,
            "unitText": "Peptides",
        },
        "author": copy.deepcopy(AUTHORS),
        "creator": copy.deepcopy(CREATOR),
        "citation": copy.deepcopy(CITATION),
        "creditText": CREDIT_TEXT,
        "conditionsOfAccess": "Open",
        "funding": copy.deepcopy(FUNDING),
        "healthCondition": [],
        "infectiousAgent": [infectious_agent],
        "license": "https://creativecommons.org/licenses/by-nc/4.0/",
        "measurementTechnique": copy.deepcopy(MEASUREMENT_TECHNIQUE),
        "variableMeasured": copy.deepcopy(VARIABLE_MEASURED),
        "topicCategory": copy.deepcopy(TOPIC_CATEGORY),
        "usageInfo": copy.deepcopy(USAGE_INFO),
        "species": copy.deepcopy(SPECIES),
        "date": today,
        "dateModified": today,
        "isAccessibleForFree": True,
        "exampleOfWork": _build_example_of_work(species),
        "isBasedOn": _build_is_based_on(species),
    }
    if version:
        record["version"] = version
    return record


def _collect_species_groups(
    peptide_limit: Optional[int] = None,
) -> dict[str, set[int]]:
    """Walk every peptide and return {normalized_species: {peptide_id, ...}}."""
    groups: dict[str, set[int]] = {}
    count = 0
    for stub in _iter_peptide_stubs(limit=peptide_limit):
        pid = stub.get("id")
        if pid is None:
            continue
        count += 1
        if count % 500 == 0:
            logger.info("Fetched %d peptide details so far...", count)
        detail = _fetch_peptide_detail(pid)
        if detail is None:
            continue
        for raw in _extract_target_species_names(detail):
            normalized = normalize_species(raw)
            if not normalized:
                continue
            groups.setdefault(normalized, set()).add(pid)
    logger.info(
        "Processed %d peptides, found %d unique target species groups",
        count, len(groups),
    )
    return groups


def parse() -> Iterator[dict[str, Any]]:
    version = _dbaasp_version()
    today = dt.datetime.now(dt.timezone.utc).date().isoformat()
    groups = _collect_species_groups()
    for species in sorted(groups.keys()):
        peptide_ids = groups[species]
        yield _build_record(species, len(peptide_ids), version, today)


if __name__ == "__main__":
    import json
    import sys

    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    version = _dbaasp_version()
    today = dt.datetime.now(dt.timezone.utc).date().isoformat()
    groups = _collect_species_groups(peptide_limit=limit)
    for species in sorted(groups.keys()):
        record = _build_record(species, len(groups[species]), version, today)
        print(json.dumps(record, ensure_ascii=False))
