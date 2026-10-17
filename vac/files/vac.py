import datetime
import logging
import math
import re
from collections import Counter

import dateutil.parser
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

VAC_URL = "https://vac.niaid.nih.gov"

# Cap on distinct values kept per output key. Values reported by the most
# adjuvants are kept.
MAX_VALUES = 100

# Longest keyword kept (Mapping rows 59, 175, 199).
MAX_KEYWORD_LENGTH = 100

# From utils/validate.py.
PLACEHOLDER_TERMS = frozenset(
    {
        "",
        "-",
        "--",
        "n/a",
        "n.a.",
        "na",
        "none",
        "null",
        "nan",
        "not available",
        "not provided",
        "not known",
        "not specified",
        "not applicable",
        "no data",
        "missing",
        "unknown",
        "normal",
        "disease",
    }
)


# ---------------------------------------------------------------------------
# Field mappings, from "VAC - Mapping.tsv", "VAC - grouping term.tsv" and
# "VAC - VAC - Data Collection.tsv".
# ---------------------------------------------------------------------------

# preclinicalHostSpecies.name -> (scientific name, NCBI taxon id).
HOST_SPECIES = {
    "Mouse": ("Mus musculus", "10090"),
    "Rat": ("Rattus norvegicus", "10116"),
    "Rabbit": ("Oryctolagus cuniculus", "9986"),
    "Ferret": ("Mustela putorius furo", "9669"),
    "Hamster": ("Mesocricetus auratus", "10036"),
    "Guinea Pig": ("Cavia porcellus", "10141"),
    "Pig": ("Sus scrofa", "9823"),
    "Human": ("Homo sapiens", "9606"),
}
HUMAN = HOST_SPECIES["Human"]

# NHP note keywords -> species. A note matching none is Primates.
PRIMATE_SPECIES = (
    (("rhesus", "mulatta"), ("Macaca mulatta", "9544")),
    (("cynomo", "fascicularis"), ("Macaca fascicularis", "9541")),
    (("baboon", "papio"), ("Papio anubis", "9555")),
)
PRIMATES = ("Primates", "9443")

# Hosts named in "Other" notes. Guinea pig is matched before pig.
HOST_WORDS = (
    (re.compile(r"\bguinea pigs?\b", re.I), "Guinea Pig"),
    (re.compile(r"\b(mouse|mice)\b", re.I), "Mouse"),
    (re.compile(r"\brats?\b", re.I), "Rat"),
    (re.compile(r"\brabbits?\b", re.I), "Rabbit"),
    (re.compile(r"\bferrets?\b", re.I), "Ferret"),
    (re.compile(r"\bhamsters?\b", re.I), "Hamster"),
    (re.compile(r"\bpigs?\b", re.I), "Pig"),
    (re.compile(r"\bhumans?\b", re.I), "Human"),
)
NHP_WORDS = re.compile(r"\b(nhps?|non-human primates?)\b", re.I)

TAXON_IDS = {
    name.casefold(): taxon_id
    for name, taxon_id in (*HOST_SPECIES.values(), *(species for _, species in PRIMATE_SPECIES), PRIMATES)
}

# Hosts whose notes name lab strains (sample.aggregateElement.associatedGenotype).
# Human notes give ages and cell lines.
STRAIN_HOSTS = frozenset(HOST_SPECIES) - {"Human"}

# preclinicalHostSpecies.name -> description wording.
HOST_DESCRIPTIONS = {"NHP": "non-human primate"}

# Host and route notes that name no subject.
NON_SUBJECT = re.compile(r"(in vitro|in silico|n/a)\b.*", re.I)

# immunizationRoute abbreviations (Mapping rows 192, 267). SC is VAC's SQ.
ROUTES = {
    "IM": "intramuscular immunization",
    "SQ": "subcutaneous immunization",
    "SC": "subcutaneous immunization",
    "IN": "intranasal immunization",
    "ID": "intradermal immunization",
    "SL": "sublingual immunization",
    "IP": "intraperitoneal immunization",
    "IV": "intravenous immunization",
    "IT": "intratracheal immunization",
}
ROUTE_ABBREVIATION = re.compile(r"\b(" + "|".join(ROUTES) + r")\b")

PHASES = {
    "Phase I": "Phase I clinical trial",
    "Phase II": "Phase II clinical trial",
    "Phase III": "Phase III clinical trial",
    "EUA": "Emergency Use Authorization",
}

# pathogen.name / clinicalPathogen.name -> healthCondition (name, url).
DISEASE_AREAS = {
    "Allergy": ("allergic disease", "http://purl.obolibrary.org/obo/MONDO_0005271"),
    "Autoimmunity": ("autoimmune disease", "http://purl.obolibrary.org/obo/MONDO_0007179"),
    "Bacteria": ("bacterial infectious disease", "http://purl.obolibrary.org/obo/MONDO_0005113"),
    "Cancer": ("cancer", "http://purl.obolibrary.org/obo/MONDO_0004992"),
    "Fungi": ("fungal infectious disease", "http://purl.obolibrary.org/obo/MONDO_0002041"),
    "Virus": ("viral infectious disease", "http://purl.obolibrary.org/obo/MONDO_0005108"),
    "Eukaryotic Pathogens": ("parasitic infectious disease", "http://purl.obolibrary.org/obo/MONDO_0005135"),
    "Haptens": ("exposure to hapten", "http://purl.obolibrary.org/obo/ECTO_9001900"),
    "Other": ("other disease", "http://purl.obolibrary.org/obo/MONDO_0000001"),
}
GROUPING_TERMS = {name.casefold(): (name, url) for name, url in DISEASE_AREAS.values()}

# Areas named in the description's "as well as" clause.
DESCRIBED_AREAS = ("Cancer", "Allergy", "Autoimmunity")

# pathogen.note / clinicalPathogen.note values naming a disease or an organism
# (Mapping row 272), casefolded. Other values go to keywords.
AREA_NOTE_TERMS = {
    "pancreatic cancer": ("healthCondition", "pancreatic cancer"),
    "alzheimer's disease": ("healthCondition", "Alzheimer disease"),
    "beta-amyloid plaques (alzheimer's disease)": ("healthCondition", "Alzheimer disease"),
    "opioid overdose": ("healthCondition", "Opioid overdose"),
    "onchocerciasis": ("healthCondition", "Onchocerciasis"),
    "brugia malayi": ("infectiousAgent", "Brugia malayi"),
    "onchocerca volvulus": ("infectiousAgent", "Onchocerca volvulus"),
    "plasmodium falciparum": ("infectiousAgent", "Plasmodium falciparum"),
    "plasmodium falciparum (fp)": ("infectiousAgent", "Plasmodium falciparum"),
    "plasmodium falciparum (pf)": ("infectiousAgent", "Plasmodium falciparum"),
    "cytomegalovirus": ("infectiousAgent", "Cytomegalovirus"),
    "ebv": ("infectiousAgent", "EBV"),
}

# pathogen.note values describing the study design (Mapping row 272).
STUDY_DESIGN_NOTE = re.compile(
    r"(in vitro|model antigen|biodistribution).*|ovalbumin|toxicology|mechanistic studies|general", re.I
)

# clinicalParentOrganism.note values naming only a trial or phase (Mapping row 66).
TRIAL_LABEL = re.compile(r"\((trial|study) [^)]*\)|phase \d+", re.I)

# pocInstitute.name spellings of one institute (Mapping row 186).
INSTITUTE_ALIASES = {
    "UCSD": "University of California, San Diego",
    "Vaxine LTD": "Vaxine Pty Ltd",
    "Molecular Express": "Molecular Express, Inc.",
    "Tulane U": "Tulane University",
    "InImmune Corp": "Inimmune Corp.",
}

NCT_ID = re.compile(r"NCT\d{8}")
# Non-NCT trial registrations, e.g. ACTRN12620001177943.
REGISTRY_ID = re.compile(r"[A-Z]{3,}\d{6,}")

# clinicalAgeCategories names: "18-50 yrs", "> 65 yrs", "0-12 mo".
AGE_BAND = re.compile(r"(?:(\d+)\s*-\s*(\d+)|>\s*(\d+))\s*(yrs|mo)")
AGE_UNITS = {"yrs": "years", "mo": "months"}

LIST_SEPARATORS = r",|;|\s+and\s+"
ROUTE_SEPARATORS = r",|;|/|\s+and\s+"
READOUT_SEPARATORS = r",|;"
PARENTHETICAL = re.compile(r"\([^()]*\)")


# ---------------------------------------------------------------------------
# Static blocks
# ---------------------------------------------------------------------------

IS_BASED_ON = {
    "@type": "Action",
    "name": "DataCollection Generation Process in the NIAID Data Ecosystem",
    "disambiguatingDescription": (
        "How this VAC ChemicalSubstance DataCollection record was generated for the NIAID Data Ecosystem."
    ),
    "description": (
        "This record aggregated all VAC ChemicalSubstance records into a single collection. Collection size and "
        "dates are derived from the underlying records; descriptive fields such as topic categories, conditions "
        "of access and usage information are manually curated."
    ),
    "actionProcess": {
        "@type": "HowTo",
        "step": [
            "1. Page through https://vac.niaid.nih.gov/api/adjuvant/search?clinicalTrial=-1&page={n}&pageSize=100 "
            "to list every adjuvant id; the response's count gives collectionSize.",
            "2. Fetch each record from https://vac.niaid.nih.gov/api/adjuvant/get?id={id} and keep those with "
            "status 'Published' and isActive true, and their active preclinical and clinical entries.",
            "3. Aggregate and de-duplicate the mapped fields across all adjuvants and their preclinical and "
            "clinical entries (organisms, host species, trial age bands, readouts, antigens, routes, receptors, "
            "points of contact, PubMed IDs, NCT IDs), keeping the 100 values of each field reported by the most "
            "adjuvants.",
            "4. Take the earliest createdDate and the latest lastUpdated across all adjuvants for dateCreated and "
            "dateModified.",
            "5. Create a single DataCollection record with the templated description and one hasPart entry per "
            "adjuvant.",
        ],
    },
}

# EDAM topics: Vaccinology, Immunology, Preclinical and clinical studies.
TOPIC_CATEGORIES = [
    {
        "@type": "DefinedTerm",
        "name": "Vaccinology",
        "identifier": "topic_3966",
        "url": "http://edamontology.org/topic_3966",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Immunology",
        "identifier": "topic_0804",
        "url": "http://edamontology.org/topic_0804",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Preclinical and clinical studies",
        "identifier": "topic_3379",
        "url": "http://edamontology.org/topic_3379",
        "inDefinedTermSet": "EDAM",
    },
]

USAGE_INFO = {
    "@type": "CreativeWork",
    "name": "How To Use VAC",
    "url": "https://www.niaid.nih.gov/research/vaccine-adjuvant-compendium-vac",
    "description": (
        "Adjuvant records in VAC list the points of contact for each adjuvant to allow vaccine researchers to "
        "contact the intellectual property (IP) holder directly."
    ),
}

CREDIT_TEXT = (
    "Data from the NIAID Vaccine Adjuvant Compendium (VAC), https://vac.niaid.nih.gov/. Authors are encouraged "
    "to cite the Vaccine Ontology ID (VO-ID) of each adjuvant in publications."
)

# nde.py types `about` as {@type, description, displayName, name, url}; the
# sheet's identifier and inDefinedTermSet are left out.
ABOUT = {
    "@type": "DefinedTerm",
    "name": "ChemicalSubstance",
    "displayName": "Chemical Substance",
    "url": "http://purl.obolibrary.org/obo/NCIT_C103205",
    "description": "For schema, consider https://schema.org/ChemicalSubstance",
}

# nde.py's exampleOfWork.encodingFormat has no displayName.
EXAMPLE_OF_WORK = {
    "@type": "CreativeWork",
    "about": {
        "@type": "DefinedTerm",
        "name": "vaccine adjuvant",
        "displayName": "Vaccine Adjuvant",
        "url": "http://purl.obolibrary.org/obo/NCIT_C2554",
        "description": "Substances used in vaccines to enhance the ability of antigens to elicit an immune response.",
        "identifier": "NCIT_C2554",
        "inDefinedTermSet": "NCIT",
    },
    # Record API, search export and per-adjuvant product sheet.
    "encodingFormat": [
        {
            "@type": "DefinedTerm",
            "name": "JSON",
            "url": "http://edamontology.org/format_3464",
            "identifier": "format_3464",
            "inDefinedTermSet": "EDAM",
        },
        {
            "@type": "DefinedTerm",
            "name": "xlsx",
            "url": "http://edamontology.org/format_3620",
            "identifier": "format_3620",
            "inDefinedTermSet": "EDAM",
        },
        {
            "@type": "DefinedTerm",
            "name": "PDF",
            "url": "http://edamontology.org/format_3508",
            "identifier": "format_3508",
            "inDefinedTermSet": "EDAM",
        },
    ],
    "schemaVersion": "https://www.niaid.nih.gov/research/vaccine-adjuvant-compendium-vac",
}

# nde.py's author has no alternateName or parentOrganization: the acronym is
# folded into the name and the parent becomes an affiliation.
BASE_AUTHOR = {
    "@type": "Organization",
    "name": "Division of Allergy, Immunology, and Transplantation (DAIT)",
    "url": "https://www.niaid.nih.gov/about/dait",
    "affiliation": {"@type": "Organization", "name": "National Institute of Allergy and Infectious Diseases"},
}

NIAID_GRANT = {
    "@type": "MonetaryGrant",
    "funder": {
        "@type": "Organization",
        "name": "National Institute of Allergy and Infectious Diseases",
        "alternateName": "NIAID",
        "parentOrganization": "NIH",
    },
}


def insert_value(d, key, value, extend=False):
    """ Insert a value into a dictionary, promoting to a list on repeats and dropping duplicates.
    """

    if extend:
        d[key] = (d[key] + " " + value).strip() if d.get(key) else value
        return

    was_list = isinstance(d.get(key), list)
    merged = list(d[key]) if was_list else ([d[key]] if key in d else [])
    for item in value if isinstance(value, list) else [value]:
        if item not in merged:
            merged.append(item)
    d[key] = merged if was_list or isinstance(value, list) or len(merged) > 1 else merged[0]

def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(str(val), ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError, OverflowError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt


def iter_vac_records():
    page = 1
    pageSize = 100
    url = f"https://vac.niaid.nih.gov/api/adjuvant/search?clinicalTrial=-1&page={page}&pageSize={pageSize}"

    data = requests.get(url).json()
    total_records = data["count"]
    total_pages = math.ceil(total_records / pageSize)
    for page in range(1, total_pages + 1):
        print(f"Fetching page {page} of {total_pages}...")
        url = f"https://vac.niaid.nih.gov/api/adjuvant/search?clinicalTrial=-1&page={page}&pageSize={pageSize}"
        data = requests.get(url).json()
        for record in data["records"]:
            print(f"Fetching record {record['id']}...")
            record = requests.get(f"https://vac.niaid.nih.gov/api/adjuvant/get?id={record['id']}")
            yield record.json()


def get_collection_size():
    """The count field of the adjuvant search."""
    resp = requests.get(
        f"{VAC_URL}/api/adjuvant/search",
        params={"clinicalTrial": -1, "page": 1, "pageSize": 1},
        timeout=120,
    )
    resp.raise_for_status()
    return resp.json()["count"]


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def _get(obj, *keys):
    """obj[key1][key2]..., or None once a level is missing."""
    for key in keys:
        if not isinstance(obj, dict):
            return None
        obj = obj.get(key)
    return obj


def _clean(value):
    """value with whitespace runs collapsed, or None unless it is a non-blank string."""
    if not isinstance(value, str):
        return None
    return " ".join(value.split()) or None


def _names(items):
    """Cleaned names of a list of {"name": ...} objects."""
    for item in items or []:
        if isinstance(item, dict) and (name := _clean(item.get("name"))):
            yield name


def _join(items):
    """'a, b and c'."""
    items = list(items)
    return items[0] if len(items) == 1 else f"{', '.join(items[:-1])} and {items[-1]}"


def is_uninformative(value):
    """True for a placeholder, a bare "Other", or text with no letters or digits."""
    text = value.strip().casefold()
    return text in PLACEHOLDER_TERMS or text == "other" or not any(ch.isalnum() for ch in text)


def split_list(text, separators):
    """Cleaned pieces of text, split on separators outside parentheses."""
    groups = []

    def hide(match):
        groups.append(match.group())
        return f"\0{len(groups) - 1}\0"

    def restore(piece):
        return re.sub(r"\0(\d+)\0", lambda match: groups[int(match.group(1))], piece)

    pieces = re.split(separators, PARENTHETICAL.sub(hide, text), flags=re.I)
    return [piece for piece in (_clean(restore(p)) for p in pieces) if piece]


class Tally:
    """Distinct values of one output key, counted by adjuvant.

    Values merge case-insensitively under their most frequent spelling.
    """

    def __init__(self):
        self.adjuvants = {}
        self.mentions = Counter()
        self.spellings = {}

    def add(self, adjuvant_id, value):
        value = _clean(value)
        if not value:
            return
        key = value.casefold()
        self.adjuvants.setdefault(key, set()).add(adjuvant_id)
        self.mentions[key] += 1
        self.spellings.setdefault(key, Counter())[value] += 1

    def __len__(self):
        return len(self.adjuvants)

    def ranked(self):
        """Every value, by number of adjuvants, then mentions, then first appearance."""
        keys = sorted(self.adjuvants, key=lambda key: (-len(self.adjuvants[key]), -self.mentions[key]))
        return [self.spellings[key].most_common(1)[0][0] for key in keys]

    def top(self, limit=MAX_VALUES):
        """The `limit` highest-ranked values, alphabetically."""
        return sorted(self.ranked()[:limit], key=str.casefold)


def tally(records, values):
    """Tally values(record) over every adjuvant."""
    counts = Tally()
    for record in records:
        for value in values(record):
            counts.add(record["id"], value)
    return counts


def is_published(record):
    return record.get("status") == "Published" and bool(record.get("isActive"))


def preclinicals(record):
    """The adjuvant's active preclinical entries."""
    return [entry for entry in record.get("preclinicals") or [] if entry.get("isActive")]


def clinicals(record):
    """The adjuvant's active clinical entries."""
    return [entry for entry in record.get("clinicals") or [] if entry.get("isActive")]


# ---------------------------------------------------------------------------
# Field values, per adjuvant
# ---------------------------------------------------------------------------


def primates(text):
    """Primate species named in an NHP note."""
    text = text.casefold()
    return [species for keywords, species in PRIMATE_SPECIES if any(keyword in text for keyword in keywords)]


def named_hosts(text):
    """Fixed-table species named in one piece of an "Other" host note."""
    found = primates(text)
    if not found and NHP_WORDS.search(text):
        found.append(PRIMATES)
    for pattern, host in HOST_WORDS:
        if pattern.search(text):
            found.append(HOST_SPECIES[host])
            break
    return found


def host_terms(entry):
    """(species names, strains) for a preclinical entry's host."""
    host = entry.get("preclinicalHostSpecies") or {}
    name = _clean(host.get("name"))
    note = _clean(host.get("note"))
    if note and (NON_SUBJECT.fullmatch(note) or is_uninformative(note)):
        note = None

    if name == "NHP":
        return [species for species, _ in primates(note or "") or [PRIMATES]], []
    if name in HOST_SPECIES:
        strains = []
        if note and name in STRAIN_HOSTS:
            strains = [
                piece
                for piece in split_list(note, LIST_SEPARATORS)
                if not (NON_SUBJECT.fullmatch(piece) or is_uninformative(piece))
            ]
        return [HOST_SPECIES[name][0]], strains

    species = [name] if name and name != "Other" else []
    for piece in split_list(note, LIST_SEPARATORS) if note else []:
        if NON_SUBJECT.fullmatch(piece) or is_uninformative(piece):
            continue
        named = [scientific_name for scientific_name, _ in named_hosts(piece)]
        species.extend(named or [_clean(PARENTHETICAL.sub("", piece))])
    return [s for s in species if s], []


def host_species(record):
    """species: preclinical hosts, plus Homo sapiens for each clinical entry."""
    for entry in preclinicals(record):
        yield from host_terms(entry)[0]
    for _ in clinicals(record):
        yield HUMAN[0]


def strains(record):
    for entry in preclinicals(record):
        yield from host_terms(entry)[1]


def disease_areas(record):
    """pathogen.name / clinicalPathogen.name of each entry."""
    for entry in preclinicals(record):
        if area := _clean(_get(entry, "pathogen", "name")):
            yield area
    for entry in clinicals(record):
        if area := _clean(_get(entry, "clinicalPathogen", "name")):
            yield area


def area_note_terms(note):
    """(output key, value) pairs for a pathogen.note or clinicalPathogen.note."""
    note = _clean(note)
    if not note or is_uninformative(note) or STUDY_DESIGN_NOTE.fullmatch(note):
        return
    for piece in [note] if note.casefold() in AREA_NOTE_TERMS else split_list(note, LIST_SEPARATORS):
        if is_uninformative(piece) or STUDY_DESIGN_NOTE.fullmatch(piece):
            continue
        yield AREA_NOTE_TERMS.get(piece.casefold(), ("keywords", piece))


def area_notes(record, key):
    """Values of output `key` from the entries' disease-area notes."""
    notes = [_get(entry, "pathogen", "note") for entry in preclinicals(record)]
    notes += [_get(entry, "clinicalPathogen", "note") for entry in clinicals(record)]
    for note in notes:
        for field, value in area_note_terms(note):
            if field == key:
                yield value


def infectious_agents(record):
    for entry in preclinicals(record):
        if name := _clean(_get(entry, "parentOrganism", "name")):
            yield name
    for entry in clinicals(record):
        if name := _clean(_get(entry, "clinicalParentOrganism", "name")):
            yield name
    yield from area_notes(record, "infectiousAgent")


def health_conditions(record):
    """Grouping terms of the disease areas, then diseases named in area notes."""
    for area in disease_areas(record):
        if area in DISEASE_AREAS:
            yield DISEASE_AREAS[area][0]
    yield from area_notes(record, "healthCondition")


def measurement_techniques(record):
    """preclinicalReadouts names; an "Other" readout's note, split."""
    for entry in preclinicals(record):
        for readout in entry.get("preclinicalReadouts") or []:
            name = _clean(readout.get("name"))
            if name and name != "Other":
                yield name
            if note := _clean(readout.get("note")):
                yield from (piece for piece in split_list(note, READOUT_SEPARATORS) if not is_uninformative(piece))


def routes(entry):
    """Immunization routes, abbreviations expanded; an "Other" route's note, split."""
    route = entry.get("immunizationRoute") or {}
    name = _clean(route.get("name"))
    if name and name != "Other":
        yield ROUTES.get(name, name)
    note = _clean(route.get("note"))
    if not note or NON_SUBJECT.fullmatch(note):
        return
    for piece in split_list(note, ROUTE_SEPARATORS):
        if NON_SUBJECT.fullmatch(piece) or is_uninformative(piece):
            continue
        abbreviations = ROUTE_ABBREVIATION.findall(piece)
        yield from (ROUTES[abbreviation] for abbreviation in abbreviations) if abbreviations else [piece]


def keyword_candidates(record):
    """Every keyword source in the Mapping, for one adjuvant and its entries."""
    for field in (
        "adjuvantCharacteristics",
        "adjuvantCombinationFormulations",
        "adjuvantFormulations",
        "adjuvantImmuneProfiles",
        "adjuvantReceptors",
    ):
        yield from _names(record.get(field))
    yield from _names(_get(record, "adjuvantChemicalNature", "values"))
    for receptor in record.get("adjuvantReceptors") or []:
        yield receptor.get("note")
    yield _get(record, "materialNature", "name")
    yield record.get("mechanismofAction")

    for entry in preclinicals(record):
        yield from routes(entry)
        yield from _names(_get(entry, "preclinicalAntigen", "values"))
        yield _get(entry, "preclinicalAntigen", "note")
        yield _get(entry, "parentOrganism", "note")
    for entry in clinicals(record):
        yield from _names(_get(entry, "clinicalAntigen", "values"))
        yield _get(entry, "clinicalAntigen", "note")
        note = _clean(_get(entry, "clinicalParentOrganism", "note"))
        if note and not TRIAL_LABEL.fullmatch(note):
            yield note
        if phase := _clean(_get(entry, "phase", "name")):
            yield PHASES.get(phase, phase)

    yield from (area for area in disease_areas(record) if area != "Other")
    yield from area_notes(record, "keywords")


def keywords(record):
    for value in keyword_candidates(record):
        value = _clean(value)
        if value and not is_uninformative(value) and len(value) <= MAX_KEYWORD_LENGTH:
            yield value


def age_bands(record):
    for entry in clinicals(record):
        yield from _names(entry.get("clinicalAgeCategories"))


def trial_ids(record):
    """clinicalTrialGovId.note without whitespace; NCT IDs upper-cased."""
    for entry in clinicals(record):
        trial_id = re.sub(r"\s", "", _get(entry, "clinicalTrialGovId", "note") or "")
        if trial_id.upper().startswith("NCT"):
            trial_id = trial_id.upper()
        if trial_id:
            yield trial_id


def pubmed_ids(record):
    for entry in preclinicals(record) + clinicals(record):
        for pmid in entry.get("pubMedIds") or []:
            pmid = str(pmid).strip()
            if pmid.isdigit():
                yield pmid
            elif pmid:
                logger.warning(f"Skipping malformed PubMed ID {pmid!r} on adjuvant {record['id']}")


# ---------------------------------------------------------------------------
# Record sections
# ---------------------------------------------------------------------------


def build_authors(records):
    """author: DAIT, then one Person per point of contact."""
    people = {}
    for record in records:
        name = _clean(record.get("pocName"))
        if not name:
            continue
        person = people.setdefault(name.casefold(), {"@type": "Person", "name": name, "role": "Point of contact"})
        if email := _clean(record.get("pocEmail")):
            insert_value(person, "email", email)
        if title := _clean(record.get("pocTitle")):
            insert_value(person, "title", title)
        if institute := _clean(_get(record, "pocInstitute", "name")):
            institute = INSTITUTE_ALIASES.get(institute, institute)
            insert_value(person, "affiliation", {"@type": "Organization", "name": institute})

    names = tally(records, lambda record: [record.get("pocName")])
    authors = [dict(BASE_AUTHOR)]
    for name in names.top(MAX_VALUES - 1):
        person = people[name.casefold()]
        person["name"] = name
        authors.append(person)
    return authors


def build_has_part(records):
    """hasPart: one CreativeWork per adjuvant, by id."""
    parts = []
    for record in sorted(records, key=lambda record: record["id"])[:MAX_VALUES]:
        part = {
            "@type": "CreativeWork",
            "identifier": str(record["id"]),
            "name": _clean(record.get("adjuvantName")),
            "url": f"{VAC_URL}/view?id={record['id']}",
        }
        if vo_id := _clean(record.get("voId")):
            part["sameAs"] = f"http://purl.obolibrary.org/obo/{vo_id}"
        parts.append({key: value for key, value in part.items() if value})
    return parts


def defined_term(name):
    """DefinedTerm for a species, agent or condition name, with its fixed identifiers."""
    if taxon_id := TAXON_IDS.get(name.casefold()):
        return {"@type": "DefinedTerm", "name": name, "identifier": taxon_id}
    if grouping := GROUPING_TERMS.get(name.casefold()):
        # The hub's terms stage gives MONDO terms this identifier form.
        grouping_name, url = grouping
        ontology, identifier = url.rsplit("/", 1)[-1].split("_", 1)
        return {
            "@type": "DefinedTerm",
            "name": grouping_name,
            "identifier": identifier,
            "inDefinedTermSet": ontology,
            "url": url,
        }
    return {"@type": "DefinedTerm", "name": name}


def age_band(name):
    """QuantitativeValue for a clinicalAgeCategories name."""
    stage = {"@type": "QuantitativeValue", "name": name}
    match = AGE_BAND.fullmatch(name)
    if not match:
        logger.warning(f"Unrecognized age category: {name!r}")
        return stage
    low, high, above, unit = match.groups()
    stage["minValue"] = int(low or above)
    if high:
        stage["maxValue"] = int(high)
    stage["unitText"] = AGE_UNITS[unit]
    return stage


def build_sample(records):
    """sample: trial age bands and lab strains, as one aggregate Sample."""
    element = {"@type": "Sample"}
    stages = [age_band(name) for name in tally(records, age_bands).top()]
    stages.sort(key=lambda stage: (stage.get("unitText") != "months", stage.get("minValue", 0)))
    if stages:
        element["developmentalStage"] = stages
    if genotypes := tally(records, strains).top():
        element["associatedGenotype"] = genotypes
    if len(element) == 1:
        return None
    return {"@type": "SampleCollection", "aggregateElement": element}


def build_trials(records):
    """(nctid, isRelatedTo) from the clinical entries' registry IDs."""
    nct_ids, registrations = [], []
    for trial_id in tally(records, trial_ids).ranked():
        if NCT_ID.fullmatch(trial_id):
            nct_ids.append(trial_id)
        elif REGISTRY_ID.fullmatch(trial_id) and not trial_id.startswith("NCT"):
            registrations.append(trial_id)
        else:
            logger.warning(f"Skipping malformed trial ID {trial_id!r}")
    related = [
        {"@type": "CreativeWork", "identifier": trial_id, "relationship": "clinical trial registration"}
        for trial_id in sorted(registrations[:MAX_VALUES])
    ]
    return sorted(nct_ids[:MAX_VALUES]), related


def build_description(records, collection_size, agents, areas):
    """The templated description from the Data Collection sheet."""
    hosts = Counter()
    for record in records:
        for entry in preclinicals(record):
            host = _clean(_get(entry, "preclinicalHostSpecies", "name"))
            if host and host != "Other":
                hosts[HOST_DESCRIPTIONS.get(host, host.lower())] += 1
    conditions = [DISEASE_AREAS[area][0] for area in areas.ranked() if area in DESCRIBED_AREAS]
    preclinical_count = sum(len(preclinicals(record)) for record in records)
    clinical_count = sum(len(clinicals(record)) for record in records)
    clinical_adjuvants = sum(1 for record in records if clinicals(record))

    in_hosts = f" in {_join(host for host, _ in hosts.most_common())}" if hosts else ""
    against = f", against {len(agents)} pathogens including {_join(agents.ranked()[:5])}" if len(agents) else ""
    as_well_as = f", as well as {_join(conditions)}" if conditions else ""
    return (
        "Vaccine adjuvant records available through the NIAID Vaccine Adjuvant Compendium (VAC). "
        f"{collection_size} adjuvants have been evaluated in {preclinical_count} preclinical studies{in_hosts}, "
        f"and {clinical_adjuvants} of them in {clinical_count} clinical studies{against}{as_well_as}. "
        "Each VAC record describes the adjuvant's receptor targets, immune profile, formulation and product grade, "
        "lists supporting publications, and names a point of contact at the intellectual property holder."
    )


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_record(records, collection_size):
    """Build the VAC DataCollection record from the published adjuvants."""
    output = {
        "@context": "http://schema.org/",
        "@type": "DataCollection",
        "_id": "vac_adjuvants",
        "url": f"{VAC_URL}/browse",
        "includedInDataCatalog": {
            "@type": "DataCatalog",
            "name": "Vaccine Adjuvant Compendium (VAC)",
            "url": f"{VAC_URL}/",
            "versionDate": datetime.date.today().isoformat(),
            "archivedAt": f"{VAC_URL}/browse",
        },
        "name": "Vaccine adjuvant records at the Vaccine Adjuvant Compendium",
        # minValue, as in pathoplexus.
        "collectionSize": {"@type": "QuantitativeValue", "minValue": collection_size, "unitText": "adjuvants"},
        "conditionsOfAccess": "Open",
        "usageInfo": USAGE_INFO,
        "creditText": CREDIT_TEXT,
        "about": ABOUT,
        "exampleOfWork": EXAMPLE_OF_WORK,
        "topicCategory": TOPIC_CATEGORIES,
        "isBasedOn": [IS_BASED_ON],
    }

    agents = tally(records, infectious_agents)
    areas = tally(records, disease_areas)
    insert_value(output, "description", build_description(records, collection_size, agents, areas), extend=True)

    created = [_to_iso_date(record["createdDate"]) for record in records if record.get("createdDate")]
    modified = [_to_iso_date(record["lastUpdated"]) for record in records if record.get("lastUpdated")]
    if created := [date for date in created if date]:
        insert_value(output, "dateCreated", min(created))
    if modified := [date for date in modified if date]:
        insert_value(output, "dateModified", max(modified))

    insert_value(output, "author", build_authors(records))
    if parts := build_has_part(records):
        insert_value(output, "hasPart", parts)

    for key, values in (
        ("species", tally(records, host_species)),
        ("infectiousAgent", agents),
        ("healthCondition", tally(records, health_conditions)),
    ):
        if names := values.top():
            insert_value(output, key, [defined_term(name) for name in names])

    if techniques := tally(records, measurement_techniques).top():
        insert_value(output, "measurementTechnique", [{"@type": "DefinedTerm", "name": name} for name in techniques])
    if words := tally(records, keywords).top():
        insert_value(output, "keywords", words)

    if sample := build_sample(records):
        insert_value(output, "sample", sample)

    nct_ids, registrations = build_trials(records)
    if nct_ids:
        insert_value(output, "nctid", nct_ids)
    if registrations:
        insert_value(output, "isRelatedTo", registrations)

    # The hub's citations stage reads pmids as one comma-separated string.
    if pmids := tally(records, pubmed_ids).top():
        insert_value(output, "pmids", ", ".join(sorted(pmids, key=int)))

    if any(_get(record, "supportedByNiad", "value") for record in records):
        insert_value(output, "funding", [NIAID_GRANT])

    return output


def parse():
    """Yield the single VAC DataCollection record."""
    collection_size = get_collection_size()
    records = [record for record in iter_vac_records() if is_published(record)]
    logger.info(f"Building the VAC DataCollection from {len(records)} published adjuvants")
    yield build_record(records, collection_size)
