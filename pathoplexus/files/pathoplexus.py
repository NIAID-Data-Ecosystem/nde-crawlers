import datetime
import logging
from functools import cache

import dateutil.parser
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

LAPIS_URL = "https://lapis.pathoplexus.org"
BACKEND_URL = "https://backend.pathoplexus.org"
BATCH_SIZE = 10000

# Cap on distinct values kept per aggregated field. value_counts orders by
# frequency, so this keeps the most common values.
MAX_VALUES = 100

# Copied verbatim from PLACEHOLDER_TERMS in
# biothings-hub/files/nde-hub/utils/validate.py. The hub strips these from
# species / infectiousAgent / healthCondition only, so every other field this
# crawler emits -- associatedPhenotype, additionalProperty, spatialCoverage --
# has to be cleaned here or the placeholders reach the index. Keep in sync.
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

# The "no useful answer" members of Loculus's own enums, which arrive as real
# values rather than as an unset field: sequencingAssayType reports OTHER, and
# hostRole / exposureSetting / purposeOfSampling each offer an Other option. A
# bare "Other" is never worth a facet value, but it is not a placeholder in the
# hub's sense, so it is kept as a separate set.
UNINFORMATIVE_TERMS = frozenset({"other", "others", "unspecified", "undetermined"})

# Superseded revisions stay in the index, so an unfiltered query roughly
# double-counts any accession that was ever revised. Revoked records are
# hidden too, which is what the website counts (andv -> 548, not 601).
CURRENT_ONLY = {"versionStatus": "LATEST_VERSION", "isRevocation": False}

# Matches the website's default sort. Descending order is only expressible
# through the POST body, so both calls below POST instead of GET.
ORDER_BY = [{"field": "sampleCollectionDate", "type": "descending"}]


# ---------------------------------------------------------------------------
# Field mappings, from "Pathoplexus - Record Mapping (1).tsv".
#
# Each table maps an output key to the Pathoplexus fields that feed it. Values
# are aggregated across the collection, split on ";", and de-duplicated. Fields
# the sheet marks "ignore" -- and those whose note says they belong on a Sample
# record rather than a DataCollection -- are left out.
# ---------------------------------------------------------------------------

# Output key -> source fields. The Elasticsearch mapping in nde.py types these
# as plain keyword/text, so they stay lists of strings.
STRING_LIST_FIELDS = {
    "alternateName": ("specimenCollectorSampleId", "submissionId"),
    "identifier": ("gcaAccession", "gisaidIsolateId"),
    "keywords": ("outbreak",),
    "sameAs": ("gcaAccession",),
}

# Output key -> (@type, source fields). nde.py types these as objects with a
# name, so each value becomes {"@type": ..., "name": value} rather than a bare
# string. `ncbiSourceDb` names the repository the record came from, which makes
# it a DataCatalog; the rest are ontology terms.
NAMED_OBJECT_FIELDS = {
    "healthCondition": ("DefinedTerm", ("hostDisease",)),
    "measurementTechnique": (
        "DefinedTerm",
        (
            "diagnosticMeasurementMethod",
            "sequencingAssayType",
            "sequencingInstrument",
            "sequencingProtocol",
        ),
    ),
    "sdPublisher": ("DataCatalog", ("ncbiSourceDb",)),
    "variableMeasured": ("DefinedTerm", ("diagnosticMeasurementUnit", "diagnosticTargetGeneName")),
}

# nde.py types spatialCoverage as {@type, geo, name, identifier}. The sheet also
# asks for administrativeType and locationType, which the mapping has no slot
# for, so they are left off and geoLocCountry/hostOriginCountry merge by name.
SPATIAL_COVERAGE_FIELDS = (
    "exposureSetting",
    "geoLocAdmin1",
    "geoLocAdmin2",
    "geoLocCity",
    "geoLocCountry",
    "geoLocSite",
    "hostOriginCountry",
    "travelHistory",
)

# sample.aggregateElement keys nde.py types as objects with a name.
SAMPLE_OBJECT_FIELDS = {
    "anatomicalStructure": ("anatomicalMaterial", "anatomicalPart"),
    "associatedPhenotype": (
        "hostHealthOutcome",
        "hostHealthState",
        "previousInfectionDisease",
        "signsAndSymptoms",
    ),
    "cellType": ("cellLine",),
    "developmentalStage": ("hostAgeBin",),
    "instrument": ("collectionDevice", "sequencingInstrument"),
    "sampleType": ("bodyProduct", "hostRole", "sampleType"),
}

# sample.aggregateElement keys nde.py types as bare keywords.
SAMPLE_KEYWORD_FIELDS = {
    "associatedGenotype": (
        "clade",
        "genotype",
        "lineage",
        "lineage_S",
        "mutationsFromOutbreakFounder",
        "outbreakLineage",
        "serotype",
        "subtype",
    ),
    "sex": ("hostGender",),
}

# sample.aggregateElement keys nde.py types as bare text.
SAMPLE_TEXT_FIELDS = {
    "collectionMethod": ("collectionMethod",),
    "identifier": ("cultureId", "specimenCollectorSampleId"),
}

# sample.aggregateElement keys nde.py types as dates. Values are passed through
# unnormalised: Pathoplexus reports partial dates ("1905", "1900-01") and the
# Elasticsearch date type accepts those, whereas _to_iso_date would invent a
# month and day for them.
SAMPLE_DATE_FIELDS = {
    "dateCollected": ("sampleCollectionDate",),
    "dateProcessed": ("sampleReceivedDate",),
}

# sample.itemListElement is typed {@type, identifier, url}.
SAMPLE_ITEM_FIELDS = ("biosampleAccession",)

# isPartOf is an object; the sheet supplies a bioproject accession for it.
IS_PART_OF_FIELDS = ("bioprojectAccession",)

AUTHOR_NAME_FIELDS = ("authorAffiliations", "sequencedByOrganization")

# Segmented organisms carry the _L/_M/_S variants, unsegmented ones the bare
# name. collect_values drops whichever the organism does not have.
HAS_PART_IDENTIFIER_FIELDS = (
    "accession",
    "accessionVersion",
    "insdcAccessionBase",
    "insdcAccessionBase_L",
    "insdcAccessionBase_M",
    "insdcAccessionBase_S",
    "insdcAccessionFull",
    "insdcAccessionFull_L",
    "insdcAccessionFull_M",
    "insdcAccessionFull_S",
)
HAS_PART_SAME_AS_FIELDS = (
    "insdcAccessionFull",
    "insdcAccessionFull_L",
    "insdcAccessionFull_M",
    "insdcAccessionFull_S",
)
HAS_PART_URL_FIELDS = ("accessionVersion",)

IS_BASED_ON_IDENTIFIER_FIELDS = (
    "assemblyReferenceGenomeAccession",
    "insdcRawReadsAccession",
)

# Pairs read from a single aggregated call so the value and its identifier stay
# matched. (name field, identifier field, name-only fallback fields)
SPECIES_PAIR = ("hostNameScientific", "hostTaxonId", ("hostNameCommon",))
INFECTIOUS_AGENT_PAIR = ("ncbiVirusName", "ncbiVirusTaxId", ("previousInfectionOrganism",))


# ---------------------------------------------------------------------------
# Static blocks, identical on every record.
# ---------------------------------------------------------------------------

IS_BASED_ON = {
    "@type": "Action",
    "name": "DataCollection Generation Process in the NIAID Data Ecosystem",
    "description": "How this Pathoplexus DataCollection record was generated for the NIAID Data Ecosystem.",
    "actionProcess": {
        "@type": "HowTo",
        "step": [
            "1. Read the organism list from the Pathoplexus backend OpenAPI spec (components.schemas.Organism).",
            "2. For each organism, query the LAPIS databaseConfig endpoint for its display name and field list.",
            "3. Query the LAPIS aggregated endpoint for the record count and dataVersion, restricted to the latest non-revoked version of each accession.",
            "4. Query the LAPIS details endpoint sorted on sampleCollectionDate in each direction to find the earliest and latest collection dates.",
            "5. Query the LAPIS aggregated endpoint once per mapped metadata field, de-duplicating values across the fields that share an output key.",
            "6. Create one DataCollection record per organism from the counts, date bounds, accessions, hosts, places, and sample attributes.",
        ],
    },
}

# EDAM topics describing every Pathoplexus collection. Labels are the official
# ontology labels from http://edamontology.org.
TOPIC_CATEGORIES = [
    {
        "@type": "DefinedTerm",
        "name": "Genomics",
        "identifier": "topic_0622",
        "url": "http://edamontology.org/topic_0622",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Infectious disease",
        "identifier": "topic_3324",
        "url": "http://edamontology.org/topic_3324",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Public health and epidemiology",
        "identifier": "topic_3305",
        "url": "http://edamontology.org/topic_3305",
        "inDefinedTermSet": "EDAM",
    },
]

USAGE_INFO = {
    "name": "Data Use Terms",
    "url": "https://pathoplexus.org/about/terms-of-use/data-use-terms",
    "@type": "CreativeWork",
}

# nde.py types `about` as {@type, description, displayName, name, url} -- it has
# no identifier or inDefinedTermSet, unlike exampleOfWork.about below.
ABOUT = {
    "@type": "DefinedTerm",
    "name": "MolecularSequence",
    "displayName": "Molecular Sequence",
    "url": "http://purl.obolibrary.org/obo/NCIT_C164396",
    "description": "Nucleotide sequence in this case. For schema, consider SequenceAnnotation https://bioschemas.org/types/SequenceAnnotation/",
}

# exampleOfWork.about and .encodingFormat are typed
# {@type, identifier, inDefinedTermSet, name, url} -- no displayName or
# description, which is the mirror image of top-level `about`.
EXAMPLE_OF_WORK = {
    "@type": "CreativeWork",
    "about": {
        "@type": "DefinedTerm",
        "name": "Nucleotide Sequence",
        "url": "http://purl.obolibrary.org/obo/NCIT_C45374",
        "identifier": "NCIT_C45374",
        "inDefinedTermSet": "NCIT",
    },
    # Both download formats Pathoplexus offers, as EDAM terms.
    "encodingFormat": [
        {
            "@type": "DefinedTerm",
            "name": "FASTA",
            "url": "http://edamontology.org/format_1929",
            "identifier": "format_1929",
            "inDefinedTermSet": "EDAM",
        },
        {
            "@type": "DefinedTerm",
            "name": "TSV",
            "url": "http://edamontology.org/format_3475",
            "identifier": "format_3475",
            "inDefinedTermSet": "EDAM",
        },
    ],
    "schemaVersion": "https://pathoplexus.org/docs/concepts/metadataformat",
}

# Organisations credited on every collection, independent of the sequence data.
# Copied per record, since insert_value appends to the list it is given.
#
# The sheet gives PHA4GE as an alternateName and the SIB as a
# parentOrganization; author in nde.py has neither, so the acronym is folded
# into the name and the parent becomes an affiliation, which is mapped.
BASE_AUTHORS = [
    {
        "@type": "Organization",
        "name": "Public Health Alliance for Genomic Epidemiology (PHA4GE)",
        "url": "https://pha4ge.org/",
    },
    {
        "@type": "Organization",
        "name": "Loculus Development Team",
        "url": "https://loculus.org/#team",
        "affiliation": {"@type": "Organization", "name": "Swiss Institute of Bioinformatics"},
    },
]


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def insert_value(d, key, value, extend=False):
    """Insert a value into a dictionary, handling existing keys by converting to lists or extending strings as needed."""

    if key in d and not extend:
        if isinstance(d[key], list):
            if isinstance(value, list):
                for item in value:
                    if item not in d[key]:
                        d[key].append(item)
            elif value not in d[key]:
                d[key].append(value)
        else:
            if isinstance(value, list):
                d[key] = [d[key]] + [v for v in value if v != d[key]]
            elif d[key] != value:
                d[key] = [d[key], value]
    elif d.get(key) and extend:
        d[key] = (d.get(key) + " " + value).strip()
    else:
        d[key] = value


def is_uninformative(value):
    """True for a placeholder or a no-answer enum member, compared case-folded.

    Matching is exact on the whole value, so only a bare "Other" is dropped --
    "Other Animal" and "Host role - Other: Mother" survive.
    """
    text = str(value).strip().lower()
    return text in PLACEHOLDER_TERMS or text in UNINFORMATIVE_TERMS


def add_split_values(target, values, separator=";", max_size=MAX_VALUES):
    """Split each string in values on separator and add the pieces to target.

    Mutates target in place and returns it, so it can be called repeatedly to
    accumulate across several lists. Pieces are stripped, and empty ones are
    dropped so a trailing separator does not leave a blank entry behind.

    Placeholders and no-answer enum members are dropped per piece, so
    "Died;unknown" contributes only "Died". Filtering here rather than in
    collect_values keeps them from consuming max_size slots that real values
    would otherwise get.

    Stops once target holds max_size pieces. value_counts orders its keys by
    frequency, so passing one of those keeps the most common values.
    """
    if isinstance(values, str):
        values = [values]
    for value in values:
        for part in str(value).split(separator):
            part = part.strip()
            if not part or is_uninformative(part) or part in target:
                continue
            if len(target) >= max_size:
                return target
            target.add(part)
    return target


def _prune(value):
    """Recurse into value for drop_empty, leaving scalars alone."""
    if isinstance(value, dict):
        return drop_empty(value)
    if isinstance(value, list):
        return [_prune(item) for item in value]
    return value


def drop_empty(record: dict) -> dict:
    """Strip keys whose value is None or an empty list/dict/string, recursively.

    add_metadata_score in the hub's utils/validate.py tests `field in doc`, not
    whether the field holds anything, so a key left at None or [] scores exactly
    like a populated one. An organism with no dated or geolocated records would
    otherwise claim dateCreated, dateModified, datePublished, species,
    spatialCoverage and temporalCoverage -- six of the twenty-one recommended
    Dataset fields.

    0 and False are values, not emptiness, so they are kept.
    """
    cleaned = {}
    for key, value in record.items():
        value = _prune(value)
        if value is None or value == [] or value == {} or value == "":
            continue
        cleaned[key] = value
    return cleaned


def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(val, ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt


def version_to_isodate(data_version, default=None):
    """Convert a LAPIS dataVersion to an ISO date, or return default.

    dataVersion is Unix epoch seconds as a string, e.g. "1786019862" ->
    "2026-08-06". Read as UTC, since that is what LAPIS reports in. Anything
    missing or unparseable falls back to default rather than raising.
    """
    try:
        moment = datetime.datetime.fromtimestamp(int(data_version), datetime.timezone.utc)
    except (TypeError, ValueError, OverflowError, OSError):
        return default
    return moment.date().isoformat()


# ---------------------------------------------------------------------------
# LAPIS access
# ---------------------------------------------------------------------------


@cache
def get_organisms():
    """Return every organism slug, read from the backend's OpenAPI spec.

    The spec declares them as an enum on components.schemas.Organism, which is
    the same list the website builds its organism picker from. Cached, since
    it only changes when Pathoplexus adds an organism.
    """
    resp = requests.get(f"{BACKEND_URL}/api-docs", timeout=120)
    resp.raise_for_status()
    return resp.json()["components"]["schemas"]["Organism"]["enum"]


@cache
def get_database_config(organism):
    """Return the organism's schema: display name, openness, field definitions."""
    resp = requests.get(f"{LAPIS_URL}/{organism}/sample/databaseConfig", timeout=120)
    resp.raise_for_status()
    return resp.json()["schema"]


@cache
def get_schema_fields(organism):
    """Field names this organism actually has.

    Required before querying, because LAPIS answers 400 Unknown field rather
    than an empty result: dengue has insdcAccessionBase while cchf has only
    insdcAccessionBase_L/_M/_S, and serotype exists on dengue alone.
    """
    return frozenset(f["name"] for f in get_database_config(organism)["metadata"])


def get_collection_size(organism, filters=None):
    """Return {"collectionSize": n, "dataVersion": ...} for an organism.

    dataVersion identifies the data snapshot the count came from. It is per
    organism and changes whenever Pathoplexus reloads that organism's data,
    so two calls returning different values are not comparable.
    """
    resp = requests.post(
        f"{LAPIS_URL}/{organism}/sample/aggregated",
        json={**CURRENT_ONLY, **(filters or {})},
        timeout=120,
    )
    resp.raise_for_status()
    payload = resp.json()
    return {
        "collectionSize": payload["data"][0]["count"],
        "dataVersion": payload["info"]["dataVersion"],
    }


def value_counts(organism, field, filters=None, top_n=None):
    """Return {value: count} for one metadata field, most common first.

    Records where the field is unset are excluded, so the result never has a
    None key. Pass {f"{field}.isNull": True} in filters to count them instead.
    """
    resp = requests.post(
        f"{LAPIS_URL}/{organism}/sample/aggregated",
        json={
            **CURRENT_ONLY,
            f"{field}.isNull": False,
            **(filters or {}),
            "fields": [field],
        },
        timeout=120,
    )
    resp.raise_for_status()
    rows = sorted(resp.json()["data"], key=lambda r: r["count"], reverse=True)
    return {r[field]: r["count"] for r in rows[:top_n]}


@cache
def field_values(organism, field, top_n=MAX_VALUES):
    """Distinct values of one field, most frequent first.

    Cached because output keys overlap: sequencingInstrument feeds both
    measurementTechnique and sample.aggregateElement.instrument, and
    gcaAccession feeds both identifier and sameAs.
    """
    return tuple(value_counts(organism, field, None, top_n))


@cache
def field_pairs(organism, field_a, field_b, top_n=MAX_VALUES):
    """Co-occurring (a, b) values.

    One aggregated call over two fields returns real combinations rather than
    two independent tallies, which is what keeps a name attached to its own
    identifier. Returns () unless the organism has both fields.
    """
    available = get_schema_fields(organism)
    if field_a not in available or field_b not in available:
        return ()
    resp = requests.post(
        f"{LAPIS_URL}/{organism}/sample/aggregated",
        json={
            **CURRENT_ONLY,
            f"{field_a}.isNull": False,
            f"{field_b}.isNull": False,
            "fields": [field_a, field_b],
        },
        timeout=120,
    )
    resp.raise_for_status()
    rows = sorted(resp.json()["data"], key=lambda r: r["count"], reverse=True)
    return tuple((r[field_a], r[field_b]) for r in rows[:top_n])


def date_range(organism, field="sampleCollectionDate", filters=None):
    """Return (earliest, latest) for a date field, ignoring records with no date.

    Dates are strings and may be year-only ("1905"), so they sort lexically.
    """
    if field not in get_schema_fields(organism):
        return (None, None)
    bounds = []
    for direction in ("ascending", "descending"):
        resp = requests.post(
            f"{LAPIS_URL}/{organism}/sample/details",
            json={
                **CURRENT_ONLY,
                **(filters or {}),
                f"{field}.isNull": False,
                "fields": [field],
                "orderBy": [{"field": field, "type": direction}],
                "limit": 1,
            },
            timeout=120,
        )
        resp.raise_for_status()
        rows = resp.json()["data"]
        bounds.append(rows[0][field] if rows else None)
    return tuple(bounds)


def collect_values(organism, fields, top_n=MAX_VALUES):
    """De-duped, split union of several metadata fields, alphabetically.

    Fields the organism does not have are skipped rather than queried.
    """
    available = get_schema_fields(organism)
    values = set()
    for field in fields:
        if field in available:
            add_split_values(values, field_values(organism, field, top_n), max_size=top_n)
    return sorted(values)


# ---------------------------------------------------------------------------
# Record sections
# ---------------------------------------------------------------------------


def build_authors(organism):
    """author: the standing organisations, then one per contributing lab."""
    authors = [dict(author) for author in BASE_AUTHORS]
    for name in collect_values(organism, AUTHOR_NAME_FIELDS):
        insert_value_list(authors, {"@type": "Organization", "name": name})
    return authors


def insert_value_list(target, value):
    """Append value to target unless an equal entry is already present."""
    if value not in target:
        target.append(value)


def build_has_part(organism):
    """hasPart: one entry per accession, with url and sameAs where they apply."""
    with_url = set(collect_values(organism, HAS_PART_URL_FIELDS))
    with_same_as = set(collect_values(organism, HAS_PART_SAME_AS_FIELDS))

    parts = []
    for accession in collect_values(organism, HAS_PART_IDENTIFIER_FIELDS):
        part = {"@type": "CreativeWork", "identifier": accession}
        if accession in with_url:
            part["url"] = (
                f"https://pathoplexus.org/{organism}/search?selectedSeq={accession}"
            )
        if accession in with_same_as:
            part["sameAs"] = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
        parts.append(part)
    return parts


def build_sample(organism):
    """One sample object; each aggregateElement key below holds the list."""
    element = {}
    for key, fields in SAMPLE_KEYWORD_FIELDS.items():
        values = collect_values(organism, fields)
        if values:
            element[key] = values
    for key, fields in SAMPLE_OBJECT_FIELDS.items():
        values = collect_values(organism, fields)
        if values:
            element[key] = [{"name": value} for value in values]

    for key, fields in SAMPLE_TEXT_FIELDS.items():
        values = collect_values(organism, fields)
        if values:
            element[key] = values
    for key, fields in SAMPLE_DATE_FIELDS.items():
        values = collect_values(organism, fields)
        if values:
            element[key] = values

    items = [
        {"@type": "BioSample", "identifier": accession,
         "url": f"https://www.ebi.ac.uk/biosamples/samples/{accession}"}
        for accession in collect_values(organism, SAMPLE_ITEM_FIELDS)
    ]

    if not element and not items:
        return None
    sample: dict = {"@type": "SampleCollection"}
    if element:
        # aggregateElement has no @type slot in nde.py, unlike its parent.
        sample["aggregateElement"] = element
    if items:
        sample["itemListElement"] = items
        sample["numberOfItems"] = {"@type": "QuantitativeValue", "value": len(items)}
    return sample


def build_spatial_coverage(organism):
    """spatialCoverage: one entry per distinct place name."""
    return [
        {"@type": "Place", "name": name}
        for name in collect_values(organism, SPATIAL_COVERAGE_FIELDS)
    ]


def build_named_pair(organism, pair):
    """Names paired with their identifier, then any name-only leftovers."""
    name_field, id_field, fallback_fields = pair

    entries = []
    seen = set()
    for name, identifier in field_pairs(organism, name_field, id_field):
        if name in seen:
            continue
        seen.add(name)
        entries.append({"@type": "DefinedTerm", "name": name, "identifier": identifier})

    for name in collect_values(organism, (name_field, *fallback_fields)):
        if name not in seen:
            seen.add(name)
            entries.append({"@type": "DefinedTerm", "name": name})
    return entries


def build_mapped_lists(organism, output):
    """Write the table-driven keys onto output, skipping any that come back empty.

    STRING_LIST_FIELDS keys are typed as keyword/text in nde.py and stay bare
    strings; NAMED_OBJECT_FIELDS keys are typed as objects and get wrapped.
    """
    for key, fields in STRING_LIST_FIELDS.items():
        values = collect_values(organism, fields)
        if values:
            insert_value(output, key, values)

    for key, (type_name, fields) in NAMED_OBJECT_FIELDS.items():
        values = collect_values(organism, fields)
        if values:
            insert_value(output, key, [{"@type": type_name, "name": value} for value in values])

    projects = collect_values(organism, IS_PART_OF_FIELDS)
    if projects:
        insert_value(
            output,
            "isPartOf",
            [{"@type": "CreativeWork", "identifier": p} for p in projects],
        )


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def get_organism_info(organism, filters=None, top_n=MAX_VALUES):
    """Collection-level facts: name, size, snapshot version, and date bounds."""
    schema = get_database_config(organism)
    date_created, date_modified = date_range(organism, filters=filters)
    published, _ = date_range(organism, "earliestReleaseDate", filters)
    collection = get_collection_size(organism, filters)
    return {
        "organism": organism,
        "name": schema["instanceName"],
        "collectionSize": collection["collectionSize"],
        # "version": collection["dataVersion"],
        "dateCreated": date_created,
        "dateModified": date_modified,
        "datePublished": published,
        "countries": value_counts(organism, "geoLocCountry", filters, top_n),
        "data_use_terms": value_counts(organism, "dataUseTerms", filters, top_n),
    }


def build_record(organism, info=None):
    """Build one NDE DataCollection record for a single organism."""
    if info is None:
        info = get_organism_info(organism)

    name = info.get("name")
    data_use_terms = info.get("data_use_terms", {})
    open_count = data_use_terms.get("OPEN", 0)
    restricted_count = data_use_terms.get("RESTRICTED", 0)
    # value_counts drops nulls but not placeholder strings, and these go
    # straight into the description sentence, so filter them here too. That can
    # empty the list -- an organism whose only geoLocCountry value is "unknown"
    # -- so the clause below is omitted rather than left dangling.
    countries = [c for c in info.get("countries", {}) if not is_uninformative(c)]
    origin = f" collected from {', '.join(countries)}" if countries else ""

    output = {
        "_id": "pathoplexus_" + organism,
        "@type": "DataCollection",
        "includedInDataCatalog": {
            "@type": "DataCatalog",
            "name": "Pathoplexus",
            "url": "https://pathoplexus.org",
            "versionDate": version_to_isodate(
                info.get("version"), datetime.date.today().isoformat()
            ),
            "archivedAt": f"https://pathoplexus.org/{organism}/search",
        },
        "name": f"{name} sequence records at Pathoplexus",
        "description": (
            f"Viral nucleotide sequence records for {name} available through "
            f"Pathoplexus. {name} viral sequences{origin}, including "
            f"{open_count} Open access sequences and {restricted_count} "
            f"embargoed sequences."
        ),
        "url": f"https://pathoplexus.org/{organism}/search",
        "dateCreated": _to_iso_date(info.get("dateCreated")),
        "dateModified": _to_iso_date(info.get("dateModified")),
        "datePublished": _to_iso_date(info.get("datePublished")),
        # nde.py types collectionSize as {@type, value, unitText, min/maxValue},
        # not a bare integer.
        "collectionSize": {
            "@type": "QuantitativeValue",
            "value": info.get("collectionSize"),
            "unitText": "sequences",
        },
        # "version": info.get("version"),
        "creativeWorkStatus": CURRENT_ONLY["versionStatus"],
        "conditionsOfAccess": "Open" if not restricted_count else "Varied",
        "license": "https://pathoplexus.org/about/terms-of-use/terms-of-service",
        "usageInfo": USAGE_INFO,
        "creditText": "To cite data from Pathoplexus, please visit https://pathoplexus.org/about/citation",
        "temporalCoverage": {
            "@type": "TemporalInterval",
            "startDate": _to_iso_date(info.get("dateCreated")),
            "endDate": _to_iso_date(info.get("dateModified")),
            "temporalType": "collection",
        },
        "about": ABOUT,
        "exampleOfWork": EXAMPLE_OF_WORK,
        "topicCategory": TOPIC_CATEGORIES,
        "isBasedOn": [IS_BASED_ON],
        "author": build_authors(organism),
        "hasPart": build_has_part(organism),
        "species": build_named_pair(organism, SPECIES_PAIR),
        "infectiousAgent": build_named_pair(organism, INFECTIOUS_AGENT_PAIR)
        or [{"@type": "DefinedTerm", "name": name}],
        "spatialCoverage": build_spatial_coverage(organism),
    }

    build_mapped_lists(organism, output)

    sample = build_sample(organism)
    if sample:
        output["sample"] = sample

    for identifier in collect_values(organism, IS_BASED_ON_IDENTIFIER_FIELDS):
        insert_value_list(output["isBasedOn"], {"@type": "CreativeWork", "identifier": identifier})

    pruned: dict = drop_empty(output)

    # Dropping the interval's empty dates leaves {@type, temporalType} behind,
    interval = pruned.get("temporalCoverage") or {}
    if not interval.get("startDate") and not interval.get("endDate"):
        pruned.pop("temporalCoverage", None)

    return pruned


def parse(organisms=None):
    """Yield one DataCollection record per organism (all 14 by default)."""
    if organisms is None:
        organisms = get_organisms()
    for organism in organisms:
        logger.info(f"Building record for {organism}...")
        yield build_record(organism)

