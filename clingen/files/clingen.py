#!/usr/bin/env python3

import copy
import csv
import datetime as dt
import io
import re
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Iterator, Optional

GENE_VALIDITY_DEFAULT_URL = (
    "https://search.clinicalgenome.org/kb/gene-validity/download"
)
VARIANT_PATH_DEFAULT_URL = (
    "http://erepo.clinicalgenome.org/evrepo/api/classifications/all"
    "?format=tabbed"
)

SOURCE_METADATA = {
    "name": "ClinGen (as DataCollection records)",
    "url": "https://clinicalgenome.org/",
    "description": (
        "The Clinical Genomics Resource (ClinGen) is an NIH-sponsored "
        "resource focused on Genes with clinical relevance. The resource is "
        "highly curated and is made available with a public domain license."
    ),
    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
}

_DDE_CURATED_BY = {
    "name": "Data Discovery Engine",
    "url": "https://discovery.biothings.io/",
    "dateModified": "2026-01-03",
}

MEASUREMENT_TECHNIQUE_DEFINED_TERMS = [
    {
        "@type": "DefinedTerm",
        "alternateName": ["Curate", "Curation"],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Curation",
        "url": "http://purl.obolibrary.org/obo/NCIT_C48292",
    },
    {
        "@type": "DefinedTerm",
        "alternateName": ["Classification", "Classified"],
        "curatedBy": _DDE_CURATED_BY,
        "inDefinedTermSet": "NCIT",
        "isCurated": True,
        "name": "Classification",
        "url": "http://purl.obolibrary.org/obo/NCIT_C25161",
    },
]

VARIABLE_MEASURED_ANNOTATION = {
    "@type": "DefinedTerm",
    "alternateName": ["Annotation"],
    "curatedBy": _DDE_CURATED_BY,
    "inDefinedTermSet": "NCIT",
    "isCurated": True,
    "name": "Annotation",
    "url": "http://purl.obolibrary.org/obo/NCIT_C44272",
}

VARIABLE_MEASURED_PATHOGENICITY = {
    "@type": "DefinedTerm",
    "alternateName": ["Pathogenicity", "pathogenicity"],
    "curatedBy": _DDE_CURATED_BY,
    "inDefinedTermSet": "NCIT",
    "isCurated": True,
    "name": "Pathogenicity",
    "url": "http://purl.obolibrary.org/obo/NCIT_C168796",
}

VARIABLE_MEASURED_GVFIC = {
    "@type": "DefinedTerm",
    "alternateName": [
        "Genetic Variation Functional Impact Classification",
        "VARIANT IMPACT CLASSIFICATION",
    ],
    "curatedBy": _DDE_CURATED_BY,
    "inDefinedTermSet": "NCIT",
    "isCurated": True,
    "name": "Genetic Variation Functional Impact Classification",
    "url": "http://purl.obolibrary.org/obo/NCIT_C181345",
}

VARIABLE_MEASURED_GENE_VALIDITY = [
    copy.deepcopy(VARIABLE_MEASURED_ANNOTATION),
    {
        "@type": "DefinedTerm",
        "name": "Inheritance Mode",
        "url": "http://purl.obolibrary.org/obo/NCIT_C45827",
        "inDefinedTermSet": "NCIT",
    },
]

VARIABLE_MEASURED_VARIANT_PATH = [
    copy.deepcopy(VARIABLE_MEASURED_PATHOGENICITY),
    copy.deepcopy(VARIABLE_MEASURED_GVFIC),
]

ENCODING_FORMAT_DEFINED_TERMS = {
    "CSV": {
        "@type": "DefinedTerm",
        "name": "CSV",
        "url": "http://edamontology.org/format_3752",
        "inDefinedTermSet": "EDAM",
    },
    "JSON": {
        "@type": "DefinedTerm",
        "name": "JSON",
        "url": "http://edamontology.org/format_3464",
        "inDefinedTermSet": "EDAM",
    },
    "XML": {
        "@type": "DefinedTerm",
        "name": "XML",
        "url": "http://edamontology.org/format_2332",
        "inDefinedTermSet": "EDAM",
    },
    "TXT": {
        "@type": "DefinedTerm",
        "name": "TXT",
        "url": "http://edamontology.org/format_2330",
        "inDefinedTermSet": "EDAM",
    },
    "SQL": {
        "@type": "DefinedTerm",
        "name": "SQL",
        "url": "http://edamontology.org/format_3788",
        "inDefinedTermSet": "EDAM",
    },
    "MS-EXCEL": {
        "@type": "DefinedTerm",
        "name": "MS-Excel",
        "url": "http://edamontology.org/format_3620",
        "inDefinedTermSet": "EDAM",
    },
}

GENE_EXAMPLE_ADDITIONAL_PROPERTY = [
    {"name": "GENE SYMBOL", "value": "AARS1", "@type": "PropertyValue"},
    {"name": "GENE ID (HGNC)", "value": "HGNC:20", "@type": "PropertyValue"},
    {
        "name": "DISEASE LABEL",
        "value": "Charcot-Marie-Tooth disease axonal type 2N",
        "@type": "PropertyValue",
    },
    {
        "name": "DISEASE ID (MONDO)",
        "value": "MONDO:0013212",
        "@type": "PropertyValue",
    },
    {"name": "MOI", "value": "AD", "@type": "PropertyValue"},
    {"name": "SOP", "value": "SOP10", "@type": "PropertyValue"},
    {
        "name": "CLASSIFICATION",
        "value": "Definitive",
        "@type": "PropertyValue",
    },
    {
        "name": "ONLINE REPORT",
        "value": (
            "https://search.clinicalgenome.org/kb/gene-validity/"
            "CGGV:assertion_92de3832-c272-4993-8586-"
            "288c6331dec2-2024-03-14T160000.000Z"
        ),
        "@type": "PropertyValue",
    },
    {
        "name": "CLASSIFICATION DATE",
        "value": "2024-03-14T16:00:00.000Z",
        "@type": "PropertyValue",
    },
    {
        "name": "GCEP",
        "value": "Charcot-Marie-Tooth Disease Gene Curation Expert Panel",
        "@type": "PropertyValue",
    },
]

VARIANT_EXAMPLE_ADDITIONAL_PROPERTY = [
    {
        "name": "#Variation",
        "value": "NM_000277.2(PAH):c.1A>G (p.Met1Val)",
        "@type": "PropertyValue",
    },
    {"name": "ClinVar Variation Id", "value": "586", "@type": "PropertyValue"},
    {
        "name": "Allele Registry Id",
        "value": "CA114360",
        "@type": "PropertyValue",
    },
    {
        "name": "HGVS Expressions",
        "value": (
            "NM_000277.2:c.1A>G, NC_000012.12:g.102917130T>C, "
            "CM000674.2:g.102917130T>C, NC_000012.11:g.103310908T>C, "
            "CM000674.1:g.103310908T>C, NC_000012.10:g.101835038T>C, "
            "NG_008690.1:g.5473A>G, NG_008690.2:g.46281A>G, "
            "NM_000277.1:c.1A>G, XM_011538422.1:c.1A>G, "
            "NM_001354304.1:c.1A>G, XM_017019370.2:c.1A>G, "
            "NM_000277.3:c.1A>G, ENST00000307000.7:c.-147A>G, "
            "ENST00000546844.1:c.1A>G, ENST00000547319.1:n.312A>G, "
            "ENST00000549111.5:n.97A>G, ENST00000551337.5:c.1A>G, "
            "ENST00000551988.5:n.90A>G, ENST00000553106.5:c.1A>G, "
            "ENST00000635500.1:n.29-4232A>G, "
            "NM_000277.2(PAH):c.1A>G (p.Met1Val)"
        ),
        "@type": "PropertyValue",
    },
    {"name": "HGNC Gene Symbol", "value": "PAH", "@type": "PropertyValue"},
    {"name": "Disease", "value": "phenylketonuria", "@type": "PropertyValue"},
    {"name": "Mondo Id", "value": "MONDO:0009861", "@type": "PropertyValue"},
    {
        "name": "Mode of Inheritance",
        "value": "Autosomal recessive inheritance",
        "@type": "PropertyValue",
    },
    {"name": "Assertion", "value": "Pathogenic", "@type": "PropertyValue"},
    {
        "name": "Applied Evidence Codes (Met)",
        "value": "PS3, PM3, PP4_Moderate, PM2",
        "@type": "PropertyValue",
    },
    {
        "name": "Applied Evidence Codes (Not Met)",
        "value": "PVS1",
        "@type": "PropertyValue",
    },
    {
        "name": "Summary of interpretation",
        "value": (
            "PAH-specific ACMG/AMP criteria applied: PM2: gnomAD MAF=0.00002; "
            "PP4_Moderate: Seen in PKU patients. BH4 disorders ruled out. "
            "(PMID:2574002); PS3: <3% (PMID:9450897). PM3: Detected in trans "
            "with known pathogenic variants. In summary this variant meets "
            "criteria to be classified as pathogenic for phenylketonuria in "
            "an autosomal recessive manner based on the ACMG/AMP criteria "
            "applied as specified by the PAH Expert Panel: "
            "(PM2, PM3, PP4_Moderate, PS3). "
            "Updated to reflect new PVS1 recommendations."
        ),
        "@type": "PropertyValue",
    },
    {
        "name": "PubMed Articles",
        "value": "9450897, 2574002, 2574002",
        "@type": "PropertyValue",
    },
    {
        "name": "Expert Panel",
        "value": "Phenylketonuria VCEP",
        "@type": "PropertyValue",
    },
    {"name": "Guideline", "value": None, "@type": "PropertyValue"},
    {"name": "Approval Date", "value": "3/23/2019", "@type": "PropertyValue"},
    {"name": "Published Date", "value": "5/10/2019", "@type": "PropertyValue"},
    {"name": "Retracted", "value": "FALSE", "@type": "PropertyValue"},
    {
        "name": "Evidence Repo Link",
        "value": (
            "https://erepo.genome.network/evrepo/ui/classification/"
            "CA114360/MONDO:0009861/006"
        ),
        "@type": "PropertyValue",
    },
    {
        "name": "Uuid",
        "value": "89f04437-ed5d-4735-8c4a-a9b1d91d10ea",
        "@type": "PropertyValue",
    },
]

GENE_HOWTO_STEPS = [
    "1. Download the Gene-Disease Validity Report at "
    "https://search.clinicalgenome.org/kb/gene-validity/download",
    "2. Group the records by MONDO Disease IDs and get counts "
    "(collectionSize) of the genes",
    "3. Generate the url to link back to the Gene-Disease Validity report "
    "and fill in the healthCondition value for each MONDO Disease ID",
    "4. Parse the earliest and latest curation dates for dateCreated and "
    "dateModified",
    "5. Parse the GCEPs for the author",
    "6. Generate the description based on the templated text",
    "7. Fill in the values from manually curated repository-level metadata "
    "for properties such as measurementTechnique, species, usageInfo, "
    "license, variableMeasured, topicCategory, isBasedOn, exampleOfWork, "
    "isAccessibleForFree",
]

VARIANT_HOWTO_STEPS = [
    "1. Download the Variant Pathogenicity Evidence Report at "
    "http://erepo.clinicalgenome.org/evrepo/api/classifications/all?"
    "format=tabbed",
    "2. Group the records by MONDO Disease IDs and get counts "
    "(collectionSize) of the variants",
    "3. Generate the url to link back to the Variant Pathogenicity Evidence "
    "report and fill in the healthCondition value for each MONDO Disease ID",
    "4. Parse the earliest and latest publish dates for dateCreated and "
    "dateModified",
    "5. Parse the ExpertPanel for the author",
    "6. Generate the description based on the templated text",
    "7. Fill in the values from manually curated repository-level metadata "
    "for properties such as measurementTechnique, species, usageInfo, "
    "license, variableMeasured, topicCategory, isBasedOn, exampleOfWork, "
    "isAccessibleForFree",
]

MODE_SETTINGS: dict[str, dict[str, Any]] = {
    "gene-validity": {
        "about": {
            "@type": "DefinedTerm",
            "description": "For schema, consider https://schema.org/Gene",
            "displayName": "Gene",
            "name": "Gene",
            "url": "http://purl.obolibrary.org/obo/NCIT_C16612",
            "inDefinedTermSet": "NCIT",
        },
        "name_suffix": "Gene-Disease Validity curation",
        "description_template": (
            "Genes associated with {healthCondition.name} "
            "({healthCondition.identifier}) and classifications on the "
            "validity "
            "of these associations by the {author.name}, a Gene Curation "
            "Expert Panel (GCEP). GCEPs from the Clinical Genomics Resource "
            "(ClinGen). These Gene-Disease Validity curations follow "
            "documented SOPs for the Gene-Disease validity classifications. "
            "In addition to Gene-Disease Validation, ClinGen also offers "
            "Variant Pathogenicity evaluations, Clinical Actionability "
            "reports "
            "for diseases, and Gene Dosage Sensitivity reports. For more "
            "details, visit: {url}"
        ),
        "url_base": "https://search.clinicalgenome.org/kb/conditions/",
        "sop": {
            "@type": "CreativeWork",
            "name": (
                "Gene-Disease Validity Curation Standard Operating "
                "Procedures"
            ),
            "url": "https://clinicalgenome.org/docs/"
            "gene-disease-validity-standard-operating-procedure/",
        },
        "howto_steps": GENE_HOWTO_STEPS,
        "variable_measured": VARIABLE_MEASURED_GENE_VALIDITY,
        "encoding_formats": ["CSV"],
        "example_additional_property": GENE_EXAMPLE_ADDITIONAL_PROPERTY,
        "topic_category": [
            "http://edamontology.org/topic_0622",
            "http://edamontology.org/topic_0199",
            "http://edamontology.org/topic_0625",
        ],
    },
    "variant-pathogenicity": {
        "about": {
            "@type": "DefinedTerm",
            "description": "Subclass of Gene",
            "displayName": "Gene Variant",
            "name": "GeneVariant",
            "url": "http://purl.obolibrary.org/obo/NCIT_C97927",
            "inDefinedTermSet": "NCIT",
        },
        "name_suffix": "Variant Pathogenicity Curation",
        "description_template": (
            "Variant Pathogenicity classifications associated with "
            "{healthCondition.name} ({healthCondition.identifier}) by the "
            "{author.name}, a Variant Curation Expert Panel (VCEP). VCEPs "
            "from the Clinical Genomics Resource (ClinGen). These Variant "
            "Pathogenicity curations follow documented SOPs. In addition to "
            "Variant Pathogenicity evaluations, ClinGen also offers "
            "Gene-Disease Validation, Clinical Actionability reports for "
            "diseases, and Gene Dosage Sensitivity reports. For more details, "
            "visit: {url}"
        ),
        "sop": {
            "@type": "CreativeWork",
            "name": (
                "Variant Pathogenicity Curation Standard Operating "
                "Procedures"
            ),
            "url": "https://clinicalgenome.org/docs/"
            "variant-curation-standard-operating-procedure/",
        },
        "howto_steps": VARIANT_HOWTO_STEPS,
        "variable_measured": VARIABLE_MEASURED_VARIANT_PATH,
        "encoding_formats": ["JSON", "XML", "CSV", "TXT", "SQL", "MS-EXCEL"],
        "example_additional_property": VARIANT_EXAMPLE_ADDITIONAL_PROPERTY,
    },
}


def _strip(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _capitalize_sentence_start(text: str) -> str:
    s = _strip(text)
    if not s:
        return s
    first = s[0]
    if first.isalpha() and first.islower():
        return first.upper() + s[1:]
    return s


def _parse_datetime(value: str) -> Optional[dt.datetime]:
    v = value.strip()
    if not v:
        return None

    try:
        if v.endswith("Z"):
            v = v[:-1] + "+00:00"
        parsed = dt.datetime.fromisoformat(v)
        if parsed.tzinfo is None:
            return parsed.replace(tzinfo=dt.timezone.utc)
        return parsed.astimezone(dt.timezone.utc)
    except Exception:
        return None


def _iso_date(value: Optional[dt.datetime]) -> Optional[str]:
    if value is None:
        return None
    return value.date().isoformat()


def _make_record_id(prefix: str, mondo_id: str) -> str:
    raw = f"{_strip(prefix)}_{_strip(mondo_id)}"
    raw = raw.replace(":", "_")
    raw = re.sub(r"[^A-Za-z0-9_]+", "_", raw)
    return re.sub(r"_+", "_", raw).strip("_")


def _apply_template_replacements(
    value: Any,
    replacements: dict[str, str],
) -> Any:
    if isinstance(value, str):
        out = value
        for key, replacement in replacements.items():
            out = out.replace(key, replacement)
        return out

    if isinstance(value, list):
        return [
            _apply_template_replacements(item, replacements)
            for item in value
        ]

    if isinstance(value, dict):
        return {
            key: _apply_template_replacements(item, replacements)
            for key, item in value.items()
        }

    return value


def _org(name: str) -> dict[str, str]:
    return {"@type": "Organization", "name": name}


def _funding_object() -> dict[str, list[dict[str, str]]]:
    return {
        "identifier": [
            {"@type": "PropertyValue", "value": "U24HG009649"},
            {"@type": "PropertyValue", "value": "U24HG009650"},
            {"@type": "PropertyValue", "value": "U24HG006834"},
        ]
    }


def _open_text_stream(
    input_path: Optional[Path],
    input_url: Optional[str],
) -> io.TextIOBase:
    if input_path is not None:
        return input_path.open("r", encoding="utf-8", newline="")

    if input_url is None:
        raise ValueError("Provide --input or --input-url")

    req = urllib.request.Request(
        input_url,
        headers={"User-Agent": "nde-clingen-parser/0.1"},
    )
    response = urllib.request.urlopen(req)
    return io.TextIOWrapper(response, encoding="utf-8", newline="")


def _iter_gene_validity_rows(
    stream: io.TextIOBase,
    limit: Optional[int],
) -> Iterator[dict[str, str]]:
    while True:
        pos = stream.tell() if stream.seekable() else None
        line = stream.readline()
        if not line:
            return

        if "GENE SYMBOL" in line and "DISEASE ID" in line:
            if pos is not None:
                stream.seek(pos)
            else:
                stream = io.StringIO(line + stream.read())
            break

    reader = csv.DictReader(stream)
    count = 0
    for row in reader:
        cleaned = {
            key: (value or "").strip()
            for key, value in row.items()
            if key is not None
        }
        if cleaned.get("GENE SYMBOL", "").startswith("++"):
            continue
        if cleaned.get("DISEASE ID (MONDO)", "").startswith("++"):
            continue

        count += 1
        if limit is not None and count > limit:
            break
        yield cleaned


def _iter_variant_rows(
    stream: io.TextIOBase,
    limit: Optional[int],
) -> Iterator[dict[str, str]]:
    reader = csv.DictReader(stream, delimiter="\t")
    for index, row in enumerate(reader, start=1):
        if limit is not None and index > limit:
            break

        yield {
            key.lstrip("#"): (value or "").strip()
            for key, value in row.items()
            if key is not None
        }


def _build_example_of_work(mode: str) -> dict[str, Any]:
    settings = MODE_SETTINGS[mode]
    encoding_formats: list[dict[str, Any]] = []
    for label in settings["encoding_formats"]:
        term = ENCODING_FORMAT_DEFINED_TERMS.get(label)
        if term is not None:
            encoding_formats.append(copy.deepcopy(term))

    return {
        "@type": "CreativeWork",
        "about": {
            "name": "Genotype/Phenotype Annotation",
            "url": "http://edamontology.org/data_0920",
            "@type": "DefinedTerm",
            "inDefinedTermSet": "NCIT",
        },
        "encodingFormat": encoding_formats,
        "additionalProperty": copy.deepcopy(
            settings["example_additional_property"]
        ),
    }


def _build_is_based_on(
    mode: str,
    *,
    about_name: str,
    condition_name: str,
    mondo_id: str,
    author_name: str,
    url: str,
) -> list[dict[str, Any]]:
    settings = MODE_SETTINGS[mode]

    action_obj = {
        "@type": "Action",
        "name": (
            "DataCollection Generation Process in the NIAID Data Ecosystem"
        ),
        "description": (
            "How this ClinGen {about.name} DataCollection Record was "
            "generated "
            "for the NIAID Data Ecosystem."
        ),
        "actionProcess": {
            "@type": "HowTo",
            "step": copy.deepcopy(settings["howto_steps"]),
        },
    }

    source_obj = {
        "@type": "nde:ResourceCatalog",
        "name": "Clinical Genomics Resource",
        "url": "https://data.niaid.nih.gov/resources?id=dde_b880e7dc2c2437e0",
    }

    replacements = {
        "{about.name}": about_name,
        "{healthCondition.name}": condition_name,
        "{healthCondition.identifier}": mondo_id,
        "{author.name}": author_name,
        "{url}": url,
    }

    return [
        _apply_template_replacements(
            copy.deepcopy(settings["sop"]),
            replacements,
        ),
        _apply_template_replacements(action_obj, replacements),
        _apply_template_replacements(source_obj, replacements),
    ]


def _build_gene_validity_record(
    mondo: str,
    group: dict[str, Any],
) -> dict[str, Any]:
    settings = MODE_SETTINGS["gene-validity"]
    disease = group["disease"] or mondo
    disease_for_name = _capitalize_sentence_start(disease)

    authors = sorted(a for a in group["authors"] if a)
    if len(authors) == 1:
        author_name = authors[0]
    elif authors:
        author_name = "; ".join(authors)
    else:
        author_name = "ClinGen"

    dates: list[dt.datetime] = group["dates"]
    created = min(dates) if dates else None
    modified = max(dates) if dates else None

    url = settings["url_base"] + mondo

    replacements = {
        "{healthCondition.name}": disease,
        "{healthCondition.identifier}": mondo,
        "{author.name}": author_name,
        "{url}": url,
    }

    record: dict[str, Any] = {
        "_id": _make_record_id("clingen_gene_validity", mondo),
        "@type": "nde:DataCollection",
        "about": copy.deepcopy(settings["about"]),
        "includedInDataCatalog": "Clinical Genomics Resource (ClinGen)",
        "name": f"{disease_for_name} {settings['name_suffix']}".strip(),
        "url": url,
        "collectionSize": len(group["genes"]),
        "author": _org(author_name),
        "creator": _org("ClinGen"),
        "healthCondition": {
            "@type": "DefinedTerm",
            "name": disease,
            "identifier": mondo,
            "inDefinedTermSet": "MONDO",
        },
        "dateCreated": _iso_date(created),
        "dateModified": _iso_date(modified),
        "date": _iso_date(modified) or _iso_date(created),
        "description": _apply_template_replacements(
            settings["description_template"],
            replacements,
        ),
        "funding": _funding_object(),
        "conditionsOfAccess": "Open",
        "license": "https://creativecommons.org/publicdomain/zero/1.0/",
        "measurementTechnique": copy.deepcopy(
            MEASUREMENT_TECHNIQUE_DEFINED_TERMS
        ),
        "species": "Homo Sapiens",
        "topicCategory": copy.deepcopy(settings["topic_category"]),
        "usageInfo": (
            "https://www.clinicalgenome.org/docs/?doc-type="
            "policies-position-statements"
        ),
        "variableMeasured": copy.deepcopy(settings["variable_measured"]),
        "isAccessibleForFree": True,
        "exampleOfWork": _build_example_of_work("gene-validity"),
        "isBasedOn": _build_is_based_on(
            "gene-validity",
            about_name=_strip(settings["about"].get("name")),
            condition_name=disease,
            mondo_id=mondo,
            author_name=author_name,
            url=url,
        ),
    }
    return record


def _build_variant_pathogenicity_record(
    mondo: str,
    group: dict[str, Any],
) -> dict[str, Any]:
    settings = MODE_SETTINGS["variant-pathogenicity"]
    disease = group["disease"] or mondo
    disease_for_name = _capitalize_sentence_start(disease)

    authors = sorted(a for a in group["authors"] if a)
    if len(authors) == 1:
        author_name = authors[0]
    elif authors:
        author_name = "; ".join(authors)
    else:
        author_name = "ClinGen"

    dates: list[dt.datetime] = group["dates"]
    created = min(dates) if dates else None
    modified = max(dates) if dates else None

    disease_q = urllib.parse.quote(disease)
    url = (
        "https://erepo.clinicalgenome.org/evrepo/ui/summary/"
        "classifications"
        f"?columns=condition&values={disease_q}"
        "&matchTypes=exact&pgSize=25&pg=1&matchMode=and"
    )

    replacements = {
        "{healthCondition.name}": disease,
        "{healthCondition.identifier}": mondo,
        "{author.name}": author_name,
        "{url}": url,
    }

    record: dict[str, Any] = {
        "_id": _make_record_id("clingen_variant_pathogenicity", mondo),
        "@type": "nde:DataCollection",
        "about": copy.deepcopy(settings["about"]),
        "includedInDataCatalog": "Clinical Genomics Resource (ClinGen)",
        "name": f"{disease_for_name} {settings['name_suffix']}".strip(),
        "url": url,
        "collectionSize": len(group["uuids"]),
        "author": _org(author_name),
        "creator": _org("ClinGen"),
        "healthCondition": {
            "@type": "DefinedTerm",
            "name": disease,
            "identifier": mondo,
            "inDefinedTermSet": "MONDO",
        },
        "dateCreated": _iso_date(created),
        "dateModified": _iso_date(modified),
        "date": _iso_date(modified) or _iso_date(created),
        "description": _apply_template_replacements(
            settings["description_template"],
            replacements,
        ),
        "funding": _funding_object(),
        "conditionsOfAccess": "Open",
        "license": "https://creativecommons.org/publicdomain/zero/1.0/",
        "measurementTechnique": copy.deepcopy(
            MEASUREMENT_TECHNIQUE_DEFINED_TERMS
        ),
        "species": "Homo Sapiens",
        "usageInfo": (
            "https://www.clinicalgenome.org/docs/?doc-type="
            "policies-position-statements"
        ),
        "variableMeasured": copy.deepcopy(settings["variable_measured"]),
        "isAccessibleForFree": True,
        "exampleOfWork": _build_example_of_work("variant-pathogenicity"),
        "isBasedOn": _build_is_based_on(
            "variant-pathogenicity",
            about_name=_strip(settings["about"].get("name")),
            condition_name=disease,
            mondo_id=mondo,
            author_name=author_name,
            url=url,
        ),
    }
    return record


def _parse_gene_validity_records(
    input_path: Optional[Path],
    input_url: Optional[str],
    limit: Optional[int],
) -> Iterator[dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}

    with _open_text_stream(input_path, input_url) as stream:
        for row in _iter_gene_validity_rows(stream, limit):
            mondo = row.get("DISEASE ID (MONDO)", "")
            if not mondo:
                continue

            disease = row.get("DISEASE LABEL", "")
            gene_symbol = row.get("GENE SYMBOL", "")
            gene_hgnc = row.get("GENE ID (HGNC)", "")
            gcep = row.get("GCEP", "")
            date_val = _parse_datetime(row.get("CLASSIFICATION DATE", ""))

            group = groups.setdefault(
                mondo,
                {
                    "disease": disease,
                    "genes": set(),
                    "authors": set(),
                    "dates": [],
                },
            )

            group["disease"] = disease or group["disease"]
            if gene_hgnc or gene_symbol:
                group["genes"].add((gene_hgnc or gene_symbol).strip())
            if gcep:
                group["authors"].add(gcep)
            if date_val is not None:
                group["dates"].append(date_val)

    for mondo in sorted(groups.keys()):
        yield _build_gene_validity_record(mondo, groups[mondo])


def _parse_variant_pathogenicity_records(
    input_path: Optional[Path],
    input_url: Optional[str],
    limit: Optional[int],
) -> Iterator[dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}

    with _open_text_stream(input_path, input_url) as stream:
        for row in _iter_variant_rows(stream, limit):
            mondo = row.get("Mondo Id", "")
            if not mondo:
                continue

            disease = row.get("Disease", "")
            uuid = row.get("Uuid", "")
            panel = row.get("Expert Panel", "")
            published = _parse_datetime(row.get("Published Date", ""))

            group = groups.setdefault(
                mondo,
                {
                    "disease": disease,
                    "uuids": set(),
                    "authors": set(),
                    "dates": [],
                },
            )

            group["disease"] = disease or group["disease"]
            if uuid:
                group["uuids"].add(uuid)
            if panel:
                group["authors"].add(panel)
            if published is not None:
                group["dates"].append(published)

    for mondo in sorted(groups.keys()):
        yield _build_variant_pathogenicity_record(mondo, groups[mondo])


def parse() -> Iterator[dict[str, Any]]:
    yield from _parse_gene_validity_records(
        None,
        GENE_VALIDITY_DEFAULT_URL,
        None,
    )
    yield from _parse_variant_pathogenicity_records(
        None,
        VARIANT_PATH_DEFAULT_URL,
        None,
    )
