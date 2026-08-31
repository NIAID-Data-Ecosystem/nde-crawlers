"""Per-document finishing touches applied to every record, by every source.

These run at the tail of the upload pipeline, after the augmentation stages:
clean the description, roll up `date`, score metadata completeness, drop
placeholder terms and validate the document before it reaches MongoDB.
"""

import datetime
import re
from typing import Dict

from config import logger
from lxml import etree, html
from scores import (
    COMPUTATIONAL_TOOL_RECOMMENDED,
    COMPUTATIONAL_TOOL_RECOMMENDED_AUGMENTED,
    COMPUTATIONAL_TOOL_REQUIRED,
    COMPUTATIONAL_TOOL_REQUIRED_AUGMENTED,
    DATA_COLLECTION_RECOMMENDED,
    DATA_COLLECTION_RECOMMENDED_AUGMENTED,
    DATA_COLLECTION_REQUIRED,
    DATA_COLLECTION_REQUIRED_AUGMENTED,
    DATASET_RECOMMENDED_AUGMENTED_FIELDS,
    DATASET_RECOMMENDED_FIELDS,
    DATASET_REQUIRED_AUGMENTED_FIELDS,
    DATASET_REQUIRED_FIELDS,
    RESOURCE_CATALOG_RECOMMENDED,
    RESOURCE_CATALOG_RECOMMENDED_AUGMENTED,
    RESOURCE_CATALOG_REQUIRED,
    RESOURCE_CATALOG_REQUIRED_AUGMENTED,
)

from .common import as_list

AUGMENTED_FLAGS = ("fromPMID", "fromGPT", "fromEXTRACT", "fromNCT")

CONDITIONS_OF_ACCESS = ["Open", "Restricted", "Closed", "Embargoed", "Varied"]
CREATIVE_WORK_STATUS = ["Bespoke", "Available", "Backordered", "Retired"]

# Placeholders sources emit as terms; nothing upstream strips them.
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

_BREAK_TAGS = re.compile(r"(<br\s*/?>|</p>)", re.IGNORECASE)

_SCORE_FIELDS = {
    "DataCollection": (
        DATA_COLLECTION_REQUIRED,
        DATA_COLLECTION_RECOMMENDED,
        DATA_COLLECTION_REQUIRED_AUGMENTED,
        DATA_COLLECTION_RECOMMENDED_AUGMENTED,
    ),
    "ComputationalTool": (
        COMPUTATIONAL_TOOL_REQUIRED,
        COMPUTATIONAL_TOOL_RECOMMENDED,
        COMPUTATIONAL_TOOL_REQUIRED_AUGMENTED,
        COMPUTATIONAL_TOOL_RECOMMENDED_AUGMENTED,
    ),
    "ResourceCatalog": (
        RESOURCE_CATALOG_REQUIRED,
        RESOURCE_CATALOG_RECOMMENDED,
        RESOURCE_CATALOG_REQUIRED_AUGMENTED,
        RESOURCE_CATALOG_RECOMMENDED_AUGMENTED,
    ),
}
_DEFAULT_SCORE_FIELDS = (
    DATASET_REQUIRED_FIELDS,
    DATASET_RECOMMENDED_FIELDS,
    DATASET_REQUIRED_AUGMENTED_FIELDS,
    DATASET_RECOMMENDED_AUGMENTED_FIELDS,
)


def clean_description(doc: Dict) -> Dict:
    """Flatten the description to a single string and strip its HTML markup."""
    description = doc.get("description")
    if not description:
        return doc

    try:
        doc["description"] = _strip_html(description)
    except etree.ParserError as e:
        # At minimum, prevent the lxml error object from escaping
        logger.warning("ParserError while processing doc %s: %s", doc.get("_id"), str(e))
    return doc


def _strip_html(description) -> str:
    """Join a string / bytes / list description into one string of plain text."""
    # Decode bytes so the regex (a str pattern) never runs against a bytes-like object
    text = " ".join(
        part.decode("utf-8", errors="replace") if isinstance(part, bytes) else part for part in as_list(description)
    )
    # Normalize line breaks before stripping the remaining HTML tags
    text = _BREAK_TAGS.sub("\n", text)
    # lxml cannot parse a str carrying an XML encoding declaration, so re-encode that case only
    if text.strip().startswith("<?xml"):
        text = text.encode("utf-8")
    return html.fromstring(text).text_content()


def drop_placeholder_terms(doc: Dict) -> Dict:
    """Strip placeholder entries from species / infectiousAgent / healthCondition."""
    for field in ("species", "infectiousAgent", "healthCondition"):
        if field not in doc:
            continue
        value = doc[field]
        if isinstance(value, dict):
            if _is_placeholder_term(value.get("name")):
                doc.pop(field, None)
        elif isinstance(value, list):
            kept = [entry for entry in value if not _is_placeholder_term(_entry_name(entry))]
            if kept:
                doc[field] = kept
            else:
                doc.pop(field, None)
    return doc


def _entry_name(entry):
    return entry.get("name") if isinstance(entry, dict) else None


def _is_placeholder_term(name):
    if not isinstance(name, str):
        return False
    return name.lower().strip() in PLACEHOLDER_TERMS


def add_date(doc: Dict) -> Dict:
    """Set `date` to the latest of date / dateCreated / dateModified / datePublished.

    `distribution.dateModified` is also considered, and is cloned up to
    `dateModified` when the document has none of its own.
    """
    distribution_dates = sorted(
        dist.get("dateModified")
        for dist in as_list(doc.get("distribution"))
        if isinstance(dist, dict) and dist.get("dateModified")
    )

    if not doc.get("dateModified") and distribution_dates:
        doc["dateModified"] = distribution_dates[-1]

    dates = [doc[field] for field in ("date", "dateCreated", "dateModified", "datePublished") if doc.get(field)]
    dates.extend(distribution_dates)
    if dates:
        dates.sort()
        doc["date"] = datetime.datetime.fromisoformat(dates[-1]).date().isoformat()

    return doc


def is_purely_augmented(field: str, field_content) -> bool:
    """True when `field_content` holds only augmented data, i.e. no original source content."""
    # Special case: includedInDataCatalog always considered "not purely augmented"
    if field == "includedInDataCatalog":
        return False
    # Strings/bools are not 'augmented' dict structures
    if isinstance(field_content, (str, bool)):
        return False

    if isinstance(field_content, dict):
        field_content = [field_content]
    if not isinstance(field_content, list):
        return False
    return all(
        isinstance(item, dict) and any(item.get(flag, False) for flag in AUGMENTED_FLAGS) for item in field_content
    )


def check_augmented_fields(doc: Dict, fields) -> list:
    """Return the subset of `fields` holding at least one augmented entry."""
    return [
        field
        for field in fields
        if any(
            isinstance(item, dict) and any(item.get(flag, False) for flag in AUGMENTED_FLAGS)
            for item in as_list(doc.get(field))
        )
    ]


def add_metadata_score(doc: Dict) -> Dict:
    """Score required/recommended field coverage into `_meta.completeness`.

    The field sets are chosen from the document's `@type`; types without a
    dedicated configuration are scored as a Dataset.
    """
    required, recommended, required_augmented, recommended_augmented = _SCORE_FIELDS.get(
        doc.get("@type"), _DEFAULT_SCORE_FIELDS
    )

    existing_required = [f for f in required if f in doc and not is_purely_augmented(f, doc[f])]
    existing_recommended = [f for f in recommended if f in doc and not is_purely_augmented(f, doc[f])]
    found_required_augmented = check_augmented_fields(doc, required_augmented)
    found_recommended_augmented = check_augmented_fields(doc, recommended_augmented)

    required_score = len(existing_required)
    recommended_score = len(existing_recommended)
    total_required = len(required)
    total_recommended = len(recommended)

    doc.setdefault("_meta", {}).update(
        {
            "required_augmented_fields": found_required_augmented,
            "recommended_augmented_fields": found_recommended_augmented,
            "required_fields": existing_required,
            "recommended_fields": existing_recommended,
            "completeness": {
                "total_score": required_score + recommended_score,
                "total_max_score": total_required + total_recommended,
                "required_score": required_score,
                "required_ratio": _ratio(required_score, total_required),
                "required_max_score": total_required,
                "recommended_score": recommended_score,
                "recommended_score_ratio": _ratio(recommended_score, total_recommended),
                "recommended_max_score": total_recommended,
                "augmented_required_ratio": _ratio(len(found_required_augmented), total_required),
                "augmented_recommended_ratio": _ratio(len(found_recommended_augmented), total_recommended),
                "total_required_augmented": len(required_augmented),
                "total_recommended_augmented": len(recommended_augmented),
            },
        }
    )
    return doc


def _ratio(score, total):
    return round(score / total, 2) if total > 0 else 0


def check_schema(doc: Dict) -> Dict:
    """Raise one error listing every NDE schema issue found in the document."""
    doc_id = doc.get("_id") if isinstance(doc, dict) else None
    issues = []

    def check(condition, message):
        if not condition:
            issues.append(message)

    if not isinstance(doc, dict):
        raise AssertionError(f"Schema validation failed [record _id={doc_id!r}]:\n- doc is not a dict")

    check(doc.get("_id"), "_id is None")
    check(doc.get("@type"), "@type is None")
    check(doc.get("url"), "url is None")

    none = object()

    def assert_types(field, expected_types, value=none):
        entries = as_list(doc.get(field) if value is none else value)
        check(
            all(isinstance(item, dict) and item.get("@type") in expected_types for item in entries),
            f"{field} needs to be of type {' or '.join(expected_types)}",
        )

    person_or_organization = ("Organization", "Person")
    for field in ("author", "creator"):
        assert_types(field, person_or_organization)

    term_fields = (
        "species",
        "infectiousAgent",
        "healthCondition",
        "measurementTechnique",
        "topicCategory",
    )
    for field in term_fields:
        assert_types(field, ("DefinedTerm",))

    citation_fields = ("citation", "citedBy", "isBasedOn", "isBasisFor", "isPartOf", "hasPart")
    for field in citation_fields:
        assert_types(field, ("ScholarlyArticle", "CreativeWork"))

    assert_types("funding", ("MonetaryGrant",))
    assert_types("sourceOrganization", ("Organization", "ResearchProject"))

    # Validate the objects nested inside the objects created by the utility stages.
    authors = as_list(doc.get("author")) + as_list(doc.get("creator"))
    for field in citation_fields:
        for citation in as_list(doc.get(field)):
            if isinstance(citation, dict):
                assert_types(f"{field}.author", person_or_organization, citation.get("author"))
                authors.extend(as_list(citation.get("author")))

    for author in authors:
        if isinstance(author, dict) and author.get("affiliation") is not None:
            assert_types("affiliation", ("Organization",), author.get("affiliation"))

    curator_types = (
        "Person",
        "Organization",
        "DataCatalog",
        "ResourceCatalog",
        "SoftwareApplication",
        "ComputationalTool",
        "ResearchProject",
    )
    for field in term_fields:
        for term in as_list(doc.get(field)):
            if isinstance(term, dict) and term.get("curatedBy") is not None:
                assert_types(f"{field}.curatedBy", curator_types, term.get("curatedBy"))

    for grant in as_list(doc.get("funding")):
        if not isinstance(grant, dict):
            continue
        if grant.get("funder") is not None:
            assert_types("funding.funder", ("Organization",), grant.get("funder"))
        for funder in as_list(grant.get("funder")):
            if isinstance(funder, dict) and funder.get("employee") is not None:
                assert_types("funding.funder.employee", ("Person",), funder.get("employee"))
        if grant.get("isBasedOn") is not None:
            assert_types("funding.isBasedOn", ("ScholarlyArticle", "CreativeWork"), grant.get("isBasedOn"))

    catalogs = doc.get("includedInDataCatalog")
    check(catalogs, "includedInDataCatalog is None")
    if catalogs:
        check(
            all(isinstance(item, dict) and item.get("archivedAt") for item in as_list(catalogs)),
            "includedInDataCatalog.archivedAt is None in one or more items",
        )

    check(doc.get("version", None) is None, "Remove version field")

    if coa := doc.get("conditionsOfAccess"):
        check(
            coa in CONDITIONS_OF_ACCESS,
            "%s is not a valid conditionsOfAccess. Allowed conditionsOfAccess: %s" % (coa, CONDITIONS_OF_ACCESS),
        )

    if doc.get("@type") == "Sample" and (cws := doc.get("creativeWorkStatus")) is not None:
        if isinstance(cws, (list, tuple, set)):
            check(cws, "creativeWorkStatus cannot be empty")
            invalid = [status for status in cws if status not in CREATIVE_WORK_STATUS]
            check(
                not invalid,
                "%s is not a valid creativeWorkStatus. Allowed creativeWorkStatus: %s"
                % (cws, CREATIVE_WORK_STATUS),
            )
        else:
            check(
                cws in CREATIVE_WORK_STATUS,
                "%s is not a valid creativeWorkStatus. Allowed creativeWorkStatus: %s"
                % (cws, CREATIVE_WORK_STATUS),
            )

    if issues:
        issue_list = "\n".join(f"- {issue}" for issue in issues)
        raise AssertionError(f"Schema validation failed [record _id={doc_id!r}]:\n{issue_list}")

    return doc
