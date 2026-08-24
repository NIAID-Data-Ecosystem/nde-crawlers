#!/usr/bin/env python3
"""PGS Catalog MONDO-trait DataCollection crawler for the NDE."""

import copy
import datetime
import hashlib
import json
import logging
import os
import re
import time
from typing import Any, Iterable, Iterator, Optional

import dateutil.parser
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger("nde-logger")


# fmt: off
def insert_value(d, key, value, extend=False):
    """ Insert a value into a dictionary, handling existing keys by converting to lists or extending strings as needed.
    """

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

def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(str(val), ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError, OverflowError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt
# fmt: on


API_BASE = "https://www.pgscatalog.org/rest"
WEB_BASE = "https://www.pgscatalog.org"
USER_AGENT = "nde-pgs-catalog-datacollection-crawler/0.1"
REQUEST_TIMEOUT = int(os.environ.get("PGS_CATALOG_REQUEST_TIMEOUT", "90"))
PAGE_SIZE = int(os.environ.get("PGS_CATALOG_PAGE_SIZE", "50"))
SCORE_CHUNK_SIZE = int(os.environ.get("PGS_CATALOG_SCORE_CHUNK_SIZE", "50"))
REQUEST_INTERVAL = float(os.environ.get("PGS_CATALOG_REQUEST_INTERVAL", "0.65"))

CONTEXT = "http://schema.org/"

AUTHORS = [
    {
        "@type": "Organization",
        "name": "European Bioinformatics Institute",
        "url": "https://www.ebi.ac.uk/",
    },
    {
        "@type": "Organization",
        "name": "University of Cambridge",
        "url": "https://www.cam.ac.uk/",
    },
]

SOURCE_ORGANIZATION = {
    "@type": "Organization",
    "name": "Polygenic Score Catalog",
    "alternateName": "PGS Catalog",
    "parentOrganization": "EMBL-EBI and University of Cambridge",
    "url": "https://www.pgscatalog.org/",
}

ABOUT = {
    "@type": "DefinedTerm",
    "name": "MedicalRiskScore",
    "displayName": "Medical Risk Score",
    "description": "For the schema, consider https://schema.org/MedicalRiskScore",
    "url": "http://purl.obolibrary.org/obo/NCIT_C15367",
}

POLYGENIC_RISK_SCORE = {
    "@type": "DefinedTerm",
    "name": "polygenic risk score",
    "displayName": "Polygenic Risk Score",
    "description": (
        "An estimate of genetic risk obtained by aggregating the effects of " "many common genomic variants."
    ),
    "identifier": "EFO:0030082",
    "url": "http://www.ebi.ac.uk/efo/EFO_0030082",
    "inDefinedTermSet": "EFO",
}

MEASUREMENT_TECHNIQUE = [
    {
        "@type": "DefinedTerm",
        "name": "computational method",
        "identifier": "MMO_0000575",
        "url": "http://purl.obolibrary.org/obo/MMO_0000575",
        "inDefinedTermSet": "MMO",
    },
    {
        "@type": "DefinedTerm",
        "name": "Genetic variation analysis",
        "identifier": "operation_3197",
        "url": "http://edamontology.org/operation_3197",
        "inDefinedTermSet": "EDAM",
    },
]

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "name": "Genetics",
        "identifier": "topic_3053",
        "url": "http://edamontology.org/topic_3053",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Genomics",
        "identifier": "topic_0622",
        "url": "http://edamontology.org/topic_0622",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Genotype and phenotype",
        "identifier": "topic_0625",
        "url": "http://edamontology.org/topic_0625",
        "inDefinedTermSet": "EDAM",
    },
]

VARIABLE_MEASURED = {
    "@type": "DefinedTerm",
    "name": "Polygenic risk score",
    "identifier": "EFO:0030082",
    "url": "http://www.ebi.ac.uk/efo/EFO_0030082",
    "inDefinedTermSet": "EFO",
}

FUNDING = [
    {
        "@type": "MonetaryGrant",
        "identifier": "1U24HG012542-01",
        "funder": {
            "@type": "Organization",
            "name": "National Human Genome Research Institute",
        },
    },
    {
        "@type": "MonetaryGrant",
        "funder": {"@type": "Organization", "name": "Health Data Research UK"},
    },
    {
        "@type": "MonetaryGrant",
        "funder": {
            "@type": "Organization",
            "name": "Baker Heart and Diabetes Institute",
        },
    },
]

CANONICAL_PUBLICATION_DATES = {
    "10.1038/s41588-024-01937-x": "2024-09-26",
}


def _is_empty(value: Any) -> bool:
    return value is None or value == "" or value == [] or value == {}


def _clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: cleaned
            for key, item in value.items()
            if not _is_empty(item) and not _is_empty(cleaned := _clean(item))
        }
    if isinstance(value, list):
        return [cleaned for item in value if not _is_empty(item) and not _is_empty(cleaned := _clean(item))]
    return value


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _unique(values: Iterable[Any]) -> list[Any]:
    unique = []
    for value in values:
        if not _is_empty(value) and value not in unique:
            unique.append(value)
    return unique


def _bool_env(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off"}


def _int_env(name: str) -> Optional[int]:
    value = os.environ.get(name)
    return int(value) if value else None


def _trait_ids_env() -> Optional[list[str]]:
    value = os.environ.get("PGS_CATALOG_TRAIT_IDS")
    if not value:
        return None
    return sorted({item for item in re.split(r"[\s,]+", value) if item})


def _get_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(
        total=6,
        connect=6,
        read=6,
        backoff_factor=2,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "HEAD"),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retries, pool_connections=4, pool_maxsize=4)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": USER_AGENT})
    session._nde_last_request = 0.0
    return session


def _request_json(
    session: requests.Session,
    url: str,
    params: Optional[dict[str, Any]] = None,
) -> Any:
    elapsed = time.monotonic() - getattr(session, "_nde_last_request", 0.0)
    if REQUEST_INTERVAL > elapsed:
        time.sleep(REQUEST_INTERVAL - elapsed)
    response = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
    session._nde_last_request = time.monotonic()
    response.raise_for_status()
    return response.json()


def _iter_paginated(
    session: requests.Session,
    url: str,
    params: Optional[dict[str, Any]] = None,
) -> Iterator[dict[str, Any]]:
    next_url = url
    next_params = params
    while next_url:
        payload = _request_json(session, next_url, next_params)
        next_params = None
        if isinstance(payload, list):
            yield from payload
            return
        for result in payload.get("results", []):
            yield result
        next_url = payload.get("next")


def _fetch_scores(
    session: requests.Session,
    score_ids: list[str],
    cache: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    missing = [score_id for score_id in score_ids if score_id not in cache]
    for start in range(0, len(missing), SCORE_CHUNK_SIZE):
        chunk = missing[start : start + SCORE_CHUNK_SIZE]
        logger.info("Fetching %s PGS score records", len(chunk))
        for score in _iter_paginated(
            session,
            f"{API_BASE}/score/all",
            {"filter_ids": ",".join(chunk), "limit": PAGE_SIZE},
        ):
            if score_id := score.get("id"):
                cache[score_id] = score
    unfetched = sorted(set(missing) - set(cache))
    for score_id in unfetched:
        logger.warning("PGS Catalog returned no score metadata for %s", score_id)
    return [cache[score_id] for score_id in score_ids if score_id in cache]


def _fetch_performance(
    session: requests.Session,
    score_ids: list[str],
    cache: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    performance = []
    for score_id in score_ids:
        if score_id not in cache:
            cache[score_id] = list(
                _iter_paginated(
                    session,
                    f"{API_BASE}/performance/search",
                    {"pgs_id": score_id, "limit": PAGE_SIZE},
                )
            )
        performance.extend(cache[score_id])
    return performance


def iter_pgs_catalog_records() -> Iterator[dict[str, Any]]:
    """Compose the API sections described by the reviewed mapping."""
    session = _get_session()
    info = _request_json(session, f"{API_BASE}/info")
    selected_traits = _trait_ids_env()
    trait_limit = _int_env("PGS_CATALOG_TRAIT_LIMIT")
    include_performance = _bool_env("PGS_CATALOG_INCLUDE_PERFORMANCE", True)
    score_cache: dict[str, dict[str, Any]] = {}
    performance_cache: dict[str, list[dict[str, Any]]] = {}

    if selected_traits:
        traits: Iterable[dict[str, Any]] = (
            _request_json(
                session,
                f"{API_BASE}/trait/{trait_id}",
                {"include_children": 1},
            )
            for trait_id in selected_traits
        )
    else:
        traits = _iter_paginated(
            session,
            f"{API_BASE}/trait/all",
            {"include_child_associated_pgs_ids": 1, "limit": PAGE_SIZE},
        )

    emitted = 0
    for trait in traits:
        trait_id = str(trait.get("id") or "")
        if not trait_id.startswith("MONDO_"):
            continue
        score_ids = sorted(
            set(_as_list(trait.get("associated_pgs_ids"))) | set(_as_list(trait.get("child_associated_pgs_ids")))
        )
        if not score_ids:
            logger.info("Skipping score-less MONDO trait %s", trait_id)
            continue

        scores = _fetch_scores(session, score_ids, score_cache)
        performance = _fetch_performance(session, score_ids, performance_cache) if include_performance else []
        yield {
            "info": info,
            "trait": trait,
            "scores": scores,
            "performance": performance,
        }
        emitted += 1
        if trait_limit is not None and emitted >= trait_limit:
            return


def _trait_url(trait_id: str) -> str:
    return f"{WEB_BASE}/trait/{trait_id}/"


def _score_url(score_id: str) -> str:
    return f"{WEB_BASE}/score/{score_id}/"


def _sample_set_url(sample_set_id: str) -> str:
    return f"{WEB_BASE}/sampleset/{sample_set_id}/"


def _property_value(name: str, value: Any, property_id: Optional[str] = None) -> Optional[dict[str, Any]]:
    if _is_empty(value):
        return None
    if isinstance(value, (dict, list)):
        value = json.dumps(value, sort_keys=True, separators=(",", ":"))
    elif isinstance(value, bool):
        value = str(value).lower()
    else:
        value = str(value)
    return _clean(
        {
            "@type": "PropertyValue",
            "name": name,
            "propertyID": property_id,
            "value": value,
        }
    )


def _fingerprint(value: Any) -> str:
    serialized = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def _human_label(value: Any) -> str:
    return str(value).replace("_", " ").strip().capitalize()


def _age_metadata(age: Any) -> tuple[Optional[dict[str, Any]], list[dict[str, Any]]]:
    if not isinstance(age, dict):
        return None, []
    estimate = age.get("estimate")
    estimate_type = age.get("estimate_type") or "participant"
    unit = age.get("unit")
    interval = age.get("interval") if isinstance(age.get("interval"), dict) else {}
    name_parts = []
    if estimate is not None:
        name_parts.append(f"{estimate_type} participant age: {estimate}{' ' + str(unit) if unit else ''}")
    if interval:
        interval_type = interval.get("type") or "interval"
        lower = interval.get("lower")
        upper = interval.get("upper")
        if lower is not None or upper is not None:
            name_parts.append(f"{interval_type}: {lower} to {upper}{' ' + str(unit) if unit else ''}")
    developmental_stage = _clean(
        {
            "@type": "QuantitativeValue",
            "name": "; ".join(name_parts) or "Participant age",
            "value": int(estimate) if isinstance(estimate, (int, float)) and float(estimate).is_integer() else None,
            "minValue": interval.get("lower"),
            "maxValue": interval.get("upper"),
            "unitText": unit,
        }
    )
    properties = []
    if estimate is not None:
        properties.append(
            _property_value(
                f"{_human_label(estimate_type)} participant age",
                f"{estimate}{' ' + str(unit) if unit else ''}",
            )
        )
    if age.get("variability") is not None:
        variability_type = age.get("variability_type") or "variability"
        properties.append(
            _property_value(
                f"Participant age {_human_label(variability_type).lower()}",
                f"{age['variability']}{' ' + str(unit) if unit else ''}",
            )
        )
    return developmental_stage, [item for item in properties if item]


def _cohort_value(cohort: dict[str, Any]) -> Optional[str]:
    short = cohort.get("name_short")
    full = cohort.get("name_full")
    others = cohort.get("name_others")
    if short and full:
        value = f"{short}: {full}"
    else:
        value = short or full
    if others:
        other_text = ", ".join(str(item) for item in others) if isinstance(others, list) else str(others)
        value = f"{value} (also: {other_text})" if value else other_text
    return value


def _sample_row_metadata(
    sample: dict[str, Any],
    fallback_phenotype: Optional[str] = None,
) -> dict[str, Any]:
    quantities = []
    for field, label in (("sample_cases", "cases"), ("sample_controls", "controls")):
        if sample.get(field) is not None:
            quantities.append(
                {
                    "@type": "QuantitativeValue",
                    "name": label,
                    "value": sample[field],
                    "unitText": "participants",
                }
            )

    developmental_stage, properties = _age_metadata(sample.get("sample_age"))
    property_fields = (
        ("sample_percent_male", "Male participants (%)", None),
        ("followup_time", "Follow-up time", None),
        ("ancestry_broad", "Broad ancestry category", None),
        ("ancestry_free", "Detailed ancestry", None),
        ("ancestry_country", "Country of recruitment", None),
        ("ancestry_additional", "Additional ancestry description", None),
        ("source_GWAS_catalog", "GWAS Catalog study", "GWAS Catalog study"),
        ("source_PMID", "Sample source PMID", "PMID"),
        ("source_DOI", "Sample source DOI", "DOI"),
        ("cohorts_additional", "Additional sample/cohort information", None),
    )
    for field, name, property_id in property_fields:
        if prop := _property_value(name, sample.get(field), property_id):
            properties.append(prop)
    for cohort in _as_list(sample.get("cohorts")):
        if isinstance(cohort, dict) and (value := _cohort_value(cohort)):
            properties.append(_property_value("Cohort", value))

    phenotype = sample.get("phenotyping_free") or fallback_phenotype
    return _clean(
        {
            "number": sample.get("sample_number"),
            "sampleQuantity": quantities,
            "developmentalStage": developmental_stage,
            "associatedPhenotype": (
                {
                    "@type": "DefinedTerm",
                    "name": phenotype,
                }
                if phenotype
                else None
            ),
            "additionalProperty": _unique(properties),
        }
    )


def _sample_collection(
    identifier: str,
    additional_type: str,
    samples: list[dict[str, Any]],
    url: Optional[str] = None,
    fallback_phenotype: Optional[str] = None,
) -> Optional[dict[str, Any]]:
    metadata = [_sample_row_metadata(sample, fallback_phenotype) for sample in samples]
    total = sum(item.get("number") or 0 for item in metadata)
    quantities: dict[str, float] = {}
    developmental_stages = []
    phenotypes = []
    properties = []
    for item in metadata:
        for quantity in _as_list(item.get("sampleQuantity")):
            name = quantity.get("name")
            if name and quantity.get("value") is not None:
                quantities[name] = quantities.get(name, 0) + quantity["value"]
        if stage := item.get("developmentalStage"):
            developmental_stages.append(stage)
        if phenotype := item.get("associatedPhenotype"):
            phenotypes.append(phenotype)
        properties.extend(_as_list(item.get("additionalProperty")))

    aggregate = _clean(
        {
            "@type": "Sample",
            "sampleQuantity": [
                {
                    "@type": "QuantitativeValue",
                    "name": name,
                    "value": value,
                    "unitText": "participants",
                }
                for name, value in sorted(quantities.items())
            ],
            "developmentalStage": _unique(developmental_stages),
            "associatedPhenotype": _unique(phenotypes),
            "additionalProperty": _unique(properties),
        }
    )
    if not aggregate or aggregate == {"@type": "Sample"}:
        aggregate = None
    if not total and not aggregate:
        return None
    return _clean(
        {
            "@type": "SampleCollection",
            "identifier": identifier,
            "url": url,
            "additionalType": additional_type,
            "numberOfItems": (
                {
                    "@type": "QuantitativeValue",
                    "value": total,
                    "unitText": "participants",
                }
                if total
                else None
            ),
            "aggregateElement": aggregate,
        }
    )


def _development_samples(scores: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    seen = set()
    for score in sorted(scores, key=lambda item: item.get("id") or ""):
        score_id = score.get("id")
        if not score_id:
            continue
        for source_field, label, suffix in (
            ("samples_variants", "PGS variant-association sample", "variants"),
            ("samples_training", "PGS score-training sample", "training"),
        ):
            for index, sample in enumerate(_as_list(score.get(source_field)), start=1):
                if not isinstance(sample, dict):
                    continue
                fingerprint = _fingerprint(sample)
                if fingerprint in seen:
                    continue
                seen.add(fingerprint)
                collection = _sample_collection(
                    f"{score_id}:{suffix}:{index}",
                    label,
                    [sample],
                    url=_score_url(score_id),
                )
                if collection:
                    output.append(collection)
    return output


def _evaluation_samples(performance: list[dict[str, Any]]) -> list[dict[str, Any]]:
    sample_sets: dict[str, dict[str, Any]] = {}
    for result in performance:
        sample_set = result.get("sampleset") or {}
        sample_set_id = sample_set.get("id")
        if not sample_set_id:
            continue
        entry = sample_sets.setdefault(
            sample_set_id,
            {"samples": [], "fingerprints": set(), "phenotypes": []},
        )
        if result.get("phenotyping_reported"):
            entry["phenotypes"].append(result["phenotyping_reported"])
        for sample in _as_list(sample_set.get("samples")):
            if not isinstance(sample, dict):
                continue
            fingerprint = _fingerprint(sample)
            if fingerprint not in entry["fingerprints"]:
                entry["fingerprints"].add(fingerprint)
                entry["samples"].append(sample)

    output = []
    for sample_set_id in sorted(sample_sets):
        entry = sample_sets[sample_set_id]
        fallback = next(iter(_unique(entry["phenotypes"])), None)
        collection = _sample_collection(
            sample_set_id,
            "PGS evaluation sample set",
            entry["samples"],
            url=_sample_set_url(sample_set_id),
            fallback_phenotype=fallback,
        )
        if collection:
            output.append(collection)
    return output


def _publication_to_citation(publication: Any) -> Optional[dict[str, Any]]:
    if not isinstance(publication, dict):
        return None
    doi = publication.get("doi")
    pmid = publication.get("PMID") or publication.get("pmid")
    title = publication.get("title") or publication.get("name")
    identifier = publication.get("id") or publication.get("identifier")
    if not any((doi, pmid, title, identifier)):
        return None
    authors_text = publication.get("authors")
    authors = []
    if publication.get("firstauthor"):
        authors.append({"@type": "Person", "name": publication["firstauthor"]})
    elif authors_text:
        names = (
            [name.strip() for name in authors_text.split(",") if name.strip()]
            if "," in authors_text and " et al" not in authors_text
            else [authors_text]
        )
        authors = [{"@type": "Person", "name": name} for name in names]
    date_published = _to_iso_date(publication.get("date_publication"))
    if not date_published and doi:
        date_published = CANONICAL_PUBLICATION_DATES.get(str(doi).lower())
    return _clean(
        {
            "@type": "ScholarlyArticle",
            "identifier": identifier,
            "name": title,
            "doi": str(doi) if doi else None,
            "pmid": str(pmid) if pmid else None,
            "journalName": publication.get("journal"),
            "datePublished": date_published,
            "author": authors,
            "url": f"https://doi.org/{doi}" if doi else None,
        }
    )


def _citation_key(citation: dict[str, Any]) -> tuple[str, str]:
    if citation.get("doi"):
        return "doi", str(citation["doi"]).lower()
    if citation.get("pmid"):
        return "pmid", str(citation["pmid"])
    if citation.get("identifier"):
        return "identifier", str(citation["identifier"])
    return "name", str(citation.get("name") or "")


def _citations(
    info: dict[str, Any],
    scores: list[dict[str, Any]],
    performance: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    publications = [info.get("citation")]
    publications.extend(_as_list(info.get("pgs_catalog_publications")))
    publications.extend(score.get("publication") for score in scores)
    publications.extend(item.get("publication") for item in performance)
    output: dict[tuple[str, str], dict[str, Any]] = {}
    for publication in publications:
        citation = _publication_to_citation(publication)
        if citation:
            key = _citation_key(citation)
            current = output.get(key)
            if not current or len(citation) > len(current):
                output[key] = citation
    return list(output.values())


def _has_parts(score_ids: list[str], scores: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_id = {score.get("id"): score for score in scores if score.get("id")}
    output = []
    for score_id in score_ids:
        score = by_id.get(score_id, {})
        output.append(
            _clean(
                {
                    "@type": "CreativeWork",
                    "additionalType": {
                        "@type": "DefinedTerm",
                        "name": "Polygenic risk score",
                        "url": "http://www.ebi.ac.uk/efo/EFO_0030082",
                    },
                    "identifier": score_id,
                    "name": score.get("name") or score_id,
                    "url": _score_url(score_id),
                    "datePublished": _to_iso_date(score.get("date_release")),
                }
            )
        )
    return output


def _distributions(scores: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    seen_urls = set()
    for score in sorted(scores, key=lambda item: item.get("id") or ""):
        score_id = score.get("id")
        if not score_id:
            continue
        date_published = _to_iso_date(score.get("date_release"))
        downloads = [
            (
                score_id,
                f"{score.get('name') or score_id} scoring file",
                score.get("ftp_scoring_file"),
            )
        ]
        harmonized = score.get("ftp_harmonized_scoring_files") or {}
        for build in sorted(harmonized):
            value = harmonized.get(build)
            content_url = value.get("positions") if isinstance(value, dict) else value
            downloads.append(
                (
                    f"{score_id}-{build}",
                    f"{score.get('name') or score_id} harmonized scoring file ({build})",
                    content_url,
                )
            )
        for identifier, name, content_url in downloads:
            if not content_url or content_url in seen_urls:
                continue
            seen_urls.add(content_url)
            output.append(
                _clean(
                    {
                        "@type": "DataDownload",
                        "identifier": identifier,
                        "name": name,
                        "contentUrl": content_url,
                        "encodingFormat": "application/gzip",
                        "datePublished": date_published,
                    }
                )
            )
    return output


def _example_of_work(
    representative_id: str,
    score: dict[str, Any],
    performance: list[dict[str, Any]],
    info: dict[str, Any],
) -> dict[str, Any]:
    properties = []
    fields = (
        ("PGS ID", representative_id),
        ("Score name", score.get("name")),
        ("Reported trait", score.get("trait_reported")),
        ("Additional trait description", score.get("trait_additional")),
        ("Development method", score.get("method_name")),
        ("Development method parameters", score.get("method_params")),
        ("Variant count", score.get("variants_number")),
        ("Variant interaction count", score.get("variants_interactions")),
        ("Original genome build", score.get("variants_genomebuild")),
        ("Weight type", score.get("weight_type")),
        ("Ancestry distribution", score.get("ancestry_distribution")),
        ("Matches publication", score.get("matches_publication")),
        ("Score-specific license or terms", score.get("license")),
        ("PGS Catalog REST API version", (info.get("rest_api") or {}).get("version")),
        ("Ensembl version", info.get("ensembl_version")),
    )
    for name, value in fields:
        if prop := _property_value(name, value):
            properties.append(prop)

    for result in performance:
        if result.get("associated_pgs_id") != representative_id:
            continue
        for name, value in (
            ("Performance ID", result.get("id")),
            ("Performance metrics", result.get("performance_metrics")),
            ("Performance covariates", result.get("covariates")),
            ("Performance comments", result.get("performance_comments")),
        ):
            if prop := _property_value(name, value):
                properties.append(prop)

    return _clean(
        {
            "@type": "CreativeWork",
            "about": copy.deepcopy(POLYGENIC_RISK_SCORE),
            "encodingFormat": [
                {
                    "@type": "DefinedTerm",
                    "name": "TSV",
                    "identifier": "format_3475",
                    "url": "http://edamontology.org/format_3475",
                    "inDefinedTermSet": "EDAM",
                }
            ],
            "schemaVersion": "PGS Catalog scoring file format 2.0",
            "potentialAction": {
                "@type": "SearchAction",
                "name": "Retrieve PGS score metadata",
                "url": f"{API_BASE}/score/{representative_id}",
            },
            "additionalProperty": _unique(properties),
        }
    )


def _is_based_on(trait_id: str, scores: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = [
        {
            "@type": "Action",
            "name": "DataCollection Generation Process in the NIAID Data Ecosystem",
            "description": (
                "This record aggregates PGS Catalog polygenic score data into a "
                "collection grouped by MONDO health-condition trait identifier. "
                "Collection size and dates come from the source; descriptive "
                "fields are manually curated."
            ),
            "actionProcess": {
                "@type": "HowTo",
                "step": [
                    "1. Retrieve MONDO traits and child-score associations from the PGS Catalog REST API.",
                    "2. Fetch and deduplicate direct and child-associated score metadata.",
                    "3. Fetch score development and evaluation sample metadata.",
                    "4. Create one DataCollection per score-bearing MONDO trait "
                    "with score members, downloads, citations, and samples.",
                ],
            },
        },
        {
            "@type": "ResourceCatalog",
            "name": "PGS Catalog",
            "url": "https://www.pgscatalog.org/",
        },
        {
            "@type": "CreativeWork",
            "name": "PGS Catalog REST API trait query",
            "url": f"{API_BASE}/trait/{trait_id}?include_children=1",
        },
    ]
    gwas_ids = set()
    for score in scores:
        for field in ("samples_variants", "samples_training"):
            for sample in _as_list(score.get(field)):
                if isinstance(sample, dict) and sample.get("source_GWAS_catalog"):
                    gwas_ids.add(str(sample["source_GWAS_catalog"]))
    output.extend(
        {
            "@type": "CreativeWork",
            "identifier": gwas_id,
            "name": f"GWAS Catalog study {gwas_id}",
            "url": f"https://www.ebi.ac.uk/gwas/studies/{gwas_id}",
        }
        for gwas_id in sorted(gwas_ids)
    )
    return output


def _transform_record(record: dict[str, Any]) -> Optional[dict[str, Any]]:
    info = record.get("info") or {}
    trait = record.get("trait") or {}
    scores = [score for score in _as_list(record.get("scores")) if isinstance(score, dict)]
    performance = [item for item in _as_list(record.get("performance")) if isinstance(item, dict)]
    trait_id = str(trait.get("id") or "")
    if not trait_id.startswith("MONDO_"):
        logger.warning("Skipping non-MONDO or identifier-less trait: %s", trait_id)
        return None

    direct_ids = sorted(set(_as_list(trait.get("associated_pgs_ids"))))
    child_ids = sorted(set(_as_list(trait.get("child_associated_pgs_ids"))))
    score_ids = sorted(set(direct_ids) | set(child_ids))
    if not score_ids:
        return None

    trait_label = trait.get("label") or trait_id
    trait_url = _trait_url(trait_id)
    release_date = _to_iso_date((info.get("latest_release") or {}).get("date"))
    output = {
        "@context": CONTEXT,
        "@type": "DataCollection",
        "_id": f"pgs_catalog_{trait_id}",
        "identifier": trait_id,
        "url": trait_url,
        "distribution": [],
        "includedInDataCatalog": {
            "@type": "DataCatalog",
            "name": "PGS Catalog",
            "url": "https://www.pgscatalog.org/",
            "versionDate": release_date or datetime.date.today().isoformat(),
            "archivedAt": trait_url,
        },
    }

    insert_value(output, "date", release_date or datetime.date.today().isoformat())
    insert_value(output, "dateModified", release_date or datetime.date.today().isoformat())
    insert_value(output, "name", f"Polygenic scores for {trait_label} from PGS Catalog")
    description = (
        f"PGS Catalog polygenic scores associated with {trait_label} "
        f"({trait_id.replace('_', ':', 1)}). This collection contains "
        f"{len(score_ids)} scores associated directly with {trait_label} or "
        "with its child traits, together with score-development and evaluation "
        "sample metadata."
    )
    insert_value(output, "description", description, extend=True)
    if trait.get("description"):
        insert_value(output, "description", str(trait["description"]), extend=True)
    insert_value(
        output,
        "collectionSize",
        {"@type": "QuantitativeValue", "value": len(score_ids), "unitText": "polygenic scores"},
    )
    insert_value(output, "about", copy.deepcopy(ABOUT))
    health_condition = _clean(
        {
            "@type": "DefinedTerm",
            "name": trait_label,
            "identifier": trait_id.replace("_", ":", 1),
            "url": trait.get("url") or f"http://purl.obolibrary.org/obo/{trait_id}",
            "inDefinedTermSet": "MONDO",
            "alternateName": _unique(_as_list(trait.get("trait_synonyms"))),
        }
    )
    insert_value(output, "healthCondition", health_condition)
    insert_value(output, "hasPart", _has_parts(score_ids, scores))
    insert_value(output, "distribution", _distributions(scores))
    samples = _development_samples(scores) + _evaluation_samples(performance)
    if samples:
        insert_value(output, "sample", samples)
    insert_value(output, "author", copy.deepcopy(AUTHORS))
    insert_value(output, "creator", copy.deepcopy(AUTHORS))
    insert_value(output, "sourceOrganization", copy.deepcopy(SOURCE_ORGANIZATION))
    citations = _citations(info, scores, performance)
    if citations:
        insert_value(output, "citation", citations)
    insert_value(output, "conditionsOfAccess", "Open")
    insert_value(output, "funding", copy.deepcopy(FUNDING))
    insert_value(output, "isAccessibleForFree", True)
    keywords = _unique(
        ["PGS Catalog", "polygenic score", "polygenic risk score", trait_label]
        + _as_list(trait.get("trait_categories"))
        + [score.get("trait_reported") for score in scores]
    )
    insert_value(output, "keywords", keywords)
    insert_value(output, "measurementTechnique", copy.deepcopy(MEASUREMENT_TECHNIQUE))
    insert_value(output, "topicCategory", copy.deepcopy(TOPIC_CATEGORY))
    insert_value(
        output,
        "usageInfo",
        {
            "@type": "CreativeWork",
            "name": "EMBL-EBI Terms of Use",
            "description": "Individual scores may have additional author-specified licenses or restrictions.",
            "url": info.get("terms_of_use") or "https://www.ebi.ac.uk/about/terms-of-use/",
        },
    )
    insert_value(output, "variableMeasured", copy.deepcopy(VARIABLE_MEASURED))

    representative_id = next((score_id for score_id in direct_ids if score_id), score_ids[0])
    representative = next(
        (score for score in scores if score.get("id") == representative_id),
        {"id": representative_id},
    )
    insert_value(
        output,
        "exampleOfWork",
        _example_of_work(representative_id, representative, performance, info),
    )
    insert_value(output, "isBasedOn", _is_based_on(trait_id, scores))
    return _clean(output)


def parse(records: Optional[Iterable[dict[str, Any]]] = None) -> Iterator[dict[str, Any]]:
    """Yield one NDE DataCollection for every score-bearing MONDO trait."""
    raw_records = records if records is not None else iter_pgs_catalog_records()
    for record in raw_records:
        output = _transform_record(record)
        if output:
            yield output
