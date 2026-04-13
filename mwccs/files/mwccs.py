#!/usr/bin/env python3
"""
MWCCS DataCollection crawler for the NIAID Data Ecosystem.

Dynamically discovers, downloads, and parses visual abstract PDFs from the
MWCCS (Multicenter AIDS Cohort Study / Women's Interagency HIV Study Combined
Cohort Study) substudy science page, producing one DataCollection record per
substudy.

Data source
-----------
PDF visual abstracts are crawled from:
    https://statepi.jhsph.edu/mwccs/substudy-science/

Each PDF follows a structured template with sections for study name,
measurement details, participant demographics, temporal coverage, and authors.
"""

import calendar
import copy
import io
import logging
import re
import time
from html.parser import HTMLParser
from typing import Any, Optional
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen

from pypdf import PdfReader

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SUBSTUDY_PAGE_URL = "https://statepi.jhsph.edu/mwccs/substudy-science/"
USER_AGENT = "nde-mwccs-crawler/0.1"
MAX_RETRIES = 3
RETRY_BACKOFF = 10

MONTHS = {
    "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
    "jul": 7, "aug": 8, "sep": 9, "sept": 9, "oct": 10, "nov": 11, "dec": 12,
}

# ---------------------------------------------------------------------------
# Curated / static metadata (from the MWCCS resource catalog record)
# ---------------------------------------------------------------------------

INCLUDED_IN_DATA_CATALOG = {
    "@type": "DataCatalog",
    "name": (
        "The Multicenter AIDS Cohort Study (MACS) / Women's Interagency "
        "HIV Study (WIHS) Combined Cohort Study (MWCCS)"
    ),
    "alternateName": "MWCCS",
    "url": "https://statepi.jhsph.edu/mwccs/",
}

ABOUT_DEFINED_TERM = {
    "@type": "DefinedTerm",
    "description": "For schema, consider https://schema.org/Patient",
    "displayName": "Patient",
    "name": "Patient",
    "url": "http://purl.obolibrary.org/obo/NCIT_C16960",
}

TOPIC_CATEGORY = [
    {
        "@type": "DefinedTerm",
        "name": "Infectious disease",
        "identifier": "topic_0804",
        "url": "http://edamontology.org/topic_0804",
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

SPECIES = {"@type": "DefinedTerm", "name": "Homo sapiens"}

FALLBACK_HEALTH_CONDITION = {
    "@type": "DefinedTerm",
    "name": "HIV Infectious Disease",
    "url": "http://purl.obolibrary.org/obo/MONDO_0005109",
    "identifier": "MONDO_0005109",
    "inDefinedTermSet": "MONDO",
}

FALLBACK_INFECTIOUS_AGENT = [
    {"@type": "DefinedTerm", "name": "HIV-1"},
    {"@type": "DefinedTerm", "name": "HIV-2"},
]

USAGE_INFO = {
    "@type": "CreativeWork",
    "name": "MWCCS DACC Resource Request Info",
    "url": (
        "https://statepi.jhsph.edu/mwccs/wp-content/uploads/2022/11/"
        "MWCCS-DACC-Resource-Request-Form_final11922.pdf"
    ),
}

IS_BASED_ON = [
    {
        "@type": "Action",
        "name": "DataCollection Generation Process in the NIAID Data Ecosystem",
        "description": (
            "How this MWCCS Study Subjects DataCollection Record was "
            "generated for the NIAID Data Ecosystem."
        ),
        "actionProcess": {
            "@type": "HowTo",
            "step": [
                (
                    "Step 1. Use a crawler to download the Substudy visual "
                    "summary PDFs from the MWCCS substudy science page."
                ),
                (
                    "Step 2: Parse the relevant information directly from "
                    "each PDF including study name, measurement techniques, "
                    "variables measured, participant demographics, temporal "
                    "coverage, and scientific leads."
                ),
                (
                    "Step 3: Calculate sampleSize information based on the "
                    "provided percentages in the PDF."
                ),
                (
                    "Step 4. Fill in manually curated fields from the MWCCS "
                    "resource catalog record such as 'usageInfo', 'license', "
                    "'topicCategory', 'conditionsOfAccess', 'about'."
                ),
            ],
        },
    },
    {
        "@type": "ResourceCatalog",
        "name": (
            "The Multicenter AIDS Cohort Study (MACS) / Women's Interagency "
            "HIV Study (WIHS) Combined Cohort Study"
        ),
        "alternateName": "MWCCS",
        "url": "https://data.niaid.nih.gov/resources?id=dde_da18bc66317fc987",
    },
]

# Funding from the MWCCS resource catalog — all grants share the same funder list
_FUNDER_LIST = [
    {"@type": "Organization", "alternateName": "NIAAA", "name": "National Institute on Alcohol Abuse and Alcoholism", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NIA", "name": "National Institute On Aging", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NIAID", "name": "National Institute Of Allergy And Infectious Diseases", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NCI", "name": "National Cancer Institute", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NIDA", "name": "National Institute On Drug Abuse", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NIDCD", "name": "National Institute on Deafness and Other Communication Disorders", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NIDCR", "name": "National Institute Of Dental & Craniofacial Research", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NIDDK", "name": "National Institute of Diabetes and Digestive and Kidney Diseases", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NICHD", "name": "Eunice Kennedy Shriver National Institute Of Child Health & Human Development", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NHLBI", "name": "National Heart, Lung, and Blood Institute", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NINDS", "name": "National Institute Of Neurological Disorders And Stroke", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NIMH", "name": "National Institute Of Mental Health", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NINR", "name": "National Institute Of Nursing Research", "parentOrganization": "National Institutes of Health"},
    {"@type": "Organization", "alternateName": "NIMHD", "name": "National Institute on Minority Health and Health Disparities", "parentOrganization": "National Institutes of Health"},
]

FUNDING = [
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146241"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146201"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146193"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146245"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146333"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146204"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146202"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146192"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146194"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146203"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146205"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146240"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146242"},
    {"@type": "MonetaryGrant", "funder": copy.deepcopy(_FUNDER_LIST), "identifier": "U01-HL146208"},
]

# Keywords that indicate a specific health condition from PDF content
HEALTH_CONDITION_KEYWORDS = {
    "heart": ["cardiovascular disease", "heart disease"],
    "cardiac": ["cardiovascular disease"],
    "echocardiogram": ["cardiovascular disease", "heart disease"],
    "liver": ["liver fibrosis", "liver steatosis"],
    "fibroscan": ["liver fibrosis", "liver steatosis"],
    "hepat": ["liver disease"],
    "cognitive": ["cognitive impairment"],
    "neurocognitive": ["cognitive impairment"],
    "memory": ["cognitive impairment"],
    "hearing": ["hearing disorder"],
    "balance": ["vestibular disorder"],
    "vertigo": ["vestibular disorder"],
    "pulmonary": ["pulmonary disease"],
    "lung": ["lung disease"],
    "spirom": ["pulmonary disease"],
    "tooth": ["oral disease", "dental disease"],
    "oral": ["oral disease"],
    "dental": ["dental disease"],
    "mental health": ["mental disorder"],
    "substance": ["substance use disorder"],
    "microbiome": ["gut microbiome alteration"],
    "stool": ["gut microbiome alteration"],
    "sleep": ["sleep disorder"],
    "insomnia": ["sleep disorder"],
}


# ---------------------------------------------------------------------------
# HTML parsing for PDF link discovery
# ---------------------------------------------------------------------------


class _PdfLinkParser(HTMLParser):
    """Collect substudy PDF links from the MWCCS substudy science page."""

    def __init__(self, base_url: str) -> None:
        super().__init__()
        self.base_url = base_url
        self.links: dict[str, str] = {}
        self._href: Optional[str] = None
        self._text: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, Optional[str]]]) -> None:
        if tag != "a":
            return
        href = dict(attrs).get("href")
        if href:
            self._href = urljoin(self.base_url, href)
            self._text = []

    def handle_data(self, data: str) -> None:
        if self._href is not None:
            self._text.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag != "a" or self._href is None:
            return
        text = " ".join(part.strip() for part in self._text if part.strip())
        if self._href.lower().endswith(".pdf") and "substudy" in text.lower():
            from pathlib import PurePosixPath
            filename = PurePosixPath(urlparse(self._href).path).name
            self.links[filename] = self._href
        self._href = None
        self._text = []


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------


def _fetch_bytes(url: str, timeout: int = 120, retries: int = MAX_RETRIES) -> bytes:
    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:
            last_exc = exc
            if attempt < retries:
                wait = RETRY_BACKOFF * (2 ** (attempt - 1))
                logger.warning(
                    "Attempt %d/%d failed for %s: %s - retrying in %ds",
                    attempt, retries, url[:120], exc, wait,
                )
                time.sleep(wait)
    raise last_exc


def _fetch_text(url: str, timeout: int = 60) -> str:
    return _fetch_bytes(url, timeout=timeout).decode("utf-8", errors="replace")


# ---------------------------------------------------------------------------
# PDF text extraction and normalization
# ---------------------------------------------------------------------------


def _read_pdf_text(pdf_bytes: bytes) -> str:
    reader = PdfReader(io.BytesIO(pdf_bytes))
    parts = [page.extract_text() or "" for page in reader.pages]
    text = "\n".join(parts)
    # Normalize special characters
    text = text.replace("\u2010", "-").replace("\u2011", "-")
    text = text.replace("\u2013", "-").replace("\u2014", "-")
    text = text.replace("\u00ae", "").replace("\uf071", " ")
    text = text.replace("\uf0a7", "\u2022")  # Windows bullet -> standard bullet
    return text


def _normalize_ws(text: str) -> str:
    return re.sub(r"[ \t]+", " ", text).strip()


def _collapse(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _clean_lines(text: str) -> list[str]:
    return [_normalize_ws(line) for line in text.splitlines() if _normalize_ws(line)]


def _slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()


# ---------------------------------------------------------------------------
# PDF section extraction (dynamic — works on any visual abstract template)
# ---------------------------------------------------------------------------


def _label_pattern(label: str) -> str:
    parts = [re.escape(part) for part in label.split()]
    return r"\s+".join(parts)


def _extract_section(text: str, label: str, next_labels: list[str], anchor_line_start: bool = False) -> str:
    lp = _label_pattern(label)
    if anchor_line_start:
        lp = r"(?:^|\n)\s*" + lp
    np = "|".join(_label_pattern(n) for n in next_labels)
    pattern = lp + r"\s*:?[\s\n]*(.*?)\s*(?=" + np + r")"
    match = re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL)
    if not match:
        return ""
    return _normalize_ws(match.group(1))


def _extract_date_phrase(text: str) -> str:
    match = re.search(
        r"\((data collect(?:ed|ion)\s+[^)]+)\)", text, flags=re.IGNORECASE,
    )
    return _normalize_ws(match.group(1)) if match else ""


def _split_bullet_items(text: str) -> list[str]:
    if not text:
        return []
    # First, split on bullet characters that may appear inline
    # The PDFs use bullet as inline bullet separators
    if "\u2022" in text:
        parts = text.split("\u2022")
        items = [_collapse(p) for p in parts if _collapse(p)]
        return items
    items: list[str] = []
    for raw_line in [line.rstrip() for line in text.splitlines() if line.strip()]:
        stripped = raw_line.strip()
        is_bullet = stripped.startswith(("*", "-", "\u2022"))
        clean = _collapse(stripped.lstrip("*-\u2022").strip())
        if not clean:
            continue
        if is_bullet or not items:
            items.append(clean)
        else:
            items[-1] = _collapse(f"{items[-1]} {clean}")
    return items


def _find_first_large_number(lines: list[str]) -> Optional[int]:
    for line in lines:
        if re.fullmatch(r"\d{3,6}", line):
            return int(line)
    return None


# ---------------------------------------------------------------------------
# Structured data extraction from PDF text
# ---------------------------------------------------------------------------


def _parse_pdf_data(pdf_text: str, pdf_url: str, filename: str) -> dict[str, Any]:
    """Extract all structured fields from a visual abstract PDF's text."""
    lines = _clean_lines(pdf_text)
    if not lines:
        return {}

    # Study name is the first line
    name = lines[0]

    # Date phrase: "(data collected Oct 2020 - Sept 2026)"
    date_phrase = _extract_date_phrase(pdf_text)

    # What is being measured
    what_measured = _extract_section(
        pdf_text,
        "What is being measured",
        ["Measures:", "Measure:", "How often:"],
    )

    # Measures (measurement techniques)
    measures_text = _extract_section(pdf_text, "Measures", ["How often:"])
    if not measures_text:
        # Use anchor_line_start to avoid matching "Measure" inside "What is being measured"
        measures_text = _extract_section(pdf_text, "Measure", ["How often:"], anchor_line_start=True)
    measures = _split_bullet_items(measures_text)

    # How often
    how_often = _extract_section(
        pdf_text,
        "How often",
        ["Who did we collect", "What can we learn"],
    )

    # Participant count — first standalone large number
    participant_count = _find_first_large_number(lines)

    # Participant description — text after the count before demographics
    participant_desc = ""
    if participant_count is not None:
        count_str = str(participant_count)
        for i, line in enumerate(lines):
            if line == count_str:
                collected = []
                for subsequent in lines[i + 1:]:
                    lowered = subsequent.lower()
                    if re.search(r"\d+\s+years old.*median age", lowered):
                        break
                    if re.fullmatch(r"\d+%", subsequent):
                        break
                    if lowered.startswith(("what can we learn", "scientific lead")):
                        break
                    collected.append(subsequent)
                participant_desc = _collapse(" ".join(collected))
                break

    # Demographics
    median_age = None
    hiv_percent = None
    women_percent = None
    for i, line in enumerate(lines):
        age_match = re.search(r"(\d+)\s+years old.*median age", line, flags=re.IGNORECASE)
        if age_match and median_age is None:
            median_age = int(age_match.group(1))
        pct_match = re.fullmatch(r"(\d+)%", line)
        if pct_match:
            next_line = lines[i + 1] if i + 1 < len(lines) else ""
            pct = int(pct_match.group(1))
            if "living with hiv" in next_line.lower() and hiv_percent is None:
                hiv_percent = pct
            elif "women" in next_line.lower() and women_percent is None:
                women_percent = pct

    # Learning points
    learning_section = _extract_section(
        pdf_text,
        "What can we learn from the collected data?",
        ["Scientific lead(s) for this study:"],
    )
    learning_points = _split_bullet_items(learning_section) if learning_section else []

    # Authors
    authors = []
    author_match = re.search(
        r"Scientific lead\(s\) for this study:\s*(.*?)\s*Learn more about MWCCS",
        pdf_text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if author_match:
        block = _normalize_ws(author_match.group(1))
        block = block.replace(", and ", ", ").replace(" and ", ", ")
        for part in [item.strip() for item in block.split(",") if item.strip()]:
            clean = re.sub(r"^Dr\.?\s*", "", part).strip()
            if clean:
                authors.append(clean)

    # Date created from PDF filename (e.g., ..._12.08.25.pdf)
    date_created = None
    date_match = re.search(r"(\d{1,2})\.(\d{1,2})\.(\d{2,4})\.pdf$", filename, flags=re.IGNORECASE)
    if date_match:
        month = int(date_match.group(1))
        day = int(date_match.group(2))
        year = int(date_match.group(3))
        if year < 100:
            year += 2000
        date_created = f"{year:04d}-{month:02d}-{day:02d}"

    # Temporal coverage from date phrase
    temporal_coverage = None
    tc_match = re.search(
        r"data collect(?:ed|ion)\s+([A-Za-z]+)\s+(\d{4})\s*-\s*([A-Za-z]+)\s+(\d{4})",
        date_phrase,
        flags=re.IGNORECASE,
    )
    if tc_match:
        start_month = MONTHS.get(tc_match.group(1)[:4].lower(), 1)
        start_year = int(tc_match.group(2))
        end_month = MONTHS.get(tc_match.group(3)[:4].lower(), 12)
        end_year = int(tc_match.group(4))
        end_day = calendar.monthrange(end_year, end_month)[1]
        temporal_coverage = {
            "startDate": f"{start_year:04d}-{start_month:02d}-01",
            "endDate": f"{end_year:04d}-{end_month:02d}-{end_day:02d}",
        }

    return {
        "name": name,
        "url": pdf_url,
        "filename": filename,
        "date_phrase": date_phrase,
        "what_measured": what_measured,
        "measures": measures,
        "how_often": how_often,
        "participant_count": participant_count,
        "participant_description": participant_desc,
        "median_age": median_age,
        "hiv_percent": hiv_percent,
        "women_percent": women_percent,
        "learning_points": learning_points,
        "authors": authors,
        "date_created": date_created,
        "temporal_coverage": temporal_coverage,
    }


# ---------------------------------------------------------------------------
# Dynamic health condition inference
# ---------------------------------------------------------------------------


def _infer_health_conditions(data: dict[str, Any]) -> list[dict[str, str]]:
    """Infer health conditions from the study name and measured variables."""
    searchable = " ".join([
        data.get("name", ""),
        data.get("what_measured", ""),
        " ".join(data.get("measures", [])),
    ]).lower()

    inferred: list[str] = []
    seen = set()
    for keyword, conditions in HEALTH_CONDITION_KEYWORDS.items():
        if keyword in searchable:
            for cond in conditions:
                if cond.lower() not in seen:
                    seen.add(cond.lower())
                    inferred.append(cond)

    terms = [{"@type": "DefinedTerm", "name": c} for c in inferred]

    # Always include the HIV fallback condition
    fallback_name = FALLBACK_HEALTH_CONDITION["name"].lower()
    if not any(t.get("name", "").lower() == fallback_name for t in terms):
        terms.append(copy.deepcopy(FALLBACK_HEALTH_CONDITION))

    return terms


# ---------------------------------------------------------------------------
# DataCollection record building
# ---------------------------------------------------------------------------


def _compute_count(total: Optional[int], percent: Optional[int]) -> Optional[int]:
    if total is None or percent is None:
        return None
    return round(total * percent / 100)


def _build_description(data: dict[str, Any]) -> str:
    parts = []
    count = data.get("participant_count")
    desc = data.get("participant_description", "")
    if count:
        parts.append(f"Data collected from {count} {desc}".strip())
    elif desc:
        parts.append(f"Data collected from {desc}")

    learning = data.get("learning_points", [])
    if learning:
        parts.append(
            "This data can help researchers learn more about "
            + "; ".join(learning)
        )

    return _collapse(". ".join(parts)) if parts else data.get("name", "")


def _build_collection_size(data: dict[str, Any]) -> list[dict[str, Any]]:
    sizes = []
    count = data.get("participant_count")
    if count:
        sizes.append({
            "@type": "QuantitativeValue",
            "minValue": count,
            "unitText": "Study Subjects",
            "unitCode": "http://purl.obolibrary.org/obo/NCIT_C41189",
        })
    return sizes


def _build_sample(data: dict[str, Any]) -> dict[str, Any]:
    sample: dict[str, Any] = {
        "@type": "Sample",
        "sex": ["Men", "Women"],
    }

    median_age = data.get("median_age")
    if median_age:
        sample["age"] = {
            "minValue": median_age,
            "unitText": "years",
            "valueType": "median",
        }

    count = data.get("participant_count")
    quantities = []

    hiv_count = _compute_count(count, data.get("hiv_percent"))
    if hiv_count is not None:
        quantities.append({
            "@type": "QuantitativeValue",
            "minValue": hiv_count,
            "unitText": "living with HIV",
        })
        without_hiv = count - hiv_count
        quantities.append({
            "@type": "QuantitativeValue",
            "minValue": without_hiv,
            "unitText": "living without HIV",
        })

    women_count = _compute_count(count, data.get("women_percent"))
    if women_count is not None:
        quantities.append({
            "@type": "QuantitativeValue",
            "minValue": women_count,
            "unitText": "Women",
        })
        men_count = count - women_count
        quantities.append({
            "@type": "QuantitativeValue",
            "minValue": men_count,
            "unitText": "Men",
        })

    if quantities:
        sample["sampleQuantity"] = quantities

    return sample


def _build_data_collection(data: dict[str, Any]) -> dict[str, Any]:
    record: dict[str, Any] = {
        "_id": f"mwccs_{_slugify(data['name'])}",
        "@type": "DataCollection",
        "name": data["name"],
        "description": _build_description(data),
        "url": data.get("url", ""),
        "about": copy.deepcopy(ABOUT_DEFINED_TERM),
        "includedInDataCatalog": [copy.deepcopy(INCLUDED_IN_DATA_CATALOG)],
        "isBasedOn": copy.deepcopy(IS_BASED_ON),
        "conditionsOfAccess": "Restricted",
        "isAccessibleForFree": False,
        "usageInfo": copy.deepcopy(USAGE_INFO),
        "funding": copy.deepcopy(FUNDING),
        "topicCategory": copy.deepcopy(TOPIC_CATEGORY),
        "species": copy.deepcopy(SPECIES),
        "infectiousAgent": copy.deepcopy(FALLBACK_INFECTIOUS_AGENT),
    }

    included_catalog = record["includedInDataCatalog"][0]
    included_catalog["archivedAt"] = data.get("url", "")


    # Date fields
    date_created = data.get("date_created")
    if date_created:
        record["dateCreated"] = date_created
        record["date"] = date_created
        record["dateModified"] = date_created

    # Measurement techniques from parsed measures
    measures = data.get("measures", [])
    if measures:
        record["measurementTechnique"] = [
            {"@type": "DefinedTerm", "name": m} for m in measures
        ]

    # Variables measured from what_measured
    what_measured = data.get("what_measured", "")
    if what_measured:
        items = _split_bullet_items(what_measured)
        if not items:
            items = [what_measured]
        record["variableMeasured"] = [
            {"@type": "DefinedTerm", "name": v} for v in items
        ]

    # Authors
    authors = data.get("authors", [])
    if authors:
        record["creator"] = []
        for name in authors:
            parts = [p for p in name.split() if p]
            author_obj = {"@type": "Person", "name": name}
            if len(parts) >= 2:
                author_obj["givenName"] = parts[0]
                author_obj["familyName"] = parts[-1]
            record["creator"].append(author_obj)

    # Collection size
    collection_size = _build_collection_size(data)
    if collection_size:
        record["collectionSize"] = collection_size

    # Sample demographics
    record["sample"] = _build_sample(data)

    # Temporal coverage
    tc = data.get("temporal_coverage")
    if tc:
        record["temporalCoverage"] = {
            "@type": "TemporalInterval",
            "startDate": tc["startDate"],
            "endDate": tc["endDate"],
            "temporalType": "collection",
        }

    # Health conditions (dynamically inferred)
    health_conditions = _infer_health_conditions(data)
    if health_conditions:
        record["healthCondition"] = health_conditions

    return record


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def parse():
    """Discover, download, and parse MWCCS substudy PDFs into DataCollection records."""
    count = 0

    # Step 1: Discover PDF links from the substudy page
    logger.info("Discovering substudy PDF links from %s", SUBSTUDY_PAGE_URL)
    html = _fetch_text(SUBSTUDY_PAGE_URL)
    parser = _PdfLinkParser(base_url=SUBSTUDY_PAGE_URL)
    parser.feed(html)
    pdf_links = parser.links

    if not pdf_links:
        logger.error("No substudy PDF links found on %s", SUBSTUDY_PAGE_URL)
        return

    logger.info("Found %d substudy PDF links", len(pdf_links))

    # Step 2: Download each PDF, parse, and yield DataCollection records
    for filename, url in sorted(pdf_links.items()):
        logger.info("Downloading and parsing %s", filename)
        try:
            pdf_bytes = _fetch_bytes(url)
            pdf_text = _read_pdf_text(pdf_bytes)
            data = _parse_pdf_data(pdf_text, url, filename)

            if not data or not data.get("name"):
                logger.warning("Could not extract study name from %s, skipping", filename)
                continue

            yield _build_data_collection(data)
            count += 1
            logger.info("Parsed %s records", count)

        except Exception as exc:
            logger.error("Failed to process %s: %s", filename, exc)
            continue

    logger.info("Finished Parsing. Total Records: %s", count)


if __name__ == "__main__":
    for record in parse():
        print(record)
