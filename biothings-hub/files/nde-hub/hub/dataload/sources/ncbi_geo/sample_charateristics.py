import json
from pathlib import Path

import regex as re
from utils.sex import _parse_mf_list, extract_sex, find_sex_number_map

try:
    from config import logger
except ImportError:
    import logging

    logger = logging.getLogger(__name__)


def insert_value(d, key, value):
    if key in d:
        if isinstance(d[key], list) and value not in d[key]:
            d[key].append(value)
        if not isinstance(d[key], list) and d[key] != value:
            d[key] = [d[key], value]
    else:
        d[key] = value


def find_age(text):
    """
    Find an age integer in a string.
    Matches patterns like:
      - "age: 60"
      - "age 60"
      - "age-60"
      - "male;age: 20"
      - "60 years old"
    Returns int age or None.
    """
    if text is None:
        return None
    s = str(text).lower()
    # common "age" followed by optional "(years)" then separator/space then number
    m = re.search(
        r"\bage\b\s*(?:\(\s*years?\s*\)\s*)?(?:\s*[:=\-]\s*|\s+)(\d{1,3})\b",
        s,
    )
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    # fallback: "<number> years old" or "<number> yrs"
    m2 = re.search(r"\b(\d{1,3})\s*(?:years?|yrs?|y)\b", s)
    if m2:
        try:
            return int(m2.group(1))
        except ValueError:
            return None
    return None


def apply_heuristics(subproperty, value):
    if value is None:
        return None
    if subproperty == "sex":
        if value.startswith("m") and not _parse_mf_list(value):
            if find_age(value):
                return "Male", {"value": find_age(value), "unitText": "year"}
        if value.startswith("f") and not _parse_mf_list(value):
            if find_age(value):
                return "Female", {"value": find_age(value), "unitText": "year"}
    if subproperty == "sex/age":
        m = re.match(r"^\s*([mf])\s*/\s*(\d{1,3})\s*$", value)
        if m:
            sex = "Male" if m.group(1) == "m" else "Female"
            try:
                age = int(m.group(2))
            except ValueError:
                age = None
            if age:
                return sex, {"value": age, "unitText": "year"}
            else:
                return sex

    if subproperty == "sex-age":
        m = re.match(r"^\s*([mf])\s*-\s*(\d{1,3})\s*$", value)
        if m:
            sex = "Male" if m.group(1) == "m" else "Female"
            try:
                age = int(m.group(2))
            except ValueError:
                age = None
            if age:
                return sex, {"value": age, "unitText": "year"}
            else:
                return sex

    if subproperty == "sex/mating type":
        # leading sex token like: "f; age:24 to 28 weeks;"
        m_sex = re.match(r"^\s*([mf])\b", value)
        if not m_sex:
            return None
        sex = "Male" if m_sex.group(1) == "m" else "Female"

        def _normalize_unit(u: str) -> str:
            u = (u or "").lower()
            if u in {"wk", "wks"} or u.startswith("week"):
                return "weeks"
            if u in {"mo", "mos"} or u.startswith("month"):
                return "months"
            if u in {"y", "yr", "yrs"} or u.startswith("year"):
                return "years"
            return u

        # specifically: has "age" and "to" (range), with weeks/months/years
        m_range = re.search(
            r"\bage\b[^0-9]*?(\d{1,3})\s*to\s*(\d{1,3})\s*" r"(weeks?|wks?|wk|months?|mos?|mo|years?|yrs?|yr|y)\b",
            value,
        )
        if m_range:
            minv = int(m_range.group(1))
            maxv = int(m_range.group(2))
            unit = _normalize_unit(m_range.group(3))
            return sex, {
                "minValue": minv,
                "maxValue": maxv,
                "unitText": unit,
            }

        else:
            return sex

    return None


def parse_sex(subproperty, value, mapping):
    """Parses the sex from the given subproperty and value.
    Returns either:
      - Tuple of (sex, developmentalStage dict)
      - sex string e.g. "Female"
      - List of sex strings e.g. ["Male", "Female"]
      - None if unable to parse
    """
    subproperty = subproperty.strip().lower()
    if isinstance(value, str):
        value = value.strip().lower()

    valid_subproperties = [
        "sex   male=1;female=2",
        "sex  (f /m)",
        "sex (0 = female, 1 = male)",
        "sex (0-female 1-male)",
        "sex (0-male 1-female)",
        "sex (0=f;1=m)",
        "sex (0=female, 1=male)",
        "sex (0=females; 1=male)",
        "sex (0=m)",
        "sex (0=male, 1=female)",
        "sex (1 - male, 2 - female)",
        "sex (1= female, 2= male)",
        "sex (1=male,2=female)",
        "sex (1=male;_2=female)",
        "sex (1=woman; 0=man)",
        "sex (clinical gender)",
        "sex (f=female, m=male)",
        "sex (female/male)",
        "sex (female=0 male=1)",
        "sex (female=0, male=1)",
        "sex (female=1, male=0)",
        "sex (fetal)",
        "sex (gender)",
        "sex (h. sapiens)",
        "sex (inferred_from_sequencing_data)",
        "sex (m, male; f, female)",
        "sex (m/f)",
        "sex (m=male, f=female)",
        "sex (m=male; f=female)",
        "sex (male.1_female.2)",
        "sex (male/female)",
        "sex (male= 1, female=2)",
        "sex (self-reported)",
        "sex (singleton or multiplexed)",
        "sex 3",
        "sex and age",
        "sex assigned_at_birth",
        "sex at_birth",
        "sex by xist counts",
        "sex chromosome complement",
        "sex chromosome ploidy",
        "sex chromosome",
        "sex chromosomes",
        "sex composition",
        "sex cv",
        "sex from y chromosome gene expression",
        "sex id",
        "sex known",
        "sex of chicken",
        "sex of child",
        "sex of mice",
        "sex of patient",
        "sex of the mice",
        "sex ontology",
        "sex predicted",
        "sex type",
        "sex_(1 = female, 0 = male)",
        "sex_age",
        "sex-age",
        "sex_chry",
        "sex_code",
        "sex_dscrp",
        "sex_nomorph",
        "sex_ontology",
        "sex_predicted",
        "sex_recipient",
        "sex, 1=m, 2=f",
        "sex,strain",
        "sex; female=1, male=2",
        "sex",
        "sex(0=female, 1=male)",
        "sex/age",
        "sex/genotype",
        "sex/germline",
        "sex/karyotype",
        "sex/mating type",
        "sex/pregnancy status",
        "sex/stage",
        "sexe",
        "sexo",
        "sexual phenotype",
    ]

    if subproperty not in valid_subproperties:
        return None

    if subproperty in ["sex", "sex/age", "sex-age", "sex/mating type"]:
        sex = apply_heuristics(subproperty, value)
        if sex:
            return sex

    if sex_map := find_sex_number_map(subproperty, ["=", "-", "."]):
        try:
            idx = int(value)
            sex = sex_map.get(idx)
            if sex:
                return sex
        except (ValueError, TypeError):
            pass

    # with open(Path(__file__).resolve().parent / "sex_mappings.json", "r") as f:
    #     mapping = json.load(f)
    return extract_sex(value, mapping)


def parse_sample_characteristics(output, value, sample_mapping, nde_mapping, sex_mapping):
    values = value if isinstance(value, list) else [value]
    for v in values:
        # Split by ":", strip whitespace, and only split on the first ":"
        parts = [p.strip() for p in v.split(":", 1)]
        if len(parts) != 2:
            logger.warning(f"Invalid sample characteristic format: {v}")
            continue

        subproperty, field_value = parts[0], parts[1]

        sex = parse_sex(subproperty, field_value, sex_mapping)
        if isinstance(sex, tuple):
            if sex[0] and isinstance(sex[0], list):
                for s in sex[0]:
                    insert_value(output, "sex", s)
            else:
                insert_value(output, "sex", sex[0])
            if sex[1]:
                insert_value(output, "developmentalStage", sex[1])
        elif sex and isinstance(sex, list):
            for s in sex:
                insert_value(output, "sex", s)
        elif sex:
            insert_value(output, "sex", sex)

        if sex:
            continue  # Skip further

        # with open(Path(__file__).resolve().parent / "nde_mapping.json", "r") as f:
        #     nde_mapping = json.load(f)
        # with open(Path(__file__).resolve().parent / "mapping_dict.json", "r") as f:
        #     mapping = json.load(f)

        subproperty = subproperty.strip().lower().replace(" ", "_")
        if subproperty in sample_mapping:
            k, v = sample_mapping[subproperty]
            if v:
                v = subproperty if v == "subproperty" else field_value
                if k in nde_mapping and nde_mapping[k][0] == "object":
                    d = {nde_mapping[k][1]: v}
                    if nde_mapping[k][1] == "sampleQuantity":
                        d["name"] = subproperty
                    if nde_mapping[k][1] == "variableMeasuered" or nde_mapping[k][1] == "anatomicalStructure":
                        d["@type"] = "DefinedTerm"
                    insert_value(output, k, d)
                    logger.info(f"Mapped sample characteristic subproperty: {subproperty} to {k} with value: {d}")
                elif k in nde_mapping and nde_mapping[k][0] == "value":
                    insert_value(output, k, v)
                    logger.info(f"Mapped sample characteristic subproperty: {subproperty} to {k} with value: {v}")
                else:
                    logger.warning(f"Unmapped nde_mapping property: {k}")
            else:
                d = {"@type": "PropertyValue", "propertyID": subproperty, "value": field_value}
                insert_value(output, "additionalProperty", d)

        else:
            logger.warning(f"Unmapped sample characteristic subproperty: {subproperty}")
