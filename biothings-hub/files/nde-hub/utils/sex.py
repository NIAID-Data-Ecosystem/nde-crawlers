from typing import Dict, Iterable, Optional, Pattern

import regex as re


def _parse_mf_list(value: str):
    """
    Accepts strings like:
      - "f;m;f"  -> ["Female", "Male"]
      - "f/f/f"  -> "Female"
      - "m=0"    -> invalid (not a pure m/f list)
      - "male;m;female" -> invalid
      - "m;m.f"  -> invalid (mixed separators)

    Rule: the entire string must be ONLY m/f tokens separated by ONE single
    non-alphanumeric, non-space symbol (consistent throughout).
    """
    s = str(value).strip()
    if not s:
        return None

    # Find all "separator" characters present (punctuation-like): not word, not whitespace
    seps = {ch for ch in s if re.match(r"[^\w\s]", ch)}
    if not seps:
        # No separator => allow a single token "m" or "f" only
        t = s.strip().lower()
        if t == "m":
            return "Male"
        if t == "f":
            return "Female"
        return None

    if len(seps) != 1:
        return None  # mixed separators, e.g. "m;m.f"

    sep = next(iter(seps))

    parts = [p.strip().lower() for p in s.split(sep)]
    if any(p == "" for p in parts):
        return None

    # Must be purely m/f tokens (no "male", "female", numbers, etc.)
    if any(p not in ("m", "f") for p in parts):
        return None

    canon = ["Male" if p == "m" else "Female" for p in parts]
    return list(set(canon))


def extract_sex(value: str, mapping=None):
    """Given a string value and optional mapping dict, return a standardized sex label or list of labels."""
    v = value.lower().strip()

    # Strict m/f list handling (and reject other separator-based formats)
    mf = _parse_mf_list(v)
    if mf is not None:
        return mf

    if mapping and v in mapping:
        return mapping[v]

    patterns = {
        "male": re.compile(r"\b(?:man|men|boy|boys|xy|male|males)\b", re.I),
        "female": re.compile(r"\b(?:woman|women|girl|girls|xx|female|females)\b", re.I),
        "intersex": re.compile(r"\b(?:hermaphrodite|hermaphrodites)\b", re.I),
    }

    labels = []
    for label, pat in patterns.items():
        if pat.search(v):
            labels.append(label.capitalize())

    if len(labels) == 1:
        return labels[0]
    if len(labels) > 1:
        return labels
    return None


def add_sex(value, output, mapping=None):
    sex_value = extract_sex(value, mapping)
    if not sex_value:
        return

    # normalize to list of labels
    sex_values = sex_value if isinstance(sex_value, list) else [sex_value]

    existing = output.get("sex")
    if existing is None:
        output["sex"] = sex_values.copy()
        return

    if not isinstance(existing, list):
        existing = [existing]

    output["sex"] = list(set(existing + sex_values))


def build_sex_number_regex(
    allowed_symbols: Optional[Iterable[str]] = None,
) -> Pattern:
    """
    Matches (either order):
      - number <symbol> (m|male|males|f|female|females)   e.g. 1=m, 2=females
      - (m|male|males|f|female|females) <symbol> number   e.g. males=0, female-2, f.2

    If allowed_symbols is None: allow ANY non-alphanumeric, non-space separator (punctuation),
    e.g. '=', '-', '.', ':', '/', etc.
    If allowed_symbols is provided: only those symbols are allowed (supports 1+ char symbols).
    """
    if allowed_symbols is None:
        sep = r"[^\w\s]+"  # any punctuation-like separator(s)
    else:
        syms = list(allowed_symbols)
        if not syms:
            raise ValueError("allowed_symbols must be None or a non-empty iterable of symbols")
        sep = r"(?:%s)" % "|".join(re.escape(s) for s in syms)

    # male/female labels (case-insensitive), including plurals
    label = r"(?:m|male(?:s)?|f|female(?:s)?)"

    left = rf"(?P<num>\d+)\s*(?P<sep>{sep})\s*(?P<label>{label})"
    right = rf"(?P<label2>{label})\s*(?P<sep2>{sep})\s*(?P<num2>\d+)"
    return re.compile(rf"\b(?:{left}|{right})\b", re.I)


def find_sex_number_map(
    text: str,
    allowed_symbols: Optional[Iterable[str]] = None,
) -> Dict[int, str]:
    """
    Example:
      text = "1=m female=2"
      -> {1: "Male", 2: "Female"}
    """
    pat = build_sex_number_regex(allowed_symbols)

    def canonical(label_raw: str) -> str:
        s = label_raw.strip().lower()
        # normalize to 'male'/'female' based on first letter
        return "Male" if s.startswith("m") else "Female"

    out = {}
    for m in pat.finditer(text):
        label_raw = m.group("label") or m.group("label2")
        num = int(m.group("num") or m.group("num2"))
        out[num] = canonical(label_raw)  # if repeated, last one wins
    return out
