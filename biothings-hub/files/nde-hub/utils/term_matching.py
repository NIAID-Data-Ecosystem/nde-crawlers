"""Conservative mention matching shared by citation and description enrichment."""

import re
from functools import lru_cache


SAFE_SHORT_MENTIONS = frozenset(
    {
        "aids",
        "covid19",
        "ebv",
        "h1n1",
        "h3n2",
        "h5n1",
        "hbv",
        "hcv",
        "hiv",
        "hpv",
        "hsv",
        "mers",
        "merscov",
        "rsv",
        "sars",
        "sarscov",
        "sarscov2",
        "tb",
        "zikv",
    }
)

_IRREGULAR_SPECIES_PLURALS = {
    "bacteria": "bacterium",
    "fungi": "fungus",
    "mice": "mouse",
}


def normalize_term_text(value):
    """Normalize a mention or ontology label for conservative comparison."""
    return " ".join(re.findall(r"[\w]+", str(value or "").casefold()))


def compact_term_text(value):
    return "".join(character for character in str(value or "").casefold() if character.isalnum())


def is_ambiguous_short_mention(name):
    """True when a short entity mention is unsafe without disambiguation."""
    compact = compact_term_text(name)
    if not compact or compact in SAFE_SHORT_MENTIONS:
        return False
    return len(compact) <= 3 or (str(name).strip().isupper() and len(compact) <= 8)


def term_labels(term):
    """Yield authoritative names without trusting the cache's input name."""
    for field in ("name", "commonName"):
        if term.get(field):
            yield term[field]

    alternate_names = term.get("alternateName")
    if alternate_names is not None:
        if not isinstance(alternate_names, list):
            alternate_names = [alternate_names]
        yield from alternate_names

    # UniProt display names commonly have the form "Human | Homo sapiens".
    if display_name := term.get("displayName"):
        yield from (part.strip() for part in str(display_name).split("|") if part.strip())


def term_expansion_mentioned(term, mention, haystacks):
    """True when context spells out a label beyond an ambiguous acronym."""
    compact_mention = compact_term_text(mention)
    return any(
        compact_term_text(label) != compact_mention and mentioned_in(label, haystacks) for label in term_labels(term)
    )


def term_matches_mention(term, mention):
    """Require a standardized label or synonym to agree with its mention."""
    normalized_mention = normalize_term_text(mention)
    compact_mention = compact_term_text(mention)
    if not normalized_mention:
        return False

    for label in term_labels(term):
        normalized_label = normalize_term_text(label)
        if normalized_label == normalized_mention:
            return True
        # Accept punctuation-only variants such as SARS-CoV-2 vs SARS-CoV2,
        # but never use substring matching.
        if compact_mention and compact_term_text(label) == compact_mention:
            return True
        if compact_mention in SAFE_SHORT_MENTIONS:
            label_tokens = normalized_label.split()
            acronym = "".join(token[0] for token in label_tokens if token)
            if compact_mention in label_tokens or compact_mention == acronym:
                return True
    return False


def _singularize_species_mention(value):
    """Normalize a simple plural species/common name without fuzzy matching."""
    words = normalize_term_text(value).split()
    if not words:
        return ""

    last_word = words[-1]
    if last_word in _IRREGULAR_SPECIES_PLURALS:
        words[-1] = _IRREGULAR_SPECIES_PLURALS[last_word]
    elif last_word.endswith("ies") and len(last_word) > 4:
        words[-1] = f"{last_word[:-3]}y"
    elif last_word.endswith(("ches", "shes", "uses", "xes", "zes")):
        words[-1] = last_word[:-2]
    elif last_word.endswith("s") and not last_word.endswith(("is", "ss", "us")):
        words[-1] = last_word[:-1]
    return " ".join(words)


def _scientific_name_abbreviation_matches(first, second):
    first_words = normalize_term_text(first).split()
    second_words = normalize_term_text(second).split()
    if len(first_words) != 2 or len(second_words) != 2 or first_words[1] != second_words[1]:
        return False
    first_genus, second_genus = first_words[0], second_words[0]
    return (
        len(first_genus) == 1
        and len(second_genus) > 1
        and first_genus == second_genus[0]
        or len(second_genus) == 1
        and len(first_genus) > 1
        and second_genus == first_genus[0]
    )


def species_term_matches_mention(term, mention):
    """Match a taxon label, allowing plurals and abbreviated genus names."""
    if term_matches_mention(term, mention):
        return True

    labels = list(term_labels(term))
    singular_mention = _singularize_species_mention(mention)
    for label in labels:
        if singular_mention and _singularize_species_mention(label) == singular_mention:
            return True

    # An initial plus species epithet is ambiguous across genera. Expand it
    # only when one authoritative label supplies the full scientific name and
    # another explicitly contains the abbreviated form. For example, UniProt
    # describes Plasmodium vivax as "malaria parasite P. vivax"; a coincidental
    # Phyllostachys vivax candidate has no such label.
    # Compare the corroborating abbreviated label after punctuation
    # normalization so both ``P. falciparum`` and ``P.falciparum`` match the
    # same authoritative synonym. The full scientific-name check above still
    # prevents the initial from expanding to an unrelated genus.
    normalized_abbreviation = normalize_term_text(mention)
    return any(_scientific_name_abbreviation_matches(mention, label) for label in labels) and any(
        mentioned_in(normalized_abbreviation, (normalize_term_text(label),)) for label in labels
    )


@lru_cache(maxsize=16_384)
def mention_pattern(name):
    words = str(name or "").strip().split()
    if not words:
        return None
    escaped = r"\s+".join(re.escape(word) for word in words)
    return re.compile(rf"(?<!\w){escaped}(?!\w)", re.IGNORECASE)


def mentioned_in(name, haystacks):
    """True only for a complete mention, never a substring of another word."""
    pattern = mention_pattern(name)
    return bool(pattern and any(pattern.search(haystack or "") for haystack in haystacks))
