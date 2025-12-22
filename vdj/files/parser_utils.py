import re
from typing import Any, Dict, Iterator, List, Sequence


def _ensure_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _add_unique_dict(target_list: List[Dict[str, Any]], seen_keys: set, data: Dict[str, Any]) -> None:
    cleaned = {k: v for k, v in data.items() if v not in (None, "", [], {}, set())}
    if not cleaned:
        return
    key = tuple(sorted((k, tuple(v) if isinstance(v, Sequence) and not isinstance(v, (str, bytes)) else v) for k, v in cleaned.items()))
    if key in seen_keys:
        return
    seen_keys.add(key)
    target_list.append(cleaned)


def _add_unique_value(target_list: List[Any], seen_values: set, value: Any) -> None:
    if value in (None, ""):
        return
    norm = value.strip() if isinstance(value, str) else value
    if norm in (None, ""):
        return
    if norm in seen_values:
        return
    seen_values.add(norm)
    target_list.append(norm)


def _iter_string_values(value: Any) -> Iterator[str]:
    if value in (None, "", [], {}, set()):
        return
    if isinstance(value, list):
        for item in value:
            yield from _iter_string_values(item)
        return
    if isinstance(value, dict):
        for key in ("label", "name", "value", "id"):
            if key in value and value[key] not in (None, ""):
                yield from _iter_string_values(value[key])
                return
        return
    if isinstance(value, str):
        parts = [part.strip() for part in re.split(r"[;,]", value) if part.strip()]
        if parts:
            for part in parts:
                yield part
            return
        yield value.strip()
        return
    yield value


def _sanitize_identifier(value: Any, fallback: str = "entry") -> str:
    """Normalize identifiers so they are filesystem/URL friendly."""
    if value in (None, ""):
        return fallback
    if not isinstance(value, str):
        value = str(value)
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value)
    cleaned = cleaned.strip("_")
    return cleaned or fallback
