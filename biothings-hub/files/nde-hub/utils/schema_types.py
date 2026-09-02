"""Resolve JSON-LD class inheritance through the Data Discovery Engine.

The DDE registry stores a transitive ``parent_classes`` path for each class, used
to test whether an emitted ``@type`` is an expected type or a descendant of one.

Exact matches need no lookup. Resolving a descendant queries the DDE; an
unreachable registry raises ``ClassRegistryUnavailable``.
"""

import json
import logging
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

logger = logging.getLogger(__name__)

DDE_CLASS_API_ROOT = "https://discovery.biothings.io/api/registry"
DDE_CLASS_API_TIMEOUT = 5

# Unprefixed NDE documents use schema.org names as well as NDE/NIAID extension
# names. Try the standard vocabulary first, then the extension registries used
# by NDE_schema.jsonld. A prefixed or full-URI type is sent only to its registry.
_DDE_CLASS_REGISTRIES = (
    "schema",
    "nde",
    "niaid",
    "bioschemas",
    "bioschemastypes",
    "bioschemastypesdrafts",
)
_DDE_PREFIX_TO_REGISTRY = {prefix: prefix for prefix in _DDE_CLASS_REGISTRIES}

# Successful lookups only; failures are not cached.
_PARENT_CACHE: dict[str, frozenset] = {}


class ClassRegistryUnavailable(RuntimeError):
    """The DDE class registry could not be reached."""


def is_type_or_descendant(actual_types, expected_types) -> bool:
    """Return whether any actual JSON-LD type is expected or inherits from it."""
    if isinstance(expected_types, str):
        expected_types = [expected_types]
    expected = {_local_type_name(value) for value in expected_types if isinstance(value, str)}
    expected.discard(None)

    if isinstance(actual_types, (list, tuple, set)):
        candidates = actual_types
    else:
        candidates = [actual_types]

    for actual_type in candidates:
        if not isinstance(actual_type, str):
            continue
        local_name = _local_type_name(actual_type)
        if local_name in expected:
            return True
        if expected.intersection(_parent_types(actual_type)):
            return True
    return False


def _parent_types(type_name: str) -> frozenset:
    """Return the type's transitive parents. Raises ClassRegistryUnavailable on outage."""
    if type_name in _PARENT_CACHE:
        return _PARENT_CACHE[type_name]

    try:
        parents = _resolve_parents(type_name)
    except (OSError, URLError, json.JSONDecodeError) as error:
        logger.error(
            "DDE class registry unreachable at %s; cannot resolve JSON-LD class %r: %s",
            DDE_CLASS_API_ROOT,
            type_name,
            error,
        )
        raise ClassRegistryUnavailable(
            f"Unable to resolve JSON-LD class {type_name!r} through the DDE class registry"
        ) from error

    _PARENT_CACHE[type_name] = parents
    return parents


def _resolve_parents(type_name: str) -> frozenset:
    """Query the DDE for a type's transitive parents. Raises if it is unreachable."""
    for registry, curie in _registry_candidates(type_name):
        class_info = _fetch_dde_class(registry, curie)
        if class_info is None:
            continue

        parents = set()
        parent_paths = class_info.get("parent_classes") or []
        if isinstance(parent_paths, str):
            parent_paths = [parent_paths]
        for path in parent_paths:
            if not isinstance(path, str):
                continue
            for parent in path.split(","):
                if local_name := _local_type_name(parent.strip()):
                    parents.add(local_name)
        return frozenset(parents)

    return frozenset()


def _fetch_dde_class(registry, curie):
    """Return one DDE class record, or ``None`` when that registry has no match."""
    encoded_curie = quote(curie, safe="")
    request = Request(
        f"{DDE_CLASS_API_ROOT}/{registry}/{encoded_curie}",
        headers={"Accept": "application/json", "User-Agent": "nde-hub-schema-validator/1.0"},
    )
    try:
        with urlopen(request, timeout=DDE_CLASS_API_TIMEOUT) as response:
            return json.load(response)
    except HTTPError as error:
        if error.code == 404:
            return None
        raise


def _registry_candidates(type_name):
    """Yield DDE registry/CURIE pairs for a compact type, CURIE, or full URI."""
    type_name = type_name.strip()
    local_name = _local_type_name(type_name)
    if not local_name:
        return ()

    prefix = _type_prefix(type_name)
    if prefix in _DDE_PREFIX_TO_REGISTRY:
        registry = _DDE_PREFIX_TO_REGISTRY[prefix]
        return ((registry, f"{prefix}:{local_name}"),)

    return tuple((registry, f"{registry}:{local_name}") for registry in _DDE_CLASS_REGISTRIES)


def _type_prefix(type_name):
    """Extract the vocabulary prefix from a CURIE or recognized DDE/schema URI."""
    if type_name.startswith(("http://schema.org/", "https://schema.org/")):
        return "schema"

    marker = "discovery.biothings.io/ns/"
    if marker in type_name:
        remainder = type_name.split(marker, 1)[1]
        return remainder.split("/", 1)[0]

    if ":" in type_name and not type_name.startswith(("http://", "https://")):
        return type_name.split(":", 1)[0]
    return None


def _local_type_name(type_name):
    """Normalize ``schema:WebPage`` and full URIs to the local name ``WebPage``."""
    if not isinstance(type_name, str) or not type_name.strip():
        return None

    type_name = type_name.strip().rstrip("/#")
    if type_name.startswith(("http://", "https://")):
        return type_name.rsplit("/", 1)[-1].rsplit("#", 1)[-1]
    if ":" in type_name:
        return type_name.split(":", 1)[1]
    return type_name
