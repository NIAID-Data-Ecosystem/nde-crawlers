# Advanced insert_value. Use instead of base_func.insert_value when a single key
# accumulates thousands of values (aggregating many sub-records into one document).
# Dedupes against a hash set instead of scanning the list: O(n) instead of O(n^2).
#
#     aggregate = {}
#     seen = {}                                  # one seen dict per output dict
#     for s in samples:
#         insert_value(aggregate, "url", s["url"], seen=seen)
#         insert_value(aggregate, "sampleType", {"@type": "DefinedTerm", "name": s["type"]}, seen=seen)
#
# One `seen` covers every key of that dict. Reusing a `seen` across two output
# dicts silently drops values. Omit `seen` for one-off fields.

import json


def _marker(value):
    """ Hashable stand-in for value; equal values give equal markers.
    """

    try:
        return json.dumps(value, sort_keys=True, default=repr)
    except TypeError:
        return repr(value)


def insert_value(d, key, value, extend=False, seen=None):
    """ Insert a value into a dictionary, promoting to a list on repeats and dropping duplicates.

    Pass ``seen``, a dict reused across calls for the same ``d``, to dedupe against a
    hash set instead of scanning the existing list. Output is identical either way.
    """

    if extend:
        d[key] = (d[key] + " " + value).strip() if d.get(key) else value
        return

    was_list = isinstance(d.get(key), list)
    if not was_list:
        merged = [d[key]] if key in d else []
    elif seen is not None and key in seen:
        merged = d[key]        # already copied on first touch; append in place
    else:
        merged = list(d[key])  # copy once: d[key] may still be owned by the caller

    incoming = value if isinstance(value, list) else [value]
    if seen is None:
        for item in incoming:
            if item not in merged:
                merged.append(item)
    else:
        if key not in seen:
            # prime from merged, not seen.setdefault: an eager default rebuilds
            # this set every call and puts the scan back
            seen[key] = {_marker(item) for item in merged}
        markers = seen[key]
        for item in incoming:
            marker = _marker(item)
            if marker not in markers:
                markers.add(marker)
                merged.append(item)

    d[key] = merged if was_list or isinstance(value, list) or len(merged) > 1 else merged[0]
