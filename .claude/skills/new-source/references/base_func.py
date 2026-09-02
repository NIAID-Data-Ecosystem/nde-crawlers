def insert_value(d, key, value, extend=False):
    """ Insert a value into a dictionary, promoting to a list on repeats and dropping duplicates.
    """

    if extend:
        d[key] = (d[key] + " " + value).strip() if d.get(key) else value
        return

    was_list = isinstance(d.get(key), list)
    merged = list(d[key]) if was_list else ([d[key]] if key in d else [])
    for item in value if isinstance(value, list) else [value]:
        if item not in merged:
            merged.append(item)
    d[key] = merged if was_list or isinstance(value, list) or len(merged) > 1 else merged[0]

def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(str(val), ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError, OverflowError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt
