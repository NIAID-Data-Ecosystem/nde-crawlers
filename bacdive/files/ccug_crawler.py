"""CCUG (Culture Collection University of Gothenburg) HTML parser.

CCUG serves strain records only as server-rendered HTML -- there is no JSON API
behind the catalogue page. Each strain lives at:

    https://www.ccug.se/strain?id=<id>

This module fetches that page and parses it into a clean ``{label: value}`` dict.
"""

import logging

import requests
from bs4 import BeautifulSoup

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

VIEW_URL = "https://www.ccug.se/strain"


def fetch_ccug_html(ccug_id, timeout=30):
    """GET the CCUG catalogue page for a strain id.

    Returns the HTML string, or None on error (logged, not raised, so one bad
    record cannot abort a crawl).
    """
    try:
        resp = requests.get(VIEW_URL, params={"id": ccug_id}, timeout=timeout)
    except requests.RequestException as e:
        logger.warning(f"CCUG id={ccug_id}: request failed: {e}")
        return None
    if resp.status_code != 200:
        logger.warning(f"CCUG id={ccug_id}: HTTP {resp.status_code}")
        return None
    return resp.text


def _clean(text):
    """Collapse all runs of whitespace (incl. newlines and nbsp) to single spaces."""
    return " ".join(text.split())


def parse_ccug_html(html):
    """Parse a CCUG catalogue page into a clean ``{label: value}`` dict.

    'Strain' and 'Name' come from the ``<h1>`` ("CCUG <n> - <organism name>"). The
    remaining fields come from the ``result_table`` blocks, where a label cell is
    marked ``class="highlight"`` and the next cell holds the value -- e.g. 'Sample
    Origin', 'Source', 'Depositor', 'Deposit Date', 'Other Collections', 'Tests OX',
    'Api ...'. Rows whose first cell is not a highlight label (the fatty-acid PEAK
    matrix) are skipped. Empty cells are kept as ''.
    """
    soup = BeautifulSoup(html, "html.parser")
    record = {}

    # Strain number + organism name: <h1> is "CCUG <n> - <organism name>".
    if h1 := soup.find("h1"):
        strain, _, name = _clean(h1.get_text()).partition(" - ")
        record["Strain"] = strain
        if name:
            record["Name"] = name

    # Detail rows: first cell is a class="highlight" label, second is the value.
    # (The fatty-acid PEAK matrix uses plain <td> cells, so it is skipped.)
    for table in soup.find_all("table", class_="result_table"):
        for row in table.find_all("tr"):
            cells = row.find_all("td")
            if len(cells) >= 2 and "highlight" in cells[0].get("class", []):
                label = _clean(cells[0].get_text()).rstrip(":").strip()
                value = _clean(cells[1].get_text())
                if label and label not in record:
                    record[label] = value

    return record


def get_ccug_data(ccug_id):
    """Fetch + parse a single CCUG strain by id.

    Returns the clean ``{label: value}`` dict, or None if the page could not be fetched.
    """
    html = fetch_ccug_html(ccug_id)
    if html is None:
        return None
    return parse_ccug_html(html)


if __name__ == "__main__":
    import json

    print(json.dumps(get_ccug_data(15571), indent=2, ensure_ascii=False))
