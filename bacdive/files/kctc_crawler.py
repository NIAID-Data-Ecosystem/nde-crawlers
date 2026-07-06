"""KCTC (Korean Collection for Type Cultures) HTML parser.

KCTC serves strain records only as server-rendered HTML -- there is no JSON API
behind the catalogue page. Each strain lives at:

    https://kctc.kribb.re.kr/en/collection/view?sn=<sn>

This module fetches that page and parses the detail table into a clean
``{label: value}`` dict.
"""

import logging

import requests
from bs4 import BeautifulSoup

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

VIEW_URL = "https://kctc.kribb.re.kr/en/collection/view"

# KCTC's server returns HTTP 400 for the default "python-requests" User-Agent,
# so send a browser-like one.
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
    )
}


def fetch_kctc_html(sn, timeout=30):
    """GET the KCTC catalogue page for a strain serial number (``sn``).

    Returns the HTML string, or None on error (logged, not raised, so one bad
    record cannot abort a crawl).
    """
    try:
        resp = requests.get(VIEW_URL, params={"sn": sn}, headers=HEADERS, timeout=timeout)
    except requests.RequestException as e:
        logger.warning(f"KCTC sn={sn}: request failed: {e}")
        return None
    if resp.status_code != 200:
        logger.warning(f"KCTC sn={sn}: HTTP {resp.status_code}")
        return None
    return resp.text


def _clean(text):
    """Collapse all runs of whitespace (incl. newlines) to single spaces."""
    return " ".join(text.split())


def parse_kctc_html(html):
    """Parse a KCTC catalogue page into a clean ``{label: value}`` dict.

    The detail table is a simple two-column layout of ``<th>label</th><td>value</td>``
    rows. Labels are kept exactly as shown on the page, e.g. 'KCTC No.', 'Name',
    'Type Strain', 'Biosafty Level' (their spelling), 'Other Collection No.',
    'Oxygen Requirement', 'Temperature', 'pH', 'KCTC Media No.', 'Price',
    'MTA Restrictions'. Empty cells are kept as ''.
    """
    soup = BeautifulSoup(html, "html.parser")
    record = {}
    for row in soup.find_all("tr"):
        th = row.find("th")
        td = row.find("td")
        if th is None or td is None:
            continue
        label = _clean(th.get_text(" "))
        value = _clean(td.get_text(" "))
        if label and label not in record:
            record[label] = value
    return record


def get_kctc_data(sn):
    """Fetch + parse a single KCTC strain by ``sn``.

    Returns the clean ``{label: value}`` dict, or None if the page could not be fetched.
    """
    html = fetch_kctc_html(sn)
    if html is None:
        return None
    return parse_kctc_html(html)


if __name__ == "__main__":
    import json

    print(json.dumps(get_kctc_data(19022), indent=2, ensure_ascii=False))
