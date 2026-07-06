"""UCCCB (University of Coimbra Culture Collection of Bacteria) HTML parser.

UCCCB serves strain records as server-rendered HTML (no JSON API is needed; the
custom plugin's admin-ajax endpoints are only for the cart). Each strain lives at:

    https://ucccb.uc.pt/strain-details/?detail=<detail>

This module fetches that page and parses it into a clean ``{label: value}`` dict.
"""

import logging

import requests
from bs4 import BeautifulSoup, NavigableString

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

VIEW_URL = "https://ucccb.uc.pt/strain-details/"


def fetch_ucccb_html(detail, timeout=30):
    """GET the UCCCB strain-details page for a ``detail`` id (e.g. "UCCCB10").

    Returns the HTML string, or None on error (logged, not raised, so one bad
    record cannot abort a crawl).
    """
    try:
        resp = requests.get(VIEW_URL, params={"detail": detail}, timeout=timeout)
    except requests.RequestException as e:
        logger.warning(f"UCCCB {detail}: request failed: {e}")
        return None
    if resp.status_code != 200:
        logger.warning(f"UCCCB {detail}: HTTP {resp.status_code}")
        return None
    return resp.text


def _clean(text):
    """Collapse all runs of whitespace (incl. newlines and nbsp) to single spaces."""
    return " ".join(text.split())


def parse_ucccb_html(html):
    """Parse a UCCCB strain-details page into a clean ``{label: value}`` dict.

    Each field is a ``<b class="full-list-unitb">Label</b>`` followed by its value,
    which is either a ``<p class="full-list-unitp">`` element or a bare text node
    between ``<br/>`` tags. Section headings (STRAIN DESIGNATION, GROWTH CONDITIONS,
    ...) are ``<h3>`` and are skipped. Empty values are kept as ''.
    """
    soup = BeautifulSoup(html, "html.parser")
    record = {}
    for b in soup.select("b.full-list-unitb"):
        label = _clean(b.get_text()).rstrip(":")
        if not label:
            continue
        value = ""
        for sib in b.next_siblings:
            name = getattr(sib, "name", None)
            if name == "p":  # value wrapped in <p class="full-list-unitp">
                value = _clean(sib.get_text())
                break
            if isinstance(sib, NavigableString) and sib.strip():  # bare text between <br/>
                value = _clean(sib)
                break
            if name in ("b", "h3", "hr"):  # reached the next label/section without a value
                break
        if label not in record:
            record[label] = value
    return record


def get_ucccb_data(detail):
    """Fetch + parse a single UCCCB strain by ``detail`` id.

    Returns the clean ``{label: value}`` dict, or None if the page could not be fetched.
    """
    html = fetch_ucccb_html(detail)
    if html is None:
        return None
    return parse_ucccb_html(html)


if __name__ == "__main__":
    import json

    print(json.dumps(get_ucccb_data("UCCCB10"), indent=2, ensure_ascii=False))
