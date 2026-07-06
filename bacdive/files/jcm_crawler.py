"""JCM (Japan Collection of Microorganisms, RIKEN BRC) HTML parser.

JCM serves strain records only as server-rendered HTML from a CGI script -- there
is no JSON API behind it. Each strain lives at:

    https://www.jcm.riken.jp/cgi-bin/jcm/jcm_number?JCM=<jcm>

This module fetches that page and parses it into a clean ``{label: value}`` dict.
"""

import logging

import requests
from bs4 import BeautifulSoup

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

VIEW_URL = "https://www.jcm.riken.jp/cgi-bin/jcm/jcm_number"


def fetch_jcm_html(jcm, timeout=30):
    """GET the JCM catalogue page for a JCM number.

    Returns the HTML string, or None on error (logged, not raised, so one bad
    record cannot abort a crawl).
    """
    try:
        resp = requests.get(VIEW_URL, params={"JCM": jcm}, timeout=timeout)
    except requests.RequestException as e:
        logger.warning(f"JCM {jcm}: request failed: {e}")
        return None
    if resp.status_code != 200:
        logger.warning(f"JCM {jcm}: HTTP {resp.status_code}")
        return None
    return resp.text


def _clean(text):
    """Collapse all runs of whitespace (incl. newlines and nbsp) to single spaces."""
    return " ".join(text.split())


def parse_jcm_html(html):
    """Parse a JCM catalogue page into a clean ``{label: value}`` dict.

    The page mixes two layouts:

    * A free-text header of "Label: value" segments -- e.g. 'Category', 'Medium',
      'Temperature', 'Rehydration fluid', 'Source', 'BacDive ID', and (on type
      strains) many more like 'G+C (mol%)', 'Phylogeny', 'Genome sequence'. Some
      are packed onto one line separated by ';'.
    * Two ``border="1"`` tables -- the delivery/use info block (Biosafety level,
      Terms and conditions, Export control, ...) and the delivery category block
      (Domestic / Overseas). A separate ``border="0"`` table only lays out the
      header, so it is skipped.

    'Name' holds the organism name shown in bold at the top. Field coverage varies
    a lot between minimal entries and fully-described type strains.
    """
    soup = BeautifulSoup(html, "html.parser")
    for tag in soup(["script", "style"]):
        tag.decompose()

    record = {}

    # Organism name: the bold species heading at the top of the entry.
    if name := soup.find("strong"):
        record["Name"] = _clean(name.get_text())

    # Header block: old-style markup separates fields with <br>/<hr>. Turn those
    # into newlines, then read "Label: value" segments (a few are packed onto one
    # line separated by ';'). First occurrence of a label wins.
    for br in soup.find_all(["br", "hr"]):
        br.replace_with("\n")
    for line in soup.get_text().split("\n"):
        for part in line.split(";"):
            if ":" in part:
                label, _, value = part.partition(":")
                label, value = _clean(label), _clean(value).rstrip(".").strip()
                if label and value and label not in record:
                    record[label] = value

    # Data tables use border="1" (label in the first cell, value in the second).
    for table in soup.find_all("table", border="1"):
        for row in table.find_all("tr"):
            cells = row.find_all("td")
            if len(cells) >= 2:
                label = _clean(cells[0].get_text())
                value = _clean(cells[1].get_text())
                if label and label not in record:
                    record[label] = value

    return record


def get_jcm_data(jcm):
    """Fetch + parse a single JCM strain by JCM number.

    Returns the clean ``{label: value}`` dict, or None if the page could not be fetched.
    """
    html = fetch_jcm_html(jcm)
    if html is None:
        return None
    return parse_jcm_html(html)


if __name__ == "__main__":
    import json

    print(json.dumps(get_jcm_data(19270), indent=2, ensure_ascii=False))
