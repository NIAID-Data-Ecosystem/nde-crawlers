"""ATCC (American Type Culture Collection) HTML parser.

ATCC has no open product API (its Insite/Optimizely /api/ endpoints require auth),
but the product page is server-rendered. Each product lives at:

    https://www.atcc.org/products/<id>

The site sits behind Cloudflare, so a browser User-Agent is sent. This module
fetches the page and parses it into a clean dict:

* the product-information fields (the simple summary table plus the detailed
  product-information accordion -- Medium, Temperature, depositor, etc.), and
* a "Permits & Restrictions" entry: a list of {title, text} for each permit.

The long "Legal disclaimers" boilerplate (Intended use / Warranty / Disclaimers)
is intentionally skipped.
"""

import json
import logging

import requests
from bs4 import BeautifulSoup

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

# ATCC is behind Cloudflare; the default "python-requests" UA gets a reduced page.
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
    )
}


def fetch_atcc_html(product_id, timeout=30):
    """GET the ATCC product page for a product id.

    Returns the HTML string, or None on error (logged, not raised, so one bad
    record cannot abort a crawl).
    """
    url = f"https://www.atcc.org/products/{product_id}"
    try:
        resp = requests.get(url, headers=HEADERS, timeout=timeout)
    except requests.RequestException as e:
        logger.warning(f"ATCC {product_id}: request failed: {e}")
        return None
    if resp.status_code != 200:
        logger.warning(f"ATCC {product_id}: HTTP {resp.status_code}")
        return None
    return resp.text


def _clean(text):
    """Collapse all runs of whitespace (incl. newlines and nbsp) to single spaces."""
    return " ".join(text.split())


def _load_product_ld(soup):
    """Return the schema.org Product JSON-LD block as a dict, or {}."""
    for s in soup.find_all("script", type="application/ld+json"):
        try:
            data = json.loads(s.string or "")
        except (ValueError, TypeError):
            continue
        if isinstance(data, dict) and data.get("@type") == "Product":
            return data
    return {}


def parse_atcc_html(html):
    """Parse an ATCC product page into a clean dict.

    Keys: 'Name' and 'ATCC Number' (from the Product JSON-LD), every
    ``product-information`` label/value pair (summary table + detail accordion),
    and 'Permits & Restrictions' -> list of {'title', 'text'}.
    """
    soup = BeautifulSoup(html, "html.parser")
    record = {}

    # Name + ATCC number from the Product JSON-LD (fall back to <h1> for the name).
    ld = _load_product_ld(soup)
    if ld.get("name"):
        record["Name"] = ld["name"]
    elif h1 := soup.find("h1"):
        record["Name"] = _clean(h1.get_text())
    if ld.get("sku"):
        record["ATCC Number"] = ld["sku"]

    # Drop the "Legal disclaimers" accordion (long boilerplate) before harvesting,
    # otherwise Intended use / Warranty / Disclaimers show up as huge fields.
    for item in soup.select(".generic-accordion__item"):
        title = item.select_one(".generic-accordion__item-title-text")
        if title and _clean(title.get_text()).lower() == "legal disclaimers":
            item.decompose()

    # All product-information dt/dd pairs (summary table + detail accordion sections).
    for dt in soup.select("dt.product-information__title"):
        dd = dt.find_next_sibling("dd")
        if dd is None:
            continue
        label, value = _clean(dt.get_text()), _clean(dd.get_text())
        if label and label not in record:
            record[label] = value

    # Permits & Restrictions: each permit has a heading + rich-text body.
    permits = []
    for permit in soup.select("section.product-permits .product-permit"):
        heading = permit.select_one(".product-permit__heading")
        body = permit.select_one(".product-permit__rtf")
        entry = {}
        if heading:
            entry["title"] = _clean(heading.get_text())
        if body:
            entry["text"] = _clean(body.get_text())
        if entry:
            permits.append(entry)
    if permits:
        record["Permits & Restrictions"] = permits

    return record


def get_atcc_data(product_id):
    """Fetch + parse a single ATCC product by id. Returns the dict, or None."""
    html = fetch_atcc_html(product_id)
    if html is None:
        return None
    return parse_atcc_html(html)


if __name__ == "__main__":
    print(json.dumps(get_atcc_data(51156), indent=2, ensure_ascii=False))
