import datetime
import logging
import re
import time

import dateutil.parser
import requests
from atcc_crawler import get_atcc_data
from ccug_crawler import get_ccug_data
from jcm_crawler import get_jcm_data
from kctc_crawler import get_kctc_data
from ucccb_crawler import get_ucccb_data

import bacdive

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

SPARQL_URL = "https://sparql.dsmz.de/api/bacdive"
STRAIN_QUERY = (
    "SELECT ?s WHERE { ?s a <https://purl.dsmz.de/schema/Strain> }"
)
BATCH_SIZE = 100    # v2 fetch endpoint cap
SLEEP = 0.05        # seconds between fetch calls

def insert_value(d, key, value, extend=False):
    """ Insert a value into a dictionary, handling existing keys by converting to lists or extending strings as needed.
    """

    if key in d and not extend:
        if isinstance(d[key], list):
            if isinstance(value, list):
                for item in value:
                    if item not in d[key]:
                        d[key].append(item)
            elif value not in d[key]:
                d[key].append(value)
        else:
            if isinstance(value, list):
                d[key] = [d[key]] + [v for v in value if v != d[key]]
            elif d[key] != value:
                d[key] = [d[key], value]
    elif d.get(key) and extend:
        d[key] = (d.get(key) + " " + value).strip()
    else:
        d[key] = value

def get_all_bacdive_ids():
    """Return a sorted list of every BacDive strain ID via SPARQL."""
    logger.info("Fetching all BacDive IDs via SPARQL...")
    resp = requests.get(
        SPARQL_URL,
        params={"query": STRAIN_QUERY},
        headers={"Accept": "application/sparql-results+json"},
        timeout=120,
    )
    resp.raise_for_status()
    bindings = resp.json()["results"]["bindings"]
    ids = []
    for b in bindings:
        uri = b["s"]["value"]
        ids.append(int(uri.rsplit("/", 1)[-1]))
    ids.sort()
    return ids


def iter_bacdive_records(ids=None, sleep=SLEEP):
    """Yield (bacdive_id, record_dict) for every ID (all of BacDive by default)."""
    if ids is None:
        ids = get_all_bacdive_ids()
    logger.info(f"Fetching BacDive records for {len(ids)} IDs...")
    client = bacdive.BacdiveClient()  # public=True -> no auth
    for i in range(0, len(ids), BATCH_SIZE):
        logger.info(f"Fetching IDs {min(i + BATCH_SIZE, len(ids))} out of {len(ids)}...")
        chunk = ids[i:i + BATCH_SIZE]
        if client.search(id=chunk) == 0:
            continue
        for entry in client.retrieve():
            bid = entry.get("General", {}).get("BacDive-ID")
            yield bid, entry
        time.sleep(sleep)



def _as_list(v):
    """Return a value as a list, or an empty list if None. If the value is already a list, return it unchanged."""
    if v is None:
        return []
    if not isinstance(v, list):
        return [v]
    return v



_LENGTH_RE = re.compile(
    r"^\s*(?P<min>\d+(?:\.\d+)?)(?:\s*-\s*(?P<max>\d+(?:\.\d+)?))?\s*(?P<unit>[µu]m|nm|mm)?\s*$"
)


def _parse_quantitative_length(val, name):
    """Parse 'X[-Y] unit' into a QuantitativeValue-shaped dict, else None."""
    if val is None:
        return None
    m = _LENGTH_RE.match(str(val))
    if not m:
        return None
    entry = {"name": name}
    mn = float(m.group("min"))
    if mx := m.group("max"):
        entry["minValue"] = mn
        entry["maxValue"] = float(mx)
    else:
        entry["value"] = mn
    if unit := m.group("unit"):
        entry["unitText"] = "um" if unit in ("µm", "um") else unit
    return entry


def _to_iso_date(val):
    if val is None:
        return None
    try:
        if isinstance(val, int):
            val = str(val)
        dt = dateutil.parser.parse(val, ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt


_TEMPERATURE_RE = re.compile(r"\s*(\d+(?:\.\d+)?)\s*(\S.*)?$")


def parse_temperature(value, minimum=None, maximum=None):
    """Parse temperature string(s) into a normalized QuantitativeValue dict.

    ``value`` is the main/optimal temperature; optional ``minimum``/``maximum`` add
    minValue/maxValue (UCCCB exposes the three separately). Sources differ in
    spacing and degree glyph ('30°C', '30 ℃', '28ºC'), so the number is pulled out
    with a regex and the unit is normalized to '°C'. Returns {} if nothing parses.
    """
    result = {}
    for key, raw in (("value", value), ("minValue", minimum), ("maxValue", maximum)):
        if raw and (m := _TEMPERATURE_RE.match(str(raw))):
            result[key] = float(m.group(1))
    if not result:
        return {}
    return {"@type": "QuantitativeValue", "unitText": "°C", **result}


def get_sample_state(dsm_number):
    """Return the DSMZ catalogue delivery forms for a DSM number, e.g.
    ['Freeze Dried', 'Active culture on request', 'DNA']. Returns [] if unavailable.

    The delivery forms are not exposed via the strains JSON API (its supplyForms field is
    empty); they live in a static per-strain price fragment on the DSMZ website.
    """
    url = (
        "https://www.dsmz.de/fileadmin/cataloguepriceinformation/"
        f"priceinformation_DSM_{dsm_number}.html"
    )
    try:
        resp = requests.get(url, timeout=20)
    except requests.RequestException as e:
        logger.warning(f"Could not fetch delivery forms for DSM {dsm_number}: {e}")
        return []
    if resp.status_code != 200:
        return []
    m = re.search(r"delivery forms and prices are:(.*?)Price Category", resp.text, re.S)
    block = m.group(1) if m else resp.text
    rows = re.findall(r'<td>([^<]+?)</td>\s*<td[^>]*>&nbsp;</td>\s*<td align="right">', block)
    return [n.strip() for n in rows if n.strip().lower() != "delivery form"]


def parse_dsm(ccn, output):
    """Parse DSMZ catalogue data for a culture collection number (e.g. 'DSM 1234') and insert into the output dict.
        Example:  https://www.dsmz.de/collection/catalogue/details/culture/DSM-14760
                  https://api.strains.dsmz.de/dsm/14760
                  https://api.bacdive.dsmz.de/v2/fetch/100
    """
    request_url = f"https://api.strains.dsmz.de/dsm/{ccn.split(' ')[1]}"
    try:
        resp = requests.get(request_url, timeout=20)
    except requests.RequestException as e:
        logger.warning(f"Could not fetch DSMZ strains data for {ccn}: {e}")
        return
    if resp.status_code != 200:
        return
    data = resp.json()
    if not _as_list(data.get("collections")):
        return
    included_in_data_catalog = {
        "@type": "DataCatalog",
        "name": _as_list(data.get("collections"))[0].get("legalName"),
        "url": "https://www.dsmz.de/",
        "versionDate": datetime.date.today().isoformat(),
        "archivedAt": f"https://www.dsmz.de/collection/catalogue/details/culture/{'-'.join(ccn.split(' '))}",
    }
    insert_value(output, "includedInDataCatalog", included_in_data_catalog)
    item_location = {
        "@type": "AdministrativeArea",
        "name": "Germany",
        "administrativeType": "Country",
    }
    insert_value(output, "itemLocation", item_location)
    if sample_state := get_sample_state(ccn.split(" ")[1]):
        insert_value(output, "sampleState", sample_state)
        for state in _as_list(output["sampleState"]):
            if "on request" in state.lower():
                insert_value(output, "creativeWorkStatus", ["Available", "Bespoke"])
                break

    sample_storage = {}

    if growth_temp := _as_list(data.get("growthConditions")):
        if growth_temp[0].get("optimalTemperature") and growth_temp[0].get("minimalTemperature") and growth_temp[0].get("maximalTemperature"):
            sample_storage["@type"] = "QuantitativeValue"
            sample_storage["value"] = growth_temp[0].get("optimalTemperature")
            sample_storage["minValue"] = growth_temp[0].get("minimalTemperature")
            sample_storage["maxValue"] = growth_temp[0].get("maximalTemperature")
            sample_storage["unitText"] = "°C"
        elif test_temp := _as_list(growth_temp[0].get("testsTemperature")):
            sample_storage["@type"] = "QuantitativeValue"
            sample_storage["minValue"] = test_temp[0].get("minimal")
            sample_storage["maxValue"] = test_temp[0].get("maximal")
            sample_storage["unitText"] = "°C"

    if sample_storage:
        insert_value(output, "sampleStorageTemperature", sample_storage)

    insert_value(output, "usageInfo", {"description": "Nagoya Protocol Restrictions"})

    if date := _as_list(data.get("origin")):
        if date[0].get("sampleDate"):
            insert_value(output, "dateCollected", _to_iso_date(date[0].get("sampleDate")))

    output["sampleAvailability"] = True

def parse_kctc(ccn, output):
    cnn = ccn.split(" ")[1]
    data = get_kctc_data(cnn)
    if not data:
        logger.warning(f"Could not fetch KCTC data for {ccn}")
        return
    included_in_data_catalog = {
        "@type": "DataCatalog",
        "name": "Korean Collection for Type Cultures (KCTC)",
        "url": "https://kctc.kribb.re.kr/",
        "versionDate": datetime.date.today().isoformat(),
        "archivedAt": f"https://kctc.kribb.re.kr/collections/view?sn={cnn}",
    }
    insert_value(output, "includedInDataCatalog", included_in_data_catalog)
    item_location = {
        "@type": "AdministrativeArea",
        "name": "South Korea",
        "administrativeType": "Country",
    }
    insert_value(output, "itemLocation", item_location)

    insert_value(output, "sampleState", data.get("Price"))

    if sample_storage := parse_temperature(data.get("Temperature")):
        insert_value(output, "sampleStorageTemperature", sample_storage)

    output["sampleAvailability"] = True

    if mta := data.get("MTA Restrictions"):
        insert_value(output, "usageInfo", {"description": mta})
    if data.get("Price"):
        insert_value(output, "creativeWorkStatus", ["Available"])

def parse_jcm(ccn, output):
    cnn = ccn.split(" ")[1]
    data = get_jcm_data(cnn)
    if not data:
        logger.warning(f"Could not fetch JCM data for {ccn}")
        return
    included_in_data_catalog = {
        "@type": "DataCatalog",
        "name": "Japan Collection of Microorganisms (JCM)",
        "url": "https://www.jcm.riken.jp/",
        "versionDate": datetime.date.today().isoformat(),
        "archivedAt": f"https://www.jcm.riken.jp/cgi-bin/jcm/jcm_number?JCM={cnn}",
    }
    insert_value(output, "includedInDataCatalog", included_in_data_catalog)

    item_location = {
        "@type": "AdministrativeArea",
        "name": "Japan",
        "administrativeType": "Country",
    }
    insert_value(output, "itemLocation", item_location)

    if sample_storage := parse_temperature(data.get("Temperature")):
        insert_value(output, "sampleStorageTemperature", sample_storage)

    delivery = []
    if d := data.get("Domestic"):
        delivery.append(f"Domestic: {d}")
    if o := data.get("Overseas"):
        delivery.append(f"Overseas: {o}")
    if delivery:
        insert_value(output, "sampleState", " ".join(delivery))

    insert_value(output, "sampleAvailability", True)

    if data.get("Domestic") or data.get("Overseas"):
        insert_value(output, "creativeWorkStatus", "Available")
    if "culture on request" in data.get("Domestic", "").lower() or "culture on request" in data.get("Overseas", "").lower():
        insert_value(output, "creativeWorkStatus", "Bespoke")

    if spatial_coverage := data.get("Source"):
        insert_value(output, "spatialCoverage", {"name": spatial_coverage})

def parse_ccug(ccn, output):
    ccn = ccn.split(" ")[1]
    data = get_ccug_data(ccn)
    if not data:
        logger.warning(f"Could not fetch CCUG data for {ccn}")
        return
    included_in_data_catalog = {
        "@type": "DataCatalog",
        "name": "Culture Collection University of Gothenburg (CCUG)",
        "url": "https://www.ccug.se/",
        "versionDate": datetime.date.today().isoformat(),
        "archivedAt": f"https://www.ccug.se/strain?id={ccn}",
    }
    insert_value(output, "includedInDataCatalog", included_in_data_catalog)

    item_location = {
        "@type": "AdministrativeArea",
        "name": "Sweden",
        "administrativeType": "Country",
    }
    insert_value(output, "itemLocation", item_location)

    insert_value(output, "sampleState", "freeze-dried (lyophilized) biomass in sealed glass ampoules")

    sample_storage = {
        "@type": "QuantitativeValue",
        "minValue": 4,
        "maxValue": 8,
        "unitText": "°C",
    }
    insert_value(output, "sampleStorageTemperature", sample_storage)
    insert_value(output, "sampleAvailability", True)
    insert_value(output, "creativeWorkStatus", "Available")
    if spatial_coverage := data.get("Sample Origin"):
        insert_value(output, "spatialCoverage", {"name": spatial_coverage})


def parse_atcc(ccn, output):
    ccn = ccn.split(" ")[1]
    data = get_atcc_data(ccn)
    if not data:
        logger.warning(f"Could not fetch ATCC data for {ccn}")
        return
    included_in_data_catalog = {
        "@type": "DataCatalog",
        "name": "American Type Culture Collection (ATCC)",
        "url": "https://www.atcc.org/",
        "versionDate": datetime.date.today().isoformat(),
        "archivedAt": f"https://www.atcc.org/products/{ccn}",
    }
    insert_value(output, "includedInDataCatalog", included_in_data_catalog)

    item_location = {
        "@type": "AdministrativeArea",
        "name": "Worldwide",
    }

    insert_value(output, "itemLocation", item_location)
    if sample_state := data.get("Product Format"):
        insert_value(output, "sampleState", sample_state)

    if sample_storage := parse_temperature(data.get("Temperature")):
        insert_value(output, "sampleStorageTemperature", sample_storage)

    insert_value(output, "sampleAvailability", True)

    # Permits & Restrictions (list of {title, text}) -> single usageInfo description
    permits = data.get("Permits & Restrictions") or []
    if usage := " ".join(
        ": ".join(part for part in (p.get("title"), p.get("text")) if part)
        for p in permits
    ):
        insert_value(output, "usageInfo", {"description": usage})

    insert_value(output, "creativeWorkStatus", "Available")

    if spatial_coverage := data.get("Geographical isolation"):
        insert_value(output, "spatialCoverage", {"name": spatial_coverage})

def parse_ucccb(ccn, output):
    data = get_ucccb_data(ccn)
    if not data:
        logger.warning(f"Could not fetch UCCCB data for {ccn}")
        return
    included_in_data_catalog = {
        "@type": "DataCatalog",
        "name": "University of Coimbra Bacteria Culture Collection (UCCCB)",
        "url": "https://ucccb.uc.pt/",
        "versionDate": datetime.date.today().isoformat(),
        "archivedAt": f"https://ucccb.uc.pt/strain-details/?detail={ccn}",
    }
    insert_value(output, "includedInDataCatalog", included_in_data_catalog)

    item_location = {
        "@type": "AdministrativeArea",
        "name": "Portugal",
        "administrativeType": "Country",
    }
    insert_value(output, "itemLocation", item_location)

    if sample_storage := parse_temperature(
        data.get("Optimal Temperature(s)"),
        data.get("Minimum Temperature"),
        data.get("Maximum Temperature"),
    ):
        if name := data.get("Preservation Procedures Used"):
            sample_storage["name"] = name
        insert_value(output, "sampleStorageTemperature", sample_storage)

    insert_value(output, "sampleAvailability", True)
    insert_value(output, "creativeWorkStatus", "Available")
    if spatial_coverage := data.get("Country"):
        insert_value(output, "spatialCoverage", {"name": spatial_coverage})
    if date := data.get("Date of isolation"):
        insert_value(output, "dateCollected", _to_iso_date(date))


def parse():
    for bid, record in iter_bacdive_records():
        general = record.get("General") or {}
        natc = record.get("Name and taxonomic classification") or {}
        isolation_root = record.get("Isolation, sampling and environmental information") or {}
        sequence_info = record.get("Sequence information") or {}
        literature = record.get("Literature") or {}
        references = record.get("Reference") or []
        morphology = record.get("Morphology") or {}

        _id = general.get("BacDive-ID")
        if not _id:
            continue
        url = f"https://bacdive.dsmz.de/strain/{_id}"

        output = {
            "@context": "http://schema.org/",
            "@type": "Sample",
            "_id": f"bacdive_{_id}",
            "identifier": str(_id),
            "url": url,
            "distribution": [{"@type": "DataDownload", "contentUrl": url}],
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "BacDive",
                "url": "https://bacdive.dsmz.de/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
            "conditionsOfAccess": "Open",
            "license": "https://creativecommons.org/licenses/by/4.0/",
            "isAccessibleForFree": False
        }

        if desc := general.get("description"):
            insert_value(output, "description", desc, extend=True)

        if kw := general.get("keywords"):
            insert_value(output, "keywords", kw)

        # NCBI tax id → infectiousAgent (strain preferred over species)
        for level in ("strain", "species"):
            match = next(
                (e for e in _as_list(general.get("NCBI tax id"))
                 if e.get("Matching level") == level),
                None,
            )
            if match and (tax_id := match.get("NCBI tax id")) is not None:
                ia = {"identifier": str(tax_id)}
                if level == "strain":
                    ia["@type"] = "DefinedTerm"
                insert_value(output, "infectiousAgent", ia)
                break

        if doi := general.get("doi"):
            insert_value(output, "doi", doi)

        for syn in _as_list((natc.get("LPSN") or {}).get("synonyms")):
            if name := syn.get("synonym"):
                insert_value(output, "infectiousAgent", {"name": name})

        # Culture-collection numbers → alternateIdentifier; DSM rows → sameAs URL
        ccns = []
        if ccn_str := literature.get("culture collection no."):
            ccns.extend(c.strip() for c in str(ccn_str).split(",") if c.strip())
        if dsm := general.get("DSM-Number"):
            if f"DSM {dsm}" not in ccns:
                ccns.append(f"DSM {dsm}")
        for ccn in dict.fromkeys(ccns):
            insert_value(output, "alternateIdentifier", ccn)
            if ccn.upper().startswith("DSM "):
                num = ccn.split(" ", 1)[1].strip()
                insert_value(
                    output,
                    "sameAs",
                    f"https://www.dsmz.de/collection/catalogue/details/culture/DSM-{num}",
                )

            if ccn.upper().startswith("DSM"):
                parse_dsm(ccn, output)
            if ccn.upper().startswith("KCTC"):
                parse_kctc(ccn, output)
            if ccn.upper().startswith("JCM"):
                parse_jcm(ccn, output)
            if ccn.upper().startswith("CCUG"):
                parse_ccug(ccn, output)
            if ccn.upper().startswith("ATCC"):
                parse_atcc(ccn, output)
            if ccn.upper().startswith("UCCCB"):
                parse_ucccb(ccn, output)

        # Isolation entries (single dict or list)
        for iso in _as_list(isolation_root.get("isolation")):
            if st := iso.get("sample type"):
                insert_value(output, "sampleType", {"name": st})

            loc = {}
            if country := iso.get("country"):
                loc["name"] = country
            if cc := iso.get("origin.country"):
                loc["identifier"] = cc
            geo = {}
            for key in ("latitude", "longitude"):
                raw = iso.get(key)
                if raw is None:
                    continue
                try:
                    geo[key] = float(raw)
                except (TypeError, ValueError):
                    pass
            if geo:
                loc["geo"] = geo
            if loc:
                insert_value(output, "locationOfOrigin", loc)

            if continent := iso.get("continent"):
                insert_value(output, "spatialCoverage", {"name": continent})
            if host := iso.get("host species"):
                insert_value(output, "species", {"name": host})
            if d := _to_iso_date(iso.get("isolation date")):
                output.setdefault("dateProcessed", d)
            if d := _to_iso_date(iso.get("sampling date")):
                output.setdefault("dateCollected", d)
            for key in ("enrichment culture composition", "isolation procedure"):
                if val := iso.get(key):
                    insert_value(output, "collectionMethod", val)

        # Isolation source categories → environmentalSystem
        for cat in _as_list(isolation_root.get("isolation source categories")):
            for _, v in cat.items():
                if v:
                    insert_value(output, "environmentalSystem", {"name": v})

        # Sequence information → isBasisFor (Genome + 16S accessions)
        for key in ("Genome sequences", "16S sequences"):
            for s in _as_list(sequence_info.get(key)):
                entry = {}
                if name := s.get("description"):
                    entry["name"] = name
                ident = s.get("INSDC accession") or s.get("accession")
                if ident:
                    entry["identifier"] = str(ident)
                if entry:
                    insert_value(output, "isBasisFor", entry)

        # Literature.literature → citedBy
        for lit in _as_list(literature.get("literature")):
            entry = {}
            if pmid := lit.get("Pubmed-ID"):
                insert_value(output, "pmids", str(pmid))
                entry["pmid"] = str(pmid)
            if d := lit.get("DOI"):
                insert_value(output, "citation", {"doi": d})
                entry["doi"] = d
            if title := lit.get("title"):
                entry["name"] = title
            if journal := lit.get("journal"):
                entry["journalName"] = journal
            if year := _to_iso_date(lit.get("year")):
                entry["datePublished"] = year
            if authors := lit.get("authors"):
                authors = [a.strip() for a in str(authors).split(",") if a.strip()]
                if authors:
                    entry["author"] = [{"name": a} for a in authors]
            if entry:
                insert_value(output, "citedBy", entry)

        if output.get("pmids"):
            if isinstance(output["pmids"], list):
                output["pmids"] = ", ".join(output["pmids"])

        # Reference → isBasedOn
        for ref in _as_list(references):
            entry = {}
            if authors := ref.get("authors"):
                entry["author"] = [{"name": authors}]
            if name := ref.get("title") or ref.get("catalogue"):
                entry["name"] = name
            if u := ref.get("doi/url"):
                entry["url"] = f"https://doi.org/{u}" if str(u).startswith("10.") else u
            if entry:
                insert_value(output, "isBasedOn", entry)

        # Morphology → associatedPhenotype (DefinedTerm-style labels only)
        for cell in _as_list(morphology.get("cell morphology")):
            gs = cell.get("gram stain")
            if gs in ("positive", "negative"):
                insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": f"gram {gs}"})
            if shape := cell.get("cell shape"):
                insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": shape})
            mot = cell.get("motility")
            if mot == "yes":
                insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": "motile"})
            elif mot == "no":
                insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": "non-motile"})
            if length := _parse_quantitative_length(cell.get("cell length"), "cell length"):
                length["@type"] = "QuantitativeValue"
                insert_value(output, "associatedPhenotype", length)
        for mc in _as_list(morphology.get("multicellular morphology")):
            if name := mc.get("complex name"):
                insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": name})
            if color := mc.get("complex color"):
                insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": color})

        if isinstance(output.get("includedInDataCatalog"), list):
            output["additionalType"] = "ExperimentalRunSample"
        else:
            output["additionalType"] = "BioSample"
        yield output
