import re
import time

import bacdive
import requests
import datetime
import dateutil.parser
import logging

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
    client = bacdive.BacdiveClient()  # public=True -> no auth
    for i in range(0, len(ids), BATCH_SIZE):
        print(f"Fetching IDs {min(i + BATCH_SIZE, len(ids))} out of {len(ids)}...")
        chunk = ids[i:i + BATCH_SIZE]
        if client.search(id=chunk) == 0:
            continue
        for entry in client.retrieve():
            bid = entry.get("General", {}).get("BacDive-ID")
            yield bid, entry
        time.sleep(sleep)



def _as_list(v) -> list[dict]:
    """Coerce BacDive's union-typed fields (dict | list[dict] | None) to list[dict]."""
    if isinstance(v, dict):
        return [v]
    if isinstance(v, list):
        return [x for x in v if isinstance(x, dict)]
    return []


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
        dt = dateutil.parser.parse(val, ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt


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

        # LPSN synonyms → infectiousAgent.alternateName
        for syn in _as_list((natc.get("LPSN") or {}).get("synonyms")):
            if name := syn.get("synonym"):
                insert_value(output, "infectiousAgent", {"alternateName": name})

        # Culture-collection numbers → alternateIdentifier; DSM rows → sameAs URL
        ccns = []
        if ccn_str := literature.get("culture collection no."):
            ccns.extend(c.strip() for c in str(ccn_str).split(",") if c.strip())
        if dsm := general.get("DSM-Number"):
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
                # insert_value(output, "pmids", str(pmid))
                entry["pmid"] = str(pmid)
            if d := lit.get("DOI"):
                # insert_value(output, "citation", d)
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

        # if output.get("pmids"):
        #     if isinstance(output["pmids"], list):
        #         output["pmids"] = ", ".join(output["pmids"])

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

        yield output
