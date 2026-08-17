import datetime
import logging
import re
import time

import bacdive
import dateutil.parser
import requests

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

def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(str(val), ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError, OverflowError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt

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


def _as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


_LENGTH_RE = re.compile(
    r"^\s*(?P<min>\d+(?:\.\d+)?)(?:\s*-\s*(?P<max>\d+(?:\.\d+)?))?\s*(?P<unit>[µu]m|nm|mm)?\s*$"
)


def _parse_quantitative_length(value, name):
    match = _LENGTH_RE.match(str(value or ""))
    if not match:
        return None
    entry = {"@type": "QuantitativeValue", "name": name}
    minimum = float(match.group("min"))
    if maximum := match.group("max"):
        entry["minValue"] = minimum
        entry["maxValue"] = float(maximum)
    else:
        entry["value"] = minimum
    if unit := match.group("unit"):
        entry["unitText"] = "um" if unit in ("µm", "um") else unit
    return entry


def parse(records=None):
    """Yield mapped records; pass `(id, record)` pairs for an offline example run."""
    records = iter_bacdive_records() if records is None else records
    for _, record in records:
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
            "isAccessibleForFree": False,
        }

        if description := general.get("description"):
            insert_value(output, "description", description, extend=True)
        if keywords := general.get("keywords"):
            insert_value(output, "keywords", keywords)
        if doi := general.get("doi"):
            insert_value(output, "doi", doi)

        taxon_entries = _as_list(general.get("NCBI tax id"))
        for level in ("strain", "species"):
            match = next((entry for entry in taxon_entries if entry.get("Matching level") == level), None)
            if match and (taxon_id := match.get("NCBI tax id")) is not None:
                insert_value(
                    output,
                    "infectiousAgent",
                    {"@type": "DefinedTerm", "identifier": str(taxon_id)},
                )
                break

        for synonym in _as_list((natc.get("LPSN") or {}).get("synonyms")):
            if name := synonym.get("synonym"):
                insert_value(output, "infectiousAgent", {"@type": "DefinedTerm", "name": name})

        culture_numbers = []
        if raw_numbers := literature.get("culture collection no."):
            culture_numbers.extend(value.strip() for value in str(raw_numbers).split(",") if value.strip())
        if dsm_number := general.get("DSM-Number"):
            dsm_identifier = f"DSM {dsm_number}"
            if dsm_identifier not in culture_numbers:
                culture_numbers.append(dsm_identifier)
        for culture_number in dict.fromkeys(culture_numbers):
            insert_value(output, "alternateIdentifier", culture_number)
            if culture_number.upper().startswith("DSM "):
                identifier = culture_number.split(" ", 1)[1].strip()
                insert_value(
                    output,
                    "sameAs",
                    f"https://www.dsmz.de/collection/catalogue/details/culture/DSM-{identifier}",
                )

        for isolation in _as_list(isolation_root.get("isolation")):
            if sample_type := isolation.get("sample type"):
                insert_value(output, "sampleType", {"@type": "DefinedTerm", "name": sample_type})

            location = {"@type": "AdministrativeArea"}
            if country := isolation.get("country"):
                location["name"] = country
                location["administrativeType"] = "country"
            if country_code := isolation.get("origin.country"):
                location["identifier"] = country_code
            geo = {"@type": "GeoCoordinates"}
            for coordinate in ("latitude", "longitude"):
                raw_value = isolation.get(coordinate)
                if raw_value is None:
                    continue
                try:
                    geo[coordinate] = float(raw_value)
                except (TypeError, ValueError):
                    logger.warning("Could not parse %s: %r", coordinate, raw_value)
            if len(geo) > 1:
                location["geo"] = geo
            if len(location) > 1:
                insert_value(output, "locationOfOrigin", location)

            if continent := isolation.get("continent"):
                insert_value(
                    output,
                    "spatialCoverage",
                    {"@type": "AdministrativeArea", "name": continent},
                )
            if host := isolation.get("host species"):
                insert_value(output, "species", {"@type": "DefinedTerm", "name": host})
            if date_processed := _to_iso_date(isolation.get("isolation date")):
                insert_value(output, "dateProcessed", date_processed)
            if date_collected := _to_iso_date(isolation.get("sampling date")):
                insert_value(output, "dateCollected", date_collected)
            for key in ("enrichment culture composition", "isolation procedure"):
                if value := isolation.get(key):
                    insert_value(output, "collectionMethod", value)

        for category in _as_list(isolation_root.get("isolation source categories")):
            for value in category.values():
                if value:
                    insert_value(
                        output,
                        "environmentalSystem",
                        {"@type": "DefinedTerm", "name": value},
                    )

        for sequence_field in ("Genome sequences", "16S sequences"):
            for sequence in _as_list(sequence_info.get(sequence_field)):
                work = {"@type": "CreativeWork"}
                if name := sequence.get("description"):
                    work["name"] = name
                identifier = sequence.get("INSDC accession") or sequence.get("accession")
                if identifier:
                    work["identifier"] = str(identifier)
                if len(work) > 1:
                    insert_value(output, "isBasisFor", work)

        for literature_entry in _as_list(literature.get("literature")):
            citation = {"@type": "ScholarlyArticle"}
            if pmid := literature_entry.get("Pubmed-ID"):
                insert_value(output, "pmids", str(pmid))
                citation["pmid"] = str(pmid)
            if doi := literature_entry.get("DOI"):
                insert_value(output, "citation", {"@type": "ScholarlyArticle", "doi": doi})
                citation["doi"] = doi
            if title := literature_entry.get("title"):
                citation["name"] = title
            if journal := literature_entry.get("journal"):
                citation["journalName"] = journal
            if published := _to_iso_date(literature_entry.get("year")):
                citation["datePublished"] = published
            if authors := literature_entry.get("authors"):
                citation["author"] = [
                    {"@type": "Person", "name": name.strip()}
                    for name in str(authors).split(",")
                    if name.strip()
                ]
            if len(citation) > 1:
                insert_value(output, "citedBy", citation)

        if isinstance(output.get("pmids"), list):
            output["pmids"] = ", ".join(output["pmids"])

        for reference in _as_list(references):
            work = {"@type": "CreativeWork"}
            if authors := reference.get("authors"):
                author_type = "Organization" if str(authors).startswith("Curators of ") else "Person"
                work["author"] = [{"@type": author_type, "name": authors}]
            if name := reference.get("title") or reference.get("catalogue"):
                work["name"] = name
            if raw_url := reference.get("doi/url"):
                work["url"] = f"https://doi.org/{raw_url}" if str(raw_url).startswith("10.") else raw_url
            if len(work) > 1:
                insert_value(output, "isBasedOn", work)

        for cell in _as_list(morphology.get("cell morphology")):
            if (gram_stain := cell.get("gram stain")) in ("positive", "negative"):
                insert_value(
                    output,
                    "associatedPhenotype",
                    {"@type": "DefinedTerm", "name": f"gram {gram_stain}"},
                )
            if shape := cell.get("cell shape"):
                insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": shape})
            if cell.get("motility") in ("yes", "no"):
                motility = "motile" if cell["motility"] == "yes" else "non-motile"
                insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": motility})
            if length := _parse_quantitative_length(cell.get("cell length"), "cell length"):
                insert_value(output, "associatedPhenotype", length)

        for multicellular in _as_list(morphology.get("multicellular morphology")):
            for key in ("complex name", "complex color"):
                if value := multicellular.get(key):
                    insert_value(output, "associatedPhenotype", {"@type": "DefinedTerm", "name": value})

        yield output
