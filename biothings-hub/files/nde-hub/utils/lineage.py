import os

import orjson
from biothings_client import get_client


def process_lineage(docs):
    mt = get_client("taxon")

    if isinstance(docs, str):
        with open(os.path.join(docs, "data.ndjson"), "rb") as f:
            records = [orjson.loads(line) for line in f]
    else:
        records = list(docs)

    # Collect all unique taxon IDs from species and infectiousAgent fields
    all_taxon_ids = set()
    for record in records:
        for field in ["species", "infectiousAgent"]:
            if field in record:
                for item in record[field]:
                    taxid = item.get("identifier")
                    if taxid and taxid.isdigit():
                        all_taxon_ids.add(int(taxid))

    # Post a batch request to get the taxon information
    taxon_info_list = mt.gettaxa(list(all_taxon_ids))

    # Lookup dictionary for taxon ID to lineage
    taxon_lineage_dict = {}
    for taxon_info in taxon_info_list:
        taxid = taxon_info.get("taxid")
        lineage = taxon_info.get("lineage", [])
        if taxid is not None:
            taxon_lineage_dict[taxid] = lineage

    # Add the lineage IDs to each record under _meta.lineage.taxon
    count = 0
    for record in records:
        count += 1
        if count % 1000 == 0:
            print(f"Processing document {count}...")

        lineage_ids = set()
        for field in ["species", "infectiousAgent"]:
            if field in record:
                for item in record[field]:
                    taxid = item.get("identifier")
                    if taxid and taxid.isdigit():
                        taxid = int(taxid)
                        lineage = taxon_lineage_dict.get(taxid, [])
                        lineage_ids.update(lineage)
        record.setdefault("_meta", {}).setdefault("lineage", {})["taxon"] = list(lineage_ids)

    return records
