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

    # Fetch taxon information for all taxon IDs
    taxon_info_list = mt.gettaxa(list(all_taxon_ids))

    taxon_lineage_dict = {}
    taxon_parent_dict = {}

    # Collect taxon IDs from the lineages
    lineage_taxon_ids = set()
    for taxon_info in taxon_info_list:
        taxid = taxon_info.get("taxid")
        lineage = taxon_info.get("lineage", [])
        parent_taxid = taxon_info.get("parent_taxid")

        if taxid is not None:
            taxon_lineage_dict[taxid] = lineage
            if parent_taxid:
                taxon_parent_dict[taxid] = parent_taxid
            else:
                # Handle root taxa with no parent
                taxon_parent_dict[taxid] = None
            lineage_taxon_ids.update(lineage)

    # Fetch taxon information for lineage taxon IDs not already fetched
    missing_taxon_ids = lineage_taxon_ids.difference(all_taxon_ids)
    if missing_taxon_ids:
        lineage_taxon_info_list = mt.gettaxa(list(missing_taxon_ids))
        for taxon_info in lineage_taxon_info_list:
            taxid = taxon_info.get("taxid")
            parent_taxid = taxon_info.get("parent_taxid")
            if taxid is not None:
                taxon_parent_dict[taxid] = parent_taxid if parent_taxid else None

    def get_lineage_entries(taxid):
        lineage_ids = taxon_lineage_dict.get(taxid, []) + [taxid]
        entries = []
        for taxon in lineage_ids:
            parent_taxon = taxon_parent_dict.get(taxon)
            entry = {"taxon": taxon}
            if parent_taxon is not None:
                entry["parent_taxon"] = parent_taxon
            entries.append(entry)
        return entries

    count = 0
    for record in records:
        count += 1
        if count % 1000 == 0:
            print(f"Processing document {count}...")

        lineage_entries_set = set()
        for field in ["species", "infectiousAgent"]:
            if field in record:
                for item in record[field]:
                    taxid = item.get("identifier")
                    if taxid and taxid.isdigit():
                        taxid = int(taxid)
                        entries = get_lineage_entries(taxid)
                        for entry in entries:
                            taxon = entry["taxon"]
                            parent_taxon = entry.get("parent_taxon")
                            lineage_entries_set.add((taxon, parent_taxon))

        lineage_entries = []
        for taxon, parent_taxon in lineage_entries_set:
            entry = {"taxon": taxon}
            if parent_taxon is not None:
                entry["parent_taxon"] = parent_taxon
            lineage_entries.append(entry)

        record.setdefault("_meta", {})["lineage"] = lineage_entries

    return records
