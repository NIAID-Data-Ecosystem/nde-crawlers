# IEDB DataCollection parser decisions

- Emit one `DataCollection` per non-empty IEDB results tab for each exact eligible source organism: Epitopes, Antigens, Assays, Receptors, and References.
- Use `_id = iedb_<tab>_<OrganismId>` and `identifier = IEDB:<TAB>:<OrganismId>` so sibling tab collections cannot merge.
- Stream the complete curation XML archive for Epitopes, Assays, and References. The archive is reference-oriented and expands to tens of gigabytes, so aggregation is disk-backed.
- Use the official `antigen_search`, `tcr_search`, and `bcr_search` Query API views for Antigens and Receptors because those tab entities do not have stable standalone records in `Curation.xsd`.
- Group only by primary epitope-structure source-organism roles. Assay `HostOrganism` and `RecipientOrganism` values are context, not collection keys.
- Resolve NCBI rank through `taxdump.tar.gz` and retain species-rank taxa plus descendants. IEDB-only nodes are retained when their `ParentTaxId` is at or below species.
- Do not populate `infectiousAgent` in the crawler. The host/source/infectious-agent distinction remains a post-first-pass review item.
- `collectionSize.unitText` matches the UI tab labels: `epitopes`, `antigens`, `assays`, `receptors`, and `references`.
- Deduplicate by the tab's stable key. `hasPart` is emitted only when the complete member set is at or below `IEDB_HAS_PART_LIMIT` (default 1000); larger collections omit `hasPart` entirely instead of publishing a truncated list.
- Omit `distribution`, top-level `version`, and `sample`. Preserve the complete XML export URL and tab query URLs in `isBasedOn`.
- The Query API `_iri_search` fields expand through taxonomy descendants. Public links and membership crosswalks use exact `source_organism_iri` or `source_organism_iris` fields instead.
