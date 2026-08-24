# PGS Catalog mapping decisions

Implemented from `PGS Mapping - Sheet1.tsv` after review:

- Prefix record IDs as `pgs_catalog_<MONDO_ID>` so MONDO-grouped collections from different sources cannot merge accidentally. The source `identifier` remains the bare MONDO ID.
- Do not emit `healthCondition.sameAs`. MONDO already provides the harmonized disease identifier needed for ontology cross-references.
- Do not emit `analyticalMethod`; that property is not available on the reviewed DataCollection profile. Method details for the deterministic representative score remain in `exampleOfWork.additionalProperty`.
- Use only the three reviewed EDAM topics in `topicCategory`. PGS Catalog trait categories are retained in `keywords`.
- Use the reviewed computational-method and genetic-variation-analysis terms in `measurementTechnique`.
- Omit collection-level `license` and top-level `version`. Per-score terms remain in the representative work metadata; the release date is recorded in `includedInDataCatalog.versionDate`.

Operational scope for the first implementation:

- Emit score-bearing MONDO traits, including traits whose scores are associated only through child traits.
- Use the complete deduplicated union of direct and child-associated score IDs and do not silently truncate large collections.
- Skip MONDO traits with neither direct nor child-associated scores.
- Keep source spellings in sample metadata unchanged.

`@context` is emitted as the schema.org URL string used by the current uploader mapping. The expanded context object shown in the review sheet is valid JSON-LD but conflicts with the uploader's text mapping.
