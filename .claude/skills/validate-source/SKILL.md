---
name: validate-source
description: Validate and repair an existing NDE crawler parser and uploader against the current NDE Elasticsearch mapping, JSON-LD schema, upload pipeline, enum rules, and nested-object @type requirements. Use when checking an existing source, fixing parser output that fails check_schema, auditing a source after shared utils change, or completing new-source scaffolding.
---

# /validate-source — validate and repair an NDE source

Invocation: `/validate-source <source-name-or-path> [--type Dataset|Sample]`

Inspect the existing source, fix every issue found, and rerun the checks. Do not scaffold a new crawler or require `/new-source`.

## Discover the source

1. Resolve the crawler directory, normally `<source>/`, and the hub package at `biothings-hub/files/nde-hub/hub/dataload/sources/<source>/`.
2. Read the crawler's `files/ndjson.py` to identify the parser module and generator it actually calls. Do not assume the parser filename.
3. Confirm where the records are built. Not every source parses in `<source>/files/`:
   - the parser may live in the hub package (`sources/<source>/parser.py`); some sources have no crawler directory;
   - a crawler directory may be dead while the live parser is hub-side — check the source's build context in `docker-compose-crawlers.yml`;
   - `uploader.py` may define schema objects itself (curated term constants, `post_process`).

   Resolve what the uploader's `load_data` imports and audit every module it reaches.
4. Read the uploader and infer `Dataset` versus `Sample` from `NDESourceUploader` or `NDESourceSampleUploader`. Use `--type` only as an explicit override; report a disagreement instead of silently changing the uploader class.
5. Locate supplied fixtures, mapping files, heuristics, and source-specific tests. If no fixture exists, perform static validation and clearly report that runtime parser validation remains unverified.
6. Preserve existing user changes and source semantics. Never change `_id` format during validation unless the user explicitly requests an identifier migration.

Do not start a full crawl merely to validate a parser. Prefer a supplied example record, a bounded iterator input, or monkeypatching the record iterator and unrelated network collectors.

## Use these sources of truth

Check the current files rather than copying rules from memory:

1. `biothings-hub/files/nde-hub/utils/validate.py` for final upload-time checks.
2. `NDESourceUploader.get_mapping()` or `NDESourceSampleUploader.get_mapping()` in `biothings-hub/files/nde-hub/hub/dataload/nde.py` for indexed fields and nested properties.
3. [NDE_schema.jsonld](https://raw.githubusercontent.com/NIAID-Data-Ecosystem/nde-schemas/main/NDE_schema.jsonld) for the semantic type of an object.
4. The source's mapping and heuristics files for intended transformations.

When these disagree, satisfy the runtime validator without violating the JSON-LD schema. If a valid schema property is missing from `get_mapping()`, add the narrow mapping property needed by the selected uploader type and explain the shared mapping change. Do not silently discard a requested source field.

## Inspect every output path

Trace the parser's generator plus every helper that inserts or returns output objects. Check conditional branches, cached objects, fallbacks used after failed lookups, and list entries—not only the main output literal.

Objects reach a field by several routes. Check each:

- a literal, `output["species"] = {"name": x}`;
- appended to a list assigned later, `items.append(obj)` then `output["funding"] = items`;
- built incrementally, `obj = {}` then `obj["name"] = x`;
- returned by a helper, `"funding": _funding_list()`;
- read from a lookup table, `TERMS[key]`, `TERMS.get(key)`, `TABLE[key].copy()`;
- spread from a constant, `{**CATALOG, "archivedAt": url}` — resolve the constant;
- taken from source data, `output["author"] = item.pop("creator")` — backfill with `setdefault("@type", ...)`.

Where one variable name serves several fields, only act when every use wants the same `@type`; otherwise inspect the site by hand.

For each emitted top-level field:

- Confirm it exists in the selected `get_mapping()`. Allowed pre-pipeline exceptions are `_id`, `@context`, `@type`, `pmids`, and `pmcs`; `_meta` is added later by the pipeline.
- Confirm every nested key exists in the mapping's `properties` or is added when the schema supports it.
- Confirm the value shape agrees with the mapping: scalar, object, or list of objects. Only add `@type` where the mapping models an object (its entry has `properties`). A JSON-LD `anyOf` can resolve to a class the index stores as a keyword; `sex` is a `keyword`.
- Preserve mapping rows that cannot be implemented confidently and report them for review.

## Enforce object `@type`

Give every emitted object the schema-correct `@type` at construction time. Do not depend on a successful API enrichment or later standardizer to repair an untyped object. Validate every entry when a field may be either one object or a list.

`check_schema` accepts each listed type or any child class of it, resolving
inheritance through the DDE class API. `WebPage` and `SoftwareApplication` are
valid where `CreativeWork` is expected. Expected parent types:

| Field or nested field | Expected `@type` parent(s) |
|---|---|
| `author`, `creator`, and a citation's `author` | `Person` or `Organization` |
| An author's `affiliation` | `Organization` |
| `species`, `infectiousAgent`, `healthCondition`, `measurementTechnique`, `topicCategory` | `DefinedTerm` |
| `citation`, `citedBy`, `isBasisFor`, `isPartOf`, `hasPart` | `ScholarlyArticle`, `CreativeWork`, or `Sample` |
| `isBasedOn` | `ScholarlyArticle`, `CreativeWork`, `Sample`, `Action`, or `ResourceCatalog` |
| `funding` | `MonetaryGrant` |
| `funding.funder` | `Organization` |
| `funding.funder.employee` | `Person` |
| `funding.isBasedOn` | `ScholarlyArticle` or `CreativeWork` |
| `sourceOrganization` | `Organization` or `ResearchProject` |
| A term's `curatedBy` | `Person`, `Organization`, `DataCatalog`, `ResourceCatalog`, `SoftwareApplication`, `ComputationalTool`, or `ResearchProject` |

Also enforce common schema objects even when `check_schema` does not yet check them:

- `distribution` → `DataDownload`
- `includedInDataCatalog` → `DataCatalog` or `ResourceCatalog`
- DOI or PMID publication stub → `ScholarlyArticle`
- `usageInfo` object → `CreativeWork`
- `locationOfOrigin` and geographic `spatialCoverage` objects → normally `AdministrativeArea`; nested `geo` → `GeoCoordinates`
- `sampleType`, `environmentalSystem`, and ontology-backed phenotype objects → `DefinedTerm`
- Numeric phenotype, quantity, or temperature objects → `QuantitativeValue`

Look up any unlisted object in the current JSON-LD schema. Do not guess from the property name.
When a parser emits a child class, keep the most specific correct `@type`; do
not replace it with its parent merely to satisfy validation.
Exact expected types need no network lookup. Resolving a child type queries the
DDE, and an outage raises `ClassRegistryUnavailable`. Report that limitation
instead of rewriting a correct child type as its parent.

## Enforce the remaining schema rules

- Require non-empty `_id`, `@type`, and `url`.
- Require `includedInDataCatalog`; require every catalog object to contain `archivedAt`.
- Do not emit `version`.
- Restrict `conditionsOfAccess` to `Open`, `Restricted`, `Closed`, `Embargoed`, or `Varied`.
- For `Sample`, restrict `creativeWorkStatus` to `Bespoke`, `Available`, `Backordered`, or `Retired`.
- Keep `description` as one string, never a list.
- Keep `pmids` and `pmcs` as comma-separated strings so the citations stage can consume them.
- Require every `infectiousAgent` object to have `@type: DefinedTerm` and at least one of `name` or `identifier`; omit empty stubs.
- Route date-like source values through the parser's date-normalization helper.
- For a `Sample` that should receive description enrichment, require `additionalType: BioSample`; ordinary Samples are intentionally ineligible.

## Fix issues

Apply fixes at the earliest constructor shared by every affected output path:

- Add `@type` when an object is created, not just before one yield.
- Backfill types on cached or externally loaded dictionaries with `setdefault` when preserving an existing valid type matters.
- Type nested authors, affiliations, funders, employees, curators, geographic objects, and work objects independently.
- Correct fallback objects as well as successful API responses.
- Update the selected Elasticsearch mapping only for schema-valid subfields that the parser now emits.
- Keep shared pipeline stages centralized. Do not duplicate standardizers, finalization, or statistics logging in the uploader.

After changing a helper that is copied into bundled examples, update those examples so future sources do not reintroduce the bug.

## Verify the repaired source

Run all checks that the available fixtures and environment support:

1. Compile every modified Python file.
2. Run the parser against at least one representative record without launching an unbounded crawl.
3. Inspect the generated document recursively and confirm that every object has its schema-correct `@type`.
4. Run the actual current `check_schema` on the representative output. It reports all detected issues in one multiline `AssertionError`; fix the complete list and rerun.
5. When dependencies, caches, and safe fixtures are available, run the representative record through `run_pipeline` to validate utility-created data and finalization. Do not make production-scale external requests solely for validation.
6. Run source-specific tests and `docker compose -f docker-compose-crawlers.yml config` when crawler integration files changed.
7. Treat `ClassRegistryUnavailable` as infrastructure, not source data. Re-run once the DDE responds; do not replace a correct child type with its parent to get past it.
8. Run `git diff --check` and review the final diff for unrelated changes.

Do not claim runtime validation passed when imports, credentials, fixtures, or external services prevented it. Continue with static checks and report the exact limitation.

After a real upload, use the pipeline's `Pipeline util stats` logs to confirm per-stage before/after counts. `records_processed` is batch-level; `records_changed` reflects tracked-field changes. `post_process` has stage-level counts only.

## Report the result

Summarize:

- The source and inferred record type.
- Every file changed and each class of issue fixed.
- Validation commands run and their outcome.
- Checks skipped because of missing fixtures, dependencies, credentials, or services.
- Mapping rows or schema conflicts that still require a user decision.
- Confirmation that `_id` behavior was preserved.
