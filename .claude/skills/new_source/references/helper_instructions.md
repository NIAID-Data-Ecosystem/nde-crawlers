# Helper instructions

Rules for choosing helper functions in the hub source `uploader.py`. The
authoritative decision table (with chain order and import paths) lives in
`SKILL.md` §7 — this file is the original notes those rules were derived from.

## General-purpose helpers

- **`load_pmid_ctfd`** — use when `pmids`, `pmcs`, or `citation.doi` are in the
  output. It also pulls species and health conditions, so when you use this
  you should also use `standardize_data` **and** `standardize_funding`.
- **`standardize_data`** — use when the output has `species` or
  `healthCondition`.
- **`standardize_funding`** — use when the output has `funding` in the parser.
- **`process_descriptions`** — use when the source has a `description`. Run
  **after** `standardize_data`.
- **`process_lineage`** — for the taxonomy browser. Currently reads all of
  `species` and `infectiousAgent`. Will be integrated into the overall utils.
- **`corrections`** — for the `sourceOrganization` property. Ginger will let
  you know when this applies. If your source is part of a NIAID program
  collection it will add the `sourceOrganization` to your record. Apply
  broadly.
- **`process_measurement_technique`** — use when the output has
  `measurementTechnique`. Requires `self.name`.
- **`add_topic_category`** — requires `self.name`. Possibly add to all sources.
  Zubair generated files for adding `topicCategory` to NDE; folder is
  `/data/nde_hub/topic_categories` on su11.
- **`nctid_helper`** — use when the parser emits an `nctid` value. The helper
  derives measurement techniques from the NCT trial info.

## Source-restricted helpers

- **`add_disambiguating_description`** — only for `immport` and `clinepidb`
  (they have a `disambiguating_description` field).
- **`load_ct_wrapper`** — `vivli`-specific.
- **`standardize_fields`** — `dde`-specific.
- **`handle_dde_docs`** — `dde`-specific.