---
name: new_source
description: Scaffold a new NDE crawler source end-to-end — the crawler container directory, the biothings-hub dumper/uploader package, and the docker-compose-crawlers.yml entry — driven by a mapping TSV/CSV (plus optional examples and heuristics).
---

# /new_source — scaffold a new NDE crawler

Invocation: `/new_source [--type Dataset|Sample]`

Default `--type` is `Dataset`. `Sample` switches the uploader base class to `NDESourceSampleUploader` and the cross-check target mapping to `NDESourceSampleUploader.get_mapping()` in [biothings-hub/files/nde-hub/hub/dataload/nde.py](biothings-hub/files/nde-hub/hub/dataload/nde.py).

## Inputs to collect from the user

Before writing any files, ask the user (with `AskUserQuestion` if not already supplied):

1. **Source name** — short, lowercase, snake_case (e.g. `bacdive`). This becomes:
   - the crawler directory name (`<name>/`)
   - the hub source package name (`biothings-hub/files/nde-hub/hub/dataload/sources/<name>/`)
   - `SRC_NAME` on the dumper and `name` on the uploader
   - the `_id` prefix (`<name>_<identifier>`)
   - the docker-compose service key (`<name>-crawler`)
2. **Input directory** — a directory containing at least:
   - **Required:** one `*.tsv` or `*.csv` mapping file (each row maps a source field path → target schema.org field, with optional notes). Example format: see [references/bacdive_example/mapping.tsv](.claude/skills/new_source/references/bacdive_example/mapping.tsv).
   - **Optional:** one or more example records (`.json`, `.ndjson`, `.csv`, `.tsv`) showing the raw input shape.
   - **Optional:** a heuristics `.csv` or `.tsv` (extra rules layered on top of the mapping). See [references/bacdive_example/heuristics.tsv](.claude/skills/new_source/references/bacdive_example/heuristics.tsv) for a representative example. When heuristics conflict with the mapping or its notes column, the heuristics file wins.
3. **How to iterate raw records** — confirm there is an existing iterator (e.g. `iter_<name>_records()`) or that the user has already provided one in a reference parser (like [references/bacdive_example/bacdive.py](.claude/skills/new_source/references/bacdive_example/bacdive.py)). The parser you generate should call it; do not invent a new fetch strategy.
4. **`--type`** — if not on the command line, default to `Dataset`. Confirm if ambiguous.

If any required input is missing, ask before proceeding. Do not invent mapping rules.

## What to create

For source name `<name>` and chosen `<Type>` (Dataset or Sample):

### 1. Crawler container directory `<name>/`

Model strictly on [bacdive/](bacdive/). Layout:

```
<name>/
├── Dockerfile
└── files/
    ├── <name>_crawler.py        # the parser
    ├── ndjson.py                # writes /data/<name>_crawled/{data.ndjson,release.txt}
    ├── requirements.txt
    ├── docker-entrypoint.sh
    ├── prod-crontab
    └── run-api-crawler.sh
```

**Naming exception:** when the source name collides with an imported Python module (as `bacdive` does with the `bacdive` PyPI package), the parser file is suffixed `_crawler` (e.g. `bacdive_crawler.py`). Otherwise the parser file is just `<name>.py`. Match whatever filename `ndjson.py` will then `import`.

Copy [bacdive/Dockerfile](bacdive/Dockerfile), [bacdive/files/docker-entrypoint.sh](bacdive/files/docker-entrypoint.sh), [bacdive/files/prod-crontab](bacdive/files/prod-crontab), [bacdive/files/run-api-crawler.sh](bacdive/files/run-api-crawler.sh), and [bacdive/files/ndjson.py](bacdive/files/ndjson.py) and adjust the source name and the parser import in `ndjson.py`. `requirements.txt` should list the third-party packages the parser actually imports.

### 2. The parser file (`<name>_crawler.py` or `<name>.py`)

Required structure, in this order:

1. **Imports** — `logging`, `datetime`, `dateutil.parser`, plus whatever the source needs.
2. **Logger** — `logger = logging.getLogger("nde-logger")`.
3. **The two boilerplate helpers, copied verbatim** from [references/base_func.py](.claude/skills/new_source/references/base_func.py):
   - `insert_value(d, key, value, extend=False)`
   - `_to_iso_date(val)`
   These must be in *every* parser. Do not modify the function bodies; do not rename them.
4. **Record iterator** — the existing `iter_<name>_records()` (or equivalent) the user pointed to. If they gave it in the input dir, copy it through. Do not invent one.
5. **Helper parsers** — small private helpers (e.g. `_as_list`, `_parse_quantitative_length`) only if the mapping needs them. Keep them minimal.
6. **`parse()` generator** — yields one dict per record. Build the boilerplate output first, then walk the mapping row-by-row.

**Boilerplate output (always present, in this shape):**

```python
output = {
    "@context": "http://schema.org/",
    "@type": "<Type>",          # "Dataset" or "Sample" per --type
    "_id": f"<name>_{_id}",
    "identifier": str(_id),
    "url": url,
    "distribution": [{"@type": "DataDownload", "contentUrl": url}],
    "includedInDataCatalog": {
        "@type": "DataCatalog",
        "name": "<HumanReadableName>",
        "url": "<catalog_root_url>",
        "versionDate": datetime.date.today().isoformat(),
        "archivedAt": url,
    },
}
```

Use [references/bacdive_example/bacdive.py](.claude/skills/new_source/references/bacdive_example/bacdive.py) and [bacdive/files/bacdive_crawler.py](bacdive/files/bacdive_crawler.py) as the canonical shape.

### 3. Driving the parser from the mapping

For each row in the mapping file (skipping rows where the `Mapping` column is empty or marked "Ignore"):

- The first column is a dotted source-field path (e.g. `General.description`, `Isolation, sampling and environmental information.isolation.country`). Walk the raw record by splitting on `.`; treat any segment ending in `[i]` as a list index (or coerce union-typed fields with a small `_as_list` helper, as bacdive does).
- The `Mapping` column is the target schema.org field. Dotted targets (e.g. `locationOfOrigin.name`) mean build a sub-object and pass it as the value.
- The `Notes` column may contain rules like "convert to URL", "Format as an NCBI Taxon DefinedTerm Object", or "set X = 'country'". Apply them. If a note conflicts with the heuristics file, the heuristics file wins.
- Always insert via `insert_value(output, target, value)` — never `output[target] = ...` directly, except for the boilerplate dict literal at the top. This is what gives correct list-merging behavior for repeated fields.
- **Description:** join to a single string and call `insert_value(output, "description", desc, extend=True)`. Never let `description` become a list. If the source has multiple description-like fields, concatenate with `extend=True` in the order they appear.
- **Dates:** route every date-like value through `_to_iso_date(...)` before inserting.
- **`infectiousAgent` minimum fields:** every emitted `infectiousAgent` object must carry at least one of `name` or `identifier`. If the mapping for a given source row sets only secondary keys (e.g. only `infectiousAgent.alternateName`), reroute the value to `infectiousAgent.name` so the object is valid. If both `name` and `identifier` would be empty after processing a record, omit the `infectiousAgent` entry entirely rather than emit a name-less stub.

### 4. Enum fields — hard rules

The upload-time validator in [biothings-hub/files/nde-hub/utils/utils.py](biothings-hub/files/nde-hub/utils/utils.py) (`check_schema`) **rejects records** with out-of-enum values for:

- `conditionsOfAccess` ∈ `{"Open", "Restricted", "Closed", "Embargoed"}`
- `creativeWorkStatus` (only when `@type == "Sample"`) ∈ `{"Bespoke", "Available", "Backordered", "Retired"}`

If the mapping wants either of these fields, map source values to the allowed enum. If you cannot map confidently, omit the field rather than emit an invalid value.

Also: `version` must not be set (`check_schema` asserts `doc.get("version") is None`).

### 5. pmids / pmcs — exception

`pmids` and `pmcs` are **not** in [nde.py](biothings-hub/files/nde-hub/hub/dataload/nde.py)'s mapping — the upload-time pmid helper handles them. If the source provides PubMed / PMC identifiers and the mapping requests `pmids` or `pmcs`:

- Emit a single **comma-separated string** (e.g. `"12345678, 23456789"`), not a list, not a list of dicts.
- Do not validate these against `get_mapping()`.

If a citation-shaped object is more appropriate (title, authors, journal), emit `citedBy` / `citation` per the schema instead.

### 6. Hub source package `biothings-hub/files/nde-hub/hub/dataload/sources/<name>/`

Three files:

- `__init__.py`:
  ```python
  from .dumper import <Name>Dumper  # noqa
  from .uploader import <Name>Uploader  # noqa
  ```
- `dumper.py` — model on [biothings-hub/files/nde-hub/hub/dataload/sources/bacdive/dumper.py](biothings-hub/files/nde-hub/hub/dataload/sources/bacdive/dumper.py). Must set `SRC_NAME = "<name>"` and a `SCHEDULE`. The `SRC_URLS` docker URI must reference `nde-crawlers-<name>-crawler`.
- `uploader.py`:
  - Subclass `NDESourceUploader` if `--type=Dataset`, `NDESourceSampleUploader` if `--type=Sample`. Import from `hub.dataload.nde`.
  - Set `name = "<name>"`.
  - **Always** decorate `load_data` with `@nde_upload_wrapper` *if* you override it.
  - **Only override `load_data`** if at least one helper from §7 is needed. If no helper applies, do not write a `load_data` method — let the base class default handle reading `data.ndjson`. The class body should then be just `name = "<name>"`.

### 7. Choosing helper functions for `load_data`

Read [references/helper_instructions.md](.claude/skills/new_source/references/helper_instructions.md) and choose helpers by these rules:

| Condition in the parser output | Helper to import & call | Notes |
|---|---|---|
| Output has `species` or `healthCondition` | `from utils.pubtator import standardize_data` | Pass `data_folder` or upstream `docs`. |
| Output has `pmids`, `pmcs`, *or* `citation.doi` | `from utils.pmid_helper import load_pmid_ctfd` | Call **before** `standardize_data`. When `load_pmid_ctfd` is used, **also include `standardize_funding` and `standardize_data`** in the chain (order per the chain block below), even if the parser output doesn't otherwise trigger those rules. |
| Output has `funding` | `from utils.funding_helper import standardize_funding` | First call in the chain — consumes `data_folder`. |
| Output has a `description` | `from utils.extract import process_descriptions` | Call **after** `standardize_data`. |
| Output has `measurementTechnique` | `from utils.measurement_technique_helper import process_measurement_technique` | Pass `self.name`. |
| Topic categories needed (Zubair's TSVs in `/data/nde_hub/topic_categories/`) | `from utils.topic_category_helper import add_topic_category` | Pass `self.name`. |
| Output has `nctid` | `from utils.nctid_helper import nctid_helper` | Derives measurement techniques from NCT trial info. |
| `disambiguating_description` field is in scope (immport / clinepidb) | `from utils.disambiguating_description import add_disambiguating_description` | Source-restricted — use only for those two. |
| Source is `vivli` | `from utils.clinical_trails_helper import load_ct_wrapper` | Source-restricted. |
| Source is `dde` | `standardize_fields` / `handle_dde_docs` | Source-restricted. |

**Chain order** when multiple apply (matches existing sources):

```
standardize_funding(data_folder)
  → load_pmid_ctfd(...)
  → standardize_data(...)
  → process_descriptions(...)
  → process_measurement_technique(..., self.name)
  → nctid_helper(...)
  → add_topic_category(..., self.name)
  → add_disambiguating_description(..., self.name)
```

Only one of the entries in the chain takes `data_folder` (the first one); the rest take the previous generator. If none of the above apply, do not write `load_data`.

Always include `from utils.utils import nde_upload_wrapper` and `@nde_upload_wrapper` on `load_data` when overridden.

### 8. docker-compose-crawlers.yml

Append a new service block, alphabetically placed if reasonable, matching the existing pattern:

```yaml
  <name>-crawler:
    build:
      context: <name>
      dockerfile: Dockerfile
    volumes:
      - data:/data
```

Use [docker-compose-crawlers.yml](docker-compose-crawlers.yml) lines around the existing bacdive entry as the exact template.

## Cross-check pass (mandatory, before reporting done)

After writing the parser, do one verification pass:

1. Open [biothings-hub/files/nde-hub/hub/dataload/nde.py](biothings-hub/files/nde-hub/hub/dataload/nde.py) and locate `NDESourceUploader.get_mapping` (line 130) if `--type=Dataset`, or `NDESourceSampleUploader.get_mapping` (line 1300) if `--type=Sample`.
2. For every top-level field your `parse()` emits, confirm it exists in that mapping. Exceptions that are allowed even though absent: `_id`, `@context`, `@type`, `pmids`, `pmcs`, and `_meta` (added later by `add_metadata_score`).
3. For every emitted field, confirm sub-keys you set match the mapping's sub-`properties`. If you've emitted `locationOfOrigin.administrativeType` but the mapping has no `administrativeType` property, either drop it or pick the correctly-named sub-key.
4. Confirm enum constraints in §4 hold.
5. Confirm `description` is a single string, not a list.
6. Confirm `pmids` / `pmcs` (if present) are single comma-separated strings.
7. If anything fails, fix the parser and re-check. Do not silently drop fields the user mapped — surface the mismatch to the user.

## Demonstration

[references/bacdive_example/](.claude/skills/new_source/references/bacdive_example/) (bundled with this skill) contains a worked example:
- `mapping.tsv` — the mapping
- `heuristics.tsv` — rules layered on top of the mapping (heuristics wins on conflict, per §"Inputs to collect")
- `example1.json` — a trimmed raw input record (one strain, one entry per section the parser touches)
- `bacdive.py` — a reference parser (matches the boilerplate-output shape this skill must emit)

The corresponding skill outputs already exist in the repo for inspection:
- [bacdive/](bacdive/) (crawler container directory)
- [biothings-hub/files/nde-hub/hub/dataload/sources/bacdive/](biothings-hub/files/nde-hub/hub/dataload/sources/bacdive/) (hub source)
- The `bacdive-crawler` entry in [docker-compose-crawlers.yml](docker-compose-crawlers.yml)

When asked to demonstrate, run the skill against the bundled `references/bacdive_example/` inputs and produce output that matches the structure of those existing files.

## Report at end

When done, post a short summary listing every file created or modified, and call out any mapping rows you could not confidently convert (so the user can review).