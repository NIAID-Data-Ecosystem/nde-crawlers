---
name: new-source
description: Scaffold a new NDE crawler source end-to-end — the crawler container directory, the biothings-hub dumper/uploader package, and the docker-compose-crawlers.yml entry — driven by a mapping TSV/CSV (plus optional examples and heuristics).
---

# /new-source — scaffold a new NDE crawler

Invocation: `/new-source [--type Dataset|Sample]`

Default `--type` is `Dataset`. `Sample` switches the uploader base class to `NDESourceSampleUploader` and the cross-check target mapping to `NDESourceSampleUploader.get_mapping()` in [biothings-hub/files/nde-hub/hub/dataload/nde.py](biothings-hub/files/nde-hub/hub/dataload/nde.py).

## Inputs to collect from the user

Before writing any files, ask the user (with `AskUserQuestion` if not already supplied):

1. **Source name** — short, lowercase, snake_case (e.g. `bacdive`). This becomes:
   - the crawler directory name (`<name>/`)
   - the hub source package name (`biothings-hub/files/nde-hub/hub/dataload/sources/<name>/`)
   - `SRC_NAME` on the dumper and `name` on the uploader
   - the `_id` prefix (`<name>_<identifier>`) **only when the identifier needs it** — see §"Deciding the `_id`: prefix vs. bare identifier"
   - the docker-compose service key (`<name>-crawler`)
2. **Input directory** — a directory containing at least:
   - **Required:** one `*.tsv` or `*.csv` mapping file (each row maps a source field path → target schema.org field, with optional notes). Example format: see [references/bacdive_example/mapping.tsv](references/bacdive_example/mapping.tsv).
   - **Optional:** one or more example records (`.json`, `.ndjson`, `.csv`, `.tsv`) showing the raw input shape.
   - **Optional:** a heuristics `.csv` or `.tsv` (extra rules layered on top of the mapping). See [references/bacdive_example/heuristics.tsv](references/bacdive_example/heuristics.tsv) for a representative example. When heuristics conflict with the mapping or its notes column, the heuristics file wins.
3. **How to iterate raw records** — confirm there is an existing iterator (e.g. `iter_<name>_records()`) or that the user has already provided one in a reference parser (like [references/bacdive_example/bacdive.py](references/bacdive_example/bacdive.py)). The parser you generate should call it; do not invent a new fetch strategy.
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
3. **The two boilerplate helpers, copied verbatim** from [references/base_func.py](references/base_func.py):
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
    "_id": f"<name>_{_id}",     # OR just str(_id) — decide per §"Deciding the _id"
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

Use [references/bacdive_example/bacdive.py](references/bacdive_example/bacdive.py) as the compact worked example and [bacdive/files/bacdive_crawler.py](bacdive/files/bacdive_crawler.py) as the full-source example.

#### Deciding the `_id`: prefix vs. bare identifier

The `_id` must be **globally unique** (records sharing an `_id` may be merged downstream). Historically we always prepended `<name>_`, but stripping that prefix later — when a mirror repo is found — breaks users' saved URLs. So **only prepend when the raw identifier actually risks a collision:**

- **Structured mix of letters + numbers** (e.g. `GSE12345`, `PRJNA398089`, `E-MTAB-1234`) → **bare identifier**: `"_id": str(_id)`.
- **All digits or all letters** → **prepend**: `"_id": f"<name>_{_id}"`, unless the scheme is well-established/well-adopted/registered with identifiers.org (then bare).
- **Ambiguous** → default to prepending and review against existing `_id`s.

Read [references/id_decision.md](references/id_decision.md) for the full decision tree and the production-change tracking rule. `identifier` always stays the raw source value (`str(_id)`) regardless of the `_id` choice. **State the chosen form and its one-line reason in the end-of-run report** so the user can override.

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

The upload-time validator in [biothings-hub/files/nde-hub/utils/validate.py](biothings-hub/files/nde-hub/utils/validate.py) (`check_schema`) **rejects records** with out-of-enum values for:

- `conditionsOfAccess` ∈ `{"Open", "Restricted", "Closed", "Embargoed", "Varied"}`
- `creativeWorkStatus` (only when `@type == "Sample"`) ∈ `{"Bespoke", "Available", "Backordered", "Retired"}`

If the mapping wants either of these fields, map source values to the allowed enum. If you cannot map confidently, omit the field rather than emit an invalid value.

Also: `version` must not be set (`check_schema` asserts `doc.get("version") is None`).

### 5. pmids / pmcs — exception

`pmids` and `pmcs` are **not** in [nde.py](biothings-hub/files/nde-hub/hub/dataload/nde.py)'s mapping — the pipeline's `citations` stage consumes them and removes them from the record. If the source provides PubMed / PMC identifiers and the mapping requests `pmids` or `pmcs`:

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
- `uploader.py` — normally just two lines of body:
  ```python
  from hub.dataload.nde import NDESourceUploader


  class <Name>Uploader(NDESourceUploader):
      name = "<name>"
  ```
  - Subclass `NDESourceUploader` if `--type=Dataset`, `NDESourceSampleUploader` if `--type=Sample`. Import from `hub.dataload.nde`.
  - Set `name = "<name>"`.
  - **Do not write a `load_data` method.** The base class already reads `data.ndjson` and runs every standardizer (see §7). Add `__metadata__` if the source needs `src_meta` / a merger.
  - Only override `load_data` when the records need custom parsing that the crawler's `data.ndjson` doesn't already give you (see §7's "When to override `load_data`").

### 7. Standardizers

`@nde_upload_wrapper` (in [biothings-hub/files/nde-hub/utils/pipeline.py](biothings-hub/files/nde-hub/utils/pipeline.py)) runs the same pipeline for every source. Each stage checks the records in front of it — and, for the curated-file stages, whether that source has a lookup file on disk — and skips itself if there is nothing to do. Uploaders do not import or call these.

| Stage | Runs when a record has | What it does |
|---|---|---|
| `citations` | `pmids`, `pmcs` or `citation.doi` | Citation + funding from NCBI E-utilities, plus PubTator species / diseases |
| `funding` | `funding` | Curated NIH grant from the funding cache; CrossRef funder names |
| `terms` | `species`, `infectiousAgent` or `healthCondition` | Standardizes them, splits hosts from infectious agents |
| `descriptions` | an eligible record type, a `description`, and missing taxonomy or `healthCondition` | Mines species + diseases out of the text via EXTRACT. Eligible types are `Dataset`, `DataCollection`, `ResourceCatalog`, or `Sample` with `additionalType: "BioSample"` |
| `measurement_technique` | `measurementTechnique` **and** `/data/nde-hub/standardizers/measurement_technique_lookup/<source>.csv` | Maps repository techniques to ontology terms |
| `nctid` | `nctid` **and** `/nvme/nde-hub/standardizers/nctid_lookup/nctid.csv` | measurementTechnique from the trial's study design |
| `topic_category` | `/data/nde-hub/topic_categories/<source>.json` | Adds curated EDAM topics |
| `disambiguating_description` | `/data/nde-hub/disambiguating_descriptions/<source>.csv` | Adds the curated summary |
| `lineage` | always | `_meta.lineage` for the taxonomy browser |

Then every record gets `sourceOrganization` corrections, `date`, `_meta.completeness`, a cleaned description, placeholder-term removal and `check_schema`.

The pipeline collects repository-wide statistics for every active stage and logs stage-level record counts plus before/after counts for each field the stage manages. Statistics are emitted even when validation stops an upload early. `records_processed` counts all records in a batch where the stage ran; use `records_changed` for the records whose tracked fields actually changed. `post_process` has stage-level counts only and no field tracking. Do not add duplicate statistics logging to a new uploader.

Two optional settings on the uploader change the pipeline:

- `post_process(self, doc)` — the source's last word on a record, after every stage but before `lineage`. Return the document, or `None` to drop it. Use it for source-specific fixups or filtering (see [bei](biothings-hub/files/nde-hub/hub/dataload/sources/bei/uploader.py), [pdb](biothings-hub/files/nde-hub/hub/dataload/sources/pdb/uploader.py), [covid_radx](biothings-hub/files/nde-hub/hub/dataload/sources/covid_radx/uploader.py)).
- `skip_stages = ("descriptions",)` — opt out of a stage. Use it when a source's volume makes a stage's per-record API calls impractical, not to express a metadata decision.

**When to override `load_data`:** only when the records don't come straight out of `data.ndjson` — a custom parser ([ncbi_geo](biothings-hub/files/nde-hub/hub/dataload/sources/ncbi_geo/gse_uploader.py)), per-file jobs ([biostudies](biothings-hub/files/nde-hub/hub/dataload/sources/biostudies/uploader.py)), or source-restricted curation ([dde](biothings-hub/files/nde-hub/hub/dataload/sources/dde/uploader.py), [vivli](biothings-hub/files/nde-hub/hub/dataload/sources/vivli/uploader.py)). Then keep `@nde_upload_wrapper` on it and yield plain dicts:

```python
from hub.dataload.nde import NDESourceUploader
from utils import nde_upload_wrapper

from .parser import my_parser


class <Name>Uploader(NDESourceUploader):
    name = "<name>"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        yield from my_parser(data_folder)
```

`utils.iter_ndjson(data_folder)` yields the crawler's records if you need them alongside custom logic. Never re-implement a pipeline stage inside `load_data`.

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

## Validate the generated source

Before reporting done, invoke and follow [`/validate-source`](../validate-source/SKILL.md) for the generated source and selected record type. Treat that validation-and-repair pass as mandatory; fix every issue it finds. Keep the verification and `@type` rules centralized there instead of duplicating them in this skill.

## Demonstration

[references/bacdive_example/](references/bacdive_example/) (bundled with this skill) contains a worked example:
- `mapping.tsv` — the mapping
- `heuristics.tsv` — rules layered on top of the mapping (heuristics wins on conflict, per §"Inputs to collect")
- `example1.json` — a trimmed raw input record (one strain, one entry per section the parser touches)
- `bacdive.py` — a reference parser (matches the boilerplate-output shape this skill must emit)

The corresponding skill outputs already exist in the repo for inspection:
- [bacdive/](bacdive/) (crawler container directory)
- [biothings-hub/files/nde-hub/hub/dataload/sources/bacdive/](biothings-hub/files/nde-hub/hub/dataload/sources/bacdive/) (hub source)
- The `bacdive-crawler` entry in [docker-compose-crawlers.yml](docker-compose-crawlers.yml)

When asked to demonstrate, run the skill against the bundled `references/bacdive_example/` inputs and produce output that matches the structure of those existing files, then run `/validate-source bacdive --type Sample`.

## Report at end

When done, post a short summary listing every file created or modified, the validation commands run, and any mapping rows you could not confidently convert. Also state the `_id` decision (prefixed `<name>_...` vs. bare `identifier`) and the one-line reason from §"Deciding the `_id`", so the user can override if they know of a collision risk. Mention that per-utility repository statistics will be available in the upload logs after the first run.
