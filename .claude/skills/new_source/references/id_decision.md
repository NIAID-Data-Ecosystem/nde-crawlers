# Deciding the `_id`: prefix vs. bare identifier

The `_id` value for each record **must be globally unique** so that URLs stay RESTful. Two records with the same `_id` may be merged downstream — which is *desirable* when the same record is mirrored across repositories (they should share an `_id`), but *harmful* when two genuinely different records collide.

Historically we prepended the repository short-name (`<name>_<identifier>`) to guarantee uniqueness. But when a new repo turns out to mirror existing records, that prefix has to be stripped so the records can merge — and that retroactive `_id` change breaks downstream users who saved the old URLs. So **only prepend the prefix when the raw identifier actually risks a collision.** Decide once, per source, using this procedure, and apply the same choice to every record.

## Procedure

1. **Does the source identifier contain a mix of letters and numbers?**
   - **Yes → is there a distinct, structured pattern** (e.g. `GSE12345`, `PRJNA398089`, `SAMN0000123`, `E-MTAB-1234`)?
     - **Yes → use the bare `identifier`** (no prefix), *unless* there is specific concern it could inadvertently merge with another source's ids.
     - **No** (letters and digits but no recognizable structure) → treat like the all-letters/all-digits case below (fall through to step 2's judgment).
   - **No → go to step 2.**
2. The identifier is **all letters or all digits** (e.g. a bare accession number `398089`, or an all-alpha slug). These are easy to confuse with another repository's ids, so ask three qualifying questions about the *resource/identifier scheme*:
   - Is the resource **over 20 years old and well-established**?
   - Is it **well-adopted within the community**?
   - Is it **registered with [identifiers.org](https://identifiers.org)**?

   Then:
   - **≥ 2 "yes" → use the bare `identifier`** (the scheme is stable and recognizable enough to stand alone).
   - **≥ 2 "no" → review case-by-case:** check the existing `_id`s already in use for similarly-structured identifiers before deciding. When in doubt, **prepend `<name>_`** — it is the safer default for ambiguous ids.
3. **Digit-only or plain-letter identifiers with no qualifying strength → prepend `<name>_`.**

## Applying the decision

- **Bare identifier:** `"_id": str(_id)` (drop the `f"<name>_{...}"`).
- **Prefixed:** `"_id": f"<name>_{_id}"`.
- Either way, `identifier` stays the raw source value (`str(_id)`).

**Surface the decision to the user** in the end-of-run report: state which form you chose and the one-line reason (which heuristic branch it fell under), so they can override if they know of a collision risk.

> **Once a source is in production, an `_id` change must be tracked.** If a prefix decision has to be revised after go-live (e.g. a duplicate repo is discovered), flag it explicitly rather than silently changing the `_id` format — downstream saved URLs depend on it.
