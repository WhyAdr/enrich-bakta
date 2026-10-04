# Refined implementation plan: InterProScan annotation enrichment

**Review date:** 2026-10-04
**Reviewed against:** `f5de296` (`origin/main`)
**Status:** design reviewed and validated; implementation remains pending

This document is the corrected implementation specification for adding
InterProScan evidence to `enrich-bakta`. It replaces the earlier Gemini draft's
non-executable pseudocode and records the checks that were actually run against
the published corpus and current package.

## Review verdict

The direction is sound: use InterProScan TSV as the runtime evidence source,
preserve the pristine Bakta bytes, promote only bounded feature-level evidence,
and keep pathway context out of GenBank qualifiers. The draft was not ready to
implement unchanged.

The main corrections are:

1. The InterProScan MD5 is TSV column 2, not column 1. The official TSV
   contract defines columns 1--13 as the required match fields, with GO and
   pathway columns 14 and 15 optional. The materialized corpus has 13-column
   unintegrated rows and 15-column integrated rows; a parser must accept those
   shapes and reject malformed shapes. See the
   [official InterProScan output-format specification](https://github.com/ebi-pf-team/interproscan-docs/blob/v5/docs/OutputFormats.rst).
2. The published evidence is under `data/<sample>/evidence/interproscan/`;
   the draft's `InterProScan-datasets/<sample>/` root directories do not exist.
   `docs/data/PUBLISHED-MANIFEST.json` is the clean-clone authority.
3. The draft's C14 eggNOG overlap value of 4,528 is not reproducible under a
   defined `(query_id, GO term)` pair key. The validated overlaps are 2,063
   pairs for C14 and 2,054 for SM. The 3,858/3,857 overlaps with existing
   Bakta GO qualifiers and the 5,549/5,475 novel pair counts are reproducible.
4. No InterProScan source module or CLI exists in the current package. The
   proposed CLI and real-data output claims are therefore acceptance targets,
   not current capabilities.
5. The draft calls `qualifier_insertion()` and `comment_insertion()` with the
   wrong signatures and constructs `CandidateDecision` objects directly. The
   current architecture expects source adapters to return insertions and
   evidence rows; `build_candidate_ledger()` then creates and validates typed
   candidate decisions centrally.
6. The existing lineage contract is
   `--translation-evidence-manifest` plus
   `--allow-imported-translations`. Do not introduce a parallel
   `--allow-unverified-lineage` flag for this feature.
7. The existing unified CLI has one `--context-report` option. InterProScan
   must extend that contract rather than add a competing
   `--interproscan-context-report` output path.
8. A wall-clock date in an automatic COMMENT marker would violate byte
   idempotency. Use the existing explicit `--merge-timestamp` policy and make
   the source hash/version marker deterministic by default.

## Evidence checked

The following checks were run without rewriting repository data:

- `python tools/validate_dataset_manifest.py docs/data/MANIFEST.json`:
  86 artifacts validated.
- `python tools/validate_dataset_manifest.py
  docs/data/PUBLISHED-MANIFEST.json`: 42 artifacts validated.
- All nine published InterProScan TSV/GFF3/JSON SHA-256 values match the
  published manifest. They are tracked through Git LFS.
- The current baseline suite passes: `146 passed`.
- The committed Python source scope passes Ruff, Ruff format, and
  `python -m mypy src`. The exact `ruff check .` command in this working
  tree also visits two untracked follow-through helper scripts; those files are
  intentionally not part of this plan or release input. A clean clone is the
  authoritative repository-wide gate.
- `gbparse` summaries/validation and the canonical raw-byte parser confirmed
  the base GenBank records and translation semantics without loading full
  flatfiles into the review context.

### Validated corpus baseline

The manifest remains authoritative for file sizes and hashes. These compact
metrics were independently recomputed by streaming the TSVs:

| Metric | BK71A | C14 | SM |
|---|---:|---:|---:|
| Target CDS features | 4,131 | 4,636 | 4,576 |
| InterProScan data rows | 32,791 | 40,234 | 39,684 |
| Unique query CDSs | 3,834 | 4,457 | 4,381 |
| CDSs without an InterProScan query | 297 | 179 | 195 |
| Rows with an InterPro accession | 22,360 | 27,420 | 27,044 |
| Rows without an InterPro accession | 10,431 | 12,814 | 12,640 |
| Unique CDS--InterPro pairs | 12,994 | 16,160 | 15,930 |
| Unique InterPro accessions | 5,415 | 6,282 | 6,254 |
| CDS--GO pairs | 6,880 | 9,407 | 9,332 |
| Unique GO terms | 1,354 | 1,725 | 1,731 |
| Pfam pairs | 5,131 | 6,219 | 6,088 |
| Pfam pairs already represented by Bakta notes | 9 | 5 | 15 |
| Novel Pfam pairs | 5,122 | 6,214 | 6,073 |
| TIGRFAM pairs | 1,454 | 1,766 | 1,758 |
| CDD pairs | 2,183 | 2,591 | 2,571 |
| Reactome pathway tokens | 3,717,389 | 4,534,168 | 4,416,062 |
| MetaCyc pathway tokens | 2,177,198 | 2,441,941 | 2,410,944 |

Identity checks found zero missing query locus tags, zero inconsistent TSV MD5s,
zero FAA MD5 mismatches for C14/SM, and zero target-GBFF translation MD5
mismatches for all three samples. The BK71A comparison is against the restored
GBFF translations and does not establish pristine Bakta input lineage.

GO overlap is defined explicitly as the same `(query_id, normalized GO term)`:

| Pair relationship | C14 | SM |
|---|---:|---:|
| InterProScan pairs already in pristine Bakta | 3,858 | 3,857 |
| Novel InterProScan pairs | 5,549 | 5,475 |
| InterProScan pairs also present in eggNOG | 2,063 | 2,054 |

The last row is a diagnostic baseline for cross-source support, not a claim
that all overlapping evidence should be emitted twice.

## Scientific and provenance policy

### Inputs and lineage

- C14 and SM use the pristine inputs
  `data/C14/bakta/C14-NMZ.gbff` / `.faa` and
  `data/SM/bakta/SM-NMZ.gbff` / `.faa`.
- BK71A uses `data/BK71A/derived/restored/BK71A-restored.gbff` only as a
  historical target. The repository has no pristine BK71A Bakta FAA and no
  bound translation-evidence manifest for that published artifact.
- A functional InterProScan merge must fail closed for a restored target unless
  `--translation-evidence-manifest` is supplied and the existing
  `--allow-imported-translations` policy is satisfied. The implementation must
  not silently treat BK71A's restored translations as pristine evidence.
- The TSV is the runtime parser input. Published GFF3 and JSON are retained for
  audit/provenance and optional one-time parity checks; production code must not
  `json.load()` the 823--997 MB JSON exports merely to obtain the version.
- `--interproscan-version` is explicit and required for a reproducible merge.
  The value is recorded in metadata, COMMENT provenance, and the sidecar. A
  repository manifest may document the version, but runtime auto-detection from
  a giant JSON file is not a supported fallback.

### Qualifier policy

| Evidence | Default output | Rule |
|---|---|---|
| InterPro accession | `/db_xref=InterPro:IPRxxxxxx` | strict accession grammar; existing exact values are `supported_existing` |
| GO term | `/db_xref=GO:ddddddd` | normalize optional source decorations; preserve raw token in evidence rows |
| Pfam | `/note=PFAM:PFxxxxx` | compare by base accession so `PFxxxxx.33` in Bakta suppresses a duplicate |
| TIGRFAM | `/note=TIGRFAM:TIGRxxxxx` | opt-in through the member-database allowlist; no free-text descriptions |
| CDD | context by default; optional normalized note policy | enable only after a dedicated CDD grammar/fixture decision |
| Reactome/MetaCyc pathways | never a feature qualifier | summarize in the context sidecar and retain the raw TSV as authoritative evidence |
| Source inference | one `/inference` per affected CDS | add only when a new InterProScan feature qualifier survives reconciliation |

The initial member-database allowlist is `Pfam,TIGRFAM`. Other member analyses
remain audit/context evidence until a specific normalized qualifier policy is
approved. Do not promote member descriptions, scores, dates, or pathway labels
to GenBank qualifiers. All outputs remain computational annotation evidence, not
proof of expression, enzyme activity, pathway completeness, or phenotype.

## Architecture aligned with the current package

### New source adapter

Create `src/enrich_bakta_lib/sources/interproscan.py` with these responsibilities:

1. `iter_interproscan_tsv(stream)` yields validated immutable hit records from a
   text stream. It must:
   - accept exactly 13, 14, or 15 columns, treating optional GO/pathway fields
     as absent when omitted;
   - validate non-empty query IDs, 32-hex MD5, positive sequence length,
     inclusive `1 <= start <= stop <= length`, and required analysis/signature
     fields;
   - normalize `-` to absent values;
   - accept bare GO IDs and official source-decorated GO tokens, normalizing to
     bare `GO:ddddddd` for promotion;
   - validate InterPro, Pfam, TIGRFAM, and optional CDD accessions with full
     matches, never prefix checks;
   - preserve row number and raw values for diagnostics; and
   - never materialize the input file as one string or list of raw rows.
2. Aggregate only the bounded per-query evidence needed for planning. Do not
   use `hits = list(iter_interproscan_tsv(...))` as in the draft. Pathway
   context must be counted or streamed into the sidecar; it must not create
   unbounded in-memory lists merely because the sidecar was requested.
3. Validate identity before planning any insertion:
   - every query maps to exactly one base CDS locus tag;
   - the TSV MD5 agrees with the normalized FAA sequence when FAA is supplied;
   - the FAA sequence agrees with the target CDS `/translation`; and
   - every repeated query has one consistent MD5 and protein length.
   Use the package's existing `parse_faa()`, `normalize_protein()`,
   `parse_genbank_bytes()`, `cds_by_locus()`, and translation-evidence loader.

The adapter should return an explicit plan object analogous to `EggnogPlan`:

```text
InterProScanPlan(
    insertions,
    evidence_rows,
    stats,
    context_report,
)
```

It must not construct `CandidateDecision` directly. For each proposed or
already-supported value, emit an `interproscan_candidate` evidence row with
`source`, query/feature identity, field, qualifier, raw and normalized values,
status, `candidate_role`, `support_class`, reason code, and emitted qualifier
summary. The central candidate-ledger builder will assign deterministic IDs,
insertion references, and cross-source support.

### Current insertion APIs to reuse

The implementation must call the existing APIs with their actual contracts:

```python
qualifier_insertion(
    base.data, feature, qualifier, value, source, source_value, order,
    candidate_role=..., evidence_class=...
)

comment_insertion(
    base, record, marker, lines, source, order
)
```

Do not reimplement qualifier formatting or byte splicing. Use deterministic
orders after the highest insertion order already planned by another source.
Use the existing `reconcile_insertions()` and `finalize_merge()` gates.

### Candidate-ledger integration

Make the smallest source-aware changes in `core/decisions.py` and
`core/merge_engine.py`:

- normalize `interproscan`/`interpro` to source ID `InterProScan`;
- recognize `interproscan_candidate` rows in `_row_source()` and
  `build_candidate_ledger()`;
- map the adapter's field names to the final qualifier/value in
  `_row_field_and_value()`;
- add `InterProScan` to unified `source_hashes` and metadata;
- map `InterProScan provenance` into the `interproscan` source family so
  provenance is removed if its functional support is removed;
- rely on the existing semantic reconciliation to collapse exact same-target,
  same-qualifier, same-value GO proposals from eggNOG and InterProScan; and
- do not add a one-off GO-specific duplicate branch unless a regression test
  proves the generic reconciliation path insufficient.

The current validator requires producer-provenance support from the same source
and evidence class. Since one InterProScan inference can support several
InterProScan evidence classes, use the recommended design: add a narrowly
scoped InterProScan provenance support-class map to the ledger builder and
validator. It must allow an InterProScan provenance candidate to support
accepted InterPro/GO/member candidates on the same target while retaining the
existing class-specific rule for Baktfold, KofamScan, and eggNOG. Do not flatten
all source evidence into one class merely to bypass this invariant, and do not
emit an unsupported provenance candidate with no valid support links.

### Context-report contract

Reuse the existing unified `--context-report` option. Preserve the existing
`enrich-bakta.eggnog-context.v1` shape for eggNOG-only runs. Define
`enrich-bakta.interproscan-context.v1` for InterProScan-only runs and a
versioned combined `enrich-bakta.context.v1` envelope when both sources are
present. The InterProScan section must include:

- source SHA-256, InterProScan version, TSV column shape, row/query counts;
- Reactome, MetaCyc, and other pathway token counts globally and per query;
- a deterministic digest or sorted bounded representation of pathway context;
- a statement that pathway tokens are not promoted to feature qualifiers; and
- no duplicate copy of the entire raw TSV when the published evidence file is
  already available.

The sidecar writer must be atomic and deterministic. No automatic timestamp may
enter its bytes.

## CLI and workflow changes

### Unified workflow

Extend the existing `enrich()` API and `enrich_bakta.py` CLI with:

```text
--interproscan PATH
--interproscan-version VERSION       required when --interproscan is used
--interproscan-member-dbs Pfam,TIGRFAM
```

Use the existing `--context-report`, `--manifest`,
`--translation-evidence-manifest`, `--allow-imported-translations`,
`--no-comment-note`, `--no-feature-provenance`, and `--merge-timestamp` options.
Do not add `--interproscan-context-report` or `--allow-unverified-lineage`.

The workflow must:

1. validate the base GenBank bytes before source planning;
2. require FAA for C14/SM InterProScan runs and validate it against both TSV
   and GBFF translations;
3. enforce the existing restored-input translation-evidence gate;
4. append the TSV to `other_inputs` so final manifests bind its bytes;
5. pass a monotonic insertion order to the adapter;
6. merge its rows into the existing reconciliation and candidate-ledger flow;
7. merge source context with any eggNOG context report; and
8. record source hash, version, parser statistics, identity counts, member-DB
   policy, lineage policy, and context-report hash in the JSON manifest.

### Standalone compatibility entry point

Retain the standalone source-workflow convention with a thin root shim
`merge_interproscan_bakta.py`. It must import the canonical package
implementation and contain no duplicate logic. Its CLI should use the same
named policy options and accept a minimal positional form documented by the
existing Kofam/eggNOG standalone entry points. The unified CLI remains the
canonical path for multi-source reconciliation.

## Implementation phases

### Phase 0: baseline and fixtures

- Keep the current 146-test baseline passing.
- Add small synthetic GBFF/FAA/TSV fixtures; do not put the 100+ MB evidence
  files into ordinary unit-test setup.
- Add a fixture with 13 columns, one with 15 columns, and one with a
  source-decorated GO token.
- Add deterministic hashes and expected counts to the test assertions, not to
  production code.

### Phase 1: validation and streaming parser

- Add field-specific value validators to `value_rules.py`.
- Implement strict row parsing, normalization, coordinate checks, and repeated
  query consistency checks.
- Test malformed column counts, MD5, numeric fields, coordinates, accessions,
  unknown query IDs, duplicate base locus tags, FAA missing/mismatch, and
  GBFF translation mismatch.

### Phase 2: planner and provenance

- Implement the aggregate/planner result object and evidence-row contract.
- Add InterPro, GO, allowlisted member signatures, inference, and COMMENT
  insertions through the existing insertion helpers.
- Test base-value suppression, versioned Pfam overlap, exact insertion order,
  unsupported member-DB rejection, and no-pathway-qualifier behavior.

### Phase 3: reconciliation and context

- Integrate source rows with `build_candidate_ledger()`.
- Add the source-family and provenance-support behavior selected above.
- Add exact duplicate collapse and cross-source `supporting_candidate_ids`
  tests for an eggNOG/InterProScan GO pair.
- Add standalone and unified context-report schema tests.

### Phase 4: workflow and compatibility

- Add unified CLI options and, if retained, the thin standalone wrapper.
- Test `--version`, missing-input errors, required version, FAA/lineage gates,
  `--context-report`, `--manifest`, and no-op reruns.
- Update README, package entry points, CHANGELOG, and release validation docs.

### Phase 5: real-data acceptance

For C14 and SM, run the TSV against the pristine Bakta GBFF/FAA and record:

- input SHA-256 values and parser row/query counts;
- zero identity/translation mismatches;
- output GenBank validation with `gbparse validate`;
- JSON manifest schema validation and candidate-ledger validation;
- absence of pathway tokens in feature qualifiers;
- deterministic context-report contents; and
- byte-identical second-run output with zero new insertions.

For BK71A, run the parser and identity checks against the restored GBFF, but
do not claim a reproducible functional enrichment acceptance run until a bound
translation-evidence manifest and the approved imported-lineage policy are
available. A failed-closed BK71A invocation is an expected acceptance case.

## Required regression tests

The implementation is not complete until tests cover all of the following:

1. 13/14/15-column parsing and rejection of other shapes.
2. MD5 parity against FAA and target GBFF translations.
3. Unknown and duplicate query identity failures.
4. Invalid InterPro, GO, Pfam, TIGRFAM, and CDD tokens.
5. Existing InterPro/GO/Pfam values become `supported_existing` without new
   qualifier lines.
6. Versioned `PFAM:PFxxxxx.33` suppresses a bare `PFxxxxx` candidate, without
   false prefix matches such as `PFxxxxx0`.
7. InterPro/GO/member additions are emitted once and are byte-idempotent.
8. Reactome/MetaCyc tokens never appear in GenBank feature qualifiers.
9. Context-only pathway records remain available in the sidecar and retain the
   source hash/row provenance.
10. Same-target GO support from eggNOG and InterProScan collapses to one output
    qualifier while both source candidates remain linked.
11. Producer provenance is linked to valid InterProScan support and is removed
    when all corresponding functional insertions are reconciled away.
12. Restored BK71A input fails without the existing translation-evidence gate
    and succeeds only with an explicitly valid lineage manifest/policy.
13. The output passes `gbparse validate`; the manifest passes
    `tests/test_manifest_schema.py` and `validate_candidate_ledger()`.

## Acceptance commands

Run these from a clean clone with Git LFS materialized. The current checkout's
untracked follow-through bundle must not be staged or used as a release input.

```powershell
$env:PYTHONDONTWRITEBYTECODE = 1
python -m pytest -q --basetemp .test-output/interproscan-test -p no:cacheprovider
python -m ruff check .
python -m ruff format --check .
python -m mypy
python tools/validate_dataset_manifest.py docs/data/MANIFEST.json
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
```

After implementation, the canonical C14 smoke run is:

```powershell
python enrich_bakta.py `
  --bakta data/C14/bakta/C14-NMZ.gbff `
  --faa data/C14/bakta/C14-NMZ.faa `
  --interproscan data/C14/evidence/interproscan/C14-NMZ.interproscan.tsv `
  --interproscan-version 5.59-91.0 `
  --output .test-output/C14-interpro-enriched.gbff `
  --manifest .test-output/C14-interpro-enriched.manifest.json `
  --context-report .test-output/C14-interpro-context.json
```

The SM run is identical with the SM paths. Validate each output with:

```powershell
gbparse validate .test-output/C14-interpro-enriched.gbff --format json
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
python -m pytest -q tests/test_manifest_schema.py tests/test_merge_pipeline.py
```

The idempotency gate must compare the output bytes from two runs and assert
zero new insertions on the second run. It must also inspect the parsed output,
not merely compare a manifest count, to prove that no pathway token entered a
feature qualifier.

## Release and handoff

This feature adds a new evidence source and CLI surface. The recommended
release target is `0.4.0`, with the exact version decision recorded before
implementation. Update the package version, CHANGELOG, README, CI matrix, and
published acceptance record together.

The plan is ready for implementation after these refinements. It does not
claim that InterProScan enrichment is already implemented. The untracked
`enrich-bakta-followthrough-bundle-2026-10-04/` directory is pre-existing local
review material and must remain outside any implementation commit unless the
owner explicitly requests otherwise.
