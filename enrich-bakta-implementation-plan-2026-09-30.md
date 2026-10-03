# Draft implementation plan: enrich-bakta reorganization and hardening

Status: draft for approval. This document is planning only: no source files,
datasets, generated artifacts, or remote branches are changed by this plan.

Date: 2026-09-30 (Asia/Jakarta)

Repository baseline:

- Remote: `https://github.com/WhyAdr/enrich-bakta.git`
- Branch: `main`, currently at `5d6754d`
- Python package version: `0.2.0`
- Current checkout path: `D:\\W\\baktfold`
- Current worktree: tracked tree clean; untracked `enrich-bakta-main-audit-2026-09-29.md`
- Audit input: `enrich-bakta-main-audit-2026-09-29.md`, Astra's 2026-09-29 audit

## Decisions already made

1. Use a domain-based package layout rather than one arbitrary folder per script.
2. Preserve the existing console commands, documented direct script invocations,
   and compatibility-facing module imports where practical.
3. Retain `normalize_baktfold.py` as a clearly isolated legacy utility.
4. Organize and publish the reproducible C14, SM, and BK71A data artifacts so
   remote agents can clone the repository and exercise the workflows.
5. Rename the local checkout from `baktfold` to `enrich-bakta` after the
   repository and documentation migration is stable.

The publication decision is subject to the data-governance gates in this plan.
"Open source" is not by itself evidence that every upstream-derived artifact
is licensed, privacy-cleared, or permitted for redistribution.

## Current script roles

| Current file | Role | Disposition |
|---|---|---|
| `merge_engine.py` | Raw GenBank model, semantic checks, byte-preserving insertion, reconciliation, atomic writes, and manifests | Supported core module |
| `enrich_bakta.py` | Canonical multi-source orchestration CLI | Supported top-level workflow |
| `graft_baktfold_additions.py` | Baktfold source adapter and planner | Supported source adapter |
| `merge_kofamscan_bakta.py` | KofamScan/KOALA source adapter and optional Baktfold merge | Supported source adapter |
| `merge_eggnog_bakta.py` | eggNOG TSV/XLSX parser, planner, merger, and context sidecar | Supported source adapter |
| `restore_bakta_translations.py` | Explicit translation-restoration workflow for selected pseudogene CDSs | Supported but safety-sensitive workflow |
| `normalize_baktfold.py` | Historical line-oriented normalizer superseded by byte-preserving insertion | Legacy; retain outside the supported package and quality gates |

The supported dependency shape is:

```text
merge_engine
├── Baktfold adapter
├── KofamScan adapter ──> Baktfold adapter
└── eggNOG adapter

unified enrichment workflow ──> all supported adapters
translation restoration      ──> eggNOG adapter + merge_engine
```

## Target repository layout

The canonical implementation will move to a real package with small public
interfaces and deep implementations behind them:

```text
enrich-bakta/
├── src/
│   └── enrich_bakta/
│       ├── __init__.py
│       ├── __main__.py
│       ├── core/
│       │   ├── __init__.py
│       │   └── merge_engine.py
│       ├── sources/
│       │   ├── __init__.py
│       │   ├── baktfold.py
│       │   ├── kofam.py
│       │   └── eggnog.py
│       ├── workflows/
│       │   ├── __init__.py
│       │   ├── enrich.py
│       │   └── restore_translations.py
│       └── legacy/
│           └── normalize_baktfold.py
├── compat/
│   └── README.md
├── tests/
├── docs/
│   ├── audits/
│   ├── architecture/
│   ├── data/
│   └── plans/
├── data/
│   ├── C14/
│   ├── SM/
│   ├── BK71A/
│   ├── fixtures/
│   ├── MANIFEST.tsv
│   └── README.md
├── .gitattributes
├── .gitignore
├── pyproject.toml
└── README.md
```

The old root-level Python names will remain as thin compatibility launchers
while the migration is validated:

```text
enrich_bakta.py                 -> enrich_bakta.workflows.enrich
graft_baktfold_additions.py    -> enrich_bakta.sources.baktfold
merge_kofamscan_bakta.py       -> enrich_bakta.sources.kofam
merge_eggnog_bakta.py          -> enrich_bakta.sources.eggnog
restore_bakta_translations.py  -> enrich_bakta.workflows.restore_translations
normalize_baktfold.py          -> enrich_bakta.legacy.normalize_baktfold
```

The shims will re-export the currently used Python symbols and call `main()`
when invoked directly. Existing `project.scripts` names remain unchanged.
Before finalizing the package name, add an isolated import-resolution test for
the coexistence of the `enrich_bakta/` package and the compatibility launcher
named `enrich_bakta.py`; if source-tree precedence is ambiguous, use a private
implementation package name behind the same public distribution and launchers
rather than sacrificing direct invocation compatibility.

## Work phases

### Phase 0: freeze the baseline and preserve evidence

1. Preserve the untracked Astra audit unchanged until it is copied into
   `docs/audits/`.
2. Record the current commit, worktree status, remote URL, Python version,
   dependency versions, CLI help/version output, and installed-wheel behavior.
3. Run the existing gates before any move:

   ```text
   python -m pytest -q
   python -m ruff check .
   python -m ruff format --check .
   python -m mypy
   python -m build
   git diff --check
   ```

4. Test the five installed commands from outside the checkout and record the
   current direct `python <script>.py` examples from the README.
5. Hash every candidate data artifact before moving it. No data publication or
   source refactor proceeds if the baseline cannot be distinguished from the
   post-migration state.

### Phase 1: reorganize the package without changing behavior

1. Move implementation modules into the target package paths.
2. Add package-qualified imports and explicit `__all__` exports at the seams.
3. Update `pyproject.toml` from flat `py-modules` to package discovery and
   point all console scripts at the canonical package implementations.
4. Keep root compatibility launchers for the old script paths and symbols.
5. Update tests to exercise canonical package imports, then add compatibility
   tests for the old imports and direct script invocations.
6. Keep `normalize_baktfold.py` outside the supported package, Ruff target,
   mypy file list, and normal release path.
7. Run the full existing gate after the move, plus an isolated wheel install
   that proves the legacy normalizer is not accidentally packaged.

No semantic or annotation-transfer policy changes belong in this phase.

### Phase 2: implement Astra's release-blocking semantic fixes

The existing byte-preservation invariant must remain separate from the semantic
model. The implementation should retain raw byte spans for splicing while using
a complete, authoritative semantic view for identity and transfer decisions.

#### A01 — complete multiline feature locations

- Accumulate all location continuation lines before the first qualifier.
- Attach complete raw features to the corresponding Biopython record/feature
  by record and feature ordinal, asserting type/count alignment.
- Compare complete location structure: part order, strand, joins, circular
  origin crossings, fuzzy boundaries, and remote references. Never reduce a
  compound location to a minimum/maximum span.
- Use the same complete location model for strict parity and eggNOG paired
  gene/CDS matching.
- Add failures for different continuation parts and successes for equivalent
  line wrapping, including `complement(join(...))`, origin crossings, and fuzzy
  locations.

#### A02 — preserve semantic qualifier values

- Keep raw spans for output, but obtain semantic qualifier values from a
  Biopython-aligned view or a field-aware decoder.
- Do not concatenate every continuation with `strip()`; free text, quoted
  sequences, translations, and inference values need distinct rules.
- Assert that newly emitted qualifier values have the intended semantic value,
  not merely that Biopython can parse the resulting bytes.
- Add tests for multiline free text, escaped quotes, sequence-like values,
  wrapped inference values, and rerun/idempotence behavior.

#### A03 — make translation restoration an explicit evidence import

- Separate `imported_translation` from independently `validated_translation`
  in the planner, manifest, and user-facing messages.
- Require an explicit restoration policy flag for the workflow; do not let a
  copied FAA sequence silently acquire ordinary FAA/GBFF identity status.
- Prefer source-run provenance where available: matching nucleotide/feature
  hashes, producer metadata, and the original FAA lineage.
- Validate supported amino-acid symbols, length plausibility, genetic-code
  metadata, codon start, exceptions, and partial/disrupted cases without
  claiming that every pseudogene should translate as an intact CDS.
- Fail closed on empty or duplicate `/translation` qualifiers. Preserve
  `/pseudo` or `/pseudogene` and record unresolved cases explicitly.
- Keep functional transfers from restored sequences behind the explicit policy;
  do not describe an imported sequence as independent genomic confirmation.
- Add tests for unrelated same-tag FAA entries, impossible lengths, invalid
  symbols, internal stops, legitimate partial/exceptional translations, empty
  translations, duplicate translations, and provenance separation.

#### A04 — require a supported eggNOG schema contract

- Add a version/schema registry for the supported v3 contract.
- Reject unknown versions, missing metadata without an explicit compatible
  schema selection, empty legends, repeated/conflicting version declarations,
  and unsupported v2 layouts.
- Permit fallback only for a documented known version and record the selected
  schema ID and source of the confidence order in the manifest.
- Apply the same contract to TSV and XLSX and add one-hot tests for all 13
  confidence positions.

#### A05 — make reconciliation statuses describe the final output

- Give candidates stable identities and distinguish `planned_status` from
  `final_status`.
- Link surviving candidates to actual insertion IDs and parsed output features.
- Reconcile functional candidates before calculating producer provenance and
  counts.
- Mark suppressed or conflict-only evidence as suppressed/conflicting rather
  than `emitted`; ensure candidate rows and final bytes agree.
- Test default skip plus both preference policies, shared support, and
  conflict-only sources end to end.

### Phase 3: implement the remaining audit hardening

#### A06 — protect authoritative paired gene values

- Validate the existing Bakta gene/CDS pair before planning Baktfold names.
- Do not add a source name to one half when the other half has an authoritative
  conflicting value.
- Detect internally inconsistent source pairs separately from cross-source
  conflicts and report both without silently repairing existing annotations.

#### A07 — centralize field-specific value validation

- Validate full and partial EC syntax according to an explicit Baktfold policy.
- Require nonempty identifiers after `afdb_v6:`, `cath:`, and `pdb:` prefixes.
- Normalize missing-name sentinels before suffix cleaning and reject empty
  results such as a suffix-only `_123` value.
- Preserve raw values and skip reasons in manifests; do not imply that syntax
  validation proves ontology membership or biochemical activity.

#### A08 — deduplicate complete structured tokens

- Replace substring matching with shared, token-boundary-aware extraction for
  GO, KO, CAZy, and related structured values.
- Distinguish recognized annotation tokens from negated or incidental prose.
- Add prefix/suffix, punctuation, numeric-suffix, and CAZy family/subfamily
  regressions across Kofam and eggNOG.

#### A09 — make identity validation consistently fail closed

- Centralize cardinality checks before building locus indexes.
- Reject duplicate identical as well as conflicting locus qualifiers where the
  interface requires uniqueness; retain cross-record duplicate rejection.
- Diagnose ambiguous or missing targeted IDs instead of silently skipping them.
- After the table header, distinguish data-shaped `#`/`##` query rows from
  recognized producer comments and footers so malformed rows cannot disappear
  without a diagnostic.

#### A10 — bind XLSX parsing to the hashed snapshot

- Parse XLSX from the captured bytes using `io.BytesIO`, not by reopening the
  mutable path.
- Close workbooks in `finally` blocks and normalize malformed archives to a
  clear `MergeError`.
- Keep formula rejection and required-header checks explicit.
- Add a path-replacement race test proving the parsed annotations and hash come
  from the same snapshot.

#### A11 — preflight and stage all output artifacts

- Reject directory destinations, file parents, and ancestor/descendant output
  relationships before replacing any existing output.
- Serialize and validate every artifact first, then promote staged artifacts.
- Preserve a prior output on preflight or staging failure.
- Document the remaining cross-file crash boundary; do not claim that multiple
  atomic renames form a single transaction.
- Add injected write/replace failure tests and retain the existing alias,
  symlink, and hard-link collision checks.

#### A12 — make manifests reconstruct decisions

- Record Baktfold candidate outcomes, parity diagnostics, policy flags, schema
  and format information, and actual feature support.
- Fix missing-paired-gene diagnostics before existing-gene handling.
- Record `clean_gene_suffix`, provenance flags, tool version in every sidecar,
  and validated/no-op restoration targets.
- Use typed JSON values where JSON is the canonical machine-readable format;
  treat display-oriented TSV as a projection rather than the reconstruction
  authority.
- Decide and record the policy for Baktfold translation mismatches instead of
  silently allowing downstream functional transfers without diagnostics.

### Phase 4: reorganize documentation and make the repository navigable

Move existing documentation without losing history or source grounding:

| Current file | Target |
|---|---|
| `bakta-enhance-integrate-review.md` | `docs/architecture/implementation-review-2026-08-25.md` |
| `enrich-bakta-audit.md` | `docs/audits/independent-audit-2026-08-25.md` |
| `enrich-bakta-main-audit-2026-09-29.md` | `docs/audits/enrich-bakta-main-audit-2026-09-29.md` |
| `integrate-eggnog-plan.md` | `docs/plans/eggnog-enrichment.md` |
| `merge-koala-baktfold-plan.md` | `docs/plans/kofam-baktfold-merge.md` |

Update all links to repository-relative paths. Remove hard-coded
`d:/W/baktfold`/`D:\\W\\baktfold` links and add a short status banner to
historical plans and audits. Keep `README.md` at the repository root as the
current user-facing entry point.

### Phase 5: organize datasets and reproducible fixtures

Treat producer outputs as data artifacts, not as Python/tool caches. The current
checkout contains approximately 0.5 GB of local data and generated state. The
largest observed individual data artifact is approximately 47 MB; the four
producer-output directories are each tens to hundreds of megabytes.

Proposed mapping:

```text
data/
├── C14/
│   ├── inputs/       # Kofam, eggNOG, timing, and other source tables
│   ├── baktfold/     # C14 Baktfold producer output
│   ├── bakta/        # C14 Bakta/enriched/restored outputs
│   └── README.md
├── SM/
│   ├── inputs/
│   ├── baktfold/
│   ├── bakta/
│   └── README.md
├── BK71A/
│   └── derived/      # restored GBFF and provenance
├── fixtures/         # small deterministic data used in normal CI
├── MANIFEST.tsv
└── README.md
```

For every published artifact, record at minimum:

- relative path and sample/cohort label;
- role: upstream input, producer output, derived output, expected output, or
  regression fixture;
- byte size and SHA-256;
- producer/tool/database versions and the command used to create it;
- input lineage and whether the file is pristine, restored, enriched, or
  historical;
- source citation, license, redistribution permission, and any restrictions;
- privacy/sensitivity review status;
- validation status and the code commit that produced the artifact.

The following remain ignored and are regenerated by each agent or CI worker:

- `.mypy_cache/`
- `.pytest_cache/`
- `.ruff_cache/`
- `__pycache__/`
- `.test-output/`, unless a deliberately selected result is copied into
  `data/fixtures/` with a manifest entry

These directories are machine state, not reproducible datasets. Publishing them
would make clones larger and less portable without helping a cloud agent test
the pipeline. The reproducible alternative is to publish the curated input and
expected-output artifacts plus the commands and versions needed to regenerate
test state.

Before staging any dataset:

1. check upstream licenses and redistribution terms;
2. check for human, clinical, identifying, or otherwise restricted metadata;
3. check whether derived annotation files inherit source restrictions;
4. validate file formats and remove accidental credentials or local paths;
5. generate and verify the manifest and checksums;
6. choose ordinary Git, Git LFS, or release-artifact storage based on the
   measured repository size and current hosting limits;
7. run the dataset validation command from a clean clone before publication.

Do not publish derived enriched outputs as scientifically validated until the
release-blocking audit work and local C14/SM reruns have passed. Raw or
historical artifacts can be published earlier only with an explicit
`historical/unverified` status in the manifest.

### Phase 6: add remote-agent validation

Add a documented, deterministic workflow that a clean VM can run:

```text
python -m pip install -e ".[dev,xlsx]"
python -m pytest -q
python -m enrich_bakta --help
python -m enrich_bakta --version
python scripts/validate_dataset_manifest.py data/MANIFEST.tsv
```

Add a separate data-validation job or manually invoked workflow for the larger
C14/SM artifacts. Keep ordinary pull-request CI focused on small fixtures and
the synthetic adversarial tests; do not make every code review regenerate the
full dataset unless runtime and hosting costs are measured and acceptable.

The data-validation workflow must verify hashes, parse expected GenBank files,
check manifest lineage, exercise representative merge commands, and report
annotation/provenance deltas rather than only feature counts.

### Phase 7: final release validation and publication

Run all of the following from a clean checkout and from an isolated installed
wheel:

- existing tests plus every new A01–A12 regression;
- Ruff check and format check;
- mypy;
- wheel/sdist build, isolated install, `pip check`, and commands from outside
  the source tree;
- `git diff --check`;
- direct legacy-compatible script invocation and all five console commands;
- C14/SM reruns with annotation, manifest, hash, and idempotence comparison;
- data manifest/checksum validation;
- review of generated files for local absolute paths, credentials, and
  unintentional temporary/cache content.

Use small commits with one reviewable concern per commit:

1. this plan and audit preservation;
2. no-op package reorganization and compatibility shims;
3. A01/A02/A09 semantic identity hardening;
4. A03 restoration evidence policy;
5. A04/A07/A08/A10 source-contract and value hardening;
6. A05/A06/A12 reconciliation and provenance ledger;
7. A11 artifact staging/preflight behavior;
8. documentation relocation and link repair;
9. curated data layout, manifest, fixtures, and validation tooling;
10. approved publication storage configuration and data commit.

Stage only the intended paths for each commit and inspect
`git diff --cached --name-only`. Do not stage caches, credentials, unrelated
working files, or the entire repository by wildcard.

After explicit publication approval, push the release commits to `origin/main`,
then verify the remote branch at the exact pushed SHA. If licensing, privacy,
hosting, or checksum review fails, keep the data local and do not substitute a
silent upload of a smaller or different artifact set.

### Phase 8: rename the local checkout

After the final local validation and before or after remote publication as
convenient, rename only the local directory:

```text
D:\\W\\baktfold  ->  D:\\W\\enrich-bakta
```

This is not a Git history operation. Reopen the renamed checkout, verify the
remote URL and `main` tracking branch, search for stale absolute paths, rerun
the installed/out-of-tree CLI checks, and keep the move reversible until the
new path is confirmed.

## Acceptance criteria

The work is complete only when all of these are true:

1. The package layout reflects core, source adapters, workflows, and legacy
   code without changing supported command names.
2. Existing direct script paths and public Python symbols continue to work or
   have explicit, tested compatibility launchers.
3. `normalize_baktfold.py` remains available for historical use but is not part
   of the supported enrichment package or quality gate.
4. Astra's A01–A12 findings are either fixed with regression coverage or
   explicitly documented as a deliberate, tested limitation with changed
   claims; A01–A05 and the restoration policy cannot remain unresolved for a
   release claiming conservative transfer and complete provenance.
5. Published data has checksums, lineage, licensing/redistribution status,
   sensitivity review, and a clean-clone validation path.
6. Generated machine caches are reproducible and ignored; curated fixtures and
   real datasets are distinguishable.
7. The local checkout is named `enrich-bakta`, tracked on `main`, and the
   eventual remote verification records the exact published commit.

## Rollback and safety

- Do not delete the legacy utility, audit records, or local data during the
  migration.
- Use Git moves for tracked files and checksum-verified moves for ignored data.
- Keep each phase in a separate commit so package migration, hardening, data
  publication, and path rename can be reverted independently.
- If a compatibility or audit regression appears, stop at the last passing
  commit rather than pushing a partially reorganized tree.
- No remote push occurs during the drafting stage.
