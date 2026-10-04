# enrich-bakta: post-push review and follow-through hardening

Prepared for Wahyu and Luna, 2026-10-04 (Asia/Jakarta).

**Reviewed main:** `cc98e4f081c2e0554a7108b3897f6679fe61cd1d`, package `0.3.0`.
The final remote recheck still resolved to this commit. The comparison baseline
was `5d6754d306e80cddcad7da5fc12b155ab44d2e85` from the previous handoff.

**Disposition:** the migration and semantic foundation are substantially improved,
but this main still needs follow-through before its scientific validation and
provenance claims are treated as closed. The attached patch is a tested local
implementation of the concrete fixes below. It is not a release or remote commit.
No remote files, branches, tags, issues, or datasets were changed.

## 1. Evidence and review scope

I inspected the five new commits, canonical modules, compatibility shims, tests,
README, CI, pinned schema fixtures, prior handoff and September 29 audit, packaging,
and dataset governance scripts. I independently reproduced behavior rather than
accepting the old findings or the new acceptance record as proof.

The curated Git LFS corpus was available this time. All 33 published artifacts
materialized and passed size/hash/lineage validation. The 77-item local inventory
is intentionally broader than the public corpus; its absent local files are not
misclassified as missing public artifacts.

| Validation | Observed result |
|---|---|
| Unmodified main suite | **85 passed** |
| Main Ruff, formatting, mypy | Passed; mypy checks 14 supported modules |
| First targeted hardening probe set on main | **22 failed**, each asserting the corrected behavior |
| Previous 19 probes on main | 16 passed, 3 failed; one failure is the intended reject-versus-skip policy and another uses superseded metadata keys, so this is not three new independent defects |
| Tested patch, including 29 new regression cases | **114 passed** |
| Patch applied to a fresh detached checkout | Applied cleanly; **114 passed**, Ruff, formatting, mypy and whitespace checks passed |
| Declared Biopython floor 1.83 | **114 passed** on Python 3.12 |
| Original build with permitted setuptools 69.5.1 | Failed at SPDX `project.license` validation |
| Revised build with setuptools 77.0.3 | Wheel and sdist built successfully |
| Final installed wheel outside checkout | Six legacy module aliases, all five commands at help/version, both module entry points, and `pip check` passed |
| C14/SM published TSV versus XLSX on main | Both XLSX files rejected: their first column is `query` |
| C14/SM TSV versus XLSX with patch | Normalized supported hit fields match exactly: **4,563 / 4,497 rows** |
| C14/SM Baktfold + Kofam on main | Succeeded and reran byte-identically with zero insertions; detailed existing-KO evidence was incomplete |
| C14/SM restoration + eggNOG with patch | 4 / 6 imports retained in the new manifests; exact accepted eggNOG GBFF hashes reproduced |
| C14/SM complete unified workflow with patch | Passed, including parent-ledger verification and all three sources |
| Patched restoration, eggNOG and unified reruns | Zero insertions, identical GBFF bytes, retained imported origins |
| BK71A | Structural parse only: 2 records, 8,640 features, 4,131 CDSs; pristine FAA/evidence absent |
| Windows execution of this new patch / GitHub CI execution | Not performed locally; the matrix is an execution gate for Luna |

The runtime and compact real-run results are included in the ZIP. Tests and data
checks validate the merge behavior; they do not independently validate upstream
enzyme specificity, Foldseek assignments, eggNOG donors, expression, or phenotype.

What is closed from the old plan: package/launcher collision, legacy exclusion,
semantic locations and free-text qualifier handling, raw CDS identity cardinality,
captured XLSX bytes, destination staging, and explicit schema selection have real
implementations and regression coverage. Keep those improvements. The remaining
problems sit mainly in translation semantics and the evidence lifecycle.

## 2. Confirmed findings and exact required behavior

### R01 — High: `validated-only` does not validate a complete CDS model

**Main location:** `workflows/restore_translations.py:51–79`.

The helper uses `coding.translate(table=...)` with `cds=False`, then compares the
generic codon translation to FAA. This has both false negatives and false positives:

- A complete table-11 `GTG AAA TAA` CDS with the correct FAA `MK` is rejected,
  because the helper translates the initiator to V.
- `ATG AAA AAA` with FAA `MKK` passes despite having no terminal stop.
- `AAA AAA TAA` with FAA `KK` passes despite its invalid complete-CDS initiator.
- `<1..9`, an explicit translation exception, and a CDS extending to `1..12`
  on a 9-bp record can all receive a `passed` genomic status.
- `/transl_table=999` escapes as `KeyError` instead of a contextual merge error.

**Tested fix:** restrict the implemented validator to complete, exact, local,
in-bounds CDS models with a coherent strand and `codon_start=1`. Use
`Seq.translate(..., cds=True)`, retaining full joined/complemented extraction.
Reject unsupported partial/offset/remote/exceptional models as unresolved under
`validated-only`. Catch expected translation errors and unknown genetic codes.

This deliberately conservative model does not declare a disrupted pseudogene
biologically impossible. An unsupported conceptual model remains unresolved; an
explicit `import-faa` operation keeps imported origin and its actual failed or
unresolved genomic diagnostic. Do not broaden the strict model by flattening
compound locations or treating fuzzy positions as exact integers.

Biopython's documented `cds=True` behavior checks the initiator, length, and
terminal stop, and handles valid alternative initiators as M. See the primary
reference in section 8.

### R02 — High: origin ledgers stop after one enrichment, and restoration reruns lose them

**Main locations:** standalone Kofam/eggNOG wrappers, `workflows/enrich.py`,
restoration's `if translations: continue`, and central finalization.

All three FAA-dependent enrichment paths load and enforce a parent restoration
manifest, but neither its structured entries nor its metadata ledger are placed
in the next manifest. The resulting GBFF still contains the restoration marker.
Using its new manifest for the following run fails with
`translation evidence manifest contains no restoration entries`; the original
manifest cannot substitute because its output hash describes the previous GBFF.

Restoration on an already-restored GBFF currently accepts the existing FAA equality
without requiring the parent manifest and writes an empty origin ledger. Passing
`validated-only` on that rerun must never silently promote an imported origin.

**Tested fix:** carry every verified origin entry, including proteins untouched by
the current source, through every subsequent JSON manifest and bind it to the new
output hash. Record the parent manifest's captured-byte SHA-256. Restoration
reruns require the parent, emit explicit no-op diagnostics, and retain each
protein's original origin and original import/validation policy. Baktfold-only
and Baktfold-containing workflows accept and preserve the same lineage, with an
explicit opt-in for imported proteins. Finalization refuses a marked base whose
origin ledger has been dropped.

The parent hash must describe the bytes actually loaded, not a second disk read.
The patch captures this in `TranslationEvidenceLedger.manifest_sha256`; a
regression replaces the file after loading and verifies that the original
snapshot hash is recorded.

### R03 — High hardening: malformed ledger origin can bypass imported-protein opt-in

**Main location:** `core/merge_engine.py:679–798`.

The reader coerces fields with `str()` and `bool()`, accepts any origin string,
and checks duplicate metadata/entry ledgers only for equal query-key sets.
Changing `imported_faa` to a typo makes the later equality check miss the opt-in
guard. A changed origin inconsistent with `restoration_policy=import-faa` is also
accepted. A string boolean is not a boolean assertion.

**Tested fix:** validate supported origin/status enums, strict booleans, digest
syntax, origin/policy/status consistency, bound output, feature identity and
encoded protein hashes. Require the duplicate metadata and entry representations
to agree in their values. Check that all restored query IDs identified by the
in-file marker remain represented, rather than trusting an independently shortened
metadata dictionary. Preserve a contextual `MergeError` for malformed values.

These are integrity and schema checks, not authentication. Someone deliberately
fabricating a mutually consistent history can fabricate its hashes too. A prior
`genomically_validated` record is an audited producer assertion; hashes do not
create independent biological truth. Revalidation or signed lineage would be a
separate feature, not an implication of this patch.

### R04 — Medium, real-data blocker: published XLSX headers no longer parse

**Main location:** `sources/eggnog.py:398–480`.

Both published C14 and SM workbooks have the first header `query`; the shared
hit builder requires `#query`. The files are not damaged or formula-driven.
The advertised optional XLSX route fails before TSV/XLSX equivalence can be
checked. The new synthetic XLSX fixture uses `#query` and misses this real export.

**Tested fix:** recognize precisely the XLSX `query` alias and normalize it to
`#query` before duplicate-column validation. A workbook containing both spellings
is ambiguous and must fail. Keep TSV's explicit `#query` header contract, captured
`BytesIO` parsing, formula rejection, and workbook closure unchanged.

### R05 — Medium: reconciliation removes substantive Kofam evidence for existing KOs

**Main locations:** `sources/kofam.py:234–276`, `reconcile_insertions()`.

The per-hit score/threshold/E-value note is tagged `KofamScan provenance`.
The reconciliation rule suppresses provenance unless a new non-provenance
annotation survived. A valid hit whose KO already exists can therefore lose
its detailed note and inference; the surviving candidate is also mislabeled
`suppressed_conflict`, although no naming conflict occurred.

**Tested fix:** model hit notes as substantive `KofamScan hit evidence`, distinct
from producer inference. An existing KO decision is `supported_existing`, checked
against affirmative base evidence. Hit notes remain eligible with or without
`--no-feature-provenance`, and can support the optional producer inference.
Reruns retain truthful support even when no new qualifiers are added.

On controlled runs with unspecified Kofam version, the patch changes C14
insertions from **7,724 → 10,522** and SM from **7,678 → 10,380**. These are
additional evidence-note/inference insertions, not thousands of new KO calls.
The patched reruns still have zero insertions and identical GBFF bytes. These
version-unspecified hashes are intentionally not compared with acceptance runs
that may have supplied a version string; output provenance depends on that flag.

### R06 — Medium: candidate decisions still misdescribe pairing and surviving values

**Main location:** `core/decisions.py`, especially 145–444; eggNOG candidate rows.

Reproduced defects:

1. An eggNOG paired name is reconstructed as one CDS candidate and a separate
   generic gene candidate, instead of one decision with both targets.
2. Two sources supporting one value have one other supporter each, but
   `len(supporting) > 1` requires at least three sources before assigning
   `shared_support`.
3. `COG_category` is not mapped to its emitted `note` qualifier. COG source
   rows can appear suppressed while a generic fallback candidate covers the
   actual insertion. Stricter existing-value validation exposed this on real data.
4. An eggNOG name blocked by a different authoritative name becomes
   `supported_existing` for the proposed name, even though that proposed value
   is absent from the output.
5. The detailed source `entries` retain `status=emitted` and rejected qualifier
   strings after conflict suppression, even when `decisions` says suppressed.
6. Validation checks each referenced insertion but does not require support on
   every declared target of a paired decision.

**Tested fix:** supply explicit exact target UIDs from the adapter, scope fallback
grouping to the locus and feature role, map COG to note, preserve both targets of
a naming decision, recognize any other shared supporter, and distinguish blocked
authoritative naming from same-value existing support. Preserve the planned
projection separately while updating source rows with final status. Validate every
accepted target and every `supported_existing` value. Cache the feature-UID map
once; rebuilding it per candidate is prohibitively expensive on real genomes.

The patch also separates functional and provenance candidate/insertion counts.
Existing total fields remain totals; consumers should not reinterpret them as
functional transfer counts.

### R07 — Medium: Baktfold gene provenance survives a suppressed gene if an EC survives

**Main location:** `reconcile_insertions()` provenance pruning.

For Baktfold name `xyz` versus eggNOG name `abc`, default skip suppresses both
names. If the Baktfold CDS also contributes an independent EC, the generic
same-locus/source survival test retains `Baktfold gene-symbol evidence:v...`.
The generic provenance-support validator also accepts this unrelated EC as
functional support for the name-evidence note.

**Tested fix:** retain Baktfold class-specific provenance only when the same
target retains evidence of that class: gene → gene-symbol note, EC → functional
note, structural xref → structure inference. Keep independent permitted EC/xref
transfers. The regression verifies both the absence of the rejected name note
and the presence of the independent EC.

### R08 — Medium: broad producer-comment recognition still hides reserved query rows

**Main location:** `sources/eggnog.py:283–348`.

`_is_recognized_producer_comment()` accepts anything beginning with a weekday's
first three characters. A full data row with query `##Monday_bad` is silently
skipped. Other data-shaped double-hash strings can overlap metadata prefixes.

**Tested fix:** reject a header-width tabular reserved-ID row before comment
classification, and full-match the documented timestamp shape. Preserve the
existing valid tabbed command-metadata fixture. A future producer with a new
comment format must get an explicit fixture/contract update rather than another
broad prefix exception.

Suffix cleanup also still turns `NA_123` into accepted gene `NA`. The patch
rechecks missing sentinels after normalization, while preserving the pre-cleaning
raw value in evidence rows.

### R09 — Medium: build-system minimum contradicts PEP 639 metadata

**Main location:** `pyproject.toml:2`.

`setuptools>=69` permits a backend that cannot parse the new SPDX string and
`project.license-files`. A controlled no-isolation build with 69.5.1 fails.
This is not visible when isolated builds happen to resolve the newest backend.

**Tested fix:** set the build minimum to `setuptools>=77.0.3`; build explicitly
with that floor and add the Biopython/backend floor job. PyPA's primary packaging
guide identifies 77.0.3 for this metadata support. The existing OS matrix still
needs to run the new patch.

## 3. What the supplied patch does not complete

The concrete regressions above are fixed and tested, but the original handoff's
goal of a fully adapter-native decision graph is not fully implemented by
`build_candidate_ledger()`. It still reconstructs candidates from insertion plans
and selected source rows. Be explicit about this remaining boundary.

Before calling the complete A05/A12 contract closed, finish these steps in a
separate commit after the tested patch:

1. **Exhaustive Baktfold outcomes.** Emit a source candidate before each decision,
   including already-present values, source/base name conflicts, unsupported
   pairs, and protein-mismatch suppression. Current Baktfold no-op runs can have
   no per-value decisions at all; metadata summaries are not an exhaustive ledger.
2. **Class-specific provenance support links.** Associate each provenance node
   with the exact accepted candidate IDs for its evidence class. The tested
   pruning prevents the reproduced false note, but generic support lists should
   not be presented as a fully normalized evidence graph.
3. **Explicit role/count schema.** Distinguish functional proposal, substantive
   hit evidence, producer provenance and context. Preserve the new separate counts;
   add schema tests rather than relying on strings containing `provenance` forever.
4. **Final semantic-output validation.** The current engine verifies all inserted
   semantic values and the new ledger checks declared target support. Extend the
   ledger's final check to the parsed output directly for existing/supporting
   relationships, and exercise all final status enums, including invalid states.
5. **Reconciliation reason per node.** A dropped candidate must explain whether
   it lost to an explicit source preference, lacked a supported pair, or failed
   protein identity. Avoid reconstructing reasons from display qualifier strings.
6. **Typed manifest schema and CI.** Add a checked-in manifest/decision schema,
   assert actual JSON against it, and give dataset byte validation and scientific
   reruns separate gates. Current `dataset-manifest` checks integrity, not biological
   output equivalence. Full clean-clone data gates require `git lfs pull`.

The README part of the tested patch narrows the exhaustive-ledger claim to match
the implementation. If full exhaustive reporting is a release requirement,
this follow-up remains a release gate, not an optional cleanup task.

The complete draft diff below has **no ellipses or pseudo-hunks** and applies to
the exact reviewed commit. The next design steps above are implementation
requirements; they are not represented as completed code in that diff.

## 4. Scientific and compatibility decisions for Luna

- Keep imported FAA equality separate from independent nucleotide-model validation.
  An import can have failed/unresolved genomic diagnostics and still be imported
  under explicit policy; do not transform it into `passed` during reruns.
- Keep `/pseudo` and `/pseudogene` markers. A restored research artifact is not an
  INSDC submission-validation claim or evidence of intact enzymatic capacity.
- Do not auto-promote the preliminary EC strings merely to make C14/SM run without
  a flag. However, ExPASy explicitly recognizes `3.5.1.n3`, `3.6.5.n1`, and
  `4.2.2.n1`. Their existing `invalid_value` label means outside this adapter's
  conservative promotion grammar, not malformed in ENZYME. A future schema should
  separate lexical/database class (`provisional`) from promotion decision.
- Do not fabricate a migration ledger for historical restored outputs whose
  original FAA and restoration provenance are absent. Recompute C14/SM from the
  pristine inputs and retain the historical files as historical. BK71A remains
  structural-only until its original evidence is supplied.
- The stricter validated-only model intentionally rejects partial/exceptional
  cases until a tested conceptual-translation model exists. The opt-in import path
  remains available and explicitly labeled.
- Source-entry `status` now reflects the final decision; the old value is retained
  as `planned_status`. Consumers that treated `entries.status` as a planner outcome
  need to switch to the explicit planned field.
- `KofamScan hit evidence` is a new insertion source label. Consumers should use
  the recorded source ID/role and qualifier, rather than exact-match the old
  display label for substantive notes.
- No version bump or tag is included. Choose the next version after review and
  synchronize `pyproject.toml`, `TOOL_VERSION`, and release notes in one commit.
- Multi-file staging remains sequential promotion. A later promotion failure can
  leave earlier artifacts published. Neither the existing nor this patch promises
  cross-file crash atomicity, concurrent-writer exclusion, or durable directory
  fsync on every platform.

## 5. Execution sequence and gates

Use a branch/worktree and preserve unrelated local changes. Fetch, verify the
baseline, then apply the single tested patch. If main has advanced, inspect the
new diff and rebase intentionally; do not apply with silent rejects.

```bash
git fetch origin main
git switch -c hardening/post-push-followthrough origin/main
git rev-parse HEAD
# Expected: cc98e4f081c2e0554a7108b3897f6679fe61cd1d
git apply --check /absolute/path/to/followthrough-hardening.patch
git apply /absolute/path/to/followthrough-hardening.patch
python -m pip install -e '.[dev,xlsx]'
python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m mypy
git diff --check
python -m build
git lfs pull
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
```

Run the bundled `run_real_checks.py` followed by `run_chain_checks.py` from the
repository root in that environment. They use the published source paths and
write derived review outputs to a sibling `followthrough-evidence/real-runs/`
directory. Do not stage generated genomic outputs or review sidecars by accident.
The scripts emit compact JSON outcomes; inspect every row as well as the exit code.

Also run a no-isolation build with setuptools **77.0.3** and the complete suite
with Biopython **1.83**, then install the wheel in a clean environment and run its
commands from a directory outside the checkout. The added CI floor job covers the
same principle. Run Windows/Linux CI after the patch and inspect each job result;
adding YAML is not evidence of a Windows pass.

Real-chain gate example (C14; substitute matching SM names):

```bash
python restore_bakta_translations.py \
  data/C14/bakta/C14-NMZ.gbff data/C14/bakta/C14-NMZ.faa \
  data/C14/evidence/C14-NMZ-query.emapper.annotations \
  .test-output/followthrough/C14-restored.gbff \
  --translation-policy import-faa \
  --manifest .test-output/followthrough/C14-restored.json

python enrich_bakta.py \
  --bakta .test-output/followthrough/C14-restored.gbff \
  --faa data/C14/bakta/C14-NMZ.faa \
  --baktfold data/C14/baktfold/C14_NMZ.gbff \
  --kofamscan data/C14/evidence/C14-KofamKOALA.txt \
  --kofamscan-version 1.3.0 \
  --eggnog data/C14/evidence/C14-NMZ-query.emapper.annotations \
  --translation-evidence-manifest .test-output/followthrough/C14-restored.json \
  --allow-imported-translations --baktfold-invalid-ec-policy skip \
  --output .test-output/followthrough/C14-unified.gbff \
  --manifest .test-output/followthrough/C14-unified.json
```

For the rerun use the **unified GBFF and its unified manifest** as the parent,
the same source evidence and policy flags, and new output paths. Require zero
insertions, exact GBFF byte identity, and four imported origin entries. SM must
retain six. Removing imported opt-in must block their FAA-dependent use.

Update `docs/data/ACCEPTANCE.md` with the exact new commit, environment, commands,
flags, hashes, insertion/evidence counts and rerun results. Its current duplicate-
run statement checks identical pristine inputs; explicitly record the additional
output-as-input idempotence chains. Keep the old acceptance record clearly dated
as historical evidence instead of silently rewriting its observations.

## 6. Controlled patched real-run results

The Baktfold+Kofam runs below omit a Kofam version flag. The unified runs explicitly
use `1.3.0`. No timestamp is supplied. These exact flags matter for hash comparison.

| Workflow | C14 | SM |
|---|---|---|
| TSV/XLSX matching normalized hits | 4,563 | 4,497 |
| Imported restoration targets | 4 | 6 |
| Baktfold+Kofam insertions | 10,522 | 10,380 |
| eggNOG qualifier candidate emissions | 29,392 | 29,086 |
| eggNOG total insertions | 33,944 | 33,567 |
| Unified total insertions | 41,500 | 41,013 |
| Each checked rerun | 0 insertions; exact byte identity | 0 insertions; exact byte identity |

```text
C14 Baktfold+Kofam: 7a2a8482346b60503cf6f03a0a654294b0ee23cc4af0ec22f8e0016295907862
SM  Baktfold+Kofam: 63bb4e91beb0035d0b18f0a4d5538429d317940a6fe6d9c4fa78f9348976360b
C14 eggNOG: 3a6c4cebcfb497a9f0c0b754ae677effa292e8af909a4249a1abea675f3affd7
SM  eggNOG: 6513327169594fb9c09b18a5ced830aca66f06ee3baa6beae60eb1fdd955878b
C14 unified: 90cdfc6897410cdc20fd4e3f83aa2b316463306831e2b30ca15566568458c117
SM  unified: ad3408c57bd5982c4e0bc9c638f290da08c7867d84f13c155ffbfe608c9fb1b6
```

The eggNOG hashes agree with the earlier acceptance record: this patch corrects
ledger semantics and lineage without changing those standalone GBFF annotations.
Kofam evidence retention intentionally changes its output; do not freeze the
historically incomplete hash as a compatibility gate.

## 7. Complete tested diff

This block is the same patch bytes as `followthrough-hardening.patch` in the ZIP.
Apply it once, as one coordinated integration unit. It changes the canonical
implementation and preserves the compatibility shims. The standalone test module
is included; it imports the existing repository fixture helpers.

```diff
diff --git a/.github/workflows/ci.yml b/.github/workflows/ci.yml
index cec29fa..517200b 100644
--- a/.github/workflows/ci.yml
+++ b/.github/workflows/ci.yml
@@ -37,6 +37,18 @@ jobs:
       - run: python -m pip install -e ".[dev,xlsx]"
       - run: python -c "import openpyxl; print(openpyxl.__version__)"
       - run: python -m pytest -q tests/test_hardening_regressions.py -k xlsx
+      - run: python -m pytest -q tests/test_followthrough.py -k xlsx
+
+  dependency-floor:
+    runs-on: ubuntu-latest
+    steps:
+      - uses: actions/checkout@v4
+      - uses: actions/setup-python@v5
+        with:
+          python-version: "3.12"
+      - run: python -m pip install -e ".[dev,xlsx]" "biopython==1.83" "setuptools==77.0.3" wheel
+      - run: python -m pytest -q
+      - run: python -m build --no-isolation
 
   package:
     runs-on: ${{ matrix.os }}
diff --git a/README.md b/README.md
index dcf92b2..80d66d8 100644
--- a/README.md
+++ b/README.md
@@ -57,8 +57,11 @@ Useful options:
 - `--no-feature-provenance` omits per-feature provenance.
 - `--merge-timestamp VALUE` adds an explicitly supplied timestamp. No timestamp
   is generated by default, so identical inputs produce identical output.
-- `--baktfold-invalid-ec-policy skip` records legacy non-grammar Baktfold EC
-  tokens as non-promoted `invalid_value` decisions. The default is `reject`.
+- `--baktfold-invalid-ec-policy skip` records tokens outside this adapter's
+  promotion grammar as non-promoted `invalid_value` decisions. The default is
+  `reject`. Preliminary ENZYME identifiers such as `3.5.1.n3` are recognized
+  database identifiers; this legacy status describes the adapter's promotion
+  policy, not their validity in ENZYME.
 - `--manifest PATH` writes TSV unless the suffix is `.json`.
 
 ## KofamScan merge
@@ -204,7 +207,10 @@ The restoration script refuses any query absent from the FAA or GBFF and any
 translationless CDS not marked `/pseudo` or `/pseudogene`. Requested CDSs that already have a
 translation must still match the FAA exactly; they are validated but not
 modified. `--translation-policy validated-only` requires an independent genomic
-translation match. `--translation-policy import-faa` records imported origin,
+translation match for a complete, local CDS with exact bounds, a supported
+genetic code, a valid initiator and terminal stop. Fuzzy, offset, remote and
+exceptional translation models remain unresolved under this policy.
+`--translation-policy import-faa` records imported origin,
 and every downstream functional merge requires both
 `--translation-evidence-manifest` and `--allow-imported-translations`. The
 script inserts only the matched FAA sequence and records the locus, feature
@@ -239,6 +245,24 @@ stored in the typed `decisions` section, while `entries` remains the detailed
 source/insertion projection. TSV is a display projection and is not accepted
 as a translation-evidence ledger.
 
+Every enrichment workflow carries a verified restoration ledger into its JSON
+manifest and binds it to the new output hash. On the next run, supply that most
+recent manifest with `--translation-evidence-manifest`. Imported proteins still
+require `--allow-imported-translations`. Baktfold supports these options too.
+Restoration reruns require the parent manifest and preserve each protein's
+original origin, even when no new translation is inserted. A TSV manifest or
+an output written without a JSON manifest cannot carry this machine-readable
+lineage to the next run.
+
+Candidate decisions distinguish source proposals from final results. A paired
+gene name has both exact targets; shared values retain all source support.
+`suppressed_authoritative_gene` records a proposal blocked by a different or
+incomplete authoritative name. Kofam hit notes are substantive profile evidence
+and remain eligible when the KO already exists, including with
+`--no-feature-provenance`. Counts distinguish functional and provenance
+insertions. Baktfold's existing/no-op and rejected proposals are not yet an
+exhaustive per-value decision ledger; their current summaries remain in metadata.
+
 These xrefs follow the project's enrichment convention. Formal submission needs
 the current [INSDC db_xref controlled vocabulary](https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/);
 Biopython parsing is not an INSDC submission-validator result. eggNOG evidence is
diff --git a/pyproject.toml b/pyproject.toml
index ba238e3..cc2800c 100644
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -1,5 +1,5 @@
 [build-system]
-requires = ["setuptools>=69"]
+requires = ["setuptools>=77.0.3"]
 build-backend = "setuptools.build_meta"
 
 [project]
diff --git a/src/enrich_bakta_lib/core/decisions.py b/src/enrich_bakta_lib/core/decisions.py
index 27e8a9c..f71e7e2 100644
--- a/src/enrich_bakta_lib/core/decisions.py
+++ b/src/enrich_bakta_lib/core/decisions.py
@@ -137,6 +137,7 @@ def _row_field_and_value(row: Mapping[str, Any]) -> tuple[str, str, str]:
         "CAZy": "db_xref",
         "PFAMs": "note",
         "eggNOG_OGs": "note",
+        "COG_category": "note",
     }.get(field, str(row.get("qualifier", field)))
     value = str(row.get("normalized_value", row.get("value", "")))
     return field, qualifier, value
@@ -251,12 +252,18 @@ def build_candidate_ledger(
                 for insertion_id in final_by_semantic.get(_semantic_key(item), [])
             }
         )
-        if emitted:
-            final_status = "emitted" if final_ids else "suppressed_conflict"
-        elif row_status in {"existing", "existing_gene"}:
+        if row_status in {"existing", "existing_gene"}:
             final_status = "supported_existing"
+        elif emitted or matches:
+            final_status = "emitted" if final_ids else "suppressed_conflict"
         else:
             final_status = row_status or "not_applicable"
+        row["planned_status"] = planned_status
+        row["planned_emitted_qualifiers"] = row.get("emitted_qualifiers", "")
+        row["final_status"] = final_status
+        row["status"] = final_status
+        if final_status not in {"emitted", "shared_support"}:
+            row["emitted_qualifiers"] = ""
         payload = {
             "entry_type": entry_type,
             "source_sha256": source_sha,
@@ -294,8 +301,8 @@ def build_candidate_ledger(
             )
         )
 
-    grouped: dict[tuple[str, str, str, str, str, str], list[Insertion]] = defaultdict(
-        list
+    grouped: dict[tuple[str, str, str, str, str, str, str, str], list[Insertion]] = (
+        defaultdict(list)
     )
     for item in planned:
         if id(item) in covered:
@@ -306,6 +313,8 @@ def build_candidate_ledger(
                 source_id,
                 _source_hash(source_id, source_hashes),
                 item.record,
+                item.locus_tag,
+                "paired" if item.qualifier == "gene" else item.feature_type,
                 item.qualifier,
                 item.value,
                 item.source_value,
@@ -315,6 +324,8 @@ def build_candidate_ledger(
         source_id,
         source_sha,
         record,
+        locus_tag,
+        feature_scope,
         qualifier,
         value,
         raw_value,
@@ -336,6 +347,8 @@ def build_candidate_ledger(
         payload = {
             "source_sha256": source_sha,
             "record": record,
+            "locus_tag": locus_tag,
+            "feature_scope": feature_scope,
             "field": qualifier,
             "raw_value": raw_value,
             "normalized_value": value,
@@ -388,7 +401,7 @@ def build_candidate_ledger(
             if other.candidate_id != candidate.candidate_id
         }
         candidate.supporting_candidate_ids = tuple(sorted(supporting))
-        if len(supporting) > 1 and candidate.final_status == "emitted":
+        if supporting and candidate.final_status == "emitted":
             candidate.final_status = "shared_support"
 
     functional_support: dict[tuple[str, str], set[str]] = defaultdict(set)
@@ -438,6 +451,20 @@ def build_candidate_ledger(
             }
         ),
         "candidate_ledger_schema": "enrich-bakta.candidate-ledger.v2",
+        "functional_candidate_count": sum(
+            not candidate.confidence_diagnostics.get("provenance", False)
+            for candidate in candidates
+        ),
+        "provenance_candidate_count": sum(
+            bool(candidate.confidence_diagnostics.get("provenance"))
+            for candidate in candidates
+        ),
+        "functional_insertion_count": sum(
+            not _is_provenance_insertion(item) for item in final
+        ),
+        "provenance_insertion_count": sum(
+            _is_provenance_insertion(item) for item in final
+        ),
     }
     return rows, counts
 
@@ -454,6 +481,11 @@ def validate_candidate_ledger(
     if len(insertion_by_id) != len(final_items):
         raise MergeError("candidate ledger validation found duplicate insertion IDs")
     feature_index = _feature_index(base) if base is not None else {}
+    feature_by_uid = (
+        {feature_uid(feature): feature for feature in base.features}
+        if base is not None
+        else {}
+    )
     valid_feature_uids = set(feature_index.values()) if base is not None else None
     ledger_rows = list(rows)
     seen: set[str] = set()
@@ -530,6 +562,42 @@ def validate_candidate_ledger(
                         f"candidate {candidate_id!r} references the wrong value"
                     )
                 referenced_final_ids.add(reference)
+            supported_targets = {
+                _insertion_feature_uid(insertion_by_id[reference], feature_index)
+                for reference in referenced
+            }
+            if base is not None:
+                supported_targets.update(
+                    target
+                    for target in targets
+                    if normalized_value in feature_by_uid[target].values(qualifier)
+                )
+                if set(targets) != supported_targets:
+                    raise MergeError(
+                        f"candidate {candidate_id!r} lacks support on every target"
+                    )
+        if status == "supported_existing" and base is not None:
+            from enrich_bakta_lib.sources.value_rules import structured_note_tokens
+
+            for target in targets:
+                feature = feature_by_uid[target]
+                values = set(feature.values(qualifier))
+                if qualifier == "db_xref":
+                    values.update(
+                        token
+                        for note in feature.values("note")
+                        for token in structured_note_tokens(note)
+                    )
+                if qualifier == "EC_number":
+                    values.update(
+                        value.removeprefix("EC:")
+                        for value in feature.values("db_xref")
+                        if value.startswith("EC:")
+                    )
+                if normalized_value not in values:
+                    raise MergeError(
+                        f"supported-existing candidate {candidate_id!r} is absent on its target"
+                    )
         if status == "suppressed_conflict" and referenced:
             raise MergeError(
                 f"suppressed candidate {candidate_id!r} still references insertions"
diff --git a/src/enrich_bakta_lib/core/merge_engine.py b/src/enrich_bakta_lib/core/merge_engine.py
index d923a10..f514d2a 100644
--- a/src/enrich_bakta_lib/core/merge_engine.py
+++ b/src/enrich_bakta_lib/core/merge_engine.py
@@ -667,9 +667,11 @@ class TranslationEvidenceLedger(dict[str, TranslationEvidence]):
         values: Mapping[str, TranslationEvidence],
         *,
         required_query_ids: Iterable[str] = (),
+        manifest_sha256: str = "",
     ) -> None:
         super().__init__(values)
         self.required_query_ids = frozenset(required_query_ids)
+        self.manifest_sha256 = manifest_sha256
 
 
 def has_translation_evidence_marker(data: bytes) -> bool:
@@ -686,7 +688,8 @@ def load_translation_evidence(
             "cannot preserve the structured evidence ledger"
         )
     try:
-        payload = json.loads(path.read_text(encoding="utf-8"))
+        snapshot = path.read_bytes()
+        payload = json.loads(snapshot.decode("utf-8"))
     except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
         raise MergeError(
             f"cannot read translation evidence manifest {path}: {exc}"
@@ -741,6 +744,60 @@ def load_translation_evidence(
         )
         if any(key not in raw for key in required):
             raise MergeError(f"translation evidence for {query_id!r} is incomplete")
+        if (
+            type(raw["artifact_sequence_match"]) is not bool
+            or raw["artifact_sequence_match"] is not True
+        ):
+            raise MergeError(
+                f"translation evidence for {query_id!r} requires a true boolean artifact match"
+            )
+        if not isinstance(raw["origin"], str) or raw["origin"] not in {
+            "imported_faa",
+            "genomically_validated",
+        }:
+            raise MergeError(
+                f"translation evidence for {query_id!r} has an unknown origin"
+            )
+        if not isinstance(raw["genomic_validation_status"], str) or raw[
+            "genomic_validation_status"
+        ] not in {
+            "passed",
+            "failed",
+            "unresolved",
+            "not_attempted",
+        }:
+            raise MergeError(
+                f"translation evidence for {query_id!r} has an invalid genomic status"
+            )
+        expected_policy = (
+            "import-faa" if raw["origin"] == "imported_faa" else "validated-only"
+        )
+        if raw["restoration_policy"] != expected_policy or (
+            raw["origin"] == "genomically_validated"
+            and raw["genomic_validation_status"] != "passed"
+        ):
+            raise MergeError(
+                f"translation evidence for {query_id!r} has inconsistent origin/policy/status"
+            )
+        for key in (
+            "protein_sha256",
+            "original_base_sha256",
+            "faa_sha256",
+            "bound_output_sha256",
+        ):
+            if (
+                not isinstance(raw[key], str)
+                or re.fullmatch(r"[0-9a-f]{64}", raw[key]) is None
+            ):
+                raise MergeError(
+                    f"translation evidence for {query_id!r} has an invalid {key}"
+                )
+        if not isinstance(raw.get("producer_lineage"), dict) or not isinstance(
+            raw["validation_reason"], str
+        ):
+            raise MergeError(
+                f"translation evidence for {query_id!r} has invalid lineage/reason fields"
+            )
         try:
             evidence = TranslationEvidence(
                 feature_uid=str(raw["feature_uid"]),
@@ -775,6 +832,14 @@ def load_translation_evidence(
             raise MergeError(
                 f"translation evidence for {query_id!r} has the wrong feature identity"
             )
+        translations = cds[query_id].values("translation")
+        if (
+            len(translations) != 1
+            or protein_sha256(translations[0]) != evidence.protein_sha256
+        ):
+            raise MergeError(
+                f"translation evidence for {query_id!r} disagrees with the encoded protein"
+            )
         result[query_id] = evidence
     if not result:
         raise MergeError(
@@ -795,7 +860,40 @@ def load_translation_evidence(
             raise MergeError(
                 "translation evidence manifest metadata and entries disagree"
             )
-    return TranslationEvidenceLedger(result, required_query_ids=required_query_ids)
+        if any(
+            metadata_evidence[query_id] != evidence.as_dict()
+            for query_id, evidence in result.items()
+        ):
+            raise MergeError(
+                "translation evidence manifest metadata and entries disagree in their values"
+            )
+    if has_translation_evidence_marker(base.data):
+        marker_ids: set[str] = set()
+        for record in base.records:
+            prefix = base.data[record.start : record.features_offset].decode("utf-8")
+            content = "\n".join(line[12:] for line in prefix.splitlines())
+            for block in content.split(TRANSLATION_EVIDENCE_MARKER)[1:]:
+                match = re.search(
+                    r"Restored CDSs:\s*(.*?)\s*Original base SHA-256:", block, re.DOTALL
+                )
+                if match is None:
+                    raise MergeError(
+                        "translation evidence marker lacks a restored-CDS list"
+                    )
+                marker_ids.update(
+                    token.strip()
+                    for token in match.group(1).replace("\n", " ").split(",")
+                )
+        if not marker_ids or not marker_ids.issubset(result):
+            raise MergeError(
+                "translation evidence manifest omits targets identified by the GBFF marker"
+            )
+        required_query_ids.update(marker_ids)
+    return TranslationEvidenceLedger(
+        result,
+        required_query_ids=required_query_ids,
+        manifest_sha256=sha256_bytes(snapshot),
+    )
 
 
 def protein_sha256(sequence: str) -> str:
@@ -1112,7 +1210,7 @@ def reconcile_insertions(
         filtered.append(item)
 
     def source_family(source: str) -> str:
-        return source.lower().removesuffix(" provenance")
+        return source.lower().removesuffix(" provenance").removesuffix(" hit evidence")
 
     surviving_functional = {
         (
@@ -1125,6 +1223,23 @@ def reconcile_insertions(
     }
     provenance_filtered: list[Insertion] = []
     for item in filtered:
+        if item.source == "Baktfold provenance" and item.qualifier != "COMMENT":
+            evidence_qualifier = (
+                "gene"
+                if item.value.startswith("Baktfold gene-symbol evidence:")
+                else "EC_number"
+                if item.value.startswith("Baktfold functional-annotation evidence:")
+                else "db_xref"
+            )
+            if not any(
+                other.source == "Baktfold"
+                and other.record == item.record
+                and other.locus_tag == item.locus_tag
+                and other.feature_type == item.feature_type
+                and other.qualifier == evidence_qualifier
+                for other in filtered
+            ):
+                continue
         if (
             item.qualifier != "COMMENT"
             and _is_provenance_source(item.source)
@@ -1405,6 +1520,12 @@ def finalize_merge(
     sidecar_path: Path | None = None,
     sidecar_payload: dict[str, Any] | None = None,
 ) -> dict[str, Any]:
+    if has_translation_evidence_marker(base_data) and not metadata.get(
+        "translation_evidence"
+    ):
+        raise MergeError(
+            "restored input requires --translation-evidence-manifest and preserved translation evidence"
+        )
     inputs = [base_path, *other_inputs]
     if paths_collide(output_path, inputs):
         raise MergeError("output path must differ from every input path")
@@ -1462,6 +1583,23 @@ def finalize_merge(
         rows = insertion_rows(applied)
         if evidence_rows is not None:
             rows = [*_bind_translation_evidence(evidence_rows, output_sha256), *rows]
+        represented = {
+            row.get("query_id")
+            for row in rows
+            if row.get("entry_type") == "translation_restoration"
+        }
+        for query_id, evidence in sorted(
+            metadata.get("translation_evidence", {}).items()
+        ):
+            if query_id not in represented:
+                rows.append(
+                    {
+                        "entry_type": "translation_restoration",
+                        "query_id": query_id,
+                        "status": "carried_forward",
+                        "translation_evidence": evidence,
+                    }
+                )
         artifacts.append((manifest_path, manifest_bytes(manifest_path, metadata, rows)))
     if sidecar_path is not None and sidecar_payload is not None:
         artifacts.append((sidecar_path, json_sidecar_bytes(sidecar_payload)))
diff --git a/src/enrich_bakta_lib/sources/baktfold.py b/src/enrich_bakta_lib/sources/baktfold.py
index a3da773..78fad4c 100644
--- a/src/enrich_bakta_lib/sources/baktfold.py
+++ b/src/enrich_bakta_lib/sources/baktfold.py
@@ -14,6 +14,7 @@ import collections
 import json
 import re
 import sys
+from collections.abc import Mapping
 from pathlib import Path
 from typing import Any
 
@@ -26,14 +27,20 @@ from enrich_bakta_lib.core.merge_engine import (
     Insertion,
     MergeError,
     RawDocument,
+    TranslationEvidence,
+    cds_by_locus,
     comment_insertion,
     finalize_merge,
+    has_translation_evidence_marker,
     index_unique_features,
+    load_translation_evidence,
+    normalize_protein,
     parse_genbank_bytes,
     qualifier_insertion,
     read_input_bytes,
     sha256_bytes,
     strict_parity_check,
+    validate_faa_gbff,
     validate_genbank_semantics,
 )
 from enrich_bakta_lib.sources.value_rules import (
@@ -129,6 +136,8 @@ def plan_baktfold_additions(
     merge_timestamp: str | None = None,
     starting_order: int = 0,
     invalid_ec_policy: str = "reject",
+    translation_evidence: Mapping[str, TranslationEvidence] | None = None,
+    allow_imported_translations: bool = False,
 ) -> tuple[list[Insertion], dict[str, Any]]:
     """Validate Baktfold parity and return the complete insertion allowlist."""
     if invalid_ec_policy not in {"reject", "skip"}:
@@ -136,6 +145,24 @@ def plan_baktfold_additions(
             f"invalid Baktfold EC policy {invalid_ec_policy!r}; expected 'reject' or 'skip'"
         )
     parity = strict_parity_check(base, baktfold)
+    if has_translation_evidence_marker(base.data) and translation_evidence is None:
+        raise MergeError(
+            "restored input requires --translation-evidence-manifest for Baktfold planning"
+        )
+    if translation_evidence:
+        cds = cds_by_locus(base)
+        encoded_proteins = {
+            query: normalize_protein(cds[query].values("translation")[0])
+            for query in translation_evidence
+        }
+        validate_faa_gbff(
+            base,
+            encoded_proteins,
+            translation_evidence,
+            source_name="Baktfold",
+            translation_evidence=translation_evidence,
+            allow_imported_translations=allow_imported_translations,
+        )
     source_index = index_unique_features(baktfold.features, label="Baktfold")
     # Also reject duplicates in the base before matching.
     index_unique_features(base.features, label="Bakta base")
@@ -349,12 +376,19 @@ def graft(
     add_feature_provenance: bool = True,
     merge_timestamp: str | None = None,
     invalid_ec_policy: str = "reject",
+    translation_evidence_manifest: Path | None = None,
+    allow_imported_translations: bool = False,
 ) -> dict[str, Any]:
     bakta_data = read_input_bytes(bakta_path, "Bakta")
     baktfold_data = read_input_bytes(baktfold_path, "Baktfold")
     validate_genbank_semantics(bakta_data, "Bakta input")
     validate_genbank_semantics(baktfold_data, "Baktfold input")
     base = parse_genbank_bytes(bakta_data, "Bakta input")
+    translation_evidence = (
+        load_translation_evidence(translation_evidence_manifest, base)
+        if translation_evidence_manifest
+        else None
+    )
     source = parse_genbank_bytes(baktfold_data, "Baktfold input")
     insertions, stats = plan_baktfold_additions(
         base,
@@ -364,6 +398,8 @@ def graft(
         add_comment_note=add_comment_note,
         merge_timestamp=merge_timestamp,
         invalid_ec_policy=invalid_ec_policy,
+        translation_evidence=translation_evidence,
+        allow_imported_translations=allow_imported_translations,
     )
     candidate_rows, candidate_counts = build_candidate_ledger(
         base,
@@ -378,7 +414,10 @@ def graft(
         base_path=bakta_path,
         base_data=bakta_data,
         output_path=output_path,
-        other_inputs=[baktfold_path],
+        other_inputs=[
+            baktfold_path,
+            *([translation_evidence_manifest] if translation_evidence_manifest else []),
+        ],
         insertions=insertions,
         manifest_path=manifest_path,
         evidence_rows=candidate_rows,
@@ -402,6 +441,19 @@ def graft(
             "invalid_ec_values": stats["invalid_ec_values"],
             "invalid_ec_value_count": stats["invalid_ec_value_count"],
             **candidate_counts,
+            "translation_evidence": {
+                query: evidence.as_dict()
+                for query, evidence in (translation_evidence or {}).items()
+            },
+            "translation_evidence_parent_sha256": getattr(
+                translation_evidence, "manifest_sha256", ""
+            ),
+            "allow_imported_translations": allow_imported_translations,
+            "policies": {
+                "add_comment_note": add_comment_note,
+                "add_feature_provenance": add_feature_provenance,
+            },
+            "parity": stats["parity"],
         },
     )
     return {**stats, **final}
@@ -426,6 +478,8 @@ def main() -> int:
     parser.add_argument("output", type=Path, help="new merged .gbff")
     parser.add_argument("--manifest", type=Path, help="TSV or .json insertion manifest")
     parser.add_argument("--no-comment-note", action="store_true")
+    parser.add_argument("--translation-evidence-manifest", type=Path)
+    parser.add_argument("--allow-imported-translations", action="store_true")
     parser.add_argument(
         "--no-feature-provenance",
         "--no-inference-provenance",
@@ -462,6 +516,8 @@ def main() -> int:
             add_feature_provenance=not args.no_feature_provenance,
             merge_timestamp=args.merge_timestamp,
             invalid_ec_policy=args.baktfold_invalid_ec_policy,
+            translation_evidence_manifest=args.translation_evidence_manifest,
+            allow_imported_translations=args.allow_imported_translations,
         )
     except MergeError as exc:
         print(f"ERROR: {exc}", file=sys.stderr)
diff --git a/src/enrich_bakta_lib/sources/eggnog.py b/src/enrich_bakta_lib/sources/eggnog.py
index 50f8f82..bd26456 100644
--- a/src/enrich_bakta_lib/sources/eggnog.py
+++ b/src/enrich_bakta_lib/sources/eggnog.py
@@ -287,7 +287,13 @@ def _is_recognized_producer_comment(line: str) -> bool:
         return True
     if body.startswith(("/", "applied filters:")):
         return True
-    return body[:3] in {"Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"}
+    return (
+        re.fullmatch(
+            r"(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) [A-Z][a-z]{2}\s+\d{1,2} \d{2}:\d{2}:\d{2} \d{4}",
+            body,
+        )
+        is not None
+    )
 
 
 def parse_eggnog_tsv(
@@ -309,6 +315,10 @@ def parse_eggnog_tsv(
         if not line:
             continue
         if line.startswith("##"):
+            if columns is not None and len(line.split("\t")) == len(columns):
+                raise MergeError(
+                    f"eggNOG row {row_number}: data row begins with reserved '#'"
+                )
             match = VERSION_RE.match(line)
             if match:
                 if hits:
@@ -431,6 +441,8 @@ def parse_eggnog_xlsx(
         columns = [
             str(cell.value) if cell.value is not None else "" for cell in header_cells
         ]
+        # Official XLSX exports use 'query'; TSV uses '#query'.
+        columns = ["#query" if column == "query" else column for column in columns]
         if len(columns) != len(set(columns)):
             raise MergeError("eggNOG XLSX header contains duplicate columns")
         if any(not column.strip() for column in columns):
@@ -867,7 +879,16 @@ def plan_eggnog_additions(
                 status = "filtered_confidence"
             elif qualifier == "gene":
                 if feature_genes or paired_genes:
-                    status = "existing_gene"
+                    status = (
+                        "existing_gene"
+                        if set(feature_genes) == {value}
+                        and set(paired_genes) == {value}
+                        else "unresolved_identity"
+                        if feature_genes
+                        and paired_genes
+                        and set(feature_genes) != set(paired_genes)
+                        else "suppressed_authoritative_gene"
+                    )
                 elif not paired:
                     status = "unpaired_gene"
                 else:
@@ -935,6 +956,14 @@ def plan_eggnog_additions(
                     "seed_ortholog": hit.seed_ortholog,
                     "e_value": hit.evalue,
                     "score": hit.score,
+                    "target_feature_uids": [
+                        feature_uid(target)
+                        for target in (
+                            (paired, feature)
+                            if qualifier == "gene" and paired
+                            else (feature,)
+                        )
+                    ],
                 }
             )
         if hit.query_id in emitted_by_query and add_feature_provenance:
@@ -1127,6 +1156,18 @@ def merge(
                 else ""
             ),
             "allow_imported_translations": allow_imported_translations,
+            "translation_evidence": {
+                query: evidence.as_dict()
+                for query, evidence in (translation_evidence or {}).items()
+            },
+            "translation_evidence_parent_sha256": getattr(
+                translation_evidence, "manifest_sha256", ""
+            ),
+            "policies": {
+                "add_comment_note": add_comment_note,
+                "add_feature_provenance": add_feature_provenance,
+                "clean_gene_suffix": clean_gene_suffix,
+            },
             **candidate_counts,
         },
     )
diff --git a/src/enrich_bakta_lib/sources/kofam.py b/src/enrich_bakta_lib/sources/kofam.py
index 7efa5fc..44dd0ad 100644
--- a/src/enrich_bakta_lib/sources/kofam.py
+++ b/src/enrich_bakta_lib/sources/kofam.py
@@ -244,7 +244,7 @@ def plan_kofam_additions(
                     feature,
                     "note",
                     note,
-                    "KofamScan provenance",
+                    "KofamScan hit evidence",
                     f"row {hit.row_number}",
                     order,
                 )
@@ -288,6 +288,7 @@ def plan_kofam_additions(
                 "e_value": hit.e_value,
                 "definition": hit.definition,
                 "ko_already_present": already_present,
+                "status": "existing" if already_present else "planned",
                 "emitted_qualifiers": " | ".join(emitted),
                 "kofam_row": hit.row_number,
             }
@@ -385,6 +386,8 @@ def merge(
             add_comment_note=add_comment_note,
             merge_timestamp=merge_timestamp,
             invalid_ec_policy=baktfold_invalid_ec_policy,
+            translation_evidence=translation_evidence,
+            allow_imported_translations=allow_imported_translations,
         )
         insertions.extend(baktfold_insertions)
         combined_stats["baktfold"] = baktfold_stats
@@ -464,6 +467,7 @@ def merge(
                 )
             ),
             "baktfold_invalid_ec_policy": baktfold_invalid_ec_policy,
+            "baktfold_parity": combined_stats.get("baktfold", {}).get("parity", {}),
             "baktfold_invalid_ec_values": combined_stats.get("baktfold", {}).get(
                 "invalid_ec_values", []
             ),
@@ -478,6 +482,17 @@ def merge(
                 else ""
             ),
             "allow_imported_translations": allow_imported_translations,
+            "translation_evidence": {
+                query: evidence.as_dict()
+                for query, evidence in (translation_evidence or {}).items()
+            },
+            "translation_evidence_parent_sha256": getattr(
+                translation_evidence, "manifest_sha256", ""
+            ),
+            "policies": {
+                "add_comment_note": add_comment_note,
+                "add_feature_provenance": add_feature_provenance,
+            },
             **candidate_counts,
         },
     )
diff --git a/src/enrich_bakta_lib/sources/value_rules.py b/src/enrich_bakta_lib/sources/value_rules.py
index d0b5c2c..461079c 100644
--- a/src/enrich_bakta_lib/sources/value_rules.py
+++ b/src/enrich_bakta_lib/sources/value_rules.py
@@ -90,6 +90,10 @@ def validate_preferred_name(value: str, *, clean_suffix: bool) -> ValueValidatio
         return ValueValidationResult(
             value, None, "invalid", "preferred name is empty after normalization"
         )
+    if candidate in MISSING_VALUES:
+        return ValueValidationResult(
+            value, None, "missing", "preferred name is missing after normalization"
+        )
     return ValueValidationResult(value, candidate, "valid", "preferred name accepted")
 
 
diff --git a/src/enrich_bakta_lib/workflows/enrich.py b/src/enrich_bakta_lib/workflows/enrich.py
index 88efc71..33709a6 100644
--- a/src/enrich_bakta_lib/workflows/enrich.py
+++ b/src/enrich_bakta_lib/workflows/enrich.py
@@ -66,18 +66,11 @@ def enrich(
         raise MergeError("--eggnog-schema requires --eggnog")
     if context_report_path is not None and eggnog_path is None:
         raise MergeError("--context-report requires --eggnog")
-    if translation_evidence_manifest is not None and not (
-        kofamscan_path or eggnog_path
-    ):
-        raise MergeError(
-            "--translation-evidence-manifest requires --kofamscan or --eggnog"
-        )
     base_data = read_input_bytes(bakta_path, "Bakta")
     validate_genbank_semantics(base_data, "Bakta input")
     base = parse_genbank_bytes(base_data, "Bakta input")
     if (
         has_translation_evidence_marker(base_data)
-        and (kofamscan_path or eggnog_path)
         and translation_evidence_manifest is None
     ):
         raise MergeError(
@@ -96,6 +89,22 @@ def enrich(
         "operation": "unified-enrichment",
         "merge_timestamp": merge_timestamp or "",
         "faa_sha256": sha256_bytes(faa_data) if faa_path else "",
+        "translation_evidence": {
+            query: evidence.as_dict()
+            for query, evidence in (translation_evidence or {}).items()
+        },
+        "translation_evidence_parent_sha256": getattr(
+            translation_evidence, "manifest_sha256", ""
+        ),
+        "translation_evidence_manifest": str(translation_evidence_manifest)
+        if translation_evidence_manifest
+        else "",
+        "allow_imported_translations": allow_imported_translations,
+        "policies": {
+            "add_comment_note": add_comment_note,
+            "add_feature_provenance": add_feature_provenance,
+            "clean_gene_suffix": clean_gene_suffix,
+        },
     }
     other_inputs: list[Path] = []
     if faa_path:
@@ -115,10 +124,13 @@ def enrich(
             merge_timestamp=merge_timestamp,
             starting_order=0,
             invalid_ec_policy=baktfold_invalid_ec_policy,
+            translation_evidence=translation_evidence,
+            allow_imported_translations=allow_imported_translations,
         )
         insertions.extend(planned)
         evidence_rows.extend(stats.get("invalid_ec_values", []))
         metadata["baktfold_sha256"] = stats["baktfold_sha256"]
+        metadata["baktfold_parity"] = stats["parity"]
         metadata["baktfold_version"] = stats["baktfold_version"]
         metadata["baktfold_version_detected"] = stats["baktfold_version_detected"]
         metadata["baktfold_translation_mismatch_policy"] = stats[
diff --git a/src/enrich_bakta_lib/workflows/restore_translations.py b/src/enrich_bakta_lib/workflows/restore_translations.py
index 1330f77..e762e99 100644
--- a/src/enrich_bakta_lib/workflows/restore_translations.py
+++ b/src/enrich_bakta_lib/workflows/restore_translations.py
@@ -9,7 +9,9 @@ import sys
 from pathlib import Path
 from typing import Any
 
+from Bio.Data.CodonTable import TranslationError
 from Bio.Seq import Seq
+from Bio.SeqFeature import ExactPosition
 
 from enrich_bakta_lib.core.merge_engine import (
     TOOL_VERSION,
@@ -22,6 +24,8 @@ from enrich_bakta_lib.core.merge_engine import (
     comment_insertion,
     feature_uid,
     finalize_merge,
+    has_translation_evidence_marker,
+    load_translation_evidence,
     normalize_protein,
     parse_faa,
     parse_genbank_bytes,
@@ -55,6 +59,23 @@ def _genomic_translation_check(
     location = feature.semantic_location
     if location is None:
         return "unresolved", "semantic CDS location is unavailable"
+    if any(
+        key in feature.qualifiers
+        for key in ("exception", "transl_except", "ribosomal_slippage")
+    ):
+        return "unresolved", "exceptional CDS translation model is unsupported"
+    if getattr(location, "operator", "join") != "join":
+        return "unresolved", "CDS location operator is unsupported"
+    record_length = len(base.records[feature.record_index].sequence)
+    for part in location.parts:
+        if part.ref is not None or part.ref_db is not None:
+            return "unresolved", "remote CDS location is unsupported"
+        if type(part.start) is not ExactPosition or type(part.end) is not ExactPosition:
+            return "unresolved", "partial or fuzzy CDS location is unsupported"
+        if not 0 <= int(part.start) < int(part.end) <= record_length:
+            return "unresolved", "CDS location lies outside the source record"
+        if part.strand not in (-1, 1) or part.strand != location.strand:
+            return "unresolved", "CDS location strand is unsupported"
     codon_start_values = feature.values("codon_start")
     transl_table_values = feature.values("transl_table")
     if len(codon_start_values) != 1 or len(transl_table_values) != 1:
@@ -66,13 +87,21 @@ def _genomic_translation_check(
         return "unresolved", "codon_start and transl_table must be integers"
     if codon_start not in (1, 2, 3) or transl_table < 1:
         return "unresolved", "codon_start or transl_table is outside supported values"
+    if codon_start != 1:
+        return "unresolved", "offset CDS translation requires a supported partial model"
     try:
         sequence = Seq(base.records[feature.record_index].sequence.decode("ascii"))
         coding = location.extract(sequence)[codon_start - 1 :]
         if len(coding) == 0 or len(coding) % 3:
             return "unresolved", "CDS sequence is incomplete after codon_start"
-        translated = normalize_protein(str(coding.translate(table=transl_table)))
-    except (UnicodeDecodeError, ValueError, TypeError) as exc:
+        translated = str(coding.translate(table=transl_table, cds=True))
+    except (
+        UnicodeDecodeError,
+        ValueError,
+        TypeError,
+        KeyError,
+        TranslationError,
+    ) as exc:
         return "unresolved", f"genomic translation is unsupported: {exc}"
     if translated != protein:
         return "failed", "FAA protein differs from the independent genomic translation"
@@ -87,6 +116,7 @@ def plan_translation_restoration(
     faa_data: bytes,
     starting_order: int = 0,
     translation_policy: str = "validated-only",
+    translation_evidence: dict[str, TranslationEvidence] | None = None,
 ) -> tuple[list[Insertion], list[dict[str, Any]], dict[str, Any]]:
     """Plan policy-controlled FAA-backed translations for pseudogene CDSs."""
     if translation_policy not in TRANSLATION_POLICIES:
@@ -112,7 +142,7 @@ def plan_translation_restoration(
     base_hash = sha256_bytes(base.data)
     faa_hash = sha256_bytes(faa_data)
     restored: list[str] = []
-    translation_evidence: dict[str, TranslationEvidence] = {}
+    origin_evidence: dict[str, TranslationEvidence] = dict(translation_evidence or {})
     translated_queries = [
         query_id
         for query_id in sorted(requested)
@@ -124,6 +154,8 @@ def plan_translation_restoration(
             proteins,
             translated_queries,
             source_name="eggNOG translation restoration",
+            translation_evidence=translation_evidence,
+            allow_imported_translations=True,
         )
     for query_id in sorted(requested):
         feature = cds[query_id]
@@ -137,6 +169,20 @@ def plan_translation_restoration(
                 f"refusing to restore {query_id!r}: empty /translation qualifier"
             )
         if translations:
+            evidence_rows.append(
+                {
+                    "entry_type": "translation_validation",
+                    "query_id": query_id,
+                    "status": "supported_existing",
+                    "artifact_sequence_match": True,
+                    "genomic_validation_status": "not_attempted",
+                    "translation_origin": (
+                        origin_evidence[query_id].origin
+                        if query_id in origin_evidence
+                        else "gbff_encoded"
+                    ),
+                }
+            )
             continue
         if (
             "pseudo" not in feature.qualifiers
@@ -187,7 +233,7 @@ def plan_translation_restoration(
         )
         order += 1
         restored.append(query_id)
-        translation_evidence[query_id] = evidence
+        origin_evidence[query_id] = evidence
         evidence_rows.append(
             {
                 "entry_type": "translation_restoration",
@@ -237,7 +283,7 @@ def plan_translation_restoration(
             "translation_policy": translation_policy,
             "translation_evidence": {
                 query_id: evidence.as_dict()
-                for query_id, evidence in translation_evidence.items()
+                for query_id, evidence in origin_evidence.items()
             },
         },
     )
@@ -253,6 +299,7 @@ def restore(
     eggnog_version: str | None = None,
     eggnog_schema: str | None = None,
     translation_policy: str = "validated-only",
+    translation_evidence_manifest: Path | None = None,
 ) -> dict[str, Any]:
     _base, base_data, _faa_data, table, insertions, evidence_rows, stats = (
         prepare_restoration(
@@ -262,13 +309,18 @@ def restore(
             eggnog_version=eggnog_version,
             eggnog_schema=eggnog_schema,
             translation_policy=translation_policy,
+            translation_evidence_manifest=translation_evidence_manifest,
         )
     )
     final = finalize_merge(
         base_path=bakta_path,
         base_data=base_data,
         output_path=output_path,
-        other_inputs=[faa_path, eggnog_path],
+        other_inputs=[
+            faa_path,
+            eggnog_path,
+            *([translation_evidence_manifest] if translation_evidence_manifest else []),
+        ],
         insertions=insertions,
         manifest_path=manifest_path,
         evidence_rows=evidence_rows,
@@ -283,6 +335,9 @@ def restore(
             "restored_translation_count": stats["restored_translation_count"],
             "translation_policy": stats["translation_policy"],
             "translation_evidence": stats["translation_evidence"],
+            "translation_evidence_parent_sha256": stats[
+                "translation_evidence_parent_sha256"
+            ],
         },
     )
     result = {**stats, **final}
@@ -304,6 +359,7 @@ def prepare_restoration(
     eggnog_version: str | None = None,
     eggnog_schema: str | None = None,
     translation_policy: str = "validated-only",
+    translation_evidence_manifest: Path | None = None,
 ) -> tuple[
     RawDocument,
     bytes,
@@ -325,14 +381,30 @@ def prepare_restoration(
         schema_id=eggnog_schema,
     )
     base = parse_genbank_bytes(base_data, "Bakta input")
+    if (
+        has_translation_evidence_marker(base_data)
+        and translation_evidence_manifest is None
+    ):
+        raise MergeError(
+            "restored input requires --translation-evidence-manifest for restoration reruns"
+        )
+    origin_evidence = (
+        load_translation_evidence(translation_evidence_manifest, base)
+        if translation_evidence_manifest
+        else None
+    )
     insertions, evidence_rows, stats = plan_translation_restoration(
         base,
         parse_faa(faa_data),
         table,
         faa_data=faa_data,
         translation_policy=translation_policy,
+        translation_evidence=origin_evidence,
     )
     stats["eggnog_sha256"] = sha256_bytes(eggnog_data)
+    stats["translation_evidence_parent_sha256"] = getattr(
+        origin_evidence, "manifest_sha256", ""
+    )
     return base, base_data, faa_data, table, insertions, evidence_rows, stats
 
 
@@ -350,6 +422,7 @@ def main() -> int:
     parser.add_argument("--manifest", type=Path)
     parser.add_argument("--eggnog-version")
     parser.add_argument("--eggnog-schema")
+    parser.add_argument("--translation-evidence-manifest", type=Path)
     parser.add_argument(
         "--translation-policy",
         choices=TRANSLATION_POLICIES,
@@ -378,6 +451,7 @@ def main() -> int:
                 eggnog_version=args.eggnog_version,
                 eggnog_schema=args.eggnog_schema,
                 translation_policy=args.translation_policy or "validated-only",
+                translation_evidence_manifest=args.translation_evidence_manifest,
             )
         else:
             assert args.output is not None
@@ -390,6 +464,7 @@ def main() -> int:
                 eggnog_version=args.eggnog_version,
                 eggnog_schema=args.eggnog_schema,
                 translation_policy=args.translation_policy or "validated-only",
+                translation_evidence_manifest=args.translation_evidence_manifest,
             )
     except MergeError as exc:
         print(f"ERROR: {exc}", file=sys.stderr)
diff --git a/tests/test_followthrough.py b/tests/test_followthrough.py
new file mode 100644
index 0000000..b9d01d7
--- /dev/null
+++ b/tests/test_followthrough.py
@@ -0,0 +1,457 @@
+from __future__ import annotations
+
+import json
+
+import pytest
+from test_merge_pipeline import eggnog_bytes, kofam_bytes, record_bytes, write
+
+from enrich_bakta import enrich
+from enrich_bakta_lib.core.decisions import validate_candidate_ledger
+from enrich_bakta_lib.core.merge_engine import (
+    MergeError,
+    load_translation_evidence,
+    parse_genbank_bytes,
+)
+from enrich_bakta_lib.sources.eggnog import merge as merge_egg
+from enrich_bakta_lib.sources.eggnog import parse_eggnog_tsv
+from enrich_bakta_lib.sources.kofam import merge as merge_ko
+from enrich_bakta_lib.sources.value_rules import validate_preferred_name
+from enrich_bakta_lib.workflows.restore_translations import restore
+
+ROW = "T_0001\tseed\t1e-4\t10\t-\tabc\t-\t-\tK00001\t-\thhhhhhhhhhhhh"
+
+
+def inputs(
+    tmp_path, *, sequence=b"atgaaataa", protein=b"MK", extra=b"", location=b"1..9"
+):
+    data = record_bytes("TEST", "T_0001").replace(
+        b'/translation="MK"', b"/pseudo" + extra
+    )
+    data = data.replace(b"atgaaataa", sequence).replace(
+        b"CDS             1..9", b"CDS             " + location
+    )
+    return (
+        write(tmp_path / "base.gbff", data),
+        write(tmp_path / "base.faa", b">T_0001\n" + protein + b"\n"),
+        write(tmp_path / "egg.tsv", eggnog_bytes(ROW)),
+    )
+
+
+def restored_inputs(tmp_path):
+    base, faa, egg = inputs(tmp_path, protein=b"MM")
+    out, manifest = tmp_path / "restored.gbff", tmp_path / "restored.json"
+    restore(
+        base, faa, egg, out, manifest_path=manifest, translation_policy="import-faa"
+    )
+    return out, faa, egg, manifest
+
+
+def test_complete_alternative_start_is_validated(tmp_path):
+    base, faa, egg = inputs(tmp_path, sequence=b"gtgaaataa")
+    restore(base, faa, egg, tmp_path / "out.gbff", translation_policy="validated-only")
+
+
+@pytest.mark.parametrize(
+    "sequence,protein", [(b"atgaaaaaa", b"MKK"), (b"aaaaaataa", b"KK")]
+)
+def test_invalid_complete_cds_is_not_validated(tmp_path, sequence, protein):
+    base, faa, egg = inputs(tmp_path, sequence=sequence, protein=protein)
+    with pytest.raises(MergeError):
+        restore(
+            base, faa, egg, tmp_path / "out.gbff", translation_policy="validated-only"
+        )
+
+
+@pytest.mark.parametrize(
+    "extra,location",
+    [
+        (b"", b"<1..9"),
+        (b'\n                     /exception="ribosomal slippage"', b"1..9"),
+        (b"", b"1..12"),
+    ],
+)
+def test_unsupported_cds_model_is_not_validated(tmp_path, extra, location):
+    base, faa, egg = inputs(tmp_path, extra=extra, location=location)
+    with pytest.raises(MergeError):
+        restore(
+            base, faa, egg, tmp_path / "out.gbff", translation_policy="validated-only"
+        )
+
+
+def test_bad_translation_table_is_contextual_error(tmp_path):
+    base, faa, egg = inputs(tmp_path)
+    base.write_bytes(
+        base.read_bytes().replace(b"/transl_table=11", b"/transl_table=999")
+    )
+    with pytest.raises(MergeError):
+        restore(
+            base, faa, egg, tmp_path / "out.gbff", translation_policy="validated-only"
+        )
+
+
+@pytest.mark.parametrize("adapter", ["kofam", "eggnog", "unified"])
+def test_restore_merge_ledger_chain_survives(tmp_path, adapter):
+    base, faa, egg, parent = restored_inputs(tmp_path)
+    out, manifest = tmp_path / "merged.gbff", tmp_path / "merged.json"
+    ko = write(tmp_path / "ko.txt", kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha"))
+    common = dict(
+        manifest_path=manifest,
+        translation_evidence_manifest=parent,
+        allow_imported_translations=True,
+    )
+    if adapter == "kofam":
+        merge_ko(base, faa, ko, out, **common)
+    elif adapter == "eggnog":
+        merge_egg(base, faa, egg, out, **common)
+    else:
+        enrich(
+            bakta_path=base, faa_path=faa, eggnog_path=egg, output_path=out, **common
+        )
+    ledger = load_translation_evidence(manifest, parse_genbank_bytes(out.read_bytes()))
+    assert ledger["T_0001"].origin == "imported_faa"
+
+
+def test_restoration_rerun_requires_parent_ledger(tmp_path):
+    base, faa, egg, _ = restored_inputs(tmp_path)
+    with pytest.raises(MergeError, match="translation-evidence-manifest"):
+        restore(
+            base,
+            faa,
+            egg,
+            tmp_path / "rerun.gbff",
+            manifest_path=tmp_path / "rerun.json",
+            translation_policy="validated-only",
+        )
+
+
+@pytest.mark.parametrize(
+    "mutation", ["unknown_origin", "inconsistent_policy", "string_boolean"]
+)
+def test_invalid_ledger_fields_rejected(tmp_path, mutation):
+    base, faa, egg, manifest = restored_inputs(tmp_path)
+    payload = json.loads(manifest.read_text())
+    evidence = next(
+        row["translation_evidence"]
+        for row in payload["entries"]
+        if row.get("entry_type") == "translation_restoration"
+    )
+    if mutation == "unknown_origin":
+        evidence["origin"] = "imported_faa_typo"
+    elif mutation == "inconsistent_policy":
+        evidence["origin"] = "genomically_validated"
+    else:
+        evidence["artifact_sequence_match"] = "false"
+    manifest.write_text(json.dumps(payload))
+    with pytest.raises(MergeError):
+        load_translation_evidence(manifest, parse_genbank_bytes(base.read_bytes()))
+
+
+def test_kofam_existing_ko_retains_hit_evidence(tmp_path):
+    base = write(
+        tmp_path / "base.gbff",
+        record_bytes("TEST", "T_0001", db_xrefs=("KEGG:K00001",)),
+    )
+    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
+    ko = write(tmp_path / "ko.txt", kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha"))
+    out, manifest = tmp_path / "out.gbff", tmp_path / "out.json"
+    merge_ko(base, faa, ko, out, manifest_path=manifest, add_comment_note=False)
+    assert b"KofamScan:K00001;threshold=1;score=2;E-value=3e-4" in out.read_bytes()
+    payload = json.loads(manifest.read_text())
+    ko_decision = next(row for row in payload["decisions"] if row["field"] == "KO")
+    assert ko_decision["final_status"] == "supported_existing"
+
+
+def test_kofam_rerun_decision_is_supported_existing(tmp_path):
+    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
+    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
+    ko = write(tmp_path / "ko.txt", kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha"))
+    first, second, manifest = (
+        tmp_path / "first.gbff",
+        tmp_path / "second.gbff",
+        tmp_path / "second.json",
+    )
+    merge_ko(base, faa, ko, first)
+    merge_ko(first, faa, ko, second, manifest_path=manifest)
+    assert first.read_bytes() == second.read_bytes()
+    row = next(
+        row
+        for row in json.loads(manifest.read_text())["decisions"]
+        if row["field"] == "KO"
+    )
+    assert row["final_status"] == "supported_existing"
+
+
+def test_gene_class_provenance_removed_when_gene_conflicts(tmp_path):
+    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
+    fold = write(
+        tmp_path / "fold.gbff",
+        record_bytes("TEST", "T_0001", gene="xyz", ec_numbers=("1.2.3.4",)),
+    )
+    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
+    egg = write(tmp_path / "egg.tsv", eggnog_bytes(ROW))
+    out = tmp_path / "out.gbff"
+    enrich(
+        bakta_path=base,
+        faa_path=faa,
+        eggnog_path=egg,
+        baktfold_path=fold,
+        output_path=out,
+        add_comment_note=False,
+    )
+    assert b"Baktfold gene-symbol evidence" not in out.read_bytes()
+    assert b'/EC_number="1.2.3.4"' in out.read_bytes()
+
+
+def test_paired_gene_is_one_candidate_with_two_targets(tmp_path):
+    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
+    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
+    egg = write(tmp_path / "egg.tsv", eggnog_bytes(ROW))
+    manifest = tmp_path / "out.json"
+    merge_egg(base, faa, egg, tmp_path / "out.gbff", manifest_path=manifest)
+    genes = [
+        row
+        for row in json.loads(manifest.read_text())["decisions"]
+        if row["qualifier"] == "gene"
+    ]
+    assert len(genes) == 1
+    assert len(genes[0]["target_feature_uids"]) == 2
+    assert len(genes[0]["insertion_ids"]) == 2
+
+
+def test_two_sources_report_shared_support(tmp_path):
+    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
+    fold = write(
+        tmp_path / "fold.gbff", record_bytes("TEST", "T_0001", ec_numbers=("1.2.3.4",))
+    )
+    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
+    egg = write(
+        tmp_path / "egg.tsv",
+        eggnog_bytes(ROW.replace("\t-\t-\tK00001\t", "\t-\t1.2.3.4\tK00001\t")),
+    )
+    manifest = tmp_path / "out.json"
+    enrich(
+        bakta_path=base,
+        faa_path=faa,
+        eggnog_path=egg,
+        baktfold_path=fold,
+        output_path=tmp_path / "out.gbff",
+        manifest_path=manifest,
+    )
+    ecs = [
+        row
+        for row in json.loads(manifest.read_text())["decisions"]
+        if row["qualifier"] == "EC_number"
+    ]
+    assert len(ecs) == 2
+    assert all(row["final_status"] == "shared_support" for row in ecs)
+
+
+def test_different_authoritative_gene_is_not_supported_existing(tmp_path):
+    base = write(
+        tmp_path / "base.gbff", record_bytes("TEST", "T_0001", gene="original")
+    )
+    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
+    egg = write(tmp_path / "egg.tsv", eggnog_bytes(ROW))
+    manifest = tmp_path / "out.json"
+    merge_egg(base, faa, egg, tmp_path / "out.gbff", manifest_path=manifest)
+    row = next(
+        row
+        for row in json.loads(manifest.read_text())["decisions"]
+        if row["field"] == "Preferred_name"
+    )
+    assert row["final_status"] != "supported_existing"
+
+
+def test_reserved_weekday_query_row_is_not_silently_skipped():
+    data = eggnog_bytes(ROW, ROW.replace("T_0001", "##Monday_bad"))
+    with pytest.raises(MergeError):
+        parse_eggnog_tsv(data)
+
+
+def test_missing_sentinel_rechecked_after_suffix_cleaning():
+    result = validate_preferred_name("NA_123", clean_suffix=True)
+    assert result.status != "valid"
+
+
+def test_cog_candidate_uses_note_qualifier(tmp_path):
+    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
+    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
+    egg = write(
+        tmp_path / "egg.tsv",
+        eggnog_bytes(ROW.replace("\t-\tabc\t", "\tCOG0001\tabc\t")),
+    )
+    manifest = tmp_path / "out.json"
+    merge_egg(base, faa, egg, tmp_path / "out.gbff", manifest_path=manifest)
+    row = next(
+        row
+        for row in json.loads(manifest.read_text())["decisions"]
+        if row["field"] == "COG_category"
+    )
+    assert row["qualifier"] == "note"
+    assert row["final_status"] == "emitted"
+
+
+def test_restoration_noop_carries_imported_origin(tmp_path):
+    base, faa, egg, parent = restored_inputs(tmp_path)
+    out, manifest = tmp_path / "rerun.gbff", tmp_path / "rerun.json"
+    stats = restore(
+        base,
+        faa,
+        egg,
+        out,
+        manifest_path=manifest,
+        translation_policy="validated-only",
+        translation_evidence_manifest=parent,
+    )
+    assert stats["restored_translation_count"] == 0
+    assert stats["insertions"] == 0
+    assert out.read_bytes() == base.read_bytes()
+    ledger = load_translation_evidence(manifest, parse_genbank_bytes(out.read_bytes()))
+    assert ledger["T_0001"].origin == "imported_faa"
+
+
+def test_merge_rerun_keeps_imported_origin_and_requires_optin(tmp_path):
+    base, faa, egg, parent = restored_inputs(tmp_path)
+    first, first_manifest = tmp_path / "first.gbff", tmp_path / "first.json"
+    second, second_manifest = tmp_path / "second.gbff", tmp_path / "second.json"
+    merge_egg(
+        base,
+        faa,
+        egg,
+        first,
+        manifest_path=first_manifest,
+        translation_evidence_manifest=parent,
+        allow_imported_translations=True,
+    )
+    with pytest.raises(MergeError, match="allow-imported-translations"):
+        merge_egg(first, faa, egg, second, translation_evidence_manifest=first_manifest)
+    merge_egg(
+        first,
+        faa,
+        egg,
+        second,
+        manifest_path=second_manifest,
+        translation_evidence_manifest=first_manifest,
+        allow_imported_translations=True,
+    )
+    assert first.read_bytes() == second.read_bytes()
+    assert (
+        load_translation_evidence(
+            second_manifest, parse_genbank_bytes(second.read_bytes())
+        )["T_0001"].origin
+        == "imported_faa"
+    )
+
+
+def test_xlsx_query_alias_and_duplicate_alias_collision(tmp_path):
+    import io
+
+    openpyxl = pytest.importorskip("openpyxl")
+    from enrich_bakta_lib.sources.eggnog import parse_eggnog_path
+
+    workbook = openpyxl.Workbook()
+    worksheet = workbook.active
+    worksheet.title = "annotations"
+    lines = eggnog_bytes(ROW).decode().splitlines()
+    worksheet.append(lines[1].replace("#query", "query").split("\t"))
+    worksheet.append(lines[2].split("\t"))
+    stream = io.BytesIO()
+    workbook.save(stream)
+    path = write(tmp_path / "egg.xlsx", stream.getvalue())
+    assert (
+        parse_eggnog_path(path, expected_version="3.0.0-beta6").hits[0].query_id
+        == "T_0001"
+    )
+    worksheet.cell(1, 2, "#query")
+    stream = io.BytesIO()
+    workbook.save(stream)
+    workbook.close()
+    path.write_bytes(stream.getvalue())
+    with pytest.raises(MergeError, match="duplicate columns"):
+        parse_eggnog_path(path, expected_version="3.0.0-beta6")
+
+
+def test_paired_candidate_cannot_claim_only_one_insertion(tmp_path):
+    from enrich_bakta_lib.core.merge_engine import (
+        feature_uid,
+        insertion_uid,
+        qualifier_insertion,
+    )
+
+    base = parse_genbank_bytes(record_bytes("TEST", "T_0001"))
+    gene, cds = base.features[1:]
+    insertion = qualifier_insertion(base.data, cds, "gene", "abc", "eggNOG", "abc", 0)
+    row = {
+        "candidate_id": "candidate:test",
+        "source_id": "eggNOG",
+        "qualifier": "gene",
+        "normalized_value": "abc",
+        "final_status": "emitted",
+        "target_feature_uids": [feature_uid(gene), feature_uid(cds)],
+        "insertion_ids": [insertion_uid(insertion)],
+        "supporting_candidate_ids": [],
+    }
+    with pytest.raises(MergeError, match="every target"):
+        validate_candidate_ledger([row], [insertion], base=base)
+
+
+def test_baktfold_retains_translation_lineage(tmp_path):
+    from enrich_bakta_lib.sources.baktfold import graft
+
+    base, _, _, parent = restored_inputs(tmp_path)
+    fold = write(
+        tmp_path / "fold.gbff",
+        base.read_bytes().replace(
+            b'/product="test protein"',
+            b'/EC_number="1.2.3.4"\n                     /product="test protein"',
+        ),
+    )
+    out, manifest = tmp_path / "fold-out.gbff", tmp_path / "fold-out.json"
+    with pytest.raises(MergeError, match="allow-imported-translations"):
+        graft(base, fold, out, translation_evidence_manifest=parent)
+    graft(
+        base,
+        fold,
+        out,
+        manifest_path=manifest,
+        translation_evidence_manifest=parent,
+        allow_imported_translations=True,
+    )
+    assert (
+        load_translation_evidence(manifest, parse_genbank_bytes(out.read_bytes()))[
+            "T_0001"
+        ].origin
+        == "imported_faa"
+    )
+
+
+def test_parent_manifest_hash_uses_the_loaded_snapshot(tmp_path, monkeypatch):
+    import enrich_bakta_lib.sources.eggnog as adapter
+    from enrich_bakta_lib.core.merge_engine import sha256_bytes
+
+    base, faa, egg, parent = restored_inputs(tmp_path)
+    expected = sha256_bytes(parent.read_bytes())
+    original = adapter.load_translation_evidence
+
+    def replace_after_loading(path, document):
+        ledger = original(path, document)
+        path.write_bytes(b"a later disk version")
+        return ledger
+
+    monkeypatch.setattr(adapter, "load_translation_evidence", replace_after_loading)
+    manifest = tmp_path / "out.json"
+    adapter.merge(
+        base,
+        faa,
+        egg,
+        tmp_path / "out.gbff",
+        manifest_path=manifest,
+        translation_evidence_manifest=parent,
+        allow_imported_translations=True,
+    )
+    assert (
+        json.loads(manifest.read_text())["metadata"][
+            "translation_evidence_parent_sha256"
+        ]
+        == expected
+    )
```

## 8. Primary references and evidence files

- [Reviewed repository commit](https://github.com/WhyAdr/enrich-bakta/tree/cc98e4f081c2e0554a7108b3897f6679fe61cd1d).
- [Biopython 1.83 translation API](https://biopython.org/docs/1.83/api/Bio.Seq.html):
  complete-CDS checks and alternative-start behavior.
- [PyPA pyproject metadata guide](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/):
  PEP 639 backend support floors.
- [Setuptools license migration](https://setuptools.pypa.io/en/stable/userguide/license_migration.html).
- [ExPASy ENZYME](https://enzyme.expasy.org/): preliminary EC convention;
  [3.5.1.n3](https://enzyme.expasy.org/EC/3.5.1.n3),
  [3.6.5.n1](https://enzyme.expasy.org/EC/3.6.5.n1),
  [4.2.2.n1](https://enzyme.expasy.org/EC/4.2.2.n1).

The ZIP contains the patch, regression tests, first baseline-failure output,
compact real/chain results, reusable real-run scripts, build-floor failure log,
final build log, environment record, and checksums. It excludes genomic datasets,
generated genomic outputs, huge manifests, environments and installed packages.
