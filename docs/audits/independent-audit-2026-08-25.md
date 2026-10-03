# Independent Audit: `WhyAdr/enrich-bakta`

**Auditor:** Super Z (independent code audit)
**Date:** 2026-08-25 (Asia/Jakarta)
**Repo state audited:** `main` @ `bbb6a1b` ("Add architectural and implementation review document"), 13 tracked files, 5,706 LOC
**Method:** Full read of every source file + planning docs; static reasoning; **empirical execution of the test suite and 6 purpose-built reproduction scripts** against Python 3.12.13 / Biopython 1.86. All "REPRODUCED" findings below were confirmed by running code, not by inspection alone.
**Intended audience:** Gemini implementing the fixes. Every finding carries a severity, an evidence trail, and (where applicable) a concrete diff.

---

## Maintainer review and refined implementation disposition

This section records the repository-grounded review performed before implementation.
It takes precedence over the original severity table, example diffs, and handoff
checklist below. The audit was useful and most reproductions identify real defects,
but several conclusions or proposed fixes do not match the current producer contract,
the real C14/SM inputs, or the compatibility requirements of this repository.

Baseline re-run on the audited commit: 21 tests passed; the repository's six-file
Ruff check, Ruff format check, and mypy check passed. Both authoritative Bakta inputs
(`C14-NMZ_bakta/C14-NMZ.gbff` and `SM-NMZ_bakta/SM-NMZ.gbff`) parse without
`BiopythonParserWarning`. The real eggNOG files declare `emapper-3.0.0-beta6`, contain
4,563 and 4,497 confidence vectors respectively, and explicitly declare this order:

```text
Preferred_name GOs EC KEGG_ko KEGG_Pathway KEGG_Module KEGG_Reaction
KEGG_rclass BRITE KEGG_TC CAZy BiGG_Reaction PFAMs
```

Neither real eggNOG table contains a partial or malformed EC token. Partial-EC support
is therefore a compatibility hardening change, not a repair needed to process the two
current datasets.

| Finding | Disposition | Refined implementation |
|---|---|---|
| F-01 | Accept | Pass the exact planned bytes to `finalize_merge` and verify that the on-disk source still has the same hash before any output is written. Remove the legacy re-read path from in-repository callers. |
| F-02 | Accept with corrected scope | Preserve syntactically valid partial four-component EC values as non-promoted manifest evidence. Continue to fail on malformed EC text. Do not generalize the grammar to undocumented two- or three-component forms. |
| F-03 | Reject the P1 premise; accept the drift guard | eggNOG-mapper 3.0.0 beta was public and its documented v3 output uses the 13-character vector. The real tables also carry an explicit field-order legend. Parse and validate that legend when present, use one ordered constant, and record the resolved contract in artifacts. Do not add a v2 `--no-confidence-column` mode: v2 uses a different database generation and provides no equivalent confidence evidence, so silently treating every value as passing would weaken the scientific gate. |
| F-04 | Accept with clarified semantics | `--no-feature-provenance` suppresses Kofam `/inference`; the score/threshold/E-value `/note` remains because it is the feature-level evidence for the emitted KO, not merely producer provenance. |
| F-05 | Accept | Unified metadata records the FAA hash and detected Baktfold version, and all manifests receive the enrich-bakta tool version centrally. |
| F-06 | Accept without an extra CLI mode | Emit one `eggnog_row` entry per input row and reference it by `row_number`; do not repeat `raw_fields` on every candidate/context entry. This keeps complete source data while making the deterministic default compact. |
| F-07 | Accept with an exact-key index | Build the paired-gene index once using the already parsed exact `(record_index, locus_tag, location)` key. Do not erase biological boundary markers or otherwise broaden location equivalence. |
| F-08 | Accept | Reject control characters and use a wrapper that cannot place `/` at the start of a continuation line. Test exact parser round-trip and idempotence. |
| F-09 | Accept the safe part | Reuse bytes already read for Bakta, FAA, and TSV eggNOG inputs. Retain the simple in-memory reversibility reconstruction; bacterial GBFF size does not justify a riskier audit rewrite. |
| F-10 | Accept | Reject reserved-marker query IDs loudly instead of silently dropping or misparsing them. |
| F-11 | Accept | Return an explicit `EggnogPlan` carrying the context report instead of placing a private payload in the stats mapping. |
| F-12 | Accept | Validate exact FAA/GBFF identity for requested CDSs that already have translations; continue to restore only translationless `/pseudogene` CDSs. Make the dry-run output positional optional only in dry-run mode. |
| F-13 | Accept selectively | Wrap non-ASCII protein hashing errors, reject empty XLSX headers, correct mapped/hit statistics, expose `--clean-gene-suffix` in the unified CLI, and report Baktfold version-detection status. Do **not** reject negative Kofam thresholds/scores: HMM scores/cutoffs are not guaranteed to be non-negative. Avoid a breaking JSON return-shape change in this pass. The proposed `CONTIG` concern is inapplicable because comment insertion occurs in the header before `FEATURES`, while GenBank `CONTIG` sequence data occurs after the feature table. |
| F-14 | Accept infrastructure; retain legacy history | Add declared project/dev dependencies, a single Ruff/pytest configuration, CI, and tool-version stamping. Keep `normalize_baktfold.py` in place with an explicit legacy warning; moving or deleting it would be a separate compatibility decision. |

Additional architecture dispositions:

- Promote `BiopythonParserWarning` to a merge error; this is compatible with both real
  authoritative Bakta inputs and closes the demonstrated warning-only corruption path.
- Reconcile the combined Baktfold/Kofam path and replace magic million-sized order
  blocks with sequential source allocation or explicit bounds.
- Do not attempt transactional multi-file promotion, directory `fsync`, schema
  finalization, or an O(1)-memory reconstruction rewrite in the correctness pass.
  Those require a separate cross-platform artifact-transaction design. JSON schemas
  should follow, rather than precede, stabilization of the manifest rows in this pass.

### Refined gated phases

1. **Correctness:** planned-byte binding, safe qualifier wrapping, parser-warning
   promotion, targeted regression tests; full test/lint/type gate; commit.
2. **Data contracts:** partial-EC evidence, confidence-legend validation and artifact
   recording, reserved-ID rejection; full gate; commit.
3. **Provenance and scale:** single-read flow, complete unified hashes/versions,
   tool-version stamping, paired-gene index, one-copy eggNOG raw rows, Kofam provenance
   semantics, uniform reconciliation; full gate plus real-input dry runs; commit.
4. **API/CLI and infrastructure:** explicit `EggnogPlan`, restoration identity check,
   dry-run UX, selected validation/stat fixes, unified clean-gene option, project config,
   CI, and README updates; full gate; commit.
5. **Final integration:** run C14 and SM enrichment into ignored workspace-local outputs,
   validate GenBank semantics, manifests/context reports, and byte-identical idempotence;
   inspect the complete committed diff and history.

---

## 1. Executive Summary

`enrich-bakta` is a genuinely well-conceived byte-preserving annotation grafting engine. Its core invariants — offset-based splicing, allowlisted qualifiers, reversibility audit, pre-write identity validation, atomic promotion — are correctly implemented for the *single-threaded, immutable-input* case, and the repo's own test suite (21 tests) passes. Determinism, CRLF preservation, and multi-source idempotence all verified clean in this audit.

However, the audit found **one P0 correctness bug that can silently corrupt output**, **two P1 hard-failures on legitimate real-world data**, and a series of P2/P3 contract, performance, and reproducibility defects:

| ID | Severity | Finding | Verified |
|----|----------|---------|----------|
| F-01 | **P0** | `finalize_merge` re-reads the base file from disk after planning — a TOCTOU window that splices planned offsets into a *different* byte stream, **passing the reversibility audit while corrupting the output** (insertion landed inside `/translation` in the repro) | ✅ executed |
| F-02 | **P1** | Partial EC numbers (`1.2.3.-`, extremely common in real eggNOG output) abort the entire merge with `MergeError` instead of being skipped and manifested | ✅ executed |
| F-03 | **P1** | No released eggNOG-mapper version (2.1.x) can be parsed at all: `annotation_confidence` is not a public emapper column; the 13-position confidence↔field mapping is an *implicit, unvalidated contract* — if the producer reorders fields, confidence gates silently apply to the wrong evidence (wrong science, no error) | ✅ executed |
| F-04 | **P2** | `--no-feature-provenance` is silently ignored by the KofamScan planner (flag plumbed to Baktfold/eggNOG only) | ✅ executed |
| F-05 | **P2** | Unified `enrich()` manifest omits `faa_sha256` and `baktfold_version` — breaks the reproducibility contract that all single-source manifests honor | ✅ executed |
| F-06 | **P2** | Manifest bloat: full `raw_fields` row embedded in *every* candidate *and* context row → ~8 KB/protein (≈40 MB manifest for a 4.6k-CDS genome) | ✅ measured |
| F-07 | **P2** | `_paired_gene` is O(hits × features): **34.3 s** for 4,600 CDS on one record; scales quadratically with genome size | ✅ measured |
| F-08 | **P2** | `format_qualifier` can wrap a value so a continuation line begins with `/`, which re-parses as a phantom qualifier (breaks dedup/idempotence, corrupts round-trip) | ✅ executed |
| F-09 | **P2** | Redundant I/O: one unified run performs 9 `read_bytes` calls (FAA ×3, base GBFF ×2, eggNOG ×2); ~3× peak memory in `apply_insertions` | ✅ counted |
| F-10 | **P3** | Query IDs starting with `#` are silently dropped from eggNOG TSV (data loss, no manifest entry); same class of issue for Kofam `*`-leading IDs | ✅ executed |
| F-11 | **P3** | `plan_eggnog_additions` smuggles the context report out via a private `_context_report` key in the public stats dict | inspection |
| F-12 | **P3** | `restore_bakta_translations` skips the FAA↔GBFF identity check for queries whose CDS already has a translation (violates the project's own identity doctrine) | inspection |
| F-13 | **P3** | Kofam parser accepts negative threshold/score; `protein_sha256` raises raw `ValueError` (not `MergeError`) on non-ASCII; XLSX accepts empty header names; `--dry-run` still demands an output path; `mapped CDSs=` comment stat counts all CDS not mapped ones; unified CLI lacks `--clean-gene-suffix`; `--json` output shapes differ per CLI | inspection |
| F-14 | **P3** | Engineering infra: no `pyproject.toml`, no CI, unpinned Biopython, no tool-version stamping in manifests, ruff scope excludes `normalize_baktfold.py`; `normalize_baktfold.py` (legacy) ships with non-atomic writes and hardcoded fabricated dates | ✅ ruff run |

**What was verified SOUND (do not churn these while patching):** reversibility-audit mechanics in the immutable-input case; qualifier allowlist; CRLF/newline detection; comment insertion placement (GenBank-spec correct: `COMMENT` is the last header key before `FEATURES`); strict parity checks; gene-conflict reconciliation policies incl. paired gene+CDS emission; 3-source unified idempotence (rerun on enriched base ⇒ 0 insertions, byte-identical); deterministic output without timestamps; atomic `os.replace` promotion; path collision guards; path traversal of manifests/sidecars.

---

## 2. Baseline Verification

```
$ python -m pytest -q
21 passed in 0.78s                       # Python 3.12.13, Biopython 1.86

$ ruff check .                            # ruff 0.16.4, repo default config
Found 8 errors.  (7× EXE001 shebang-not-executable — Windows checkout artifact;
                  1× I001 unsorted-imports in normalize_baktfold.py)

$ ruff format --check .
2 files would be reformatted              # normalize_baktfold.py + a fenced code
                                          # block inside the review .md — both
                                          # outside the README's ruff target list

$ mypy --ignore-missing-imports <all 7 modules>
Success: no issues found in 7 source files
```

Note: the README's validation commands lint only 5 of 7 Python files; `normalize_baktfold.py` and `restore_baktfold.py` escape the gate. No `ruff.toml`/`pyproject.toml` exists, so "passing ruff" depends on each contributor invoking the same explicit file list.

---

## 3. P0 — F-01: Base-file TOCTOU produces silently corrupted, "audited" output

### Mechanism

Every caller does:

```python
base_data = bakta_path.read_bytes()          # read #1 — planning
base = parse_genbank_bytes(base_data, ...)   # offsets derived from read #1
insertions, ... = plan_..._additions(base, ...)
final = finalize_merge(base_path=bakta_path, ...)
```

and `finalize_merge` (`merge_engine.py:944`) then does:

```python
base_data = base_path.read_bytes()           # read #2 — splicing
merged, applied = apply_insertions(base_data, insertions, ...)
```

Insertion offsets are byte offsets into **read #1**, but they are spliced into **read #2**. The reversibility audit reconstructs by stripping the recorded payloads from the spliced stream and compares against **read #2** — so the audit is *self-consistent by construction* and cannot detect that the offsets no longer point where the planner intended. The only backstop is `validate_genbank_semantics`, which is a Biopython parse that emits **warnings, not exceptions**, for many classes of corruption.

### Empirical proof (scripts/audit_repros2.py, condensed)

A 4-byte in-place growth of the base file between planning and finalize (e.g., `/product="test protein"` → `/product="test protein XX"`) produced:

```
merge SUCCEEDED, self_check=True, output_sha256=ac09571cc61af333...
Biopython reparse OK. CDS qualifiers: {'locus_tag': ['T_0001'],
  'product': ['test protein XX'],
  'translation': ['M/EC_number="1.2.3.4'],      <-- CORRUPTION
  'inference': ['DESCRIPTION:similar to AA sequence:eggNOG:seed.1']}
EC_number insertion context: '.../translation="M                     /EC_number="1.2.3.4"\n...'
```

The `/EC_number="1.2.3.4"` insertion was spliced **into the middle of the `/translation` value**. The merge reported success; `self_check: True`; the manifest's `base_sha256` describes the *second* read while the planner validated the first. Realistic triggers: a pipeline writing the GBFF concurrently (Bakta still finishing, NFS caching, an upstream step rewriting in place, or the same run invoked twice against a partially promoted file).

### Fix (Diff A)

`merge_engine.py` — accept the validated bytes; never re-read:

```diff
--- a/merge_engine.py
+++ b/merge_engine.py
@@ def finalize_merge(
 def finalize_merge(
     *,
     base_path: Path,
+    base_data: bytes | None = None,
     output_path: Path,
     other_inputs: Iterable[Path],
     insertions: list[Insertion],
@@
-    base_data = base_path.read_bytes()
+    if base_data is None:
+        # Legacy convenience path for direct API callers. CLI code paths
+    # MUST pass the bytes they parsed, closing the plan/splice TOCTOU window:
+    # insertion offsets are only meaningful for the exact stream that was
+    # planned against.
+        base_data = base_path.read_bytes()
+    if sha256_bytes(base_data) != sha256_file(base_path):
+        raise MergeError(
+            "base file changed on disk since planning; refusing to splice "
+            "planned offsets into a different byte stream (re-run the merge)"
+        )
     merged, applied = apply_insertions(
         base_data, insertions, extra_allowed_qualifiers=extra_allowed_qualifiers
     )
```

The disk-hash comparison is a cheap (streamed, 1 MB chunks via existing `sha256_file`) guard that turns silent corruption into a loud error if someone passes stale bytes. Update all five call sites:

```diff
--- a/graft_baktfold_additions.py
+++ b/graft_baktfold_additions.py
@@     final = finalize_merge(
         base_path=bakta_path,
+        base_data=bakta_data,
         output_path=output_path,

--- a/merge_kofamscan_bakta.py
+++ b/merge_kofamscan_bakta.py
@@     final = finalize_merge(
         base_path=bakta_path,
+        base_data=base_data,
         output_path=output_path,

--- a/merge_eggnog_bakta.py
+++ b/merge_eggnog_bakta.py
@@     final = finalize_merge(
         base_path=bakta_path,
+        base_data=base_data,
         output_path=output_path,

--- a/enrich_bakta.py
+++ b/enrich_bakta.py
@@     return finalize_merge(
         base_path=bakta_path,
+        base_data=base_data,
         output_path=output_path,

--- a/restore_bakta_translations.py
+++ b/restore_baktfold_translations.py-equivalent (restore_bakta_translations.py)
@@ def restore(...):
-    _base, _faa_data, table, insertions, evidence_rows, stats = prepare_restoration(...)
+    _base, base_data, _faa_data, table, insertions, evidence_rows, stats = (
+        prepare_restoration(...)
+    )
     final = finalize_merge(
         base_path=bakta_path,
+        base_data=base_data,
         output_path=output_path,
@@ def prepare_restoration(...) -> tuple[...]:
-    ) -> tuple[RawDocument, bytes, EggnogTable, list[Insertion], list[dict[str, Any]], dict[str, Any]]:
+    ) -> tuple[RawDocument, bytes, bytes, EggnogTable, list[Insertion], list[dict[str, Any]], dict[str, Any]]:
@@     return base, faa_data, table, insertions, evidence_rows, stats
+    return base, base_data, faa_data, table, insertions, evidence_rows, stats
```

This also removes one of the two base-file reads (see F-09) and guarantees `manifest.base_sha256` describes the bytes the planner actually validated.

**Regression test to add (T-01):**

```python
def test_finalize_splices_the_bytes_that_were_planned(tmp_path, monkeypatch):
    # Plan against pristine base; mutate the file before finalize; the merge
    # must either use the passed bytes or abort loudly — never splice planned
    # offsets into a different stream.
    ...
```

---

## 4. P1 — F-02: Partial EC numbers abort the whole eggNOG merge

### Mechanism

`merge_eggnog_bakta.py`:

```python
EC_RE = re.compile(r"\d+\.\d+\.\d+\.\d+\Z")   # fully specified only
...
ec=_tokens(values["EC"].strip(), field="EC", row_number=row_number,
            pattern=EC_RE, prefix="ec:"),
```

`_tokens` raises `MergeError` on any token that fails the pattern. Real eggNOG-mapper output routinely contains partial ECs (`1.2.3.-`, `2.7.-.-`, `4.2.1`); a **single** such token anywhere in the table kills the entire run — even for proteins whose rows are clean. This also contradicts the README's own contract ("Feature-level fields are … **fully specified EC**" — i.e., partials should be *not promoted*, not *fatal*) and is inconsistent with the Baktfold path, where `EC:1.2.3.-` is happily grafted (see `tests/test_merge_pipeline.py:124`, source db_xrefs include `"EC:1.2.3.-"`).

```
REPRODUCED -> MergeError: eggNOG row 3: invalid EC value '1.2.3.-'
```

### Fix (Diff B)

Tolerate partial ECs as a distinct evidence class: parse, manifest, never emit:

```diff
--- a/merge_eggnog_bakta.py
+++ b/merge_eggnog_bakta.py
@@ EC_RE = re.compile(r"\d+\.\d+\.\d+\.\d+\Z")
+EC_PARTIAL_RE = re.compile(
+    r"\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\Z"      # 1.2.3.-   (4 components)
+    r"|\d+\.(?:\d+|-)\.(?:\d+|-)\Z"                # 1.2.-     (3 components)
+    r"|\d+\.(?:\d+|-)\Z"                           # 1.-       (2 components)
+)
@@ class EggnogHit:
     ec: tuple[str, ...]
+    ec_partial: tuple[str, ...] = ()
@@ def _build_hit(values, row_number):
+    ec_full, ec_partial = _split_ec(values["EC"].strip(), row_number)
     return EggnogHit(
         ...
-        ec=_tokens(values["EC"].strip(), field="EC", row_number=row_number,
-                   pattern=EC_RE, prefix="ec:"),
+        ec=ec_full,
+        ec_partial=ec_partial,
         ...
     )
+
+def _split_ec(value: str, row_number: int) -> tuple[tuple[str, ...], tuple[str, ...]]:
+    """Fully specified ECs are promotable; partial ECs are manifested but
+    never emitted; anything else is a hard error."""
+    if value.strip() in MISSING_ANNOTATION_VALUES:
+        return (), ()
+    full: list[str] = []
+    partial: list[str] = []
+    for raw in value.split(","):
+        token = raw.strip()
+        if token.lower().startswith("ec:"):
+            token = token[3:]
+        if not token:
+            continue
+        if EC_RE.fullmatch(token):
+            if token not in full:
+                full.append(token)
+        elif EC_PARTIAL_RE.fullmatch(token):
+            if token not in partial:
+                partial.append(token)
+        else:
+            raise MergeError(
+                f"eggNOG row {row_number}: invalid EC value {raw!r}"
+            )
+    return tuple(full), tuple(partial)
```

And in `plan_eggnog_additions`, emit one manifest row per partial with `"status": "skipped_partial_ec"` (mirroring `filtered_confidence`), so the audit trail preserves the evidence without aborting:

```diff
@@     for hit in table.hits:
         feature = cds[hit.query_id]
+        for partial in hit.ec_partial:
+            evidence.append({
+                "entry_type": "eggnog_candidate",
+                "evidence_class": "direct_annotation",
+                "row_number": hit.row_number,
+                "query_id": hit.query_id,
+                "record": feature.record_id,
+                "feature_type": feature.feature_type,
+                "locus_tag": feature.locus_tag or "",
+                "field": "EC",
+                "raw_value": partial,
+                "source_token": partial,
+                "normalized_value": partial,
+                "confidence_code": hit.confidence[CONFIDENCE_FIELD_INDEX["EC"]],
+                "status": "skipped_partial_ec",
+                "emitted_qualifiers": "",
+                "seed_ortholog": hit.seed_ortholog,
+                "e_value": hit.evalue,
+                "score": hit.score,
+            })
```

Also add `stats["partial_ec_skipped"] = sum(len(h.ec_partial) for h in table.hits)`.

**Test to add (T-02):** a table containing `1.2.3.-` must merge cleanly, emit no `/EC_number` for the partial, and record a `skipped_partial_ec` manifest row.

---

## 5. P1 — F-03: The eggNOG confidence-column contract is implicit, unvalidated, and incompatible with every released emapper

### 5.1 No released eggNOG-mapper emits `annotation_confidence`

```
REPRODUCED -> MergeError: eggNOG table is missing required column 'annotation_confidence'
```

 fed with a standard `emapper-2.1.12` annotations file (21 canonical columns). Every public emapper release (2.1.x line) lacks this column, so the tool cannot ingest any stock producer output today. The repo's planning docs reference "eggNOG-mapper 3.0.0-beta6" and real workspace TSV/XLSX pairs carrying a 13-code confidence string — i.e., the contract is bound to a **non-public, custom or pre-release producer**. Nothing in the repo documents that producer, its provenance, or how to regenerate it. If these inputs were produced by a bespoke fork, reproducibility of the whole pipeline now depends on an undocumented external binary.

### 5.2 The 13-position ↔ field mapping is a silent-failure contract

```python
CONFIDENCE_FIELD_INDEX = {"Preferred_name": 0, "GOs": 1, "EC": 2, "KEGG_ko": 3,
    "KEGG_Pathway": 4, "KEGG_Module": 5, "KEGG_Reaction": 6, "KEGG_rclass": 7,
    "BRITE": 8, "KEGG_TC": 9, "CAZy": 10, "BiGG_Reaction": 11, "PFAMs": 12}
```

`_validate_confidence` only checks "13 chars from {l,m,h,-}". Nothing ties position *i* of that string to the field the constant claims. If the producer reorders scored fields (or drops one), `--min-eggnog-confidence medium` silently gates the **wrong evidence class** — no parse error, wrong annotation decisions, confident manifest rows with mislabeled `confidence_code`. This is the highest "wrong science" risk in the codebase.

### 5.3 Fix (Diff C): make the contract explicit, validated, and recorded

1. **Declare the contract** as a single ordered tuple used everywhere (index built from it, not hand-written):

```diff
--- a/merge_eggnog_bakta.py
+++ b/merge_eggnog_bakta.py
@@-CONFIDENCE_FIELD_INDEX = {
@@-    "Preferred_name": 0, "GOs": 1, "EC": 2, "KEGG_ko": 3, "KEGG_Pathway": 4,
@@-    "KEGG_Module": 5, "KEGG_Reaction": 6, "KEGG_rclass": 7, "BRITE": 8,
@@-    "KEGG_TC": 9, "CAZy": 10, "BiGG_Reaction": 11, "PFAMs": 12,
@@-}
+CONFIDENCE_SCORED_FIELDS: tuple[str, ...] = (
+    # ORDER IS THE CONTRACT with the annotation_confidence column emitted by
+    # emapper 3.0.0-beta6+ (custom build). Every scored field occupies exactly
+    # one position; '-' means "not scored for this row".
+    "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
+    "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
+    "BiGG_Reaction", "PFAMs",
+)
+CONFIDENCE_FIELD_INDEX = {name: i for i, name in enumerate(CONFIDENCE_SCORED_FIELDS)}
+assert len(CONFIDENCE_FIELD_INDEX) == len(CONFIDENCE_SCORED_FIELDS) == 13
```

2. **Record the contract in every artifact** so downstream consumers can detect drift:

```diff
@@     stats.update({
@@         ...
+        "confidence_field_order": list(CONFIDENCE_SCORED_FIELDS),
+        "confidence_source": "emapper annotation_confidence column (13 codes)",
@@         "_context_report": {
@@             "metadata": {
+                "confidence_field_order": list(CONFIDENCE_SCORED_FIELDS),
```

3. **Cross-check against the actual header** whenever the scored columns are present, so a reordered producer fails loudly instead of silently:

```diff
@@ def parse_eggnog_tsv(...):
     ...
     if columns is None:
         raise MergeError("eggNOG TSV is missing the #query header")
+    scored_present = [c for c in columns if c in CONFIDENCE_FIELD_INDEX]
+    header_order = [c for c in columns if c in CONFIDENCE_FIELD_INDEX]
+    expected_prefix = [c for c in CONFIDENCE_SCORED_FIELDS if c in header_order]
+    if header_order != expected_prefix:
+        raise MergeError(
+            "eggNOG scored columns appear in an order that conflicts with the "
+            "annotation_confidence position contract "
+            f"(header: {header_order}); refusing to apply confidence codes to "
+            "the wrong fields"
+        )
```

(Note: emapper's canonical column order is *not* the confidence order — `Preferred_name` precedes `GOs`/`EC`/`KEGG_ko`, but `KEGG_Pathway`…`PFAMs` follow `eggNOG_OGs`-adjacent columns. The check above must therefore assert *relative order compatibility*, not strict equality; simplest correct form: verify the header's relative order of scored fields is a subsequence of `CONFIDENCE_SCORED_FIELDS`. If the producer's header order genuinely differs from the confidence order, the correct check is subsequence-validity, and a mismatch is exactly the drift you want to catch.)

4. **v2 compatibility mode** so stock emapper 2.1.x files are usable (confidence un-gated, everything manifested as `not_scored`):

```diff
+parser.add_argument(
+    "--no-confidence-column",
+    action="store_true",
+    help="accept tables without annotation_confidence; all fields treated as "
+         "not_scored and min-confidence gating is disabled",
+)
```

with `_build_hit` making `annotation_confidence` conditionally required and `_confidence_passes` returning `True` when the column was absent. Keep the default strict — but document the producer requirement in the README ("which emapper build emits this column, where to get it").

---

## 6. P2 Findings (contract violations, performance, parser hazards)

### 6.1 F-04 — `--no-feature-provenance` ignored by the KofamScan planner

```
per-hit /note present despite --no-feature-provenance: True
per-CDS /inference present despite --no-feature-provenance: True
```

`plan_kofam_additions` has no `add_feature_provenance` parameter; `merge()` forwards the flag only to `plan_baktfold_additions`. The README states the flag "omits per-feature provenance" without scoping it to Baktfold.

**Semantics decision needed.** The Kofam `/note` carries threshold/score/E-value — that is *evidence*, not mere provenance — while `/inference="profile:KofamScan[:version]"` is pure provenance. Proposed resolution (matches README wording and the graft module's precedent):

```diff
--- a/merge_kofamscan_bakta.py
+++ b/merge_kofamscan_bakta.py
@@ def plan_kofam_additions(
@@     kofamscan_version: str | None,
@@     add_comment_note: bool,
+    add_feature_provenance: bool = True,
@@     merge_timestamp: str | None,
@@     starting_order: int = 0,
 ):
@@-        inference_values = planned[(feature.start, "inference")]
-        if hit.query_id not in inferred_queries and inference not in inference_values:
+        inference_values = planned[(feature.start, "inference")]
+        if (
+            add_feature_provenance
+            and hit.query_id not in inferred_queries
+            and inference not in inference_values
+        ):
             insertions.append(...)

@@ def merge(...):
@@     kofam_insertions, evidence_rows, kofam_stats = plan_kofam_additions(
@@         base, proteins, hits,
@@         kofam_data=kofam_data,
@@         faa_data=faa_data,
@@         kofamscan_version=kofamscan_version,
+        add_feature_provenance=add_feature_provenance,
@@         add_comment_note=add_comment_note,
```

If the maintainer prefers suppressing the score notes too, add a separate `--no-kofam-hit-notes` rather than overloading the provenance flag — the score note is the only record of *why* a KO was assigned and should not be silently droppable under a "provenance" name. Also mirror the parameter in `enrich_bakta.enrich()`.

### 6.2 F-05 — Unified manifest drops `faa_sha256` and `baktfold_version`

```
unified manifest metadata keys: [base_sha256, eggnog_sha256, eggnog_version,
  exact_duplicates_collapsed, gene_conflict_policy, gene_conflicts,
  input_insertions, insertions, kofam_sha256, kofamscan_version,
  merge_timestamp, min_eggnog_confidence, operation, output_sha256,
  reconciled_insertions, records]
faa_sha256 present: False
```

Every single-source manifest records the FAA hash; the unified manifest — the one most likely to be archived as the canonical provenance record — cannot prove which FAA produced the KO/eggNOG evidence. Fix together with the F-09 single-read refactor:

```diff
--- a/enrich_bakta.py
+++ b/enrich_bakta.py
@@-    base_data = bakta_path.read_bytes()
+    base_data = bakta_path.read_bytes()
+    faa_data = faa_path.read_bytes() if faa_path else b""
@@-    proteins = parse_faa(faa_path.read_bytes()) if faa_path else None
+    proteins = parse_faa(faa_data) if faa_path else None
@@     metadata: dict[str, Any] = {
         "operation": "unified-enrichment",
         "merge_timestamp": merge_timestamp or "",
+        "faa_sha256": sha256_bytes(faa_data) if faa_path else "",
     }
@@         metadata["baktfold_sha256"] = stats["baktfold_sha256"]
+        metadata["baktfold_version"] = stats["baktfold_version"]
@@-            faa_data=faa_path.read_bytes() if faa_path else b"",
+            faa_data=faa_data,
@@-            faa_data=faa_path.read_bytes() if faa_path else b"",
+            faa_data=faa_data,
@@         "eggnog_version": stats["eggnog_version"],
+        "tool_version": __version__,
```

(Add a module-level `__version__ = "0.2.0"` — see §9.3. The last two hunks remove the F-09 double reads.)

### 6.3 F-06 — Manifest bloat from per-row `raw_fields` duplication

Measured: 200 proteins × ~7 candidates → **1.62 MB manifest (8,090 bytes/protein)**. At C14/SM scale (~4.6k CDS), expect ~35–40 MB JSON per run, most of it the same `raw_fields` dict repeated for every candidate and every context field of the same query (the `seed_ortholog` key alone appeared thousands of times). Proposal — emit raw rows once, reference by row number:

```diff
@@ def plan_eggnog_additions(...):
+    for hit in table.hits:
+        evidence.append({
+            "entry_type": "eggnog_row",
+            "row_number": hit.row_number,
+            "query_id": hit.query_id,
+            "raw_fields": dict(hit.raw_fields),
+        })
     for hit in table.hits:
         ...
         evidence.append({
             ...
-            "raw_fields": dict(hit.raw_fields),
+            # raw fields are emitted once per row as eggnog_row entries;
+            # join on row_number
```

Apply the same removal in `_collect_context`'s manifest rows. Keep `seed_ortholog`, `e_value`, `score` inline (they're small and hot). Estimated size reduction ≈ 70–85% at genome scale. Gate behind `--manifest-raw-rows {once,always,never}` (default `once`) if backward compatibility matters.

### 6.4 F-07 — `_paired_gene` is O(hits × features): 34.3 s at 4,600 CDS

```
features: 9200, CDS hits simulated: 4600
_paired_gene over 4600 hits: 34.34s (7465 us/hit)
```

`_paired_gene` linearly rescans every feature of the record for *every* eggNOG row. With the repo's own C14 data (9,644 features, 4,563 rows) this alone burns ~30–40 s per merge, and it is pure waste. Pre-index:

```diff
--- a/merge_eggnog_bakta.py
+++ b/merge_eggnog_bakta.py
@@-def _paired_gene(base: RawDocument, cds: RawFeature) -> RawFeature | None:
-    candidates = [
-        feature
-        for feature in base.records[cds.record_index].features
-        if feature.feature_type == "gene"
-        and feature.locus_tag == cds.locus_tag
-        and feature.location == cds.location
-    ]
-    if len(candidates) != 1:
-        return None
-    return candidates[0]
+def _canon_location(location: str) -> str:
+    """Compare locations modulo whitespace so '<1..456' vs '1..456' style
+    formatting drift does not silently break gene/CDS pairing."""
+    return re.sub(r"\s+", "", location)
+
+
+def _index_paired_genes(base: RawDocument) -> dict[tuple[int, str, str], RawFeature | None]:
+    """One pass. Value is None when the key is ambiguous (mirrors the old
+    'len(candidates) != 1 -> None' semantics)."""
+    index: dict[tuple[int, str, str], RawFeature | None] = {}
+    for record in base.records:
+        for feature in record.features:
+            if feature.feature_type != "gene" or not feature.locus_tag:
+                continue
+            key = (record.index, feature.locus_tag, _canon_location(feature.location))
+            if key in index:
+                index[key] = None  # ambiguous
+            else:
+                index[key] = feature
+    return index
+
+
+def _paired_gene(
+    index: dict[tuple[int, str, str], RawFeature | None], cds: RawFeature
+) -> RawFeature | None:
+    return index.get((cds.record_index, cds.locus_tag, _canon_location(cds.location)))
@@ def plan_eggnog_additions(...):
     cds, validation = validate_faa_gbff(...)
+    paired_genes = _index_paired_genes(base)
@@-        paired = _paired_gene(base, feature)
+        paired = _paired_gene(paired_genes, feature)
```

Note the added `_canon_location` also fixes a robustness wart: raw string comparison of locations breaks on whitespace or `<`/`>` boundary formatting differences between the `gene` and `CDS` lines, silently downgrading `Preferred_name` to `unpaired_gene`. Expected speedup: O(N·F) → O(N+F), ~34 s → <0.1 s at the benchmarked scale.

### 6.5 F-08 — Qualifier wrapping can start a continuation line with `/`

```
                     /note="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                     /bogus_key=oopscccccccccc"
re-parsed qualifiers: {'note': [...], 'bogus_key': ['oopscccccccccc"']}
=> wrapped value is re-parsed as TWO qualifiers
```

`format_qualifier` wraps at 59 columns with `break_long_words=True`. If the character at the wrap boundary is `/`, the continuation line begins with `                     /…`, which both this repo's `_parse_qualifiers` and Biopython re-parse as a **new qualifier key**. Consequences: (a) `_existing_values`/idempotence checks miss the value on re-merge → duplicate insertion; (b) the merged file's qualifier set is silently wrong; (c) `validate_genbank_semantics` passes because the file is still parseable. Long `eggNOG_OGs` tokens (matched by the permissive `ANNOTATION_TOKEN_RE = r"[^\s,]+\Z"`, which admits `/`) are the realistic trigger. Fix — never break immediately before a `/`:

```diff
--- a/merge_engine.py
+++ b/merge_engine.py
@@ def format_qualifier(key: str, value: str, newline: bytes) -> bytes:
     if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
         raise MergeError(f"invalid qualifier name: {key!r}")
+    if any(character in value for character in "\r\n\t"):
+        raise MergeError(
+            f"qualifier /{key} value contains control characters: {value[:48]!r}"
+        )
     escaped = value.replace('"', '""')
     payload = f'/{key}="{escaped}"'
     width = 80 - len(QUALIFIER_INDENT)
     chunks = (
         [payload]
         if key == "inference"
-        else textwrap.wrap(
-            payload, width=width, break_long_words=True,
-            break_on_hyphens=False, replace_whitespace=False,
-            drop_whitespace=True,
-        )
-        or [payload]
+        else _wrap_qualifier_payload(payload, width)
     )
+
+def _wrap_qualifier_payload(payload: str, width: int) -> list[str]:
+    """Fixed-width wrap that never lets a continuation line begin with '/'.
+
+    A '/' at column 0 of a continuation line re-parses as a new qualifier
+    key (both in this engine and in Biopython), corrupting the qualifier
+    set and defeating idempotence checks. Shrink the cut by one whenever
+    the char that would start the next line is '/'.
+    """
+    if len(payload) <= width:
+        return [payload]
+    chunks: list[str] = []
+    rest = payload
+    while len(rest) > width:
+        cut = width
+        while cut > 1 and rest[cut] == "/":
+            cut -= 1
+        chunks.append(rest[:cut])
+        rest = rest[cut:]
+    if rest:
+        chunks.append(rest)
+    return chunks
```

The control-character guard closes the related hole that a `\n`/`\t` inside a value with `replace_whitespace=False` would produce structurally invalid GenBank. **Compatibility note:** for space-free payloads (all values this codebase currently emits), fixed-width chunking is byte-identical to the old `textwrap.wrap(break_long_words=True)` output except at `/` boundaries — i.e., only previously-corrupt outputs change. Add regression test T-03 asserting every wrapped payload re-parses to exactly one qualifier.

### 6.6 F-09 — Redundant I/O and memory profile

One unified run performs 9 `read_bytes` calls (`base.faa` ×3, `base.gbff` ×2, `ann.tsv` ×2, plus one each for baktfold/kofam) and holds `output` + `reconstructed` bytearrays + base bytes simultaneously (~3× peak) in `apply_insertions`. The Diff A/E changes remove the duplicate base read and both duplicate FAA reads; `parse_eggnog_path` should also accept pre-read bytes for TSV to remove the duplicate eggNOG read:

```diff
--- a/merge_eggnog_bakta.py
+++ b/merge_eggnog_bakta.py
@@ def parse_eggnog_path(path: Path, *, expected_version: str | None = None,
+                        data: bytes | None = None) -> EggnogTable:
     if path.suffix.lower() == ".xlsx":
         return parse_eggnog_xlsx(path, expected_version=expected_version)
-    return parse_eggnog_tsv(path.read_bytes(), expected_version=expected_version)
+    return parse_eggnog_tsv(data if data is not None else path.read_bytes(),
+                            expected_version=expected_version)
```

For very large inputs, `apply_insertions` could drop the `reconstructed` bytearray and verify reversibility by re-walking `applied` offsets against `base_data` directly (O(1) extra memory); keep the current version if simplicity is preferred — 3× a bacterial GBFF is fine, but document the ceiling for eukaryotic-scale inputs.

### 6.7 F-10 — Query IDs starting with `#` silently vanish

```
parsed hits: ['T_0001']
'#weird' silently dropped: True
```

`parse_eggnog_tsv` treats every `#`-leading line as a comment *after* the header has been seen, so a legitimate query id `#weird` is dropped without any manifest entry — an audit-trail hole (the FAA contains the protein, eggNOG annotated it, the merge pretends it never existed). Kofam has the mirror-image problem: a query id beginning with `*` would be consumed as the hit marker. Fix:

```diff
--- a/merge_eggnog_bakta.py
+++ b/merge_eggnog_bakta.py
@@         if line.startswith("#"):
+            if columns is not None and "\t" in line:
+                raise MergeError(
+                    f"eggNOG row {row_number}: data row begins with '#'; "
+                    "query IDs must not start with '#'"
+                )
             candidate = line.split("\t")
```

```diff
--- a/merge_kofamscan_bakta.py
+++ b/merge_kofamscan_bakta.py
@@         fields = line[1:].strip().split(maxsplit=5)
         if len(fields) != 6:
             raise MergeError(...)
         query_id, ko, threshold, score, e_value, definition = fields
+        if query_id.startswith(("*", "#")):
+            raise MergeError(
+                f"Kofam row {row_number}: query ID {query_id!r} starts with a "
+                "reserved marker ('*' or '#')"
+            )
```

---

## 7. P3 Findings (robustness, consistency, UX)

These are individually small but collectively erode the project's "strict validation, loud failures" doctrine. Each includes a minimal fix.

### 7.1 F-11 — Context report smuggled through the stats dict

`plan_eggnog_additions` returns `(insertions, evidence, stats)` where `stats["_context_report"]` is a full sidecar payload that **both** callers must remember to `.pop()`. Any third caller (or a future `--json` dump of stats) will leak the entire report into the manifest/stats output. Replace with an explicit result object:

```diff
--- a/merge_eggnog_bakta.py
+++ b/merge_eggnog_bakta.py
+from dataclasses import dataclass, field
+
+@dataclass(frozen=True)
+class EggnogPlan:
+    insertions: list[Insertion]
+    evidence_rows: list[dict[str, Any]]
+    stats: dict[str, Any]
+    context_report: dict[str, Any]
+
@@ def plan_eggnog_additions(...) -> EggnogPlan:
@@-    stats.update({... "_context_report": {...}})
+    context_report = {...}
+    return EggnogPlan(insertions, evidence, stats, context_report)
```

Callers (`merge_eggnog_bakta.merge`, `enrich_bakta.enrich`, `restore` does not use it) destructure explicitly; the awkward `context_report if context_report_path else None` / undefined-variable dance in `enrich_bakta.py:119,145` disappears.

### 7.2 F-12 — `restore_bakta_translations` skips the identity doctrine for translated CDSs

`plan_translation_restoration` checks query existence in FAA/GBFF but never that the FAA sequence equals an existing `/translation`. A corrupted FAA (or a stale one from a re-annotation) sails through `restore` and then fails — or worse, is *not even checked* — downstream. Since the whole point of the restore→merge workflow is exact FAA/GBFF identity, restore should run the same `validate_faa_gbff` gate for queries it does *not* restore:

```diff
@@ def plan_translation_restoration(...):
@@     if missing_faa or missing_cds:
         raise MergeError(...)
+    # Queries whose CDS already carries a translation must still satisfy the
+    # FAA/GBFF identity contract — otherwise 'restore' certifies a file that
+    # the subsequent merge will reject (or, if skipped, never verifies).
+    translated_queries = [
+        q for q in sorted(requested)
+        if any(v.strip() for v in cds[q].values("translation"))
+    ]
+    if translated_queries:
+        validate_faa_gbff(
+            base, proteins, translated_queries,
+            source_name="eggNOG (translation restoration)",
+        )
```

Additionally, `--dry-run` currently *requires* the positional `output` path even though nothing is written (`python restore_bakta_translations.py BAKTA.gbff BAKTA.faa query.tsv RESTORED.gbff --dry-run`). Change `output` to `nargs="?"` (error only when not `--dry-run`).

### 7.3 F-13 — Small validation gaps

1. **Kofam accepts negative threshold/score** (only the E-value sign is checked):
```diff
@@-        _decimal(threshold, "threshold", row_number)
-        _decimal(score, "score", row_number)
+        if _decimal(threshold, "threshold", row_number) < 0:
+            raise MergeError(f"Kofam row {row_number}: negative threshold {threshold!r}")
+        if _decimal(score, "score", row_number) < 0:
+            raise MergeError(f"Kofam row {row_number}: negative score {score!r}")
```

2. **`protein_sha256` leaks a raw `ValueError`** on non-ASCII sequences (crash with traceback + exit 1 instead of a clean `MergeError` exit 2):
```diff
--- a/merge_engine.py
+++ b/merge_engine.py
@@ def protein_sha256(sequence: str) -> str:
-    return sha256_bytes(normalize_protein(sequence).encode("ascii"))
+    normalized = normalize_protein(sequence)
+    try:
+        return sha256_bytes(normalized.encode("ascii"))
+    except UnicodeEncodeError as exc:
+        raise MergeError(
+            f"protein sequence contains non-ASCII characters: {normalized[:32]!r}"
+        ) from exc
```

3. **XLSX empty header names accepted** (a trailing empty header cell creates a nameless column that flows into `raw_fields`); truncated data rows then surface as a confusing "must contain 13 l/m/h/- codes" error:
```diff
@@     columns = [str(cell.value) if cell.value is not None else "" for cell in header_cells]
     if len(columns) != len(set(columns)):
         raise MergeError("eggNOG XLSX header contains duplicate columns")
+    if any(not name.strip() for name in columns):
+        raise MergeError("eggNOG XLSX header contains empty column names")
```

4. **Misleading stat**: the eggNOG COMMENT line writes `mapped CDSs={len(cds)}` — that is the *total* CDS count, not mapped queries. Use `len({hit.query_id for hit in table.hits})` or rename to `CDSs=`. Same class of issue: `stats["hit_cds"]` is `len(emitted_by_query)` (CDSs that emitted something), while the kofam planner's `hit_cds` means CDSs with ≥1 hit — same key, two meanings across modules.

5. **Unified CLI missing `--clean-gene-suffix`** (single-source eggNOG CLI has it). Add the flag to `enrich_bakta.py` and thread through to `plan_eggnog_additions`; record `clean_gene_suffix: bool` in unified metadata.

6. **Inconsistent `--json` shapes**: `graft()` returns flat stats; kofam/eggnog `merge()` return `{"kofam": ..., **final}` / `{"eggnog": ..., **final}`. Pick one convention (suggest always `{"<source>": ..., **final}`) and document it.

7. **`newline_for_offset` fallback** returns the *first* newline style found anywhere in the document when the insertion point has no preceding newline — for mixed-EOL files this can emit the wrong style locally. Rare, but trivially improved by scanning backwards from `offset` (bounded) instead of from position 0.

8. **`comment_insertion` + `CONTIG`**: records carrying a `CONTIG` header key (WGS-style) would get COMMENT continuation lines inserted after CONTIG, mis-associating them. Bakta complete-genome output has no CONTIG, so this is latent; guard by refusing comment insertion when `re.search(rb"(?m)^CONTIG\b", record_prefix)` and falling back to a `NOTE`-style manifest-only provenance.

9. **`detect_baktfold_version`** silently returns `"unknown"`, degrading all Baktfold provenance strings, when the version regex misses. Record `baktfold_version_detected: bool` in stats so silent provenance loss is visible in the manifest.

### 7.4 F-14 — Engineering infrastructure

- **No `pyproject.toml`.** Dependencies (biopython, pytest, ruff; optional openpyxl) are unpinned and undeclared; tests only run via `python -m pytest` from the repo root. The repo's own review doc acknowledges this; nothing has landed.
- **Biopython is unpinned**, yet `validate_genbank_semantics` is the only backstop between a splice error and a corrupted file, and Biopython's GenBank parser behavior (what warns vs. raises) shifts across minor versions. Pin a floor (and ideally a tested ceiling), e.g. `biopython>=1.83,<1.90`.
- **No CI.** The README's validation commands are folklore; nothing enforces them per-commit.
- **No tool version stamping.** Manifests record `operation` but not which version of enrich-bakta produced them — an audit gap for a tool whose entire raison d'être is auditability.
- **Ruff scope mismatch**: README lints 5 files; `normalize_baktfold.py` (and `restore_bakta_translations.py`) escape both `ruff check` and `ruff format --check`. The 7× `EXE001` findings also indicate the shebang--bearing entry points lost their executable bits (Windows checkout) — either set the bit via `git update-index --chmod=+x` or drop the shebangs and route through `python -m`.
- **`normalize_baktfold.py` (legacy)** violates every standard the rest of the repo sets: non-atomic write, no input validation, no collision checks, platform-dependent default encoding (`open("r")` with no `encoding=`), fabricated hardcoded dates (`--date 08-JUL-2026`, `--annotation-date "07/08/2026, 11:50:32"` — the latter embeds a specific run's timestamp as a *default*, which is provenance fabrication), and a docstring claiming "streaming" while buffering the whole output in memory. It is excluded from the byte-preserving path and should be quarantined (`legacy/` directory, README warning) or deleted.

---

## 8. Architecture & Workflow Review (non-blocking observations)

1. **Two ordering strategies coexist.** `enrich_bakta.py` allocates magic order blocks (baktfold `0`, kofam `1_000_000`, eggnog `2_000_000`); `merge_kofamscan_bakta.merge()` computes `max(order)+1`. Both are deterministic today, but the magic blocks silently assume no source ever plans >1M insertions. Replace `order: int` with `order: tuple[int, int]` (source sequence, within-source counter) or at least assert the invariant (`assert max_baktfold_order < 1_000_000`).
2. **The one-pass Baktfold+Kofam path skips `reconcile_insertions`** and is safe only because those two sources cannot emit identical `(qualifier, value)` pairs by construction (Baktfold: gene/EC_number/structural db_xrefs; Kofam: KEGG db_xrefs/score notes/inference). That invariant is nowhere asserted or documented. Either call `reconcile_insertions` in `merge_kofamscan_bakta.merge()` too (cheap, uniform, future-proof), or add a comment + test pinning the invariant.
3. **`validate_faa_gbff` runs twice in unified mode** (once per source) rebuilding `cds_by_locus` and re-comparing translations. Acceptable at bacterial scale; worth caching if eukaryotic support is ever desired.
4. **`validate_genbank_semantics` costs a full Biopython parse** per input/output (3+ parses per unified run) and — as demonstrated in §3 — Biopython *warnings* do not fail the merge. Harden it by capturing warnings:
```python
import warnings
from Bio import BiopythonParserWarning
with warnings.catch_warnings():
    warnings.simplefilter("error", BiopythonParserWarning)
    records = list(SeqIO.parse(io.StringIO(text), "genbank"))
```
This single change would have turned the §3 corruption demo into a hard failure (the malformed-LOCUS and unescaped-quote warnings both fired there).
5. **Output-before-manifest write ordering**: `finalize_merge` promotes the output, then writes manifest/sidecar. A manifest write failure leaves an enriched GBFF with no manifest — an audit-trail orphan. Consider staging all artifacts and promoting manifest first (its content is fully determined before the output write), or at least documenting the partial-failure semantics.
6. **`atomic_write` lacks a directory fsync** — after `os.replace`, the rename itself may not be durable across a crash on some filesystems. One-line hardening:
```python
dir_fd = os.open(path.parent, os.O_RDONLY)
try:
    os.fsync(dir_fd)
finally:
    os.close(dir_fd)
```
7. **Whole-file in-memory model**: the engine holds base + output + reconstruction (~3× input). Fine for bacterial GBFFs (tens of MB); document the ceiling, or adopt the §6.6 O(1)-memory reconstruction check, before anyone points this at a 3-GB eukaryotic file.
8. **`enrich()` error messaging** for the common pseudogene case: when `validate_faa_gbff` reports `missing_or_ambiguous_translation`, the error should point users at `restore_bakta_translations.py` — today it's a bare JSON blob with no remediation hint.

---

## 9. Proposed Additional Schemas & Logics

### 9.1 JSON Schema for `enrich-bakta.merge-manifest.v1`

Commit as `schemas/merge-manifest.v1.schema.json`; validate in a test that every writer produces conforming manifests. Sketch:

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://github.com/WhyAdr/enrich-bakta/schemas/merge-manifest.v1.schema.json",
  "title": "enrich-bakta merge manifest v1",
  "type": "object",
  "required": ["schema", "metadata", "entries"],
  "properties": {
    "schema": {"const": "enrich-bakta.merge-manifest.v1"},
    "metadata": {
      "type": "object",
      "required": ["operation", "base_sha256", "output_sha256", "insertions",
                   "records", "merge_timestamp", "tool_version"],
      "properties": {
        "operation": {"enum": ["baktfold-graft", "kofam-merge",
                                "baktfold-kofam-one-pass", "eggnog-merge",
                                "unified-enrichment",
                                "pseudogene-translation-restoration"]},
        "base_sha256": {"pattern": "^[0-9a-f]{64}$"},
        "output_sha256": {"pattern": "^[0-9a-f]{64}$"},
        "faa_sha256": {"pattern": "^[0-9a-f]{64}$"},
        "baktfold_sha256": {"pattern": "^[0-9a-f]{64}$"},
        "kofam_sha256": {"pattern": "^[0-9a-f]{64}$"},
        "eggnog_sha256": {"pattern": "^[0-9a-f]{64}$"},
        "insertions": {"type": "integer", "minimum": 0},
        "records": {"type": "integer", "minimum": 1},
        "merge_timestamp": {"type": "string"},
        "tool_version": {"type": "string"},
        "gene_conflict_policy": {"enum": ["skip", "prefer-eggnog", "prefer-baktfold"]},
        "min_eggnog_confidence": {"enum": ["low", "medium", "high"]},
        "confidence_field_order": {"type": "array", "items": {"type": "string"}}
      }
    },
    "entries": {
      "type": "array",
      "items": {
        "type": "object",
        "required": ["entry_type"],
        "oneOf": [
          {"properties": {"entry_type": {"const": "insertion"},
                          "qualifier": {"enum": ["gene", "EC_number", "db_xref",
                                                  "note", "inference", "COMMENT",
                                                  "translation"]}}},
          {"properties": {"entry_type": {"const": "kofam_hit"}}},
          {"properties": {"entry_type": {"const": "eggnog_candidate"},
                          "status": {"enum": ["emitted", "existing", "existing_gene",
                                               "unpaired_gene", "filtered_confidence",
                                               "skipped_partial_ec"]}}},
          {"properties": {"entry_type": {"const": "eggnog_context"},
                          "status": {"const": "sidecar_only"}}},
          {"properties": {"entry_type": {"const": "eggnog_row"}}},
          {"properties": {"entry_type": {"const": "reconciliation"}}},
          {"properties": {"entry_type": {"const": "translation_restoration"}}}
        ]
      }
    }
  }
}
```

### 9.2 JSON Schema for `enrich-bakta.eggnog-context.v1`

Same pattern: `schema` const, `metadata` (operation, versions, hashes, `confidence_field_order`), `entries[]` with `query_id`/`locus_tag`/`context{field → {raw, values[], confidence_code, confidence_status}}`.

### 9.3 Tool-version stamping (new logic)

```python
# merge_engine.py
__version__ = "0.2.0"
TOOL_IDENTIFIER = f"enrich-bakta {__version__}"
```

Emit `metadata["tool_version"] = __version__` in `finalize_merge` (single choke point) and `--version` on every CLI. This makes manifests self-describing across future format changes (`merge-manifest.v2`, …), which the current design has no way to signal.

### 9.4 Packaging & CI (new logic)

```toml
# pyproject.toml
[project]
name = "enrich-bakta"
version = "0.2.0"
requires-python = ">=3.10"
dependencies = ["biopython>=1.83,<1.90"]

[project.optional-dependencies]
xlsx = ["openpyxl>=3.1"]
dev = ["pytest>=8", "ruff>=0.6", "mypy>=1.10"]

[tool.ruff]
line-length = 100

[tool.ruff.lint]
select = ["E", "F", "W", "I", "UP", "B", "SIM", "EXE"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

```yaml
# .github/workflows/ci.yml
name: ci
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    strategy: {matrix: {python: ["3.10", "3.12"]}}
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: {python-version: "${{ matrix.python }}"}
      - run: pip install -e .[dev]
      - run: ruff check . && ruff format --check .
      - run: mypy --ignore-missing-imports .
      - run: python -m pytest -q
```

---

## 10. Test-Suite Gaps (ordered by value)

| ID | Proposed test | Catches |
|----|--------------|---------|
| T-01 | `test_finalize_uses_planned_bytes_or_aborts` — mutate base file between planning and finalize | F-01 (P0) |
| T-02 | `test_partial_ec_skipped_not_fatal` — row with `1.2.3.-` merges; manifest row `skipped_partial_ec`; no `/EC_number` | F-02 |
| T-03 | `test_wrapped_qualifier_reparses_as_single_value` — property: for every emitted insertion, re-parsing the output yields exactly the planned qualifier/value | F-08 + future wrap regressions |
| T-04 | `test_unified_manifest_records_all_input_hashes` — faa/baktfold/kofam/eggnog hashes + versions present in unified metadata | F-05 |
| T-05 | `test_kofam_provenance_flag` — `add_feature_provenance=False` suppresses `/inference` (and, per chosen semantics, notes) | F-04 |
| T-06 | `test_hash_and_star_query_ids_rejected` | F-10 |
| T-07 | `test_biopython_parser_warnings_fail_validation` — a GBFF that triggers `BiopythonParserWarning` must raise | hardens the §3 backstop |
| T-08 | `test_unified_three_source_idempotence` — enrich(enrich(x)) == enrich(x), 0 insertions on rerun | (verified manually in this audit; pin it) |
| T-09 | `test_manifest_conforms_to_schema` — validate writer output against `schemas/*.schema.json` (needs `jsonschema` dev dep) | schema drift |
| T-10 | `test_manifest_size_budget` — N=200 fixture must stay under, say, 1 MB | F-06 |
| T-11 | `test_paired_gene_index_equivalence` — random shuffled features; index path must equal naive scan | F-07 refactor |
| T-12 | `test_xlsx_empty_header_rejected` | F-13.3 |

The suite is otherwise in good shape: it covers parity failures, CRLF, allowlist rejection, idempotence (2-source), confidence gating, reconciliation, and restoration constraints.

---

## 11. Prioritized Handoff Checklist for Gemini

**Phase 1 — correctness (do first, independently testable)**
1. F-01 / Diff A: `finalize_merge(base_data=...)` + 5 call sites + T-01.
2. F-08 / Diff G: boundary-safe qualifier wrapping + control-char guard + T-03.
3. §8.4: promote `BiopythonParserWarning` to errors inside `validate_genbank_semantics` + T-07.

**Phase 2 — data compatibility**
4. F-02 / Diff B: partial-EC tolerance + T-02.
5. F-03 / Diff C: explicit confidence contract, manifest recording, header order check, `--no-confidence-column` v2 mode; document the exact emapper producer in the README.
6. F-10 / Diff: reserved-marker query IDs + T-06.

**Phase 3 — contracts & performance**
7. F-05 + F-09: single-read refactor; unified metadata completeness (`faa_sha256`, `baktfold_version`, `tool_version`) + T-04.
8. F-07 / Diff F: paired-gene index + T-11.
9. F-06: `eggnog_row` raw-field dedup (+ optional `--manifest-raw-rows`) + T-10.
10. F-04 / Diff D: kofam provenance flag semantics + T-05.

**Phase 4 — hygiene & infrastructure**
11. F-11: `EggnogPlan` dataclass (removes `_context_report`).
12. F-12/F-13 bundle: restore identity gate, `--dry-run` output arg, negative-score checks, `protein_sha256` wrap, XLSX header check, stat fixes, `--clean-gene-suffix` parity, JSON shape unification.
13. F-14: `pyproject.toml`, pin biopython, CI workflow, ruff scope = all files, quarantine/delete `normalize_baktfold.py`.
14. Commit `schemas/*.schema.json` + T-09; add `__version__` stamping (§9.3).

**Explicitly out of scope for this pass** (flagged, deliberately not diffed): GFF3 output support; batch/cohort wrappers; parallelizing per-record parsing. All three are worthwhile but should land after Phase 1–2 stabilize the engine contracts.

---

## Appendix A — Reproduction Artifacts

All repro scripts used in this audit are preserved under `/home/z/my-project/scripts/`:

| Script | Demonstrates |
|--------|--------------|
| `audit_repros.py` | F-02 (partial EC), F-03 (v2 header), F-04 (provenance flag), F-05 (missing hashes), F-06 (manifest size), F-08 (wrap `/`), F-10 (`#` query) |
| `audit_repros2.py` | F-01 precise corruption (`translation` polluted while `self_check=True`) and exact wrap-boundary repro |
| `audit_read_counts.py` | F-09 — 9 `read_bytes` calls per unified run |
| `audit_perf.py` | F-07 — 34.34 s `_paired_gene` at 4,600 CDS × 9,200 features |
| `audit_idempotence.py` | 3-source unified idempotence (pass1==pass2==pass3, rerun insertions=0) |
| `audit_xlsx.py` | F-13.3 — empty header names / sparse-row behavior |

Environment: Python 3.12.13, Biopython 1.86, ruff 0.16.4, mypy (latest), Linux x86_64.
