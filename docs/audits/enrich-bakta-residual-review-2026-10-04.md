# enrich-bakta: residual review and execution handoff

Prepared for Wahyu, 2026-10-04 (Asia/Jakarta).

Reviewed remote main: `84fdc49e69648f2563dd8de05e696d84c40f558f`.
Comparison baseline: `cc98e4f081c2e0554a7108b3897f6679fe61cd1d`.

The nine earlier fixes landed. The subsequent Baktfold candidate, explicit-role,
class-specific support, parsed-output validation, and checked-in schema work
substantially covers the design follow-through from the previous handoff.
This review found three reproducible residual manifest defects and one remaining
release-validation gap. No remote changes were made.

## Confirmed residual defects

| ID | Reproduction on current main | Required behavior |
|---|---|---|
| F01 | Baktfold `xyz` conflicts with eggNOG `abc`. Under default skip both decisions are `suppressed_conflict` but both say `reason_code: inserted`. Under either explicit preference the losing candidate also says `inserted`. | Final rejection must override the planner reason: `gene_conflict_no_preference` or `gene_conflict_source_preference`. Preserve the planner reason separately. |
| F02 | Both sources propose EC `1.2.3.4`. Decisions become `shared_support`, while corresponding source `entries.status` and `entries.final_status` stay `emitted`. | Synchronize source projections after shared-support reconciliation; link each projection to its final `candidate_id` and final reason. |
| F03 | An authoritative base gene rejection reports zero suppressed candidates despite one rejected node. A protein mismatch reports zero despite three rejected nodes. A source-pair conflict reports zero despite two rejected nodes. | `suppressed_candidate_count` must include every `suppressed_*` final status. Confidence filtering, skipped partial ECs, invalid values, and unresolved identities remain distinct categories. |

These are audit/manifest correctness defects. The controlled probes produced the
expected GBFF annotations. The draft changes only `core/decisions.py` and adds seven
regression cases; it changes neither annotation policy nor insertion generation.

## Verified evidence

- Current main: **133 tests passed** with Python 3.12 and Biopython 1.88;
  Ruff, formatting, mypy (14 source modules), and whitespace checks passed.
- Seven new regression cases fail against current main; the complete drafted
  implementation passes **140 tests** and the same static checks.
- Five controlled GBFF outputs (three conflict policies, an authoritative-name
  rejection, and a protein mismatch) are byte-identical before and after the draft.
- The complete embedded diff applies cleanly to the reviewed main.
- [GitHub CI run 37171609464](https://github.com/WhyAdr/enrich-bakta/actions/runs/37171609464)
  passed all nine active jobs: Linux/Windows Python 3.10 and 3.12 tests,
  Linux/Windows package jobs, XLSX, dependency floor, and manifest schema.
  Dataset bytes and scientific reruns were **skipped** on this push.
- No `workflow_dispatch` run was present when checked. The checked-in Windows
  acceptance record documents successful local C14/SM runs at `e5def78`; this
  review did not independently repeat those large genomic workflows.

## Remaining release-validation work

1. Strengthen `tools/run_followthrough_real_checks.py` and
   `tools/run_followthrough_chain_checks.py` before treating their success as a
   golden scientific regression gate. They currently check TSV/XLSX equality,
   byte-identical reruns, several origin checks, and opt-in enforcement, but
   mostly **report** first-pass output hashes/counts. They do not compare those
   first-pass results with the recorded accepted baseline. An altered but
   idempotent annotation result can therefore pass.
2. Check in a baseline JSON with the accepted input hashes, producer versions,
   policy flags, tool version, expected outputs, and counts from
   `docs/data/ACCEPTANCE.md`. Compare actual first-pass results against it, exit
   nonzero on divergence, and require an explicit reviewed baseline update for
   intentional changes. If tool-version comments change bytes, review and update
   the corresponding byte baseline deliberately or compare a defined semantic
   annotation baseline; do not silently ignore arbitrary hash drift.
3. Assert exact imported query IDs and `origin=imported_faa` on **every** downstream
   and no-op ledger, including eggNOG and unified reruns. Some current checks only
   return the ledger length, and some only check its count. Assert source
   projections, final reasons, and suppression totals on generated manifests;
   structural JSON Schema validation alone cannot enforce those relationships.
4. Run the CI workflow manually on the final revision so `dataset-bytes` and
   `scientific-reruns` execute. Inspect both reports and retain their evidence.
   Ensure release/tag decisions actually require this scientific validation;
   ordinary push CI currently excludes these jobs.
5. After these gates pass, choose the next release version and synchronize
   `pyproject.toml`, `TOOL_VERSION`, changelog, and release notes. Current hardening
   remains under `Unreleased` while the runtime version is still `0.3.0`.

The golden-baseline and broader real-ledger assertions above are execution
requirements, **not completed code in the embedded patch**. BK71A remains
structural-only until its pristine FAA/evidence inputs are available. Complete
entry-type/translation-ledger schema definitions would be useful further
hardening: the current schema types decisions strictly but accepts arbitrary
objects in `entries` and arbitrary additional metadata.

## Apply and verify the tested manifest patch

Use a clean branch or worktree; preserve unrelated work. If main has advanced,
review the new commits before adapting this diff.

```bash
git fetch origin main
git switch -c hardening/residual-manifest origin/main
git rev-parse HEAD
# Reviewed baseline: 84fdc49e69648f2563dd8de05e696d84c40f558f
```

Extract the complete patch below from this document:

```bash
python - /absolute/path/to/enrich-bakta-residual-review-2026-10-04.md /absolute/path/to/residual-manifest.patch <<'PY'
from pathlib import Path
import sys
text = Path(sys.argv[1]).read_text(encoding="utf-8")
start = text.index("```diff\n") + len("```diff\n")
end = text.index("\n```", start)
Path(sys.argv[2]).write_text(text[start:end] + "\n", encoding="utf-8")
PY
git apply --check /absolute/path/to/residual-manifest.patch
git apply /absolute/path/to/residual-manifest.patch
python -m pip install -e '.[dev,xlsx]'
python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m mypy
git diff --check
```

Expected suite result at this baseline: **140 passed**. Run remote CI after
committing; this draft itself has not run on Windows. Manifest consumers should
expect corrected suppression totals and final source-entry status/reason fields.
The source entries gain `candidate_id` and `planned_reason_code`; the manifest
schema already allows these entry properties.

## Complete tested patch

```diff
diff --git a/src/enrich_bakta_lib/core/decisions.py b/src/enrich_bakta_lib/core/decisions.py
index e6c7d9e..89fe5cf 100644
--- a/src/enrich_bakta_lib/core/decisions.py
+++ b/src/enrich_bakta_lib/core/decisions.py
@@ -184,6 +184,18 @@ def _reason_code(
     source_id: str,
     gene_conflict_policy: str,
 ) -> str:
+    # Reconciliation can change an adapter's inserted proposal into a rejection.
+    # Resolve that outcome before consulting its planned reason.
+    if final_status == "suppressed_conflict" and qualifier == "gene":
+        preferred = {
+            "prefer-eggnog": "eggNOG",
+            "prefer-baktfold": "Baktfold",
+        }.get(gene_conflict_policy)
+        return (
+            "gene_conflict_source_preference"
+            if preferred is not None and source_id != preferred
+            else "gene_conflict_no_preference"
+        )
     explicit = row.get("reason_code")
     if isinstance(explicit, str) and explicit:
         return explicit
@@ -199,16 +211,6 @@ def _reason_code(
         return "protein_identity_mismatch"
     if final_status == "suppressed_unsupported_pair":
         return "unsupported_pair"
-    if final_status == "suppressed_conflict" and qualifier == "gene":
-        preferred = {
-            "prefer-eggnog": "eggNOG",
-            "prefer-baktfold": "Baktfold",
-        }.get(gene_conflict_policy)
-        return (
-            "gene_conflict_source_preference"
-            if preferred is not None and source_id != preferred
-            else "gene_conflict_no_preference"
-        )
     return {
         "filtered_confidence": "confidence_below_threshold",
         "skipped_partial_ec": "partial_ec_not_promoted",
@@ -279,6 +281,7 @@ def build_candidate_ledger(
         )
     covered: set[int] = set()
     candidates: list[CandidateDecision] = []
+    source_rows_by_candidate: dict[str, dict[str, Any]] = {}
 
     def matching_insertions(
         source_id: str,
@@ -380,6 +383,7 @@ def build_candidate_ledger(
             gene_conflict_policy=gene_conflict_policy,
         )
         row["planned_status"] = planned_status
+        row["planned_reason_code"] = row.get("reason_code", "")
         row["planned_emitted_qualifiers"] = row.get("emitted_qualifiers", "")
         row["final_status"] = final_status
         row["status"] = final_status
@@ -395,9 +399,11 @@ def build_candidate_ledger(
             "normalized_value": normalized,
             "target_feature_uids": sorted(target_uids),
         }
+        candidate_id = stable_id("candidate", payload)
+        source_rows_by_candidate[candidate_id] = row
         candidates.append(
             CandidateDecision(
-                candidate_id=stable_id("candidate", payload),
+                candidate_id=candidate_id,
                 source_id=source_id,
                 source_sha256=source_sha,
                 target_feature_uids=tuple(sorted(target_uids)),
@@ -563,6 +569,14 @@ def build_candidate_ledger(
             }
             candidate.supporting_candidate_ids = tuple(sorted(supporting))
 
+    # Publish source projections only after shared-support reconciliation is done.
+    for candidate in candidates:
+        source_row = source_rows_by_candidate.get(candidate.candidate_id)
+        if source_row is not None:
+            source_row["candidate_id"] = candidate.candidate_id
+            source_row["status"] = candidate.final_status
+            source_row["final_status"] = candidate.final_status
+            source_row["reason_code"] = candidate.reason_code
     candidates.sort(key=lambda item: item.candidate_id)
     rows = [candidate.as_dict() for candidate in candidates]
     counts = {
@@ -574,7 +588,9 @@ def build_candidate_ledger(
             candidate.final_status in EMITTED_STATUSES for candidate in candidates
         ),
         "suppressed_candidate_count": sum(
-            candidate.final_status == "suppressed_conflict" for candidate in candidates
+            candidate.final_status is not None
+            and candidate.final_status.startswith("suppressed_")
+            for candidate in candidates
         ),
         "candidate_insertion_count": len(final),
         "candidate_emitting_feature_count": len(
diff --git a/tests/test_residual_manifest.py b/tests/test_residual_manifest.py
new file mode 100644
index 0000000..604c08c
--- /dev/null
+++ b/tests/test_residual_manifest.py
@@ -0,0 +1,118 @@
+from __future__ import annotations
+
+import json
+
+import pytest
+from test_merge_pipeline import eggnog_bytes, record_bytes, write
+
+from enrich_bakta_lib.workflows.enrich import enrich
+
+
+def run_conflict(tmp_path, policy):
+    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
+    fold = write(
+        tmp_path / "fold.gbff",
+        record_bytes("TEST", "T_0001", gene="xyz", ec_numbers=("1.2.3.4",)),
+    )
+    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
+    egg = write(
+        tmp_path / "egg.tsv",
+        eggnog_bytes(
+            "T_0001\tseed\t1e-4\t10\t-\tabc\t-\t1.2.3.4\tK00001\t-\thhhhhhhhhhhhh"
+        ),
+    )
+    manifest = tmp_path / "out.json"
+    enrich(
+        bakta_path=base,
+        faa_path=faa,
+        eggnog_path=egg,
+        baktfold_path=fold,
+        output_path=tmp_path / "out.gbff",
+        manifest_path=manifest,
+        gene_conflict_policy=policy,
+    )
+    return json.loads(manifest.read_text())
+
+
+@pytest.mark.parametrize("policy", ["skip", "prefer-eggnog", "prefer-baktfold"])
+def test_gene_reconciliation_reason_overrides_planned_insertion(tmp_path, policy):
+    payload = run_conflict(tmp_path, policy)
+    rejected = [
+        row
+        for row in payload["decisions"]
+        if row["qualifier"] == "gene" and row["final_status"] == "suppressed_conflict"
+    ]
+    assert len(rejected) == (2 if policy == "skip" else 1)
+    expected_reason = (
+        "gene_conflict_no_preference"
+        if policy == "skip"
+        else "gene_conflict_source_preference"
+    )
+    assert {row["reason_code"] for row in rejected} == {expected_reason}
+    by_id = {row["candidate_id"]: row for row in rejected}
+    projections = [
+        row for row in payload["entries"] if row.get("candidate_id") in by_id
+    ]
+    assert len(projections) == len(rejected)
+    assert all(row["reason_code"] == expected_reason for row in projections)
+    assert all(row["planned_reason_code"] == "inserted" for row in projections)
+
+
+def test_shared_support_source_projection_matches_final_decision(tmp_path):
+    payload = run_conflict(tmp_path, "skip")
+    shared = {
+        row["candidate_id"]: row
+        for row in payload["decisions"]
+        if row["qualifier"] == "EC_number" and row["final_status"] == "shared_support"
+    }
+    assert len(shared) == 2
+    projections = [
+        row
+        for row in payload["entries"]
+        if row.get("entry_type") in {"baktfold_candidate", "eggnog_candidate"}
+        and row.get("field") == "EC"
+    ]
+    assert len(projections) == 2
+    for row in projections:
+        decision = shared[row["candidate_id"]]
+        assert row["status"] == row["final_status"] == decision["final_status"]
+        assert row["reason_code"] == decision["reason_code"] == "shared_support"
+        assert row["planned_reason_code"] == "inserted"
+
+
+@pytest.mark.parametrize("case", ["authoritative", "protein-mismatch", "pair-conflict"])
+def test_suppression_total_includes_each_rejection_kind(tmp_path, case):
+    base_data = record_bytes(
+        "TEST", "T_0001", gene="original" if case == "authoritative" else None
+    )
+    fold_data = record_bytes(
+        "TEST", "T_0001", gene="xyz", ec_numbers=("1.2.3.4",), db_xrefs=("pdb:1ABC",)
+    )
+    if case == "protein-mismatch":
+        fold_data = fold_data.replace(b'/translation="MK"', b'/translation="MM"')
+    if case == "pair-conflict":
+        fold_data = fold_data.replace(b'/gene="xyz"', b'/gene="other"', 1)
+    base = write(tmp_path / "base.gbff", base_data)
+    fold = write(tmp_path / "fold.gbff", fold_data)
+    manifest = tmp_path / "out.json"
+    enrich(
+        bakta_path=base,
+        baktfold_path=fold,
+        output_path=tmp_path / "out.gbff",
+        manifest_path=manifest,
+    )
+    payload = json.loads(manifest.read_text())
+    rejected = [
+        row
+        for row in payload["decisions"]
+        if row["final_status"].startswith("suppressed_")
+    ]
+    assert (
+        len(rejected)
+        == {
+            "authoritative": 1,
+            "protein-mismatch": 3,
+            "pair-conflict": 2,
+        }[case]
+    )
+    assert payload["metadata"]["suppressed_candidate_count"] == len(rejected)
```
