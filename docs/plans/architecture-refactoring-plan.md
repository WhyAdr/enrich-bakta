# Architectural Refactoring Implementation Plan

- **Date:** 2026-11-25
- **Baseline:** `0.4.1` (post InterProScan audit-hardening fixes; 220 unit tests passing; 42 published LFS artifacts validated)
- **Scope:** Internal modularity, workflow deduplication, adapter contract standardization, and GenBank parsing performance
- **Execution recipient:** Next Gemini session (implementation handoff)
- **Explicit non-goals (deferred):** `--summary-report`, InterProScan TSV/JSON/GFF3 cross-audit tool, `enrich-bakta diff`, and new evidence adapters are deferred to [deferred-features-notes.md](deferred-features-notes.md). Do not bundle new user-facing features or scientific policy changes into this refactor.

---

## 1. Goals and Non-Negotiable Invariants

This refactor addresses four structural friction points accumulated across `0.3.0`–`0.4.1`:

1. **Decompose `core/merge_engine.py`** (~1,730 lines mixing six subsystems) into small, single-responsibility modules under `src/enrich_bakta_lib/core/` with a 100% backwards-compatible re-export facade.
2. **Eliminate redundant Biopython GenBank parses** (currently **6 full parses** of the base/output `.gbff` per merge, reducing to **2**) for a major speedup on real C14/SM genomes with zero behavioral change.
3. **Standardize source adapters** onto a uniform `MergeContext`, `SourcePlan`, and `SourceAdapter` contract across `baktfold`, `kofam`, `eggnog`, and `interproscan`.
4. **Unify standalone and multi-source merge orchestration** into a shared pipeline helper (`workflows/pipeline.py`) to eliminate duplicated snapshot capture, lineage binding, reconciliation, and metadata boilerplate.

### Invariants That Must Not Change

- **Byte-identical GBFF outputs:** Every C14 and SM output hash in [followthrough-scientific-baseline.json](../data/followthrough-scientific-baseline.json) and [interproscan-scientific-baseline.json](../data/interproscan-scientific-baseline.json) must remain identical.
- **Public API and shim identity:** All root compatibility shims (`merge_engine.py`, `enrich_bakta.py`, `graft_baktfold_additions.py`, `merge_kofamscan_bakta.py`, `merge_eggnog_bakta.py`, `merge_interproscan_bakta.py`, `restore_bakta_translations.py`) and existing imports from `enrich_bakta_lib.core.merge_engine` must remain valid (`tests/test_compatibility.py` must pass unmodified).
- **Schema and manifest stability:** `enrich-bakta.merge-manifest.v2`, `enrich-bakta.candidate-ledger.v3`, `enrich-bakta.interproscan-context.v1`, `enrich-bakta.eggnog-context.v1`, and `enrich-bakta.context.v1` must validate without weakening any schema rule or changing `operation` names.
- **Single-snapshot integrity:** FAA, Bakta GBFF, and parent translation-evidence manifest bytes/digests must continue to be captured once and checked against on-disk mutation before promotion.

---

## 2. Workstream 1 — Decompose `core/merge_engine.py`

### Target Module Layout & Acyclic Dependency Graph

```text
src/enrich_bakta_lib/core/
  types.py          # Data models, constants, MergeError, hashing & insertion identity
  genbank.py        # Raw + semantic GenBank parsing, feature indexing, parity, FAA parsing
  lineage.py        # Translation evidence ledger loading/binding, legacy COMMENT detection
  splicer.py        # Qualifier/COMMENT byte formatting, splicing, reversibility & semantic verify
  reconcile.py      # Overlap deduplication, gene conflict policy, provenance pruning
  artifacts.py      # Collision checks, manifest/sidecar serialization, atomic staging, finalize_merge
  decisions.py      # (existing) Candidate ledger construction and semantic validation
  merge_engine.py   # Thin backwards-compatible re-export facade (__all__)
```

Dependency direction (strictly acyclic):

```text
types.py
  ├── genbank.py
  │     ├── lineage.py
  │     └── splicer.py
  ├── reconcile.py
  └── artifacts.py (depends on types, genbank, lineage, splicer, decisions)
        └── merge_engine.py (facade re-exporting all public & tested symbols)
```

### Exact Symbol Ownership Map

| Target module | Symbols moved from `core/merge_engine.py` |
|---|---|
| `core/types.py` | `TOOL_VERSION`, `MergeError`, `CANDIDATE_ROLES`, `RawFeature`, `RawRecord`, `RawDocument`, `Insertion`, `AppliedInsertion`, `TranslationEvidence`, `TranslationEvidenceLedger`, `sha256_bytes`, `sha256_file`, `read_input_bytes`, `insertion_uid`, `insertion_role`, `insertion_evidence_class` |
| `core/genbank.py` | `_FEATURE_RE`, `_QUALIFIER_RE`, `_LOCUS_LENGTH_RE`, `_split_raw_lines`, `_decode`, `parse_genbank_bytes`, `_semantic_records`, `validate_genbank_semantics`, `feature_key`, `feature_uid`, `index_unique_features`, `strict_parity_check`, `normalize_protein`, `parse_faa`, `cds_by_locus`, `protein_sha256` |
| `core/lineage.py` | `TRANSLATION_EVIDENCE_MARKER`, `LEGACY_RESTORATION_COMMENT_RE`, `_extract_comment_blocks`, `has_translation_evidence_marker`, `record_has_legacy_restoration_comment`, `has_legacy_restoration_comment`, `validate_faa_gbff`, `load_translation_evidence`, `_bind_translation_evidence`, `_bind_metadata_translation_evidence` |
| `core/splicer.py` | `QUALIFIER_INDENT`, `newline_for_offset`, `format_qualifier`, `qualifier_insertion`, `comment_insertion`, `apply_insertions`, `verify_insertion_semantics`, `insertion_rows` |
| `core/reconcile.py` | `reconcile_insertions` |
| `core/artifacts.py` | `paths_collide`, `atomic_write`, `manifest_bytes`, `_manifest_cell`, `write_manifest`, `json_sidecar_bytes`, `write_json_sidecar`, `preflight_artifact_paths`, `stage_artifacts`, `finalize_merge` |
| `core/merge_engine.py` | Re-exports every symbol above via explicit imports + `__all__` |

### Rough Pseudocode — `core/merge_engine.py` Facade & `TOOL_VERSION` Single Source

```python
# src/enrich_bakta_lib/core/types.py
TOOL_VERSION = "0.4.1"  # single source of truth; re-exported by merge_engine.py
...

# src/enrich_bakta_lib/core/merge_engine.py
"""Backwards-compatible facade over focused core merge primitives."""

from __future__ import annotations

from enrich_bakta_lib.core.artifacts import (
    atomic_write,
    finalize_merge,
    json_sidecar_bytes,
    manifest_bytes,
    paths_collide,
    preflight_artifact_paths,
    stage_artifacts,
    write_json_sidecar,
    write_manifest,
)
from enrich_bakta_lib.core.genbank import (
    cds_by_locus,
    feature_key,
    feature_uid,
    index_unique_features,
    normalize_protein,
    parse_faa,
    parse_genbank_bytes,
    protein_sha256,
    strict_parity_check,
    validate_genbank_semantics,
)
from enrich_bakta_lib.core.lineage import (
    LEGACY_RESTORATION_COMMENT_RE,
    TRANSLATION_EVIDENCE_MARKER,
    has_legacy_restoration_comment,
    has_translation_evidence_marker,
    load_translation_evidence,
    record_has_legacy_restoration_comment,
    validate_faa_gbff,
)
from enrich_bakta_lib.core.reconcile import reconcile_insertions
from enrich_bakta_lib.core.splicer import (
    QUALIFIER_INDENT,
    apply_insertions,
    comment_insertion,
    format_qualifier,
    insertion_rows,
    newline_for_offset,
    qualifier_insertion,
    verify_insertion_semantics,
)
from enrich_bakta_lib.core.types import (
    CANDIDATE_ROLES,
    TOOL_VERSION,
    AppliedInsertion,
    Insertion,
    MergeError,
    RawDocument,
    RawFeature,
    RawRecord,
    TranslationEvidence,
    TranslationEvidenceLedger,
    insertion_evidence_class,
    insertion_role,
    insertion_uid,
    read_input_bytes,
    sha256_bytes,
    sha256_file,
)

__all__ = [ ... ]  # explicit list of all exported symbols
```

> **Caution on `monkeypatch` in tests:** Check whether any test patches attributes directly on `enrich_bakta_lib.core.merge_engine` (e.g., `monkeypatch.setattr(merge_engine, "os.replace", ...)` or `TOOL_VERSION`). Grep `tests/` for `monkeypatch` and `patch(` before moving functions so internal references don't silently bypass a patched target.

---

## 3. Workstream 2 — Eliminate Redundant Biopython GenBank Parses (`6 → 2`)

### Problem Diagnosis

Every enrichment run currently invokes Biopython's `SeqIO.parse(..., "genbank")` **six times** on the base/merged `.gbff`:

1. Caller calls `validate_genbank_semantics(base_data)` → `_semantic_records(base_data)`
2. Caller calls `parse_genbank_bytes(base_data)` → `_semantic_records(base_data)` *(exact duplicate of #1!)*
3. `finalize_merge()` calls `validate_genbank_semantics(merged)` → `_semantic_records(merged)`
4. `finalize_merge()` calls `parse_genbank_bytes(merged)` → `_semantic_records(merged)` *(exact duplicate of #3!)*
5. `finalize_merge()` calls `parse_genbank_bytes(base_data)` for `verify_insertion_semantics(...)` *(exact duplicate of #2!)*
6. `finalize_merge()` calls `parse_genbank_bytes(base_data)` for `validate_candidate_ledger(...)` *(exact duplicate of #2!)*

(And Baktfold adds two more parses: `validate_genbank_semantics(baktfold_data)` followed immediately by `parse_genbank_bytes(baktfold_data)`.)

Because `parse_genbank_bytes()` **already calls `_semantic_records()`** with `BiopythonParserWarning` promoted to `MergeError`, every preceding `validate_genbank_semantics()` call is redundant, and `finalize_merge()` can accept the already-parsed `RawDocument`.

### Rough Pseudocode — Reusable `RawDocument` Parsing

```python
# In core/artifacts.py: finalize_merge(...)
def finalize_merge(
    *,
    base_path: Path,
    base_data: bytes,
    base_doc: RawDocument | None = None,  # NEW optional parameter
    output_path: Path,
    other_inputs: Iterable[Path],
    insertions: list[Insertion],
    manifest_path: Path | None,
    metadata: dict[str, Any],
    evidence_rows: list[dict[str, Any]] | None = None,
    extra_allowed_qualifiers: Iterable[str] = (),
    sidecar_path: Path | None = None,
    sidecar_payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    ...
    if base_doc is not None and base_doc.data != base_data:
        raise MergeError("pre-parsed base_doc bytes do not match base_data snapshot")

    base_sha256 = sha256_bytes(base_data)
    if sha256_file(base_path) != base_sha256:
        raise MergeError(
            "base file changed on disk since planning; re-run the merge before writing"
        )

    merged, applied = apply_insertions(
        base_data, insertions, extra_allowed_qualifiers=extra_allowed_qualifiers
    )
    # Single parse of merged bytes (includes Biopython warning-as-error semantic validation)
    parsed_output = parse_genbank_bytes(merged, "merged output")
    # Reuse pre-parsed base_doc if provided; otherwise parse once
    parsed_base = (
        base_doc
        if base_doc is not None
        else parse_genbank_bytes(base_data, "Bakta input")
    )
    verify_insertion_semantics(parsed_base, parsed_output, applied)
    if evidence_rows is not None:
        candidate_rows = [
            row
            for row in evidence_rows
            if row.get("entry_type") == "candidate_decision"
        ]
        if candidate_rows:
            validate_candidate_ledger(
                candidate_rows,
                insertions,
                base=parsed_base,
                output=parsed_output,
            )
    ...
```

> **Verification edge case:** Check if any test specifically asserts the exact error message prefix when an invalid GenBank file fails `validate_genbank_semantics()` vs `parse_genbank_bytes()`. Note that `parse_genbank_bytes()` checks raw `LOCUS` / `FEATURES` / `ORIGIN` structure *before* calling `_semantic_records()`, whereas `validate_genbank_semantics()` called `_semantic_records()` first. If a test expects a Biopython error message on a malformed GenBank input that has no `FEATURES` table, ensure error tests still pass (or keep the raw-vs-semantic check ordering aligned with tests).

---

## 4. Workstream 3 — Uniform `MergeContext`, `SourcePlan`, and `SourceAdapter` Contract

### Problem Diagnosis

Today the four adapters expose inconsistent planning return types and metadata conventions:

- `sources/baktfold.py:plan_baktfold_additions(...)` → returns `(list[Insertion], list[dict], dict)`
- `sources/kofam.py:plan_kofam_additions(...)` → returns `(list[Insertion], list[dict], dict)`
- `sources/eggnog.py:plan_eggnog_additions(...)` → returns `EggnogPlan(insertions, evidence_rows, stats, context_report)`
- `sources/interproscan.py:plan_interproscan(...)` → returns `InterProScanPlan(insertions, evidence_rows, stats, context_report, source_sha256, version, layout, member_dbs)`

### Rough Pseudocode — `sources/base.py`

```python
# src/enrich_bakta_lib/sources/base.py
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

from enrich_bakta_lib.core.types import Insertion, RawDocument, TranslationEvidence


@dataclass(frozen=True)
class MergeContext:
    """Single-snapshot input bundle shared across every adapter in a merge."""

    bakta_path: Path
    base_data: bytes
    base_doc: RawDocument
    base_sha256: str
    faa_path: Path | None = None
    faa_data: bytes = b""
    faa_sha256: str = ""
    proteins: Mapping[str, str] | None = None
    translation_evidence_manifest: Path | None = None
    translation_evidence: Mapping[str, TranslationEvidence] | None = None
    allow_imported_translations: bool = False
    add_comment_note: bool = True
    add_feature_provenance: bool = True
    merge_timestamp: str | None = None
    clean_gene_suffix: bool = False


@dataclass(frozen=True)
class SourcePlan:
    """Normalized output returned by every source adapter."""

    source_id: str  # "Baktfold" | "KofamScan" | "eggNOG" | "InterProScan"
    source_sha256: str
    insertions: list[Insertion]
    evidence_rows: list[dict[str, Any]]
    stats: dict[str, Any]
    manifest_metadata: dict[str, Any]
    context_report: dict[str, Any] | None = None
    context_key: str | None = None  # "eggnog" | "interproscan" | None
    input_paths: tuple[Path, ...] = ()


class SourceAdapter(Protocol):
    source_id: str

    def plan(self, ctx: MergeContext, *, starting_order: int = 0) -> SourcePlan:
        """Validate source inputs against ctx and return a deterministic SourcePlan."""
        ...
```

### Backwards-Compatibility Strategy for Adapter Functions

Keep every existing public planner function (`plan_baktfold_additions`, `plan_kofam_additions`, `plan_eggnog_additions`, `plan_interproscan`) with its current signature and return type so external callers and unit tests keep passing. Implement thin adapter classes (`BaktfoldAdapter`, `KofamAdapter`, `EggnogAdapter`, `InterProScanAdapter`) that wrap those planners and return a `SourcePlan`.

```python
# Example adapter wrapper in sources/interproscan.py
@dataclass(frozen=True)
class InterProScanAdapter:
    interproscan_path: Path
    version: str
    layout: str = DEFAULT_TSV_LAYOUT
    member_dbs: str | Sequence[str] = DEFAULT_MEMBER_DBS
    source_id: str = "InterProScan"

    def plan(self, ctx: MergeContext, *, starting_order: int = 0) -> SourcePlan:
        if ctx.faa_path is None:
            raise MergeError("--faa is required with --interproscan")
        parsed_dbs = parse_member_dbs(self.member_dbs)
        ips_plan = plan_interproscan(
            ctx.base_doc,
            ctx.faa_path,
            self.interproscan_path,
            version=self.version,
            layout=self.layout,
            member_dbs=parsed_dbs,
            add_comment_note=ctx.add_comment_note,
            add_feature_provenance=ctx.add_feature_provenance,
            merge_timestamp=ctx.merge_timestamp,
            translation_evidence=ctx.translation_evidence,
            allow_imported_translations=ctx.allow_imported_translations,
            faa_data=ctx.faa_data,
            faa_hash=ctx.faa_sha256,
            proteins=ctx.proteins,
            starting_order=starting_order,  # add starting_order parameter!
        )
        return SourcePlan(
            source_id="InterProScan",
            source_sha256=ips_plan.source_sha256,
            insertions=ips_plan.insertions,
            evidence_rows=ips_plan.evidence_rows,
            stats=ips_plan.stats,
            manifest_metadata={
                "interproscan_sha256": ips_plan.source_sha256,
                "interproscan_version": ips_plan.version,
                "interproscan_tsv_layout": ips_plan.layout,
                "member_dbs": list(ips_plan.member_dbs),
                "interproscan_member_dbs": list(ips_plan.member_dbs),
            },
            context_report=ips_plan.context_report,
            context_key="interproscan",
            input_paths=(self.interproscan_path,),
        )
```

> **Subtle detail to check:** Notice `plan_interproscan()` currently starts `order = 0` unconditionally instead of accepting `starting_order`. In `workflows/enrich.py`, passing `starting_order = next_order` keeps insertion order strictly monotonic across sources, though verify whether changing `order` in `plan_interproscan()` alters C14/SM four-source GBFF insertion ordering/SHA-256! **Check C14/SM four-source GBFF hashes carefully before changing InterProScan's insertion `order` offsets.** If four-source `output_sha256` depends on `order=0` vs `starting_order`, preserve exact `order` behavior!

---

## 5. Workstream 4 — Unify Standalone and Unified Merge Orchestration

### Problem Diagnosis

Currently, five separate entrypoints replicate the same orchestration logic with subtle divergences:

- `sources/baktfold.py:graft()` (no `reconcile_insertions()`, sets `operation="baktfold-graft"`)
- `sources/kofam.py:merge()` (supports optional `--baktfold`, sets `operation="kofam-merge"` or `"baktfold-kofam-one-pass"`)
- `sources/eggnog.py:merge()` (sets `operation="eggnog-merge"`, returns `{"eggnog": stats, **final}`)
- `sources/interproscan.py:merge()` (sets `operation="interproscan-merge"`, returns `final`)
- `workflows/enrich.py:enrich()` (sets `operation="unified-enrichment"`, builds `enrich-bakta.context.v1` when both eggNOG and InterProScan contexts exist)

### Rough Pseudocode — `workflows/pipeline.py`

```python
# src/enrich_bakta_lib/workflows/pipeline.py
from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.artifacts import finalize_merge
from enrich_bakta_lib.core.decisions import (
    build_candidate_ledger,
    validate_candidate_ledger,
)
from enrich_bakta_lib.core.genbank import parse_faa, parse_genbank_bytes
from enrich_bakta_lib.core.lineage import (
    has_translation_evidence_marker,
    load_translation_evidence,
)
from enrich_bakta_lib.core.reconcile import reconcile_insertions
from enrich_bakta_lib.core.types import Insertion, MergeError, read_input_bytes, sha256_bytes
from enrich_bakta_lib.sources.base import MergeContext, SourceAdapter, SourcePlan


def capture_merge_context(
    *,
    bakta_path: Path,
    faa_path: Path | None = None,
    translation_evidence_manifest: Path | None = None,
    allow_imported_translations: bool = False,
    add_comment_note: bool = True,
    add_feature_provenance: bool = True,
    merge_timestamp: str | None = None,
    clean_gene_suffix: bool = False,
    require_manifest_message: str = (
        "restored input requires --translation-evidence-manifest for functional validation"
    ),
) -> MergeContext:
    """Capture single-snapshot bytes, hashes, parsed GenBank, FAA, and bound lineage once."""
    base_data = read_input_bytes(bakta_path, "Bakta")
    base_doc = parse_genbank_bytes(base_data, "Bakta input")
    if (
        has_translation_evidence_marker(base_data)
        and translation_evidence_manifest is None
    ):
        raise MergeError(require_manifest_message)
    translation_evidence = (
        load_translation_evidence(translation_evidence_manifest, base_doc)
        if translation_evidence_manifest is not None
        else None
    )
    faa_data = read_input_bytes(faa_path, "FAA") if faa_path is not None else b""
    faa_sha256 = sha256_bytes(faa_data) if faa_path is not None else ""
    proteins = parse_faa(faa_data) if faa_path is not None else None
    return MergeContext(
        bakta_path=bakta_path,
        base_data=base_data,
        base_doc=base_doc,
        base_sha256=sha256_bytes(base_data),
        faa_path=faa_path,
        faa_data=faa_data,
        faa_sha256=faa_sha256,
        proteins=proteins,
        translation_evidence_manifest=translation_evidence_manifest,
        translation_evidence=translation_evidence,
        allow_imported_translations=allow_imported_translations,
        add_comment_note=add_comment_note,
        add_feature_provenance=add_feature_provenance,
        merge_timestamp=merge_timestamp,
        clean_gene_suffix=clean_gene_suffix,
    )


def assemble_context_sidecar(
    plans: Sequence[SourcePlan], *, operation: str
) -> dict[str, Any] | None:
    """Build single-source or combined enrich-bakta.context.v1 sidecar payload."""
    contexts = {
        plan.context_key: plan.context_report
        for plan in plans
        if plan.context_report is not None and plan.context_key is not None
    }
    if not contexts:
        return None
    for payload in contexts.values():
        payload["metadata"]["operation"] = operation
    if "eggnog" in contexts and "interproscan" in contexts:
        return {
            "schema": "enrich-bakta.context.v1",
            "metadata": {
                "operation": operation,
                "sources": ["eggNOG", "InterProScan"],
            },
            "eggnog": contexts["eggnog"],
            "interproscan": contexts["interproscan"],
        }
    return next(iter(contexts.values()))


def execute_merge_pipeline(
    ctx: MergeContext,
    adapters: Sequence[SourceAdapter],
    *,
    operation: str,
    output_path: Path,
    manifest_path: Path | None = None,
    context_report_path: Path | None = None,
    gene_conflict_policy: str = "skip",
    reconcile: bool = True,
    extra_metadata: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], list[SourcePlan]]:
    """Run adapters, reconcile, build/validate candidate ledger, and finalize merge."""
    plans: list[SourcePlan] = []
    insertions: list[Insertion] = []
    evidence_rows: list[dict[str, Any]] = []
    other_inputs: list[Path] = []
    if ctx.faa_path is not None:
        other_inputs.append(ctx.faa_path)
    if ctx.translation_evidence_manifest is not None:
        other_inputs.append(ctx.translation_evidence_manifest)

    metadata: dict[str, Any] = {
        "operation": operation,
        "merge_timestamp": ctx.merge_timestamp or "",
        "faa_sha256": ctx.faa_sha256,
        "translation_evidence": {
            query: ev.as_dict()
            for query, ev in (ctx.translation_evidence or {}).items()
        },
        "translation_evidence_parent_sha256": getattr(
            ctx.translation_evidence, "manifest_sha256", ""
        ),
        "translation_evidence_manifest": (
            str(ctx.translation_evidence_manifest)
            if ctx.translation_evidence_manifest is not None
            else ""
        ),
        "allow_imported_translations": ctx.allow_imported_translations,
        "policies": {
            "add_comment_note": ctx.add_comment_note,
            "add_feature_provenance": ctx.add_feature_provenance,
        },
    }

    for adapter in adapters:
        next_order = max((item.order for item in insertions), default=-1) + 1
        plan = adapter.plan(ctx, starting_order=next_order)
        plans.append(plan)
        insertions.extend(plan.insertions)
        evidence_rows.extend(plan.evidence_rows)
        metadata.update(plan.manifest_metadata)
        other_inputs.extend(plan.input_paths)

    planned_insertions = list(insertions)
    if reconcile:
        insertions, reconciliation_rows, reconciliation = reconcile_insertions(
            insertions, gene_conflict_policy=gene_conflict_policy
        )
        evidence_rows.extend(reconciliation_rows)
        metadata.update(reconciliation)

    source_hashes = {plan.source_id: plan.source_sha256 for plan in plans}
    candidate_rows, candidate_counts = build_candidate_ledger(
        ctx.base_doc,
        planned_insertions,
        insertions,
        evidence_rows,
        source_hashes=source_hashes,
        gene_conflict_policy=gene_conflict_policy,
    )
    evidence_rows.extend(candidate_rows)
    validate_candidate_ledger(candidate_rows, insertions, base=ctx.base_doc)
    metadata.update(candidate_counts)
    if extra_metadata:
        metadata.update(extra_metadata)

    sidecar_payload = (
        assemble_context_sidecar(plans, operation=operation)
        if context_report_path is not None
        else None
    )
    final = finalize_merge(
        base_path=ctx.bakta_path,
        base_data=ctx.base_data,
        base_doc=ctx.base_doc,
        output_path=output_path,
        other_inputs=other_inputs,
        insertions=insertions,
        manifest_path=manifest_path,
        evidence_rows=evidence_rows,
        metadata=metadata,
        sidecar_path=context_report_path,
        sidecar_payload=sidecar_payload,
    )
    return final, plans
```

> **Preserve return-dict shapes per standalone entrypoint:**
> - `baktfold.graft()` returns `{**baktfold_stats, **final}` (with `candidate_counts` merged into `stats`).
> - `kofam.merge()` returns `{**combined_stats, **final}` (`"kofam"`, optional `"baktfold"`, `"candidate_ledger"`).
> - `eggnog.merge()` returns `{"eggnog": eggnog_stats, **final}` (with `candidate_counts` merged into `eggnog_stats`).
> - `interproscan.merge()` returns `final`.
> - `workflows/enrich.py:enrich()` returns `final`.
> Keep existing manifest metadata keys (e.g., `clean_gene_suffix` in `policies` for eggNOG/unified, `source_hashes` in InterProScan) intact so no existing test or manifest consumer breaks.

---

## 6. Phased Execution Plan

Execute in four small, independently testable commits so regressions are immediately isolated:

| Phase | Scope | Exit criteria |
|---|---|---|
| **Phase 1** | Split `core/merge_engine.py` into `types.py`, `genbank.py`, `lineage.py`, `splicer.py`, `reconcile.py`, `artifacts.py`; keep `core/merge_engine.py` facade | All 220+ unit tests, `test_compatibility.py`, `ruff`, `mypy`, and wheel build pass with zero edits to callers |
| **Phase 2** | Add `base_doc: RawDocument | None = None` to `finalize_merge()`; remove redundant `validate_genbank_semantics()` + duplicate `parse_genbank_bytes()` calls | All unit tests pass; add unit test asserting `SeqIO.parse` is called exactly twice (1x base, 1x merged) during a single merge; verify byte-identical C14/SM outputs and record wall-time delta |
| **Phase 3** | Add `sources/base.py` (`MergeContext`, `SourcePlan`, `SourceAdapter`) and thin adapter classes in `baktfold.py`, `kofam.py`, `eggnog.py`, `interproscan.py` | Unit tests pass; existing direct planner functions (`plan_*`) remain backwards-compatible |
| **Phase 4** | Add `workflows/pipeline.py` (`capture_merge_context`, `assemble_context_sidecar`, `execute_merge_pipeline`); refactor `workflows/enrich.py` and standalone `merge()`/`graft()` entrypoints | All unit tests, schema tests, compatibility tests, and all real-data scientific gates pass with byte-identical C14/SM outputs and identical manifest/context schemas |

---

## 7. Required Verification Commands

Run on a clean Python environment with `.[dev,xlsx]` installed (on Windows sandboxed terminals, pass `--basetemp .test-output/pytest-tmp` to pytest):

```bash
python -m ruff format --check .
python -m ruff check .
python -m mypy
git diff --check
python -m pytest -q -p no:cacheprovider --basetemp .test-output/pytest-tmp
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
python tools/run_followthrough_real_checks.py
python tools/run_followthrough_chain_checks.py
python tools/run_interproscan_real_checks.py
python -m build
```

For wheel verification, inspect the built `.whl` / `.tar.gz` to confirm all new `enrich_bakta_lib/core/*.py`, `enrich_bakta_lib/sources/base.py`, and `enrich_bakta_lib/workflows/pipeline.py` modules are packaged, install the wheel into an isolated venv outside the checkout, and run `enrich-bakta --version` and `enrich-bakta-interproscan --version`.

---

## 8. Handoff Prompt for Next Session

> Execute the architectural refactoring plan in `docs/plans/architecture-refactoring-plan.md` phase by phase (Phases 1–4). Keep every commit focused on a single phase, preserve 100% backwards compatibility on `enrich_bakta_lib.core.merge_engine` and all root launcher shims, and keep every scientific baseline output hash in `docs/data/followthrough-scientific-baseline.json` and `docs/data/interproscan-scientific-baseline.json` byte-identical. Do not implement deferred features from `docs/plans/deferred-features-notes.md`. Run the full verification commands (`ruff`, `mypy`, `pytest`, `git diff --check`, dataset manifest validation, real-data gates, and wheel check) and report wall-time improvements from eliminating redundant GenBank parses before committing.
