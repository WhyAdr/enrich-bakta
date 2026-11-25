"""Shared single-snapshot merge orchestration pipeline."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
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
from enrich_bakta_lib.core.types import (
    Insertion,
    MergeError,
    RawDocument,
    TranslationEvidence,
    read_input_bytes,
    sha256_bytes,
)
from enrich_bakta_lib.sources.base import MergeContext, SourceAdapter, SourcePlan


@dataclass(frozen=True)
class PipelineResult:
    """Finalized artifacts, source plans, and candidate ledger counts."""

    final: dict[str, Any]
    plans: tuple[SourcePlan, ...]
    candidate_counts: dict[str, Any]


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
    translation_evidence_loader: Callable[
        [Path, RawDocument], dict[str, TranslationEvidence]
    ] = load_translation_evidence,
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
        translation_evidence_loader(translation_evidence_manifest, base_doc)
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
    summary_report_path: Path | None = None,
    gene_conflict_policy: str = "skip",
    reconcile: bool = True,
    include_clean_gene_suffix_policy: bool = False,
    extra_metadata: Mapping[str, Any] | None = None,
    metadata_override: (
        Callable[[dict[str, Any], list[SourcePlan], dict[str, Any]], dict[str, Any]]
        | None
    ) = None,
) -> PipelineResult:
    """Run adapters, reconcile, build/validate candidate ledger, and finalize merge."""
    plans: list[SourcePlan] = []
    insertions: list[Insertion] = []
    evidence_rows: list[dict[str, Any]] = []
    other_inputs: list[Path] = []
    if ctx.faa_path is not None:
        other_inputs.append(ctx.faa_path)
    if ctx.translation_evidence_manifest is not None:
        other_inputs.append(ctx.translation_evidence_manifest)

    policies: dict[str, Any] = {
        "add_comment_note": ctx.add_comment_note,
        "add_feature_provenance": ctx.add_feature_provenance,
    }
    if include_clean_gene_suffix_policy:
        policies["clean_gene_suffix"] = ctx.clean_gene_suffix

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
        "policies": policies,
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
    reconciliation: dict[str, Any] = {}
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
    if metadata_override is not None:
        metadata = metadata_override(metadata, plans, candidate_counts)

    sidecar_payload = (
        assemble_context_sidecar(plans, operation=operation)
        if context_report_path is not None or summary_report_path is not None
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
        summary_report_path=summary_report_path,
    )
    return PipelineResult(
        final=final, plans=tuple(plans), candidate_counts=candidate_counts
    )
