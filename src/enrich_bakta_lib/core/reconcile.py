"""Overlap deduplication, gene conflict policy, and provenance pruning."""

from __future__ import annotations

import collections
from collections.abc import Iterable
from typing import Any

from enrich_bakta_lib.core.types import (
    Insertion,
    MergeError,
    insertion_evidence_class,
    insertion_role,
)


def reconcile_insertions(
    insertions: Iterable[Insertion], *, gene_conflict_policy: str = "skip"
) -> tuple[list[Insertion], list[dict[str, Any]], dict[str, Any]]:
    """Collapse exact overlaps and resolve paired gene-symbol disagreements."""
    if gene_conflict_policy not in {"skip", "prefer-eggnog", "prefer-baktfold"}:
        raise MergeError(f"invalid gene conflict policy {gene_conflict_policy!r}")
    ordered = sorted(
        insertions,
        key=lambda item: (
            item.record,
            item.locus_tag,
            item.feature_type,
            item.qualifier,
            item.value,
            item.order,
            item.source,
        ),
    )
    rows: list[dict[str, Any]] = []
    gene_values: dict[tuple[str, str], set[str]] = collections.defaultdict(set)
    for item in ordered:
        if item.qualifier == "gene":
            gene_values[(item.record, item.locus_tag)].add(item.value)
    conflicts = {key: values for key, values in gene_values.items() if len(values) > 1}
    preferred_source = {
        "prefer-eggnog": "eggnog",
        "prefer-baktfold": "baktfold",
    }.get(gene_conflict_policy)
    filtered: list[Insertion] = []
    for item in ordered:
        conflict_values = conflicts.get((item.record, item.locus_tag))
        if (
            item.qualifier == "gene"
            and conflict_values
            and (preferred_source is None or item.source.lower() != preferred_source)
        ):
            continue
        filtered.append(item)

    def source_family(item: Insertion) -> str:
        if item.source.startswith("Baktfold"):
            return "baktfold"
        if item.source.startswith("KofamScan"):
            return "kofamscan"
        if item.source.startswith("eggNOG"):
            return "eggnog"
        return item.source.lower()

    surviving_functional = {
        (
            item.record,
            item.locus_tag,
            source_family(item),
        )
        for item in filtered
        if insertion_role(item) in {"functional_proposal", "substantive_evidence"}
    }
    provenance_filtered: list[Insertion] = []
    for item in filtered:
        if item.source == "Baktfold provenance" and item.qualifier != "COMMENT":
            if not any(
                other.source == "Baktfold"
                and other.record == item.record
                and other.locus_tag == item.locus_tag
                and other.feature_type == item.feature_type
                and insertion_evidence_class(other) == insertion_evidence_class(item)
                for other in filtered
            ):
                continue
        if (
            item.qualifier != "COMMENT"
            and insertion_role(item) == "producer_provenance"
            and (
                item.record,
                item.locus_tag,
                source_family(item),
            )
            not in surviving_functional
        ):
            rows.append(
                {
                    "entry_type": "reconciliation",
                    "record": item.record,
                    "feature_type": item.feature_type,
                    "locus_tag": item.locus_tag,
                    "qualifier": item.qualifier,
                    "value": item.value,
                    "source": item.source,
                    "status": "provenance_without_surviving_annotation_suppressed",
                    "supporting_sources": item.source,
                }
            )
            continue
        provenance_filtered.append(item)
    filtered = provenance_filtered
    unique: dict[tuple[str, str, str, str, str], Insertion] = {}
    sources: dict[tuple[str, str, str, str, str], set[str]] = collections.defaultdict(
        set
    )
    for item in filtered:
        key = (
            item.record,
            item.feature_type,
            item.locus_tag,
            item.qualifier,
            item.value,
        )
        sources[key].add(item.source)
        unique.setdefault(key, item)
    for key, supported_by in sorted(sources.items()):
        if len(supported_by) > 1:
            rows.append(
                {
                    "entry_type": "reconciliation",
                    "record": key[0],
                    "feature_type": key[1],
                    "locus_tag": key[2],
                    "qualifier": key[3],
                    "value": key[4],
                    "status": "exact_duplicate_collapsed",
                    "supporting_sources": " | ".join(sorted(supported_by)),
                }
            )
    for (record, locus_tag), values in sorted(conflicts.items()):
        rows.append(
            {
                "entry_type": "reconciliation",
                "record": record,
                "locus_tag": locus_tag,
                "qualifier": "gene",
                "value": " | ".join(sorted(values)),
                "status": "gene_conflict_skipped"
                if preferred_source is None
                else gene_conflict_policy,
                "supporting_sources": " | ".join(
                    sorted(
                        {
                            item.source
                            for item in ordered
                            if item.record == record
                            and item.locus_tag == locus_tag
                            and item.qualifier == "gene"
                        }
                    )
                ),
            }
        )
    reconciled = sorted(unique.values(), key=lambda item: (item.offset, item.order))
    return (
        reconciled,
        rows,
        {
            "input_insertions": len(ordered),
            "reconciled_insertions": len(reconciled),
            "exact_duplicates_collapsed": len(filtered) - len(reconciled),
            "gene_conflicts": len(conflicts),
            "gene_conflict_policy": gene_conflict_policy,
        },
    )
