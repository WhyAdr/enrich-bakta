#!/usr/bin/env python3
"""Restore only translationless pseudogene CDSs from a matched Bakta FAA."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from merge_eggnog_bakta import EggnogTable, parse_eggnog_path
from merge_engine import (
    Insertion,
    MergeError,
    RawDocument,
    cds_by_locus,
    finalize_merge,
    parse_faa,
    parse_genbank_bytes,
    protein_sha256,
    qualifier_insertion,
    sha256_bytes,
    validate_genbank_semantics,
)


def plan_translation_restoration(
    base: RawDocument,
    proteins: dict[str, str],
    table: EggnogTable,
    *,
    faa_data: bytes,
    starting_order: int = 0,
) -> tuple[list[Insertion], list[dict[str, Any]], dict[str, Any]]:
    """Plan only FAA-backed translations for eggNOG-referenced pseudogene CDSs."""
    cds = cds_by_locus(base)
    insertions: list[Insertion] = []
    evidence_rows: list[dict[str, Any]] = []
    order = starting_order
    requested = {hit.query_id for hit in table.hits}
    missing_faa = sorted(requested - proteins.keys())
    missing_cds = sorted(requested - cds.keys())
    if missing_faa or missing_cds:
        raise MergeError(
            "eggNOG queries do not map to restoration inputs: "
            + json.dumps(
                {
                    "missing_from_faa": missing_faa[:25],
                    "missing_from_gbff": missing_cds[:25],
                },
                indent=2,
            )
        )
    restored: list[str] = []
    for query_id in sorted(requested):
        feature = cds[query_id]
        translations = [
            value for value in feature.values("translation") if value.strip()
        ]
        if translations:
            continue
        if not feature.values("pseudogene"):
            raise MergeError(
                f"refusing to restore {query_id!r}: translationless CDS is not marked /pseudogene"
            )
        sequence = proteins[query_id]
        insertions.append(
            qualifier_insertion(
                base.data,
                feature,
                "translation",
                sequence,
                "Bakta FAA translation restoration",
                protein_sha256(sequence),
                order,
            )
        )
        order += 1
        restored.append(query_id)
        evidence_rows.append(
            {
                "entry_type": "translation_restoration",
                "query_id": query_id,
                "record": feature.record_id,
                "feature_type": feature.feature_type,
                "locus_tag": feature.locus_tag or "",
                "source": "matched Bakta FAA",
                "protein_sha256": protein_sha256(sequence),
                "status": "restored_pseudogene_translation",
            }
        )
    return (
        insertions,
        evidence_rows,
        {
            "requested_eggnog_queries": len(requested),
            "restored_translation_count": len(restored),
            "restored_locus_tags": restored,
            "faa_sha256": sha256_bytes(faa_data),
            "planned_insertions": len(insertions),
        },
    )


def restore(
    bakta_path: Path,
    faa_path: Path,
    eggnog_path: Path,
    output_path: Path,
    *,
    manifest_path: Path | None = None,
    eggnog_version: str | None = None,
) -> dict[str, Any]:
    _base, _faa_data, table, insertions, evidence_rows, stats = prepare_restoration(
        bakta_path, faa_path, eggnog_path, eggnog_version=eggnog_version
    )
    final = finalize_merge(
        base_path=bakta_path,
        output_path=output_path,
        other_inputs=[faa_path, eggnog_path],
        insertions=insertions,
        manifest_path=manifest_path,
        evidence_rows=evidence_rows,
        extra_allowed_qualifiers=("translation",),
        metadata={
            "operation": "pseudogene-translation-restoration",
            "faa_sha256": stats["faa_sha256"],
            "eggnog_version": table.version or "",
            "restored_translation_count": stats["restored_translation_count"],
        },
    )
    return {**stats, **final}


def prepare_restoration(
    bakta_path: Path,
    faa_path: Path,
    eggnog_path: Path,
    *,
    eggnog_version: str | None = None,
) -> tuple[
    RawDocument,
    bytes,
    EggnogTable,
    list[Insertion],
    list[dict[str, Any]],
    dict[str, Any],
]:
    """Validate restoration inputs and return a no-write insertion plan."""
    base_data = bakta_path.read_bytes()
    faa_data = faa_path.read_bytes()
    validate_genbank_semantics(base_data, "Bakta input")
    table = parse_eggnog_path(eggnog_path, expected_version=eggnog_version)
    base = parse_genbank_bytes(base_data, "Bakta input")
    insertions, evidence_rows, stats = plan_translation_restoration(
        base, parse_faa(faa_data), table, faa_data=faa_data
    )
    return base, faa_data, table, insertions, evidence_rows, stats


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Restore only translationless eggNOG-referenced pseudogene CDSs from a matched Bakta FAA."
    )
    parser.add_argument("bakta", type=Path)
    parser.add_argument("faa", type=Path)
    parser.add_argument("eggnog", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--eggnog-version")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    for label, path in (
        ("Bakta", args.bakta),
        ("FAA", args.faa),
        ("eggNOG", args.eggnog),
    ):
        if not path.is_file():
            parser.error(f"{label} input not found: {path}")
    try:
        if args.dry_run:
            *_ignored, stats = prepare_restoration(
                args.bakta,
                args.faa,
                args.eggnog,
                eggnog_version=args.eggnog_version,
            )
        else:
            stats = restore(
                args.bakta,
                args.faa,
                args.eggnog,
                args.output,
                manifest_path=args.manifest,
                eggnog_version=args.eggnog_version,
            )
    except MergeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(stats, indent=2, sort_keys=True))
    elif args.dry_run:
        print(
            f"Validated {stats['restored_translation_count']} pseudogene translation restorations; no output written."
        )
    else:
        print(
            f"Restored {stats['restored_translation_count']} pseudogene translations."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
