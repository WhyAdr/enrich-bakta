#!/usr/bin/env python3
"""Graft genuine Baktfold additions onto a pristine Bakta GenBank file.

The Bakta file remains the authoritative byte stream. Baktfold contributes
only missing gene symbols, EC qualifiers, and exact-prefix structural xrefs.
All inputs are validated before a temporary output is written and atomically
promoted.
"""

from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from pathlib import Path
from typing import Any

from merge_engine import (
    Insertion,
    MergeError,
    RawDocument,
    comment_insertion,
    finalize_merge,
    index_unique_features,
    parse_genbank_bytes,
    qualifier_insertion,
    sha256_bytes,
    strict_parity_check,
    validate_genbank_semantics,
)

STRUCTURAL_PREFIXES = {"afdb_v6", "cath", "pdb"}
COMMENT_MARKER = "##enrich-bakta:Baktfold:v1##"


def detect_baktfold_version(data: bytes) -> str:
    match = re.search(rb"Software:\s*v?([0-9]+(?:\.[0-9]+)*)", data[:20000])
    return match.group(1).decode("ascii") if match else "unknown"


def _structural_xrefs(values: list[str]) -> set[str]:
    result: set[str] = set()
    for value in values:
        prefix, separator, _identifier = value.partition(":")
        if separator and prefix in STRUCTURAL_PREFIXES:
            result.add(value)
    return result


def _ec_values(feature_values: dict[str, list[str]]) -> set[str]:
    result = {
        value.strip() for value in feature_values.get("EC_number", []) if value.strip()
    }
    result.update(
        value[3:].strip()
        for value in feature_values.get("db_xref", [])
        if value.startswith("EC:") and value[3:].strip()
    )
    return result


def plan_baktfold_additions(
    base: RawDocument,
    baktfold: RawDocument,
    *,
    baktfold_data: bytes,
    add_feature_provenance: bool = True,
    add_comment_note: bool = True,
    merge_timestamp: str | None = None,
    starting_order: int = 0,
) -> tuple[list[Insertion], dict[str, Any]]:
    """Validate Baktfold parity and return the complete insertion allowlist."""
    parity = strict_parity_check(base, baktfold)
    source_index = index_unique_features(baktfold.features, label="Baktfold")
    # Also reject duplicates in the base before matching.
    index_unique_features(base.features, label="Bakta base")
    version = detect_baktfold_version(baktfold_data)
    source_hash = sha256_bytes(baktfold_data)
    insertions: list[Insertion] = []
    order = starting_order
    stats: dict[str, Any] = {
        "features_processed": len(base.features),
        "features_matched": 0,
        "features_with_any_addition": 0,
        "gene_added": 0,
        "ec_added": 0,
        "xref_added": 0,
        "gene_by_type": collections.Counter(),
        "ec_by_type": collections.Counter(),
        "xref_by_type": collections.Counter(),
        "baktfold_version": version,
        "baktfold_sha256": source_hash,
        "parity": parity,
    }

    for feature in base.features:
        if not feature.locus_tag:
            continue
        key = (feature.record_id, feature.feature_type, feature.locus_tag)
        source_feature = source_index.get(key)
        if source_feature is None:
            continue
        stats["features_matched"] += 1
        additions_by_class: dict[str, list[tuple[str, str, str]]] = {
            "structural": [],
            "ec": [],
            "gene": [],
        }

        base_xrefs = set(feature.values("db_xref"))
        source_structural = _structural_xrefs(source_feature.values("db_xref"))
        for value in sorted(source_structural - base_xrefs):
            additions_by_class["structural"].append(("db_xref", value, value))

        base_ecs = _ec_values(feature.qualifiers)
        source_ecs = _ec_values(source_feature.qualifiers)
        for value in sorted(source_ecs - base_ecs):
            additions_by_class["ec"].append(("EC_number", value, f"EC:{value}"))

        base_genes = [
            value.strip() for value in feature.values("gene") if value.strip()
        ]
        source_genes = list(
            dict.fromkeys(
                value.strip()
                for value in source_feature.values("gene")
                if value.strip()
            )
        )
        if not base_genes and len(source_genes) > 1:
            raise MergeError(
                f"Baktfold feature {key!r} has ambiguous gene symbols: {source_genes!r}"
            )
        if not base_genes and source_genes:
            additions_by_class["gene"].append(
                ("gene", source_genes[0], source_genes[0])
            )

        if not any(additions_by_class.values()):
            continue
        stats["features_with_any_addition"] += 1
        existing_notes = set(feature.values("note"))
        existing_inferences = set(feature.values("inference"))

        for evidence_class in ("structural", "ec", "gene"):
            additions = additions_by_class[evidence_class]
            for qualifier, value, source_value in additions:
                insertions.append(
                    qualifier_insertion(
                        base.data,
                        feature,
                        qualifier,
                        value,
                        "Baktfold",
                        source_value,
                        order,
                    )
                )
                order += 1
                if evidence_class == "structural":
                    stats["xref_added"] += 1
                    stats["xref_by_type"][feature.feature_type] += 1
                elif evidence_class == "ec":
                    stats["ec_added"] += 1
                    stats["ec_by_type"][feature.feature_type] += 1
                else:
                    stats["gene_added"] += 1
                    stats["gene_by_type"][feature.feature_type] += 1

            if not additions or not add_feature_provenance:
                continue
            if evidence_class == "structural":
                qualifier = "inference"
                provenance = f"protein structure similarity:Baktfold Foldseek:{version}"
                exists = provenance in existing_inferences
            elif evidence_class == "gene":
                qualifier = "note"
                provenance = f"Baktfold gene-symbol evidence:v{version}"
                exists = provenance in existing_notes
            else:
                qualifier = "note"
                provenance = f"Baktfold functional-annotation evidence:v{version}"
                exists = provenance in existing_notes
            if not exists:
                insertions.append(
                    qualifier_insertion(
                        base.data,
                        feature,
                        qualifier,
                        provenance,
                        "Baktfold provenance",
                        evidence_class,
                        order,
                    )
                )
                order += 1

    if add_comment_note:
        lines = [
            f"Source Baktfold v{version}; sha256={source_hash[:16]}",
            (
                f"Added gene={stats['gene_added']}, EC={stats['ec_added']}, "
                f"structural-xref={stats['xref_added']}."
            ),
            "Bakta bytes are otherwise unchanged; see the merge manifest.",
        ]
        if merge_timestamp:
            lines.append(f"Merge timestamp: {merge_timestamp}")
        for record in base.records:
            insertion = comment_insertion(
                base,
                record,
                f"{COMMENT_MARKER}:{source_hash[:16]}",
                lines,
                "Baktfold provenance",
                order,
            )
            if insertion is not None:
                insertions.append(insertion)
                order += 1

    stats["planned_insertions"] = len(insertions)
    return insertions, stats


def graft(
    bakta_path: Path,
    baktfold_path: Path,
    output_path: Path,
    *,
    manifest_path: Path | None = None,
    add_comment_note: bool = True,
    add_feature_provenance: bool = True,
    merge_timestamp: str | None = None,
) -> dict[str, Any]:
    bakta_data = bakta_path.read_bytes()
    baktfold_data = baktfold_path.read_bytes()
    validate_genbank_semantics(bakta_data, "Bakta input")
    validate_genbank_semantics(baktfold_data, "Baktfold input")
    base = parse_genbank_bytes(bakta_data, "Bakta input")
    source = parse_genbank_bytes(baktfold_data, "Baktfold input")
    insertions, stats = plan_baktfold_additions(
        base,
        source,
        baktfold_data=baktfold_data,
        add_feature_provenance=add_feature_provenance,
        add_comment_note=add_comment_note,
        merge_timestamp=merge_timestamp,
    )
    final = finalize_merge(
        base_path=bakta_path,
        output_path=output_path,
        other_inputs=[baktfold_path],
        insertions=insertions,
        manifest_path=manifest_path,
        metadata={
            "operation": "baktfold-graft",
            "baktfold_sha256": stats["baktfold_sha256"],
            "baktfold_version": stats["baktfold_version"],
            "gene_added": stats["gene_added"],
            "ec_added": stats["ec_added"],
            "xref_added": stats["xref_added"],
            "merge_timestamp": merge_timestamp or "",
        },
    )
    return {**stats, **final}


def _serializable_stats(stats: dict[str, Any]) -> dict[str, Any]:
    return {
        key: dict(value) if isinstance(value, collections.Counter) else value
        for key, value in stats.items()
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Graft strictly validated Baktfold additions onto pristine Bakta GBFF bytes."
    )
    parser.add_argument("bakta", type=Path, help="authoritative Bakta .gbff")
    parser.add_argument("baktfold", type=Path, help="matching Baktfold .gbff")
    parser.add_argument("output", type=Path, help="new merged .gbff")
    parser.add_argument("--manifest", type=Path, help="TSV or .json insertion manifest")
    parser.add_argument("--no-comment-note", action="store_true")
    parser.add_argument(
        "--no-feature-provenance",
        "--no-inference-provenance",
        action="store_true",
        help="omit per-feature Baktfold provenance (legacy alias retained)",
    )
    parser.add_argument(
        "--merge-timestamp",
        help="optional explicit timestamp; omitted by default for deterministic output",
    )
    parser.add_argument("--json", action="store_true", help="print summary as JSON")
    parser.add_argument("--force", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.force:
        parser.error(
            "--force is no longer supported; failed identity checks cannot be bypassed"
        )
    for label, path in (("Bakta", args.bakta), ("Baktfold", args.baktfold)):
        if not path.is_file():
            parser.error(f"{label} input not found: {path}")
    try:
        stats = graft(
            args.bakta,
            args.baktfold,
            args.output,
            manifest_path=args.manifest,
            add_comment_note=not args.no_comment_note,
            add_feature_provenance=not args.no_feature_provenance,
            merge_timestamp=args.merge_timestamp,
        )
    except MergeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    summary = _serializable_stats(stats)
    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print("Baktfold graft passed strict parity and byte-preservation audits.")
        print(f"  matched features: {stats['features_matched']:,}")
        print(f"  genes added:      {stats['gene_added']:,}")
        print(f"  ECs added:        {stats['ec_added']:,}")
        print(f"  structural xrefs: {stats['xref_added']:,}")
        print(f"  output SHA-256:   {stats['output_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
