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
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.decisions import (
    build_candidate_ledger,
    validate_candidate_ledger,
)
from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    Insertion,
    MergeError,
    RawDocument,
    TranslationEvidence,
    cds_by_locus,
    comment_insertion,
    finalize_merge,
    has_translation_evidence_marker,
    index_unique_features,
    load_translation_evidence,
    normalize_protein,
    parse_genbank_bytes,
    qualifier_insertion,
    read_input_bytes,
    sha256_bytes,
    strict_parity_check,
    validate_faa_gbff,
    validate_genbank_semantics,
)
from enrich_bakta_lib.sources.value_rules import (
    structural_xref,
    validate_ec_value,
)

STRUCTURAL_PREFIXES = {"afdb_v6", "cath", "pdb"}
COMMENT_MARKER = "##enrich-bakta:Baktfold:v1##"
TRANSLATION_MISMATCH_POLICY = "suppress"


def detect_baktfold_version(data: bytes) -> str:
    match = re.search(rb"Software:\s*v?([0-9]+(?:\.[0-9]+)*)", data[:20000])
    return match.group(1).decode("ascii") if match else "unknown"


def _structural_xrefs(values: list[str]) -> set[str]:
    result: set[str] = set()
    for value in values:
        parsed = structural_xref(value, source="Baktfold")
        if parsed is not None:
            result.add(parsed)
    return result


def _ec_values(
    feature: Any,
    *,
    role: str,
    invalid_policy: str,
    invalid_values: list[dict[str, Any]],
) -> set[str]:
    result: set[str] = set()
    values = [
        *feature.values("EC_number"),
        *[value for value in feature.values("db_xref") if value.startswith("EC:")],
    ]
    for raw in values:
        validation = validate_ec_value(raw, source="Baktfold")
        if validation.status == "missing":
            continue
        if validation.status == "valid" and validation.normalized is not None:
            result.add(validation.normalized)
            continue
        if invalid_policy == "reject":
            raise MergeError(validation.reason)
        invalid_values.append(
            {
                "entry_type": "baktfold_invalid",
                "source": "Baktfold",
                "source_role": role,
                "record": feature.record_id,
                "feature_type": feature.feature_type,
                "locus_tag": feature.locus_tag or "",
                "target_location": feature.location_key,
                "field": "EC",
                "raw_value": validation.raw,
                "normalized_value": validation.normalized or "",
                "status": "invalid_value",
                "validation_status": validation.status,
                "reason": validation.reason,
                "emitted_qualifiers": "",
            }
        )
    return result


def _paired_name_key(feature: Any) -> tuple[int, str, str] | None:
    if feature.feature_type not in {"gene", "CDS"} or not feature.locus_tag:
        return None
    return feature.record_index, feature.locus_tag, feature.location_key


def _paired_names(document: RawDocument) -> dict[tuple[int, str, str], set[str]]:
    names: dict[tuple[int, str, str], set[str]] = collections.defaultdict(set)
    for feature in document.features:
        key = _paired_name_key(feature)
        if key is not None:
            names[key].update(
                value.strip() for value in feature.values("gene") if value.strip()
            )
    return names


def plan_baktfold_additions(
    base: RawDocument,
    baktfold: RawDocument,
    *,
    baktfold_data: bytes,
    add_feature_provenance: bool = True,
    add_comment_note: bool = True,
    merge_timestamp: str | None = None,
    starting_order: int = 0,
    invalid_ec_policy: str = "reject",
    translation_evidence: Mapping[str, TranslationEvidence] | None = None,
    allow_imported_translations: bool = False,
) -> tuple[list[Insertion], dict[str, Any]]:
    """Validate Baktfold parity and return the complete insertion allowlist."""
    if invalid_ec_policy not in {"reject", "skip"}:
        raise MergeError(
            f"invalid Baktfold EC policy {invalid_ec_policy!r}; expected 'reject' or 'skip'"
        )
    parity = strict_parity_check(base, baktfold)
    if has_translation_evidence_marker(base.data) and translation_evidence is None:
        raise MergeError(
            "restored input requires --translation-evidence-manifest for Baktfold planning"
        )
    if translation_evidence:
        cds = cds_by_locus(base)
        encoded_proteins = {
            query: normalize_protein(cds[query].values("translation")[0])
            for query in translation_evidence
        }
        validate_faa_gbff(
            base,
            encoded_proteins,
            translation_evidence,
            source_name="Baktfold",
            translation_evidence=translation_evidence,
            allow_imported_translations=allow_imported_translations,
        )
    source_index = index_unique_features(baktfold.features, label="Baktfold")
    # Also reject duplicates in the base before matching.
    index_unique_features(base.features, label="Bakta base")
    version = detect_baktfold_version(baktfold_data)
    source_hash = sha256_bytes(baktfold_data)
    base_pair_names = _paired_names(base)
    source_pair_names = _paired_names(baktfold)
    mismatch_loci = set(parity.get("translation_mismatches", []))
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
        "baktfold_version_detected": version != "unknown",
        "baktfold_sha256": source_hash,
        "parity": parity,
        "translation_mismatch_policy": TRANSLATION_MISMATCH_POLICY,
        "translation_mismatch_loci": sorted(mismatch_loci),
        "translation_mismatch_suppressed_features": 0,
        "base_pair_name_conflicts": [],
        "source_pair_name_conflicts": [],
        "baktfold_invalid_ec_policy": invalid_ec_policy,
        "invalid_ec_values": [],
    }

    for base_key, names in sorted(base_pair_names.items()):
        if len(names) > 1:
            stats["base_pair_name_conflicts"].append(
                {"pair": list(base_key), "names": sorted(names)}
            )
    for source_key, names in sorted(source_pair_names.items()):
        if len(names) > 1:
            stats["source_pair_name_conflicts"].append(
                {"pair": list(source_key), "names": sorted(names)}
            )

    for feature in base.features:
        if not feature.locus_tag:
            continue
        key = (feature.record_id, feature.feature_type, feature.locus_tag)
        source_feature = source_index.get(key)
        if source_feature is None:
            continue
        stats["features_matched"] += 1
        if feature.locus_tag in mismatch_loci and feature.feature_type in {
            "gene",
            "CDS",
        }:
            stats["translation_mismatch_suppressed_features"] += 1
            continue
        additions_by_class: dict[str, list[tuple[str, str, str]]] = {
            "structural": [],
            "ec": [],
            "gene": [],
        }

        base_xrefs = set(feature.values("db_xref"))
        source_structural = _structural_xrefs(source_feature.values("db_xref"))
        for value in sorted(source_structural - base_xrefs):
            additions_by_class["structural"].append(("db_xref", value, value))

        base_ecs = _ec_values(
            feature,
            role="base",
            invalid_policy=invalid_ec_policy,
            invalid_values=stats["invalid_ec_values"],
        )
        source_ecs = _ec_values(
            source_feature,
            role="source",
            invalid_policy=invalid_ec_policy,
            invalid_values=stats["invalid_ec_values"],
        )
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
        pair_key = _paired_name_key(feature)
        pair_base_names = base_pair_names.get(pair_key, set()) if pair_key else set()
        pair_source_names = (
            source_pair_names.get(pair_key, set()) if pair_key else set()
        )
        pair_name_blocked = len(pair_base_names) > 1 or len(pair_source_names) > 1
        if not base_genes and len(source_genes) > 1:
            raise MergeError(
                f"Baktfold feature {key!r} has ambiguous gene symbols: {source_genes!r}"
            )
        if (
            not base_genes
            and not pair_base_names
            and not pair_name_blocked
            and (source_genes or pair_source_names)
        ):
            source_name = (
                source_genes[0] if source_genes else sorted(pair_source_names)[0]
            )
            additions_by_class["gene"].append(("gene", source_name, source_name))

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
    stats["invalid_ec_value_count"] = len(stats["invalid_ec_values"])
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
    invalid_ec_policy: str = "reject",
    translation_evidence_manifest: Path | None = None,
    allow_imported_translations: bool = False,
) -> dict[str, Any]:
    bakta_data = read_input_bytes(bakta_path, "Bakta")
    baktfold_data = read_input_bytes(baktfold_path, "Baktfold")
    validate_genbank_semantics(bakta_data, "Bakta input")
    validate_genbank_semantics(baktfold_data, "Baktfold input")
    base = parse_genbank_bytes(bakta_data, "Bakta input")
    translation_evidence = (
        load_translation_evidence(translation_evidence_manifest, base)
        if translation_evidence_manifest
        else None
    )
    source = parse_genbank_bytes(baktfold_data, "Baktfold input")
    insertions, stats = plan_baktfold_additions(
        base,
        source,
        baktfold_data=baktfold_data,
        add_feature_provenance=add_feature_provenance,
        add_comment_note=add_comment_note,
        merge_timestamp=merge_timestamp,
        invalid_ec_policy=invalid_ec_policy,
        translation_evidence=translation_evidence,
        allow_imported_translations=allow_imported_translations,
    )
    candidate_rows, candidate_counts = build_candidate_ledger(
        base,
        insertions,
        insertions,
        stats["invalid_ec_values"],
        source_hashes={"Baktfold": stats["baktfold_sha256"]},
    )
    stats.update(candidate_counts)
    validate_candidate_ledger(candidate_rows, insertions, base=base)
    final = finalize_merge(
        base_path=bakta_path,
        base_data=bakta_data,
        output_path=output_path,
        other_inputs=[
            baktfold_path,
            *([translation_evidence_manifest] if translation_evidence_manifest else []),
        ],
        insertions=insertions,
        manifest_path=manifest_path,
        evidence_rows=candidate_rows,
        metadata={
            "operation": "baktfold-graft",
            "baktfold_sha256": stats["baktfold_sha256"],
            "baktfold_version": stats["baktfold_version"],
            "baktfold_version_detected": stats["baktfold_version_detected"],
            "gene_added": stats["gene_added"],
            "ec_added": stats["ec_added"],
            "xref_added": stats["xref_added"],
            "merge_timestamp": merge_timestamp or "",
            "translation_mismatch_policy": stats["translation_mismatch_policy"],
            "translation_mismatch_loci": stats["translation_mismatch_loci"],
            "translation_mismatch_suppressed_features": stats[
                "translation_mismatch_suppressed_features"
            ],
            "base_pair_name_conflicts": stats["base_pair_name_conflicts"],
            "source_pair_name_conflicts": stats["source_pair_name_conflicts"],
            "baktfold_invalid_ec_policy": stats["baktfold_invalid_ec_policy"],
            "invalid_ec_values": stats["invalid_ec_values"],
            "invalid_ec_value_count": stats["invalid_ec_value_count"],
            **candidate_counts,
            "translation_evidence": {
                query: evidence.as_dict()
                for query, evidence in (translation_evidence or {}).items()
            },
            "translation_evidence_parent_sha256": getattr(
                translation_evidence, "manifest_sha256", ""
            ),
            "allow_imported_translations": allow_imported_translations,
            "policies": {
                "add_comment_note": add_comment_note,
                "add_feature_provenance": add_feature_provenance,
            },
            "parity": stats["parity"],
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
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument("bakta", type=Path, help="authoritative Bakta .gbff")
    parser.add_argument("baktfold", type=Path, help="matching Baktfold .gbff")
    parser.add_argument("output", type=Path, help="new merged .gbff")
    parser.add_argument("--manifest", type=Path, help="TSV or .json insertion manifest")
    parser.add_argument("--no-comment-note", action="store_true")
    parser.add_argument("--translation-evidence-manifest", type=Path)
    parser.add_argument("--allow-imported-translations", action="store_true")
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
    parser.add_argument(
        "--baktfold-invalid-ec-policy",
        choices=("reject", "skip"),
        default="reject",
        help="reject non-grammar EC tokens or record them as non-promoted values",
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
            invalid_ec_policy=args.baktfold_invalid_ec_policy,
            translation_evidence_manifest=args.translation_evidence_manifest,
            allow_imported_translations=args.allow_imported_translations,
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
