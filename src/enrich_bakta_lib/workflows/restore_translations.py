#!/usr/bin/env python3
"""Restore only translationless pseudogene CDSs from a matched Bakta FAA."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition

from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    TRANSLATION_EVIDENCE_MARKER,
    Insertion,
    MergeError,
    RawDocument,
    TranslationEvidence,
    cds_by_locus,
    comment_insertion,
    feature_uid,
    finalize_merge,
    has_translation_evidence_marker,
    load_translation_evidence,
    normalize_protein,
    parse_faa,
    parse_genbank_bytes,
    protein_sha256,
    qualifier_insertion,
    read_input_bytes,
    sha256_bytes,
    validate_faa_gbff,
    validate_genbank_semantics,
)
from enrich_bakta_lib.sources.eggnog import EggnogTable, parse_eggnog_path

TRANSLATION_POLICIES = ("validated-only", "import-faa")
VALID_PROTEIN_SYMBOLS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def _validate_protein_artifact(sequence: str, query_id: str) -> str:
    normalized = normalize_protein(sequence)
    invalid = sorted(set(normalized) - VALID_PROTEIN_SYMBOLS)
    if not normalized or invalid:
        raise MergeError(
            f"FAA protein {query_id!r} contains invalid symbols: {invalid!r}"
        )
    return normalized


def _genomic_translation_check(
    base: RawDocument, feature: Any, protein: str
) -> tuple[str, str]:
    """Compare FAA protein to genomic translation where qualifiers support it."""
    location = feature.semantic_location
    if location is None:
        return "unresolved", "semantic CDS location is unavailable"
    if any(
        key in feature.qualifiers
        for key in ("exception", "transl_except", "ribosomal_slippage")
    ):
        return "unresolved", "exceptional CDS translation model is unsupported"
    if getattr(location, "operator", "join") != "join":
        return "unresolved", "CDS location operator is unsupported"
    record_length = len(base.records[feature.record_index].sequence)
    for part in location.parts:
        if part.ref is not None or part.ref_db is not None:
            return "unresolved", "remote CDS location is unsupported"
        if type(part.start) is not ExactPosition or type(part.end) is not ExactPosition:
            return "unresolved", "partial or fuzzy CDS location is unsupported"
        if not 0 <= int(part.start) < int(part.end) <= record_length:
            return "unresolved", "CDS location lies outside the source record"
        if part.strand not in (-1, 1) or part.strand != location.strand:
            return "unresolved", "CDS location strand is unsupported"
    codon_start_values = feature.values("codon_start")
    transl_table_values = feature.values("transl_table")
    if len(codon_start_values) != 1 or len(transl_table_values) != 1:
        return "unresolved", "codon_start and transl_table must each be single-valued"
    try:
        codon_start = int(codon_start_values[0])
        transl_table = int(transl_table_values[0])
    except ValueError:
        return "unresolved", "codon_start and transl_table must be integers"
    if codon_start not in (1, 2, 3) or transl_table < 1:
        return "unresolved", "codon_start or transl_table is outside supported values"
    if codon_start != 1:
        return "unresolved", "offset CDS translation requires a supported partial model"
    try:
        sequence = Seq(base.records[feature.record_index].sequence.decode("ascii"))
        coding = location.extract(sequence)[codon_start - 1 :]
        if len(coding) == 0 or len(coding) % 3:
            return "unresolved", "CDS sequence is incomplete after codon_start"
        translated = str(coding.translate(table=transl_table, cds=True))
    except (
        UnicodeDecodeError,
        ValueError,
        TypeError,
        KeyError,
        TranslationError,
    ) as exc:
        return "unresolved", f"genomic translation is unsupported: {exc}"
    if translated != protein:
        return "failed", "FAA protein differs from the independent genomic translation"
    return "passed", "FAA protein matches the independent genomic translation"


def plan_translation_restoration(
    base: RawDocument,
    proteins: dict[str, str],
    table: EggnogTable,
    *,
    faa_data: bytes,
    starting_order: int = 0,
    translation_policy: str = "validated-only",
    translation_evidence: dict[str, TranslationEvidence] | None = None,
) -> tuple[list[Insertion], list[dict[str, Any]], dict[str, Any]]:
    """Plan policy-controlled FAA-backed translations for pseudogene CDSs."""
    if translation_policy not in TRANSLATION_POLICIES:
        raise MergeError(f"invalid translation policy {translation_policy!r}")
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
    base_hash = sha256_bytes(base.data)
    faa_hash = sha256_bytes(faa_data)
    restored: list[str] = []
    origin_evidence: dict[str, TranslationEvidence] = dict(translation_evidence or {})
    translated_queries = [
        query_id
        for query_id in sorted(requested)
        if any(value.strip() for value in cds[query_id].values("translation"))
    ]
    if translated_queries:
        validate_faa_gbff(
            base,
            proteins,
            translated_queries,
            source_name="eggNOG translation restoration",
            translation_evidence=translation_evidence,
            allow_imported_translations=True,
        )
    for query_id in sorted(requested):
        feature = cds[query_id]
        translations = feature.values("translation")
        if len(translations) > 1:
            raise MergeError(
                f"refusing to restore {query_id!r}: duplicate /translation qualifiers"
            )
        if translations and not translations[0].strip():
            raise MergeError(
                f"refusing to restore {query_id!r}: empty /translation qualifier"
            )
        if translations:
            evidence_rows.append(
                {
                    "entry_type": "translation_validation",
                    "query_id": query_id,
                    "status": "supported_existing",
                    "artifact_sequence_match": True,
                    "genomic_validation_status": "not_attempted",
                    "translation_origin": (
                        origin_evidence[query_id].origin
                        if query_id in origin_evidence
                        else "gbff_encoded"
                    ),
                }
            )
            continue
        if (
            "pseudo" not in feature.qualifiers
            and "pseudogene" not in feature.qualifiers
        ):
            raise MergeError(
                f"refusing to restore {query_id!r}: translationless CDS is not marked /pseudo or /pseudogene"
            )
        sequence = _validate_protein_artifact(proteins[query_id], query_id)
        genomic_status, validation_reason = _genomic_translation_check(
            base, feature, sequence
        )
        if translation_policy == "validated-only" and genomic_status != "passed":
            raise MergeError(
                f"refusing to restore {query_id!r} under validated-only policy: "
                f"{validation_reason}"
            )
        origin = (
            "genomically_validated"
            if translation_policy == "validated-only"
            else "imported_faa"
        )
        evidence = TranslationEvidence(
            feature_uid=feature_uid(feature),
            protein_sha256=protein_sha256(sequence),
            origin=origin,
            artifact_sequence_match=True,
            genomic_validation_status=genomic_status,
            validation_reason=validation_reason,
            original_base_sha256=base_hash,
            faa_sha256=faa_hash,
            restoration_policy=translation_policy,
            producer_lineage={
                "eggnog_version": table.version or "",
                "eggnog_schema": table.schema_id,
            },
        )
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
        origin_evidence[query_id] = evidence
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
                "translation_policy": translation_policy,
                "translation_evidence": evidence.as_dict(),
            }
        )
    if restored:
        records = sorted({cds[query_id].record_index for query_id in restored})
        for record_index in records:
            marker = (
                f"{TRANSLATION_EVIDENCE_MARKER}:{base_hash[:16]}:{faa_hash[:16]}"
                f":{translation_policy}"
            )
            insertion = comment_insertion(
                base,
                base.records[record_index],
                marker,
                [
                    f"Restored CDSs: {', '.join(query_id for query_id in restored if cds[query_id].record_index == record_index)}",
                    f"Original base SHA-256: {base_hash}",
                    f"FAA SHA-256: {faa_hash}",
                    f"Translation policy: {translation_policy}",
                ],
                "translation restoration provenance",
                order,
            )
            if insertion is not None:
                insertions.append(insertion)
                order += 1
    return (
        insertions,
        evidence_rows,
        {
            "requested_eggnog_queries": len(requested),
            "restored_translation_count": len(restored),
            "restored_locus_tags": restored,
            "faa_sha256": sha256_bytes(faa_data),
            "planned_insertions": len(insertions),
            "translation_policy": translation_policy,
            "translation_evidence": {
                query_id: evidence.as_dict()
                for query_id, evidence in origin_evidence.items()
            },
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
    eggnog_schema: str | None = None,
    translation_policy: str = "validated-only",
    translation_evidence_manifest: Path | None = None,
) -> dict[str, Any]:
    _base, base_data, _faa_data, table, insertions, evidence_rows, stats = (
        prepare_restoration(
            bakta_path,
            faa_path,
            eggnog_path,
            eggnog_version=eggnog_version,
            eggnog_schema=eggnog_schema,
            translation_policy=translation_policy,
            translation_evidence_manifest=translation_evidence_manifest,
        )
    )
    final = finalize_merge(
        base_path=bakta_path,
        base_data=base_data,
        output_path=output_path,
        other_inputs=[
            faa_path,
            eggnog_path,
            *([translation_evidence_manifest] if translation_evidence_manifest else []),
        ],
        insertions=insertions,
        manifest_path=manifest_path,
        evidence_rows=evidence_rows,
        extra_allowed_qualifiers=("translation",),
        metadata={
            "operation": "pseudogene-translation-restoration",
            "faa_sha256": stats["faa_sha256"],
            "eggnog_sha256": stats["eggnog_sha256"],
            "eggnog_version": table.version or "",
            "eggnog_schema": table.schema_id,
            "schema_selection_source": table.schema_selection_source,
            "restored_translation_count": stats["restored_translation_count"],
            "translation_policy": stats["translation_policy"],
            "translation_evidence": stats["translation_evidence"],
            "translation_evidence_parent_sha256": stats[
                "translation_evidence_parent_sha256"
            ],
        },
    )
    result = {**stats, **final}
    result["translation_evidence"] = {
        query_id: {
            **evidence,
            "bound_output_sha256": final["output_sha256"],
        }
        for query_id, evidence in stats["translation_evidence"].items()
    }
    return result


def prepare_restoration(
    bakta_path: Path,
    faa_path: Path,
    eggnog_path: Path,
    *,
    eggnog_version: str | None = None,
    eggnog_schema: str | None = None,
    translation_policy: str = "validated-only",
    translation_evidence_manifest: Path | None = None,
) -> tuple[
    RawDocument,
    bytes,
    bytes,
    EggnogTable,
    list[Insertion],
    list[dict[str, Any]],
    dict[str, Any],
]:
    """Validate restoration inputs and return a no-write insertion plan."""
    base_data = read_input_bytes(bakta_path, "Bakta")
    faa_data = read_input_bytes(faa_path, "FAA")
    eggnog_data = read_input_bytes(eggnog_path, "eggNOG")
    validate_genbank_semantics(base_data, "Bakta input")
    table = parse_eggnog_path(
        eggnog_path,
        expected_version=eggnog_version,
        data=eggnog_data,
        schema_id=eggnog_schema,
    )
    base = parse_genbank_bytes(base_data, "Bakta input")
    if (
        has_translation_evidence_marker(base_data)
        and translation_evidence_manifest is None
    ):
        raise MergeError(
            "restored input requires --translation-evidence-manifest for restoration reruns"
        )
    origin_evidence = (
        load_translation_evidence(translation_evidence_manifest, base)
        if translation_evidence_manifest
        else None
    )
    insertions, evidence_rows, stats = plan_translation_restoration(
        base,
        parse_faa(faa_data),
        table,
        faa_data=faa_data,
        translation_policy=translation_policy,
        translation_evidence=origin_evidence,
    )
    stats["eggnog_sha256"] = sha256_bytes(eggnog_data)
    stats["translation_evidence_parent_sha256"] = getattr(
        origin_evidence, "manifest_sha256", ""
    )
    return base, base_data, faa_data, table, insertions, evidence_rows, stats


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Restore only translationless eggNOG-referenced pseudogene CDSs from a matched Bakta FAA."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument("bakta", type=Path)
    parser.add_argument("faa", type=Path)
    parser.add_argument("eggnog", type=Path)
    parser.add_argument("output", type=Path, nargs="?")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--eggnog-version")
    parser.add_argument("--eggnog-schema")
    parser.add_argument("--translation-evidence-manifest", type=Path)
    parser.add_argument(
        "--translation-policy",
        choices=TRANSLATION_POLICIES,
        help="validated-only requires genomic translation agreement; import-faa requires downstream opt-in",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    if not args.dry_run and args.output is None:
        parser.error("output is required unless --dry-run is used")
    if not args.dry_run and args.translation_policy is None:
        parser.error("--translation-policy is required for restoration writes")
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
                eggnog_schema=args.eggnog_schema,
                translation_policy=args.translation_policy or "validated-only",
                translation_evidence_manifest=args.translation_evidence_manifest,
            )
        else:
            assert args.output is not None
            stats = restore(
                args.bakta,
                args.faa,
                args.eggnog,
                args.output,
                manifest_path=args.manifest,
                eggnog_version=args.eggnog_version,
                eggnog_schema=args.eggnog_schema,
                translation_policy=args.translation_policy or "validated-only",
                translation_evidence_manifest=args.translation_evidence_manifest,
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
