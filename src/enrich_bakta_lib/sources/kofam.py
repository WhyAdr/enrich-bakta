#!/usr/bin/env python3
"""Merge KofamScan/KOALA hits, optionally with Baktfold, onto Bakta GBFF."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    Insertion,
    MergeError,
    RawDocument,
    RawFeature,
    TranslationEvidence,
    comment_insertion,
    qualifier_insertion,
    read_input_bytes,
    sha256_bytes,
)
from enrich_bakta_lib.core.merge_engine import (
    validate_faa_gbff as _validate_faa_gbff,
)
from enrich_bakta_lib.sources.baktfold import BaktfoldAdapter
from enrich_bakta_lib.sources.base import MergeContext, SourceAdapter, SourcePlan
from enrich_bakta_lib.sources.value_rules import structured_note_tokens

KO_RE = re.compile(r"K\d{5}\Z")
KEGG_RE = re.compile(r"(?:^|(?<=[^A-Za-z0-9]))KEGG:(K\d{5})\b")
COMMENT_MARKER = "##enrich-bakta:KofamScan:v1##"


@dataclass(frozen=True)
class KofamHit:
    row_number: int
    query_id: str
    ko: str
    threshold: str
    score: str
    e_value: str
    definition: str


def _decimal(value: str, field: str, row_number: int) -> Decimal:
    try:
        parsed = Decimal(value)
    except InvalidOperation as exc:
        raise MergeError(f"Kofam row {row_number}: invalid {field} {value!r}") from exc
    if not parsed.is_finite():
        raise MergeError(f"Kofam row {row_number}: non-finite {field} {value!r}")
    return parsed


def parse_kofam_table(data: bytes) -> list[KofamHit]:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise MergeError(f"Kofam table is not UTF-8 at byte {exc.start}") from exc
    header_seen = False
    hits: list[KofamHit] = []
    for row_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            normalized = line[1:].lower()
            if all(
                name in normalized
                for name in (
                    "gene name",
                    "ko",
                    "thrshld",
                    "score",
                    "e-value",
                    "definition",
                )
            ):
                header_seen = True
            continue
        if not line.startswith("*"):
            # Standard Kofam detail output may contain below-threshold rows;
            # only leading-* hits are assignments.
            continue
        fields = line[1:].strip().split(maxsplit=5)
        if len(fields) != 6:
            raise MergeError(
                f"Kofam row {row_number}: expected 6 hit fields, found {len(fields)}"
            )
        query_id, ko, threshold, score, e_value, definition = fields
        if query_id.startswith(("*", "#")):
            raise MergeError(
                f"Kofam row {row_number}: query ID {query_id!r} starts with a reserved marker"
            )
        if not KO_RE.fullmatch(ko):
            raise MergeError(f"Kofam row {row_number}: invalid KO identifier {ko!r}")
        _decimal(threshold, "threshold", row_number)
        _decimal(score, "score", row_number)
        evalue = _decimal(e_value, "E-value", row_number)
        if evalue < 0:
            raise MergeError(f"Kofam row {row_number}: negative E-value {e_value!r}")
        hits.append(
            KofamHit(
                row_number=row_number,
                query_id=query_id,
                ko=ko,
                threshold=threshold,
                score=score,
                e_value=e_value,
                definition=definition,
            )
        )
    if not header_seen:
        raise MergeError(
            "Kofam table is missing the expected header-driven column declaration"
        )
    if not hits:
        raise MergeError("Kofam table contains no leading-* hit rows")
    return hits


def validate_faa_gbff(
    base: RawDocument,
    proteins: Mapping[str, str],
    hits: list[KofamHit],
    *,
    translation_evidence: Mapping[str, TranslationEvidence] | None = None,
    allow_imported_translations: bool = False,
) -> tuple[dict[str, RawFeature], dict[str, Any]]:
    return _validate_faa_gbff(
        base,
        proteins,
        (hit.query_id for hit in hits),
        source_name="Kofam",
        translation_evidence=translation_evidence,
        allow_imported_translations=allow_imported_translations,
    )


def _existing_kos(feature: RawFeature) -> set[str]:
    result: set[str] = set()
    for value in feature.values("db_xref"):
        match = re.fullmatch(r"KEGG:(K\d{5})", value)
        if match:
            result.add(match.group(1))
    for value in feature.values("note"):
        result.update(
            token.removeprefix("KEGG:").strip()
            for token in structured_note_tokens(value)
            if token.startswith("KEGG:")
        )
    return result


def plan_kofam_additions(
    base: RawDocument,
    proteins: Mapping[str, str],
    hits: list[KofamHit],
    *,
    kofam_data: bytes,
    faa_data: bytes,
    kofamscan_version: str | None,
    add_comment_note: bool,
    add_feature_provenance: bool = True,
    merge_timestamp: str | None,
    starting_order: int = 0,
    translation_evidence: Mapping[str, TranslationEvidence] | None = None,
    allow_imported_translations: bool = False,
) -> tuple[list[Insertion], list[dict[str, Any]], dict[str, Any]]:
    cds, validation = validate_faa_gbff(
        base,
        proteins,
        hits,
        translation_evidence=translation_evidence,
        allow_imported_translations=allow_imported_translations,
    )
    inference = "profile:KofamScan" + (
        f":{kofamscan_version}" if kofamscan_version else ""
    )
    insertions: list[Insertion] = []
    evidence_rows: list[dict[str, Any]] = []
    order = starting_order
    existing_by_query = {
        query: _existing_kos(feature) for query, feature in cds.items()
    }
    planned: dict[tuple[int, str], set[str]] = {}
    for feature in cds.values():
        for qualifier in ("db_xref", "note", "inference"):
            planned[(feature.start, qualifier)] = set(feature.values(qualifier))

    new_pairs = 0
    existing_pairs = 0
    counted_pairs: set[tuple[str, str]] = set()
    inferred_queries: set[str] = set()
    for hit in hits:
        feature = cds[hit.query_id]
        already_present = hit.ko in existing_by_query[hit.query_id]
        pair = (hit.query_id, hit.ko)
        if pair not in counted_pairs:
            if already_present:
                existing_pairs += 1
            else:
                new_pairs += 1
            counted_pairs.add(pair)
        emitted: list[str] = []

        xref = f"KEGG:{hit.ko}"
        xref_values = planned[(feature.start, "db_xref")]
        if not already_present and xref not in xref_values:
            insertions.append(
                qualifier_insertion(
                    base.data,
                    feature,
                    "db_xref",
                    xref,
                    "KofamScan",
                    hit.ko,
                    order,
                    candidate_role="functional_proposal",
                    evidence_class="ko",
                )
            )
            order += 1
            xref_values.add(xref)
            emitted.append(f'/db_xref="{xref}"')

        note = (
            f"KofamScan:{hit.ko};threshold={hit.threshold};score={hit.score};"
            f"E-value={hit.e_value}"
        )
        note_values = planned[(feature.start, "note")]
        if note not in note_values:
            insertions.append(
                qualifier_insertion(
                    base.data,
                    feature,
                    "note",
                    note,
                    "KofamScan hit evidence",
                    f"row {hit.row_number}",
                    order,
                    candidate_role="substantive_evidence",
                    evidence_class="ko",
                )
            )
            order += 1
            note_values.add(note)
            emitted.append(f'/note="{note}"')

        inference_values = planned[(feature.start, "inference")]
        if (
            add_feature_provenance
            and hit.query_id not in inferred_queries
            and inference not in inference_values
        ):
            insertions.append(
                qualifier_insertion(
                    base.data,
                    feature,
                    "inference",
                    inference,
                    "KofamScan provenance",
                    kofamscan_version or "unspecified version",
                    order,
                    candidate_role="producer_provenance",
                    evidence_class="ko",
                )
            )
            order += 1
            inference_values.add(inference)
            emitted.append(f'/inference="{inference}"')
        inferred_queries.add(hit.query_id)

        evidence_rows.append(
            {
                "entry_type": "kofam_hit",
                "query_id": hit.query_id,
                "record": feature.record_id,
                "feature_type": "CDS",
                "locus_tag": feature.locus_tag or "",
                "ko": hit.ko,
                "threshold": hit.threshold,
                "score": hit.score,
                "e_value": hit.e_value,
                "definition": hit.definition,
                "ko_already_present": already_present,
                "status": "existing" if already_present else "planned",
                "candidate_role": "functional_proposal",
                "support_class": "ko",
                "reason_code": (
                    "value_already_present" if already_present else "inserted"
                ),
                "emitted_qualifiers": " | ".join(emitted),
                "kofam_row": hit.row_number,
            }
        )

    kofam_hash = sha256_bytes(kofam_data)
    faa_hash = sha256_bytes(faa_data)
    if add_comment_note:
        version = kofamscan_version or "unspecified"
        lines = [
            f"Source KofamScan {version}; table sha256={kofam_hash[:16]}",
            f"FAA sha256={faa_hash[:16]}; hits={len(hits)}, hit-CDSs={len(inferred_queries)}.",
            "KO calls are profile evidence, not proof of activity or phenotype.",
        ]
        if merge_timestamp:
            lines.append(f"Merge timestamp: {merge_timestamp}")
        for record in base.records:
            insertion = comment_insertion(
                base,
                record,
                f"{COMMENT_MARKER}:{kofam_hash[:12]}:{faa_hash[:12]}",
                lines,
                "KofamScan provenance",
                order,
            )
            if insertion is not None:
                insertions.append(insertion)
                order += 1

    stats = {
        **validation,
        "hits": len(hits),
        "hit_cds": len({hit.query_id for hit in hits}),
        "new_gene_ko_pairs": new_pairs,
        "already_represented_pairs": existing_pairs,
        "kofam_sha256": kofam_hash,
        "faa_sha256": faa_hash,
        "kofamscan_version": kofamscan_version or "",
        "feature_provenance": add_feature_provenance,
        "planned_insertions": len(insertions),
    }
    return insertions, evidence_rows, stats


@dataclass(frozen=True)
class KofamAdapter:
    """Adapter wrapping KofamScan table parsing and insertion planning."""

    kofam_path: Path
    kofamscan_version: str | None = None
    source_id: str = "KofamScan"

    def plan(self, ctx: MergeContext, *, starting_order: int = 0) -> SourcePlan:
        if ctx.faa_path is None:
            raise MergeError("--faa is required with --kofamscan")
        kofam_data = read_input_bytes(self.kofam_path, "Kofam")
        hits = parse_kofam_table(kofam_data)
        insertions, evidence_rows, stats = plan_kofam_additions(
            ctx.base_doc,
            ctx.proteins or {},
            hits,
            kofam_data=kofam_data,
            faa_data=ctx.faa_data,
            kofamscan_version=self.kofamscan_version,
            add_comment_note=ctx.add_comment_note,
            add_feature_provenance=ctx.add_feature_provenance,
            merge_timestamp=ctx.merge_timestamp,
            starting_order=starting_order,
            translation_evidence=ctx.translation_evidence,
            allow_imported_translations=ctx.allow_imported_translations,
        )
        return SourcePlan(
            source_id=self.source_id,
            source_sha256=stats["kofam_sha256"],
            insertions=insertions,
            evidence_rows=evidence_rows,
            stats=stats,
            manifest_metadata={
                "kofam_sha256": stats["kofam_sha256"],
                "kofamscan_version": stats["kofamscan_version"],
            },
            input_paths=(self.kofam_path,),
        )


def merge(
    bakta_path: Path,
    faa_path: Path,
    kofam_path: Path,
    output_path: Path,
    *,
    baktfold_path: Path | None = None,
    manifest_path: Path | None = None,
    kofamscan_version: str | None = None,
    add_comment_note: bool = True,
    add_feature_provenance: bool = True,
    merge_timestamp: str | None = None,
    translation_evidence_manifest: Path | None = None,
    allow_imported_translations: bool = False,
    baktfold_invalid_ec_policy: str = "reject",
) -> dict[str, Any]:
    from enrich_bakta_lib.workflows.pipeline import (
        capture_merge_context,
        execute_merge_pipeline,
    )

    ctx = capture_merge_context(
        bakta_path=bakta_path,
        faa_path=faa_path,
        translation_evidence_manifest=translation_evidence_manifest,
        allow_imported_translations=allow_imported_translations,
        add_comment_note=add_comment_note,
        add_feature_provenance=add_feature_provenance,
        merge_timestamp=merge_timestamp,
        require_manifest_message=(
            "restored input requires --translation-evidence-manifest for Kofam validation"
        ),
    )
    adapters: list[SourceAdapter] = []
    if baktfold_path is not None:
        adapters.append(
            BaktfoldAdapter(
                baktfold_path=baktfold_path,
                invalid_ec_policy=baktfold_invalid_ec_policy,
            )
        )
    adapters.append(
        KofamAdapter(
            kofam_path=kofam_path,
            kofamscan_version=kofamscan_version,
        )
    )
    baktfold_defaults: dict[str, Any] = {
        "baktfold_sha256": "",
        "baktfold_version": "",
        "baktfold_version_detected": False,
        "baktfold_translation_mismatch_policy": "",
        "baktfold_translation_mismatch_loci": [],
        "baktfold_translation_mismatch_suppressed_features": 0,
        "baktfold_invalid_ec_policy": baktfold_invalid_ec_policy,
        "baktfold_parity": {},
        "baktfold_invalid_ec_values": [],
        "baktfold_invalid_ec_value_count": 0,
    }

    def _kofam_metadata(
        meta: dict[str, Any],
        _plans: list[SourcePlan],
        _candidate_counts: dict[str, Any],
    ) -> dict[str, Any]:
        return {**baktfold_defaults, **meta}

    result = execute_merge_pipeline(
        ctx,
        adapters,
        operation="baktfold-kofam-one-pass" if baktfold_path else "kofam-merge",
        output_path=output_path,
        manifest_path=manifest_path,
        reconcile=True,
        metadata_override=_kofam_metadata,
    )
    combined_stats: dict[str, Any] = {}
    if baktfold_path is not None:
        combined_stats["baktfold"] = result.plans[0].stats
    combined_stats["candidate_ledger"] = result.candidate_counts
    combined_stats["kofam"] = result.plans[-1].stats
    return {**combined_stats, **result.final}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Merge validated KofamScan hits onto Bakta GBFF, optionally with Baktfold in one pass."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument("bakta", type=Path)
    parser.add_argument("faa", type=Path)
    parser.add_argument("kofam", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--baktfold",
        type=Path,
        help="matching Baktfold GBFF for one-pass combined merge",
    )
    parser.add_argument(
        "--manifest", type=Path, help="TSV or .json evidence/insertion manifest"
    )
    parser.add_argument("--kofamscan-version")
    parser.add_argument("--no-comment-note", action="store_true")
    parser.add_argument("--no-feature-provenance", action="store_true")
    parser.add_argument("--merge-timestamp")
    parser.add_argument("--translation-evidence-manifest", type=Path)
    parser.add_argument("--allow-imported-translations", action="store_true")
    parser.add_argument(
        "--baktfold-invalid-ec-policy",
        choices=("reject", "skip"),
        default="reject",
        help="reject non-grammar Baktfold EC tokens or record them as non-promoted",
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    for label, path in (
        ("Bakta", args.bakta),
        ("FAA", args.faa),
        ("Kofam", args.kofam),
        ("Baktfold", args.baktfold),
        ("Translation evidence", args.translation_evidence_manifest),
    ):
        if path is not None and not path.is_file():
            parser.error(f"{label} input not found: {path}")
    try:
        stats = merge(
            args.bakta,
            args.faa,
            args.kofam,
            args.output,
            baktfold_path=args.baktfold,
            manifest_path=args.manifest,
            kofamscan_version=args.kofamscan_version,
            add_comment_note=not args.no_comment_note,
            add_feature_provenance=not args.no_feature_provenance,
            merge_timestamp=args.merge_timestamp,
            translation_evidence_manifest=args.translation_evidence_manifest,
            allow_imported_translations=args.allow_imported_translations,
            baktfold_invalid_ec_policy=args.baktfold_invalid_ec_policy,
        )
    except MergeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(stats, indent=2, sort_keys=True, default=dict))
    else:
        kofam = stats["kofam"]
        print("KofamScan merge passed FAA/GBFF identity and byte-preservation audits.")
        print(f"  hits:                  {kofam['hits']:,}")
        print(f"  hit CDSs:              {kofam['hit_cds']:,}")
        print(f"  new gene-KO pairs:     {kofam['new_gene_ko_pairs']:,}")
        print(f"  already represented:   {kofam['already_represented_pairs']:,}")
        print(f"  output SHA-256:         {stats['output_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
