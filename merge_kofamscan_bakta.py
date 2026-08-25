#!/usr/bin/env python3
"""Merge KofamScan/KOALA hits, optionally with Baktfold, onto Bakta GBFF."""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

from graft_baktfold_additions import plan_baktfold_additions
from merge_engine import (
    Insertion,
    MergeError,
    RawDocument,
    RawFeature,
    comment_insertion,
    finalize_merge,
    parse_faa,
    parse_genbank_bytes,
    qualifier_insertion,
    sha256_bytes,
    validate_genbank_semantics,
)
from merge_engine import (
    validate_faa_gbff as _validate_faa_gbff,
)

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
    base: RawDocument, proteins: dict[str, str], hits: list[KofamHit]
) -> tuple[dict[str, RawFeature], dict[str, Any]]:
    return _validate_faa_gbff(
        base, proteins, (hit.query_id for hit in hits), source_name="Kofam"
    )


def _existing_kos(feature: RawFeature) -> set[str]:
    result: set[str] = set()
    for value in feature.values("db_xref"):
        match = re.fullmatch(r"KEGG:(K\d{5})", value)
        if match:
            result.add(match.group(1))
    for value in feature.values("note"):
        result.update(match.group(1) for match in KEGG_RE.finditer(value))
    return result


def plan_kofam_additions(
    base: RawDocument,
    proteins: dict[str, str],
    hits: list[KofamHit],
    *,
    kofam_data: bytes,
    faa_data: bytes,
    kofamscan_version: str | None,
    add_comment_note: bool,
    merge_timestamp: str | None,
    starting_order: int = 0,
) -> tuple[list[Insertion], list[dict[str, Any]], dict[str, Any]]:
    cds, validation = validate_faa_gbff(base, proteins, hits)
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
                    base.data, feature, "db_xref", xref, "KofamScan", hit.ko, order
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
                    "KofamScan provenance",
                    f"row {hit.row_number}",
                    order,
                )
            )
            order += 1
            note_values.add(note)
            emitted.append(f'/note="{note}"')

        inference_values = planned[(feature.start, "inference")]
        if hit.query_id not in inferred_queries and inference not in inference_values:
            insertions.append(
                qualifier_insertion(
                    base.data,
                    feature,
                    "inference",
                    inference,
                    "KofamScan provenance",
                    kofamscan_version or "unspecified version",
                    order,
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
        "planned_insertions": len(insertions),
    }
    return insertions, evidence_rows, stats


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
) -> dict[str, Any]:
    base_data = bakta_path.read_bytes()
    faa_data = faa_path.read_bytes()
    kofam_data = kofam_path.read_bytes()
    validate_genbank_semantics(base_data, "Bakta input")
    base = parse_genbank_bytes(base_data, "Bakta input")
    proteins = parse_faa(faa_data)
    hits = parse_kofam_table(kofam_data)

    insertions: list[Insertion] = []
    combined_stats: dict[str, Any] = {}
    other_inputs = [faa_path, kofam_path]
    if baktfold_path is not None:
        baktfold_data = baktfold_path.read_bytes()
        validate_genbank_semantics(baktfold_data, "Baktfold input")
        baktfold = parse_genbank_bytes(baktfold_data, "Baktfold input")
        baktfold_insertions, baktfold_stats = plan_baktfold_additions(
            base,
            baktfold,
            baktfold_data=baktfold_data,
            add_feature_provenance=add_feature_provenance,
            add_comment_note=add_comment_note,
            merge_timestamp=merge_timestamp,
        )
        insertions.extend(baktfold_insertions)
        combined_stats["baktfold"] = baktfold_stats
        other_inputs.append(baktfold_path)

    next_order = max((insertion.order for insertion in insertions), default=-1) + 1
    kofam_insertions, evidence_rows, kofam_stats = plan_kofam_additions(
        base,
        proteins,
        hits,
        kofam_data=kofam_data,
        faa_data=faa_data,
        kofamscan_version=kofamscan_version,
        add_comment_note=add_comment_note,
        merge_timestamp=merge_timestamp,
        starting_order=next_order,
    )
    insertions.extend(kofam_insertions)
    combined_stats["kofam"] = kofam_stats
    final = finalize_merge(
        base_path=bakta_path,
        base_data=base_data,
        output_path=output_path,
        other_inputs=other_inputs,
        insertions=insertions,
        manifest_path=manifest_path,
        evidence_rows=evidence_rows,
        metadata={
            "operation": "baktfold-kofam-one-pass" if baktfold_path else "kofam-merge",
            "faa_sha256": kofam_stats["faa_sha256"],
            "kofam_sha256": kofam_stats["kofam_sha256"],
            "kofamscan_version": kofam_stats["kofamscan_version"],
            "baktfold_sha256": (
                combined_stats.get("baktfold", {}).get("baktfold_sha256", "")
            ),
            "merge_timestamp": merge_timestamp or "",
        },
    )
    return {**combined_stats, **final}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Merge validated KofamScan hits onto Bakta GBFF, optionally with Baktfold in one pass."
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
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    for label, path in (
        ("Bakta", args.bakta),
        ("FAA", args.faa),
        ("Kofam", args.kofam),
        ("Baktfold", args.baktfold),
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
