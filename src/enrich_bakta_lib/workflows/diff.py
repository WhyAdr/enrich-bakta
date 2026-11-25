"""Compare annotation qualifiers or inspect a merge manifest's decisions."""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from collections import Counter, defaultdict
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.artifacts import json_sidecar_bytes, stage_artifacts
from enrich_bakta_lib.core.genbank import parse_genbank_bytes, strict_parity_check
from enrich_bakta_lib.core.types import (
    TOOL_VERSION,
    MergeError,
    RawDocument,
    RawRecord,
    read_input_bytes,
    sha256_bytes,
)

_ADDED_STATUSES = {"emitted", "shared_support"}
LocusLocation = tuple[str, str]


def _comment_lines(document: RawDocument, record: RawRecord) -> list[str]:
    """Read top-level COMMENT lines, including their continuations."""
    comments = []
    in_comment = False
    for line in document.data[
        record.start : record.features_offset or record.end
    ].splitlines():
        field = line[:12].strip()
        if field:
            in_comment = field == b"COMMENT"
        if in_comment and line[12:].strip():
            comments.append(line[12:].strip().decode(document.text_encoding))
    return comments


def _value_delta(before: list[str], after: list[str]) -> tuple[list[str], list[str]]:
    left, right = Counter(before), Counter(after)
    return sorted((right - left).elements()), sorted((left - right).elements())


def diff_gbff_files(
    base_path: Path,
    enriched_path: Path,
    *,
    loci: set[str] | None = None,
    qualifiers: set[str] | None = None,
) -> dict[str, Any]:
    base_bytes = read_input_bytes(base_path, "base GBFF")
    enriched_bytes = read_input_bytes(enriched_path, "enriched GBFF")
    base = parse_genbank_bytes(base_bytes, "base GBFF")
    enriched = parse_genbank_bytes(enriched_bytes, "enriched GBFF")
    strict_parity_check(base, enriched)

    deltas = []
    comments = []
    for left, right in zip(base.records, enriched.records, strict=True):
        if not loci and (not qualifiers or "COMMENT" in qualifiers):
            added_comments, removed_comments = _value_delta(
                _comment_lines(base, left), _comment_lines(enriched, right)
            )
            if added_comments or removed_comments:
                comments.append(
                    {
                        "record_id": left.record_id,
                        "record_index": left.index,
                        "added": added_comments,
                        "removed": removed_comments,
                    }
                )
        for before, after in zip(left.features, right.features, strict=True):
            if loci and not loci.intersection(before.locus_tags):
                continue
            keys = set(before.qualifiers) | set(after.qualifiers)
            if qualifiers:
                keys &= qualifiers
            added, removed = {}, {}
            for key in sorted(keys):
                increases, decreases = _value_delta(
                    before.values(key), after.values(key)
                )
                if increases:
                    added[key] = increases
                if decreases:
                    removed[key] = decreases
            if added or removed:
                deltas.append(
                    {
                        "record_id": left.record_id,
                        "record_index": left.index,
                        "locus_tag": before.locus_tag or "",
                        "feature_type": before.feature_type,
                        "location": before.location_key,
                        "added": added,
                        "removed": removed,
                    }
                )
    return {
        "schema": "enrich-bakta.diff.v1",
        "mode": "gbff",
        "base_sha256": sha256_bytes(base_bytes),
        "enriched_sha256": sha256_bytes(enriched_bytes),
        "changed_features": len(deltas),
        "changed_records_comments": len(comments),
        "deltas": deltas,
        "comment_deltas": comments,
    }


def _manifest_objects(payload: dict[str, Any], field: str) -> list[dict[str, Any]]:
    rows = payload.get(field)
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise MergeError(f"merge manifest {field!r} must be an array of objects")
    return rows


def _strings(row: dict[str, Any], field: str) -> tuple[str, ...]:
    items = row.get(field, [])
    if not isinstance(items, list) or any(not isinstance(item, str) for item in items):
        raise MergeError(f"merge manifest {field!r} must be an array of strings")
    return tuple(sorted(set(items)))


def _entry_locations(
    entries: list[dict[str, Any]],
) -> tuple[
    dict[str, set[LocusLocation]], dict[str, LocusLocation], dict[str, LocusLocation]
]:
    features: dict[str, set[LocusLocation]] = defaultdict(set)
    insertions: dict[str, LocusLocation] = {}
    candidates: dict[str, LocusLocation] = {}
    for entry in entries:
        record = entry.get("record", "")
        locus = entry.get("locus_tag") or entry.get("query_id")
        if not isinstance(record, str) or not isinstance(locus, str) or not locus:
            continue
        location = (record, locus)
        if isinstance(entry.get("candidate_id"), str):
            candidates[entry["candidate_id"]] = location
        if not record:
            continue
        for uid in _strings(entry, "target_feature_uids"):
            features[uid].add(location)
        if isinstance(entry.get("feature_uid"), str):
            features[entry["feature_uid"]].add(location)
        if isinstance(entry.get("insertion_id"), str):
            insertions[entry["insertion_id"]] = location
    return features, insertions, candidates


def diff_manifest(
    manifest_path: Path,
    *,
    loci: set[str] | None = None,
    sources: set[str] | None = None,
    qualifiers: set[str] | None = None,
    statuses: set[str] | None = None,
) -> dict[str, Any]:
    try:
        payload = json.loads(read_input_bytes(manifest_path, "merge manifest"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise MergeError(f"invalid merge manifest JSON: {exc}") from exc
    if (
        not isinstance(payload, dict)
        or payload.get("schema") != "enrich-bakta.merge-manifest.v2"
    ):
        raise MergeError("diff requires enrich-bakta.merge-manifest.v2 JSON")
    features, insertions, candidates = _entry_locations(
        _manifest_objects(payload, "entries")
    )
    decisions = _manifest_objects(payload, "decisions")

    # Insertion entries also locate decisions whose source evidence rows lack UIDs.
    for row in decisions:
        locations = {
            insertions[item]
            for item in _strings(row, "insertion_ids")
            if item in insertions
        }
        if len(locations) == 1:
            for uid in _strings(row, "target_feature_uids"):
                features[uid].update(locations)

    groups: dict[tuple[str, str, tuple[str, ...]], list[dict[str, Any]]] = defaultdict(
        list
    )
    added_candidates = set()
    for row in decisions:
        candidate, source, qualifier, status, value, reason = (
            row.get(key)
            for key in (
                "candidate_id",
                "source_id",
                "qualifier",
                "final_status",
                "normalized_value",
                "reason_code",
            )
        )
        if not all(
            isinstance(item, str)
            for item in (candidate, source, qualifier, status, value, reason)
        ):
            raise MergeError(
                "merge manifest decision is missing a required string field"
            )
        if (
            (sources and source not in sources)
            or (qualifiers and qualifier not in qualifiers)
            or (statuses and status not in statuses)
        ):
            continue
        uids = _strings(row, "target_feature_uids")
        locations = {location for uid in uids for location in features.get(uid, ())}
        if not locations:
            locations = {candidates.get(str(candidate), ("", ""))}
        diagnostics = row.get("confidence_diagnostics", {})
        if not isinstance(diagnostics, dict):
            raise MergeError("merge manifest confidence_diagnostics must be an object")
        decision = {
            "candidate_id": candidate,
            "source_id": source,
            "qualifier": qualifier,
            "normalized_value": value,
            "final_status": status,
            "reason_code": reason,
            "supporting_candidate_ids": _strings(row, "supporting_candidate_ids"),
            "support_witness": diagnostics.get("support_witness"),
        }
        for record, locus in locations:
            if loci and locus not in loci:
                continue
            groups[record, locus, uids].append(decision)
            if status in _ADDED_STATUSES:
                added_candidates.add(candidate)
    return {
        "schema": "enrich-bakta.diff.v1",
        "mode": "manifest",
        "added_decisions": len(added_candidates),
        "groups": [
            {
                "record_id": record,
                "locus_tag": locus,
                "target_feature_uids": uids,
                "decisions": sorted(
                    group,
                    key=lambda row: (
                        row["qualifier"],
                        row["source_id"],
                        row["candidate_id"],
                    ),
                ),
            }
            for (record, locus, uids), group in sorted(groups.items())
        ],
    }


_COLUMNS = (
    "record_id",
    "locus_tag",
    "feature_type",
    "location",
    "target_feature_uids",
    "source_id",
    "qualifier",
    "status",
    "value",
    "candidate_id",
    "reason_code",
    "supporting_candidate_ids",
    "support_witness",
)


def _flat_rows(report: dict[str, Any]) -> Iterator[dict[str, Any]]:
    if report["mode"] == "manifest":
        for group in report["groups"]:
            for row in group["decisions"]:
                yield {
                    **{
                        key: group[key]
                        for key in ("record_id", "locus_tag", "target_feature_uids")
                    },
                    "source_id": row["source_id"],
                    "qualifier": row["qualifier"],
                    "status": row["final_status"],
                    "value": row["normalized_value"],
                    "candidate_id": row["candidate_id"],
                    "reason_code": row["reason_code"],
                    "supporting_candidate_ids": row["supporting_candidate_ids"],
                    "support_witness": row["support_witness"],
                }
        return
    for row in report["deltas"]:
        for status in ("added", "removed"):
            for qualifier, values in row[status].items():
                for value in values:
                    yield {
                        **{
                            key: row[key]
                            for key in (
                                "record_id",
                                "locus_tag",
                                "feature_type",
                                "location",
                            )
                        },
                        "status": status,
                        "qualifier": qualifier,
                        "value": value,
                    }
    for row in report["comment_deltas"]:
        for status in ("added", "removed"):
            for value in row[status]:
                yield {
                    "record_id": row["record_id"],
                    "status": status,
                    "qualifier": "COMMENT",
                    "value": value,
                }


def _render(report: dict[str, Any], output_format: str) -> bytes:
    if output_format == "json":
        return json_sidecar_bytes(report)
    rows = _flat_rows(report)
    if output_format == "tsv":
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, _COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, ensure_ascii=False, sort_keys=True)
                    if isinstance(value, (dict, tuple, list))
                    else value
                    for key, value in row.items()
                }
            )
        return buffer.getvalue().encode("utf-8")

    output = []
    previous = None
    for row in rows:
        heading = tuple(row.get(key) for key in _COLUMNS[:5])
        if heading != previous:
            output.append(
                " / ".join(
                    [
                        str(heading[0] or "(unknown record)"),
                        str(heading[1] or "(no locus)"),
                    ]
                    + [str(part) for part in heading[2:] if part]
                )
            )
            previous = heading
        details = " ".join(
            f"{key}={json.dumps(row[key], ensure_ascii=False, sort_keys=True)}"
            for key in (
                "source_id",
                "candidate_id",
                "reason_code",
                "supporting_candidate_ids",
                "support_witness",
            )
            if row.get(key) is not None
        )
        value = json.dumps(row["value"], ensure_ascii=False)
        output.append(
            f"  {row['status']} {row['qualifier']}={value} {details}".rstrip()
        )
    return (
        "\n".join(output or ["No matching annotations or decisions."]) + "\n"
    ).encode("utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--manifest", type=Path)
    mode.add_argument("--base", type=Path)
    parser.add_argument("--enriched", type=Path)
    parser.add_argument("--locus", action="append")
    parser.add_argument("--source", action="append")
    parser.add_argument("--qualifier", action="append")
    parser.add_argument("--status", action="append")
    parser.add_argument("--format", choices=("text", "json", "tsv"), default="text")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fail-on-diff", action="store_true")
    args = parser.parse_args(argv)
    if args.manifest and args.enriched:
        parser.error("--enriched is only valid with --base")
    if args.base and (args.enriched is None or args.source or args.status):
        parser.error(
            "--base requires --enriched and does not accept --source or --status"
        )
    filters = {"loci": set(args.locus or []), "qualifiers": set(args.qualifier or [])}
    try:
        if args.manifest:
            inputs = [args.manifest]
            report = diff_manifest(
                args.manifest,
                **filters,
                sources=set(args.source or []),
                statuses=set(args.status or []),
            )
            changed = report["added_decisions"]
        else:
            inputs = [args.base, args.enriched]
            report = diff_gbff_files(args.base, args.enriched, **filters)
            changed = report["changed_features"] or report["changed_records_comments"]
        rendered = _render(report, args.format)
        if args.output:
            stage_artifacts([(args.output, rendered)], inputs)
        else:
            print(rendered.decode("utf-8"), end="")
    except (MergeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return int(bool(args.fail_on_diff and changed))


if __name__ == "__main__":
    raise SystemExit(main())
