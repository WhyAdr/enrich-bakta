"""Deterministic views of finalized enrichment manifests and context reports."""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.artifacts import preflight_artifact_paths, stage_artifacts
from enrich_bakta_lib.core.types import (
    TOOL_VERSION,
    MergeError,
    read_input_bytes,
    sha256_bytes,
)

QUALIFIERS = ("gene", "EC_number", "db_xref", "note", "inference")
ACCEPTED = ("emitted", "shared_support", "supported_existing")
SOURCES = (
    ("Baktfold", "baktfold", "baktfold_sha256"),
    ("KofamScan", "kofamscan", "kofam_sha256"),
    ("eggNOG", "eggnog", "eggnog_sha256"),
    ("InterProScan", "interproscan", "interproscan_sha256"),
)
OVERLAP_CAVEAT = (
    "Shared support across bioinformatic tools may stem from common upstream reference "
    "databases and does not constitute independent experimental confirmation."
)
IMPORT_CAVEAT = (
    "Imported FAA translations reflect explicit protein-artifact restoration, "
    "not independent genomic CDS validation."
)


@dataclass(frozen=True)
class SummaryMetrics:
    metadata: dict[str, Any]
    insertion_counts: Counter[str]
    by_source: dict[str, Counter[str]]
    by_qualifier_source: dict[tuple[str, str], Counter[str]]
    cross_source_overlaps: Counter[tuple[str, str]]
    gene_suppressions: Counter[str]
    suppressed_samples: list[dict[str, Any]]
    anti_fam_hits: list[dict[str, Any]]
    translation_evidence: dict[str, dict[str, Any]]
    interproscan_context_summary: dict[str, Any] | None
    eggnog_context_summary: dict[str, Any] | None


def resolve_report_format(path: Path, explicit: str = "auto") -> str:
    if explicit in {"md", "html"}:
        return explicit
    if explicit != "auto":
        raise MergeError(f"unsupported summary report format {explicit!r}")
    suffix = path.suffix.lower()
    if suffix in {".md", ".markdown"}:
        return "md"
    if suffix in {".html", ".htm"}:
        return "html"
    raise MergeError(
        f"unsupported summary report extension {path.suffix!r}; use .md or .html"
    )


def _object(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise MergeError(f"{label} must be a JSON object")
    return value


def _rows(value: object, label: str) -> list[dict[str, Any]]:
    if not isinstance(value, list) or any(not isinstance(row, dict) for row in value):
        raise MergeError(f"{label} must be an array of JSON objects")
    return value


def _contexts(
    payload: dict[str, Any] | None,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    if payload is None:
        return None, None
    schema = payload.get("schema")
    ips_schema = "enrich-bakta.interproscan-context.v1"
    egg_schema = "enrich-bakta.eggnog-context.v1"
    if schema == ips_schema:
        return payload, None
    if schema == egg_schema:
        return None, payload
    if schema == "enrich-bakta.context.v1":
        ips = _object(payload.get("interproscan"), "InterProScan context")
        egg = _object(payload.get("eggnog"), "eggNOG context")
        if ips.get("schema") == ips_schema and egg.get("schema") == egg_schema:
            return ips, egg
        raise MergeError("combined context contains unsupported source schemas")
    raise MergeError(f"unsupported context report schema: {schema!r}")


def build_summary_metrics(
    manifest: dict[str, Any],
    context: dict[str, Any] | None = None,
) -> SummaryMetrics:
    if manifest.get("schema") != "enrich-bakta.merge-manifest.v2":
        raise MergeError("summary report requires enrich-bakta.merge-manifest.v2 JSON")
    metadata = _object(manifest.get("metadata"), "manifest metadata")
    entries = _rows(manifest.get("entries"), "manifest entries")
    decisions = _rows(manifest.get("decisions"), "manifest decisions")
    by_source: dict[str, Counter[str]] = defaultdict(Counter)
    by_qualifier_source: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    overlap_groups: dict[tuple[tuple[str, ...], str, str], set[str]] = defaultdict(set)
    gene_suppressions: Counter[str] = Counter()
    suppressed = []
    for decision in decisions:
        source, qualifier, status = (
            decision.get("source_id"),
            decision.get("qualifier"),
            decision.get("final_status"),
        )
        if (
            not isinstance(source, str)
            or not isinstance(qualifier, str)
            or not isinstance(status, str)
            or not (source and qualifier and status)
        ):
            raise MergeError(
                "manifest decision requires source_id, qualifier, and final_status"
            )
        by_source[source][status] += 1
        by_source[source]["total"] += 1
        if qualifier in QUALIFIERS:
            by_qualifier_source[(qualifier, source)][status] += 1
        if status == "shared_support":
            targets = decision.get("target_feature_uids")
            normalized = decision.get("normalized_value")
            if (
                not isinstance(targets, list)
                or any(not isinstance(uid, str) for uid in targets)
                or not isinstance(normalized, str)
            ):
                raise MergeError(
                    "shared-support decision requires target feature UIDs and normalized value"
                )
            overlap_groups[(tuple(sorted(targets)), qualifier, normalized)].add(source)
        if qualifier == "gene" and status in {
            "suppressed_conflict",
            "suppressed_authoritative_gene",
        }:
            gene_suppressions[status] += 1
            suppressed.append(decision)
    overlaps: Counter[tuple[str, str]] = Counter()
    for (_, qualifier, _), sources in overlap_groups.items():
        if len(sources) > 1:
            overlaps[(qualifier, " + ".join(sorted(sources)))] += 1

    # Decisions carry feature UIDs. Their source entries carry the display locus.
    suppressed_ids = {
        str(row["candidate_id"]) for row in suppressed if "candidate_id" in row
    }
    locations = {}
    for entry in entries:
        candidate_id = entry.get("candidate_id")
        if isinstance(candidate_id, str) and candidate_id in suppressed_ids:
            locus = entry.get("locus_tag") or entry.get("query_id")
            if locus:
                locations[candidate_id] = (entry.get("record", ""), locus)
    for index, row in enumerate(suppressed):
        if (location := locations.get(str(row.get("candidate_id", "")))) is not None:
            suppressed[index] = {**row, "record": location[0], "locus_tag": location[1]}

    insertion_counts: Counter[str] = Counter()
    for entry in entries:
        if entry.get("entry_type") != "insertion":
            continue
        if entry.get("qualifier") == "COMMENT":
            insertion_counts["Comment"] += 1
        elif entry.get("candidate_role") == "producer_provenance":
            insertion_counts["Provenance"] += 1
        else:
            insertion_counts["Functional and substantive evidence"] += 1

    ips, egg = _contexts(context)
    anti_fam: list[dict[str, Any]] = []
    ips_summary = None
    if ips is not None:
        meta = _object(ips.get("metadata"), "InterProScan context metadata")
        qc = _object(meta.get("qc_flags", {}), "InterProScan QC flags")
        anti_fam = _rows(qc.get("anti_fam", []), "AntiFam alerts")
        pathway = _object(
            meta.get("pathway_counters", {}), "InterProScan pathway counters"
        )
        ips_summary = {
            "pathway_total": pathway.get("total_occurrences", 0),
            "pathway_by_database": _object(
                pathway.get("by_database", {}), "pathway databases"
            ),
            **{
                key: meta.get(key, "")
                for key in ("support_digest", "match_evidence_digest", "pathway_digest")
            },
            "truncated_queries": sum(
                row.get("context_signatures_truncated") is True
                for row in _rows(ips.get("entries"), "InterProScan context entries")
            ),
        }

    translation = {
        locus: _object(value, "translation evidence entry")
        for locus, value in _object(
            metadata.get("translation_evidence", {}), "translation evidence"
        ).items()
    }
    return SummaryMetrics(
        metadata=metadata,
        insertion_counts=insertion_counts,
        by_source=dict(by_source),
        by_qualifier_source=dict(by_qualifier_source),
        cross_source_overlaps=overlaps,
        gene_suppressions=gene_suppressions,
        suppressed_samples=sorted(
            suppressed,
            key=lambda row: tuple(
                str(row.get(key, ""))
                for key in ("record", "locus_tag", "source_id", "candidate_id")
            ),
        )[:50],
        anti_fam_hits=sorted(
            anti_fam,
            key=lambda row: tuple(
                str(row.get(key, ""))
                for key in ("query_id", "signature_accession", "start", "end")
            ),
        ),
        translation_evidence=translation,
        interproscan_context_summary=ips_summary,
        eggnog_context_summary=_object(egg.get("metadata"), "eggNOG context metadata")
        if egg is not None
        else None,
    )


@dataclass(frozen=True)
class _Section:
    title: str
    headers: tuple[str, ...]
    rows: list[tuple[Any, ...]]
    note: str = ""


def _sections(metrics: SummaryMetrics) -> list[_Section]:
    meta = metrics.metadata
    overview = [
        (label, meta.get(key, ""))
        for label, key in (
            ("Tool version", "tool_version"),
            ("Operation", "operation"),
            ("Base SHA-256", "base_sha256"),
            ("Output SHA-256", "output_sha256"),
            ("Records", "records"),
            ("Insertions", "insertions"),
        )
    ]
    overview.extend(
        (name, metrics.insertion_counts[name])
        for name in ("Functional and substantive evidence", "Provenance", "Comment")
    )
    overview.append(
        ("Merge timestamp", meta.get("merge_timestamp") or "omitted (deterministic)")
    )
    source_hashes = _object(meta.get("source_hashes", {}), "source hashes")
    source_bindings = [
        (label, meta.get(f"{key}_version") or "unspecified", digest)
        for label, key, hash_key in SOURCES
        if (digest := meta.get(hash_key) or source_hashes.get(label) or "")
    ]

    statuses = sorted(
        {
            status
            for counts in metrics.by_source.values()
            for status in counts
            if status.startswith("suppressed_")
        }
    )
    columns = (
        *ACCEPTED,
        *statuses,
        "filtered_confidence",
        "skipped_partial_ec",
        "invalid_value",
        "total",
    )
    results = [
        _Section("Run overview", ("Metric", "Value"), overview),
        _Section(
            "Source provenance", ("Source", "Version", "SHA-256"), source_bindings
        ),
        _Section(
            "Per-source annotation gain",
            ("Source", *columns),
            [
                (source, *(counts[status] for status in columns))
                for source, counts in sorted(metrics.by_source.items())
            ],
        ),
        _Section(
            "Per-qualifier breakdown",
            ("Qualifier", "Source", *ACCEPTED, "suppressed"),
            [
                (
                    qualifier,
                    source,
                    *(counts[status] for status in ACCEPTED),
                    sum(
                        count
                        for status, count in counts.items()
                        if status.startswith("suppressed_")
                    ),
                )
                for (qualifier, source), counts in sorted(
                    metrics.by_qualifier_source.items()
                )
            ],
        ),
        _Section(
            "Cross-source overlap",
            ("Qualifier", "Sources", "Unique annotations"),
            [
                (qualifier, sources, count)
                for (qualifier, sources), count in sorted(
                    metrics.cross_source_overlaps.items()
                )
            ],
            OVERLAP_CAVEAT,
        ),
        _Section(
            "QC, conflicts and suppressions",
            ("Metric", "Value"),
            [
                ("Gene conflict policy", meta.get("gene_conflict_policy", "")),
                ("Gene conflicts", metrics.gene_suppressions["suppressed_conflict"]),
                (
                    "Authoritative gene collisions",
                    metrics.gene_suppressions["suppressed_authoritative_gene"],
                ),
                (
                    "Baktfold invalid EC policy",
                    meta.get("baktfold_invalid_ec_policy", ""),
                ),
                (
                    "Skipped partial ECs",
                    sum(
                        counts["skipped_partial_ec"]
                        for counts in metrics.by_source.values()
                    ),
                ),
                (
                    "Invalid EC tokens",
                    sum(
                        counts["invalid_value"]
                        for (
                            qualifier,
                            _,
                        ), counts in metrics.by_qualifier_source.items()
                        if qualifier == "EC_number"
                    ),
                ),
                (
                    "AntiFam alerts in context",
                    len(metrics.anti_fam_hits)
                    if metrics.interproscan_context_summary is not None
                    else "context unavailable",
                ),
            ],
        ),
        _Section(
            "Gene suppression sample (up to 50)",
            ("Record", "Locus", "Source", "Status", "Reason"),
            [
                tuple(
                    row.get(key, "")
                    for key in (
                        "record",
                        "locus_tag",
                        "source_id",
                        "final_status",
                        "reason_code",
                    )
                )
                for row in metrics.suppressed_samples
            ],
        ),
        _Section(
            "AntiFam QC alerts",
            ("Query", "Signature", "Description", "Score", "Coordinates"),
            [
                (
                    row.get("query_id"),
                    row.get("signature_accession"),
                    row.get("signature_description"),
                    row.get("score"),
                    f"{row.get('start', '')}..{row.get('end', '')}",
                )
                for row in metrics.anti_fam_hits
            ],
        ),
        _Section(
            "Translation lineage audit",
            (
                "Locus",
                "Origin",
                "Restoration policy",
                "Genomic validation",
                "Bound output SHA-256",
            ),
            [
                (
                    locus,
                    *(
                        row.get(key)
                        for key in (
                            "origin",
                            "restoration_policy",
                            "genomic_validation_status",
                            "bound_output_sha256",
                        )
                    ),
                )
                for locus, row in sorted(metrics.translation_evidence.items())
            ],
            IMPORT_CAVEAT
            if any(
                row.get("origin") == "imported_faa"
                for row in metrics.translation_evidence.values()
            )
            else "",
        ),
        _Section(
            "Translation evidence settings",
            ("Metric", "Value"),
            [
                (
                    "Parent manifest SHA-256",
                    meta.get("translation_evidence_parent_sha256", ""),
                ),
                (
                    "Allow imported translations",
                    meta.get("allow_imported_translations", ""),
                ),
            ],
        ),
    ]
    if (ips := metrics.interproscan_context_summary) is not None:
        results.append(
            _Section(
                "InterProScan context",
                ("Metric", "Value"),
                [
                    ("Pathway occurrences", ips["pathway_total"]),
                    *(
                        (f"Pathway occurrences ({database})", count)
                        for database, count in sorted(
                            ips["pathway_by_database"].items()
                        )
                    ),
                    ("Support digest", ips["support_digest"]),
                    ("Match evidence digest", ips["match_evidence_digest"]),
                    ("Pathway digest", ips["pathway_digest"]),
                    (
                        "Queries with truncated signature lists",
                        ips["truncated_queries"],
                    ),
                ],
                "Reactome taxon labels describe database hits and do not establish bacterial pathway completeness.",
            )
        )
    if (egg := metrics.eggnog_context_summary) is not None:
        results.append(
            _Section(
                "eggNOG context",
                ("Metric", "Value"),
                [
                    ("Context hits", egg.get("context_hits", 0)),
                    ("Context values", egg.get("context_values", 0)),
                ],
            )
        )
    return results


def _display(value: object) -> str:
    return "" if value is None else str(value)


def _markdown_cell(value: object) -> str:
    escaped = html.escape(" ".join(_display(value).splitlines()), quote=False)
    return re.sub(r"([\\`*_\[\]])", r"\\\1", escaped).replace("|", "&#124;")


def render_markdown(metrics: SummaryMetrics) -> str:
    lines = ["# Enrichment summary"]
    for section in _sections(metrics):
        lines.extend(("", f"## {section.title}", ""))
        if section.note:
            lines.extend((section.note, ""))
        lines.append("| " + " | ".join(map(_markdown_cell, section.headers)) + " |")
        lines.append("| " + " | ".join(["---"] * len(section.headers)) + " |")
        lines.extend(
            "| " + " | ".join(map(_markdown_cell, row)) + " |" for row in section.rows
        )
    return "\n".join(lines) + "\n"


def render_html(metrics: SummaryMetrics) -> str:
    lines = [
        '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
        "<title>Enrichment summary</title>",
        "<style>body{font:1rem system-ui,sans-serif;max-width:90rem;margin:2rem auto;padding:0 1rem;color:#16212b}"
        "table{border-collapse:collapse;margin:1rem 0;width:100%;overflow-wrap:anywhere}"
        "th,td{border:1px solid #cbd5e1;padding:.45rem;text-align:left;vertical-align:top}"
        "th{background:#e8eef5}section{overflow-x:auto}</style><body><h1>Enrichment summary</h1>",
    ]
    for section in _sections(metrics):
        lines.append(f"<section><h2>{html.escape(section.title)}</h2>")
        if section.note:
            lines.append(f"<p>{html.escape(section.note)}</p>")
        headings = "".join(
            f'<th scope="col">{html.escape(heading)}</th>'
            for heading in section.headers
        )
        lines.append(f"<table><thead><tr>{headings}</tr></thead><tbody>")
        for row in section.rows:
            lines.append(
                "<tr>"
                + "".join(f"<td>{html.escape(_display(value))}</td>" for value in row)
                + "</tr>"
            )
        lines.append("</tbody></table></section>")
    lines.append("</body></html>")
    return "\n".join(lines) + "\n"


def summary_report_bytes(
    manifest: dict[str, Any],
    context: dict[str, Any] | None,
    report_format: str,
) -> bytes:
    metrics = build_summary_metrics(manifest, context)
    renderer = render_markdown if report_format == "md" else render_html
    return renderer(metrics).encode("utf-8")


def _read_json(data: bytes, label: str) -> dict[str, Any]:
    try:
        return _object(json.loads(data), label)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise MergeError(f"invalid {label} JSON: {exc}") from exc


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render a report from an existing enrichment manifest."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--context-report", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--format", choices=("auto", "md", "html"), default="auto")
    args = parser.parse_args()
    inputs = [args.manifest]
    if args.context_report is not None:
        inputs.append(args.context_report)
    try:
        preflight_artifact_paths([args.output], inputs)
        report_format = resolve_report_format(args.output, args.format)
        manifest = _read_json(read_input_bytes(args.manifest, "manifest"), "manifest")
        context = None
        if args.context_report is not None:
            context_data = read_input_bytes(args.context_report, "context report")
            metadata = _object(manifest.get("metadata"), "manifest metadata")
            if "context_report_sha256" in metadata and metadata[
                "context_report_sha256"
            ] != sha256_bytes(context_data):
                raise MergeError("context report SHA-256 does not match the manifest")
            context = _read_json(context_data, "context report")
        stage_artifacts(
            [(args.output, summary_report_bytes(manifest, context, report_format))],
            inputs,
        )
    except (MergeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
