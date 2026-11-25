#!/usr/bin/env python3
"""Verify the installed audit, report, and diff commands on the published C14 data."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

from jsonschema import Draft202012Validator


def run_checks(root: Path) -> None:
    output = root / ".test-output/c14-cli-smoke"
    output.mkdir(parents=True, exist_ok=True)
    gates = root / ".test-output/interproscan-gates"
    prefix = "data/C14/evidence/interproscan/C14-NMZ.interproscan"
    exports = {kind: f"{prefix}.{kind}" for kind in ("tsv", "json", "gff3")}
    version = "5.59-91.0"
    audit_path = output / "C14-audit.json"
    subprocess.run(
        [
            "enrich-bakta-audit-interproscan",
            "--tsv",
            exports["tsv"],
            "--json",
            exports["json"],
            "--gff3",
            exports["gff3"],
            "--interproscan-version",
            version,
            "--strict",
            "--output",
            str(audit_path),
            "--markdown",
            str(output / "C14-audit.md"),
        ],
        cwd=root,
        check=True,
    )
    audit = json.loads(audit_path.read_bytes())
    schema = json.loads(
        (root / "schemas/interproscan-audit.v1.schema.json").read_bytes()
    )
    Draft202012Validator(schema).validate(audit)
    published = json.loads((root / "docs/data/PUBLISHED-MANIFEST.json").read_bytes())
    hashes = {row["path"]: row["sha256"] for row in published["artifacts"]}
    for kind, relative in exports.items():
        assert audit["inputs"][kind]["sha256"] == hashes[relative], kind
    assert audit["status"] == "passed"
    assert audit["tsv_summary"]["occurrences_total"] == 40234
    assert audit["json_audit"]["protein_objects"] == 4634
    assert audit["json_audit"]["total_xrefs"] == 4636
    assert audit["json_audit"]["multi_xref_proteins"]["count"] == 2
    assert audit["json_audit"]["md5_mismatches"] == 0
    for audit_mode in ("json_audit", "gff3_audit"):
        assert audit[audit_mode]["missing_tsv_occurrences"]["count"] == 0, audit_mode
        assert audit[audit_mode]["extra_occurrences_by_analysis"] == {"pirsr": 1944}, (
            audit_mode
        )

    manifest = gates / "C14-interpro-enriched.manifest.json"
    context = gates / "C14-interpro-context.json"
    enriched = gates / "C14-interpro-enriched.gbff"
    base = root / "data/C14/bakta/C14-NMZ.gbff"
    summary = output / "C14-summary.html"
    subprocess.run(
        [
            "enrich-bakta-report",
            "--manifest",
            str(manifest),
            "--context-report",
            str(context),
            "--output",
            str(summary),
        ],
        cwd=root,
        check=True,
    )
    html = summary.read_text(encoding="utf-8")
    expected = json.loads(
        (root / "docs/data/interproscan-scientific-baseline.json").read_bytes()
    )
    gbff_hash = expected["samples"]["C14"]["expected"]["interproscan_only"][
        "output_sha256"
    ]
    for marker in ("<!doctype html>", version, "Pathway occurrences", gbff_hash):
        assert marker in html, marker

    with TemporaryDirectory(dir=output, prefix="diff-") as directory:
        scratch = Path(directory).resolve()
        assert scratch.parent == output.resolve()
        diff_path = scratch / "diff.json"
        subprocess.run(
            [
                "enrich-bakta-diff",
                "--base",
                str(base),
                "--enriched",
                str(enriched),
                "--format",
                "json",
                "--output",
                str(diff_path),
            ],
            cwd=root,
            check=True,
        )
        gbff_diff = json.loads(diff_path.read_bytes())
        assert gbff_diff["changed_features"] == 4390
        assert gbff_diff["changed_records_comments"] == 1

        subprocess.run(
            [
                "enrich-bakta-diff",
                "--manifest",
                str(manifest),
                "--format",
                "json",
                "--output",
                str(diff_path),
            ],
            cwd=root,
            check=True,
        )
        manifest_diff = json.loads(diff_path.read_bytes())
        assert len(manifest_diff["groups"]) == 4390
        assert manifest_diff["added_decisions"] == 34079
        assert all(group["locus_tag"] for group in manifest_diff["groups"])

    receipt = {
        "status": "passed",
        "published_input_hashes_verified": len(exports),
        "report_output_sha256": gbff_hash,
        "gbff_changed_features": gbff_diff["changed_features"],
        "gbff_changed_record_comments": gbff_diff["changed_records_comments"],
        "manifest_locus_groups": len(manifest_diff["groups"]),
        "manifest_added_decisions": manifest_diff["added_decisions"],
    }
    (output / "C14-smoke.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    run_checks(Path(__file__).resolve().parents[1])
