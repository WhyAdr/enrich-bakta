from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest
from test_merge_pipeline import eggnog_bytes, record_bytes, write

from enrich_bakta_lib.core.types import MergeError, sha256_bytes
from enrich_bakta_lib.workflows.enrich import enrich
from enrich_bakta_lib.workflows.enrich import main as enrich_main
from enrich_bakta_lib.workflows.summary_report import (
    build_summary_metrics,
    main,
    render_html,
)


@pytest.fixture
def merge_inputs(tmp_path: Path) -> dict:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    egg = write(
        tmp_path / "egg.tsv",
        eggnog_bytes(
            "T_0001\tseed\t1e-4\t10\tC\tabc\tGO:0000001\t1.2.3.4\tK00001\t-\thhhhhhhhhhhhh"
        ),
    )
    md5 = hashlib.md5(b"MK").hexdigest()
    ips = write(
        tmp_path / "ips.tsv",
        (
            f"T_0001\t{md5}\t2\tPfam\tPF00001\tdesc\t1\t2\t1.0\tT\t25-08-2026\tIPR000001\tipr desc\tGO:0000001\tReactome:R-XYZ-123\n"
            f"T_0001\t{md5}\t2\tAntiFam\tANF00001\t<script>alert(1)</script>\t1\t2\t1.0\tT\t25-08-2026\t-\t-\n"
        ).encode(),
    )
    return dict(
        bakta_path=base,
        faa_path=faa,
        eggnog_path=egg,
        interproscan_path=ips,
        interproscan_version="5.59-91.0",
        output_path=tmp_path / "output.gbff",
        manifest_path=tmp_path / "manifest.json",
        context_report_path=tmp_path / "context.json",
    )


@pytest.mark.parametrize("suffix", [".markdown", ".htm"])
def test_inline_bytes_are_deterministic_and_existing_artifacts_identical(
    merge_inputs: dict, tmp_path: Path, suffix: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    stats = enrich(**merge_inputs)
    paths = [
        merge_inputs[key]
        for key in ("output_path", "manifest_path", "context_report_path")
    ]
    expected = [path.read_bytes() for path in paths]
    report = tmp_path / f"summary{suffix}"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "enrich-bakta",
            "--bakta",
            str(merge_inputs["bakta_path"]),
            "--faa",
            str(merge_inputs["faa_path"]),
            "--eggnog",
            str(merge_inputs["eggnog_path"]),
            "--interproscan",
            str(merge_inputs["interproscan_path"]),
            "--interproscan-version",
            merge_inputs["interproscan_version"],
            "--output",
            str(paths[0]),
            "--manifest",
            str(paths[1]),
            "--context-report",
            str(paths[2]),
            "--summary-report",
            str(report),
        ],
    )
    assert enrich_main() == 0
    assert [path.read_bytes() for path in paths] == expected
    report_bytes = report.read_bytes()
    assert enrich(**merge_inputs, summary_report_path=report) == stats
    assert report.read_bytes() == report_bytes
    for text in (
        "Cross-source overlap",
        "InterProScan + eggNOG",
        "eggNOG context",
        "InterProScan context",
        "does not constitute independent experimental confirmation",
        "omitted (deterministic)",
        stats["output_sha256"],
    ):
        assert text.encode() in report_bytes
    assert b"<script>" not in report_bytes
    assert b"&lt;script&gt;alert(1)&lt;/script&gt;" in report_bytes


def test_inline_without_manifest_or_context_sidecar(
    merge_inputs: dict, tmp_path: Path
) -> None:
    report = tmp_path / "summary.md"
    result = enrich(
        **{**merge_inputs, "manifest_path": None, "context_report_path": None},
        summary_report_path=report,
    )
    assert result["output_sha256"] in report.read_text(encoding="utf-8")
    assert "InterProScan context" in report.read_text(encoding="utf-8")
    assert not merge_inputs["manifest_path"].exists()
    assert not merge_inputs["context_report_path"].exists()


@pytest.mark.parametrize("context_kind", ["combined", "interproscan", "eggnog"])
def test_offline_reports(
    merge_inputs: dict,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    context_kind: str,
) -> None:
    if context_kind != "combined":
        drop_source = "eggnog" if context_kind == "interproscan" else "interproscan"
        merge_inputs[f"{drop_source}_path"] = None
        if drop_source == "interproscan":
            merge_inputs["interproscan_version"] = None
    inline = tmp_path / "inline.html"
    enrich(**merge_inputs, summary_report_path=inline)
    offline = tmp_path / "standalone.report"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "enrich-bakta-report",
            "--manifest",
            str(merge_inputs["manifest_path"]),
            "--context-report",
            str(merge_inputs["context_report_path"]),
            "--output",
            str(offline),
            "--format",
            "html",
        ],
    )
    assert main() == 0
    assert inline.read_bytes() == offline.read_bytes()
    assert (
        b"eggNOG context" in offline.read_bytes()
        if context_kind != "interproscan"
        else b"eggNOG context" not in offline.read_bytes()
    )


def test_offline_rejects_wrong_context_hash_and_input_collision(
    merge_inputs: dict,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    enrich(**merge_inputs)
    manifest, context = (
        merge_inputs[key] for key in ("manifest_path", "context_report_path")
    )
    report = tmp_path / "summary.md"
    context_data = context.read_bytes()
    context.write_bytes(context_data + b" ")
    argv = [
        "enrich-bakta-report",
        "--manifest",
        str(manifest),
        "--context-report",
        str(context),
        "--output",
    ]
    monkeypatch.setattr(sys, "argv", [*argv, str(report)])
    assert main() == 2
    assert "SHA-256 does not match" in capsys.readouterr().err
    assert not report.exists()
    context.write_bytes(context_data)
    monkeypatch.setattr(sys, "argv", [*argv, str(context)])
    assert main() == 2
    assert "aliases" in capsys.readouterr().err
    assert context.read_bytes() == context_data


@pytest.mark.parametrize(
    "key", ["bakta_path", "output_path", "manifest_path", "context_report_path"]
)
def test_inline_refuses_collisions_before_writing(merge_inputs: dict, key: str) -> None:
    path = merge_inputs[key]
    if key == "bakta_path":
        expected = sha256_bytes(path.read_bytes())
    with pytest.raises(MergeError, match="aliases"):
        enrich(**merge_inputs, summary_report_path=path)
    assert not merge_inputs["output_path"].exists()
    if key == "bakta_path":
        assert sha256_bytes(path.read_bytes()) == expected


def test_rejects_unsupported_inline_extension(
    merge_inputs: dict, tmp_path: Path
) -> None:
    with pytest.raises(MergeError, match="unsupported summary report extension"):
        enrich(**merge_inputs, summary_report_path=tmp_path / "summary.txt")
    assert not merge_inputs["manifest_path"].exists()


def test_summary_counts_unique_overlap_and_reports_translation_caveat() -> None:
    shared = dict(
        final_status="shared_support",
        qualifier="db_xref",
        normalized_value="GO:0000001",
        target_feature_uids=["feature:1"],
        source_id="InterProScan",
    )
    manifest = {
        "schema": "enrich-bakta.merge-manifest.v2",
        "metadata": {
            "translation_evidence": {
                "T_0001": {
                    "origin": "imported_faa",
                    "restoration_policy": "import-faa",
                    "genomic_validation_status": "unsupported",
                    "bound_output_sha256": "a" * 64,
                }
            },
            "translation_evidence_parent_sha256": "b" * 64,
            "allow_imported_translations": True,
        },
        "entries": [
            {
                "candidate_id": "candidate:gene",
                "record": "TEST",
                "locus_tag": "T_0002",
            },
            {
                "entry_type": "insertion",
                "qualifier": "inference",
                "candidate_role": "producer_provenance",
            },
            {"entry_type": "insertion", "qualifier": "COMMENT"},
            {
                "entry_type": "insertion",
                "qualifier": "db_xref",
                "candidate_role": "functional_proposal",
            },
        ],
        "decisions": [
            shared,
            shared,
            {**shared, "source_id": "eggNOG"},
            {
                "candidate_id": "candidate:gene",
                "source_id": "Baktfold",
                "qualifier": "gene",
                "final_status": "suppressed_authoritative_gene",
            },
        ],
    }
    metrics = build_summary_metrics(manifest)
    assert metrics.cross_source_overlaps == {("db_xref", "InterProScan + eggNOG"): 1}
    assert sum(metrics.insertion_counts.values()) == 3
    assert metrics.by_source["InterProScan"]["shared_support"] == 2
    assert metrics.suppressed_samples[0]["record"] == "TEST"
    assert metrics.suppressed_samples[0]["locus_tag"] == "T_0002"
    report = render_html(metrics)
    for text in (
        "T_0001",
        "imported_faa",
        "import-faa",
        "unsupported",
        "a" * 64,
        "b" * 64,
        "not independent genomic CDS validation",
    ):
        assert text in report
    with pytest.raises(MergeError, match="unsupported context report schema"):
        build_summary_metrics(manifest, {"schema": "unexpected"})
    with pytest.raises(MergeError, match="requires enrich-bakta.merge-manifest.v2"):
        build_summary_metrics(json.loads("{}"))
