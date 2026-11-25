from __future__ import annotations

import csv
import io
import json
from pathlib import Path
from typing import Any

import pytest
from test_merge_pipeline import eggnog_bytes, record_bytes, write

from enrich_bakta_lib.core.genbank import feature_uid, parse_genbank_bytes
from enrich_bakta_lib.core.types import MergeError, sha256_bytes
from enrich_bakta_lib.workflows import cds_export
from enrich_bakta_lib.workflows.cds_export import export_cds_tsv, main
from enrich_bakta_lib.workflows.restore_translations import restore


def _rows(data: bytes) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(data.decode("utf-8")), delimiter="\t"))


def _json(path: Path, payload: dict[str, Any]) -> Path:
    return write(path, json.dumps(payload).encode("utf-8"))


def _context(sha256: str) -> dict[str, Any]:
    return {
        "schema": "enrich-bakta.interproscan-context.v1",
        "metadata": {
            "output_sha256": sha256,
            "qc_flags": {
                "anti_fam": [
                    {"query_id": "SAME", "signature_accession": "ANF00001"},
                    {"query_id": "HIT", "signature_accession": "ANF00002"},
                ]
            },
        },
        "entries": [{"query_id": "REVIEWED"}],
    }


def _manifest(sha256: str, *, context_hash: str | None = None) -> dict[str, Any]:
    metadata = {"output_sha256": sha256}
    if context_hash is not None:
        metadata["context_report_sha256"] = context_hash
    return {
        "schema": "enrich-bakta.merge-manifest.v2",
        "metadata": metadata,
        "entries": [],
        "decisions": [],
    }


def test_export_preserves_order_multivalues_and_utf8(tmp_path: Path) -> None:
    first = record_bytes(
        "Z",
        "SAME",
        gene="abcA",
        notes=(
            "PFAM:PF00001.2",
            "TIGRFAM:TIGR00001",
            "KofamScan:K00001;threshold=1;score=2;E-value=3e-4",
            "negative remark about GO:0000002",
        ),
        db_xrefs=("GO:0000001", "InterPro:IPR000001", "GO:0000001"),
        ec_numbers=("1.2.3.4", "1.2.3.-"),
    ).replace(
        b'/product="test protein"',
        '/product="test café"\n                     /inference="protein motif:InterProScan:5.59-91.0"'.encode(),
    )
    gbff = write(tmp_path / "input.gbff", first + record_bytes("A", "SAME"))
    rendered = export_cds_tsv(gbff)
    assert rendered == export_cds_tsv(gbff)
    rows = _rows(rendered)
    assert [row["record_id"] for row in rows] == ["Z", "A"]
    assert [row["record_index"] for row in rows] == ["0", "1"]
    assert [row["feature_index"] for row in rows] == ["2", "2"]
    assert {row["gbff_sha256"] for row in rows} == {sha256_bytes(gbff.read_bytes())}
    assert rows[0]["feature_uid"] != rows[1]["feature_uid"]
    expected = {
        "product": ["test café"],
        "gene": ["abcA"],
        "ec": ["1.2.3.4", "1.2.3.-"],
        "go": ["GO:0000001", "GO:0000001"],
        "ko": ["KEGG:K00001"],
        "interpro": ["InterPro:IPR000001"],
        "pfam": ["PFAM:PF00001.2"],
        "tigrfam": ["TIGRFAM:TIGR00001"],
        "inference": ["protein motif:InterProScan:5.59-91.0"],
        "accepted_decisions": None,
    }
    for key, value in expected.items():
        assert json.loads(rows[0][key]) == value
    assert rows[0]["antifam_status"] == "unavailable"
    assert rows[0]["translation_present"] == "true"
    assert rows[0]["translation_origin"] == "unavailable"
    assert b'"MK"' not in rendered
    assert b"atgaaataa" not in rendered


def test_manifest_joins_exact_uid_and_context_marks_coverage(tmp_path: Path) -> None:
    labels = ("SAME", "SAME", "HIT", "REVIEWED", "MISSING")
    gbff = write(
        tmp_path / "input.gbff",
        b"".join(record_bytes(f"R{i}", locus) for i, locus in enumerate(labels)),
    )
    sha256 = sha256_bytes(gbff.read_bytes())
    ips = _context(sha256)
    context = _json(
        tmp_path / "context.json",
        {
            "schema": "enrich-bakta.context.v1",
            "metadata": {"output_sha256": sha256},
            "interproscan": ips,
            "eggnog": {"schema": "enrich-bakta.eggnog-context.v1", "metadata": {}},
        },
    )
    manifest_payload = _manifest(
        sha256, context_hash=sha256_bytes(context.read_bytes())
    )
    raw = parse_genbank_bytes(gbff.read_bytes())
    uid = feature_uid(raw.records[1].features[2])
    decision = {
        "candidate_id": "candidate:trusted",
        "source_id": "eggNOG",
        "source_sha256": "1" * 64,
        "candidate_role": "functional_proposal",
        "qualifier": "db_xref",
        "normalized_value": "GO:0000001",
        "final_status": "supported_existing",
    }
    manifest_payload["decisions"] = [
        {**decision, "target_feature_uids": [uid]},
        {
            **decision,
            "target_feature_uids": [uid],
            "final_status": "suppressed_conflict",
        },
        {
            **decision,
            "target_feature_uids": [uid],
            "qualifier": "translation",
            "normalized_value": "MK",
        },
    ]
    manifest = _json(tmp_path / "manifest.json", manifest_payload)

    rendered = export_cds_tsv(gbff, manifest_path=manifest, context_path=context)
    rows = _rows(rendered)
    assert [row["antifam_status"] for row in rows] == [
        "ambiguous_locus",
        "ambiguous_locus",
        "hit",
        "not_reported",
        "unavailable",
    ]
    assert json.loads(rows[2]["antifam_accessions"]) == ["ANF00002"]
    assert json.loads(rows[0]["accepted_decisions"]) == []
    assert json.loads(rows[1]["accepted_decisions"]) == [decision]
    assert json.loads(rows[1]["go"]) == []
    assert rows[1]["translation_origin"] == "not_recorded"
    assert b'"MK"' not in rendered


def test_bound_inputs_and_output_alias_fail_closed(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    gbff = write(tmp_path / "input.gbff", record_bytes("R", "CDS"))
    sha256 = sha256_bytes(gbff.read_bytes())
    context = _json(tmp_path / "context.json", _context(sha256))
    context_expected = sha256_bytes(context.read_bytes())
    manifest = _json(tmp_path / "manifest.json", _manifest(sha256))
    with pytest.raises(MergeError, match="manifest link"):
        export_cds_tsv(gbff, manifest_path=manifest, context_path=context)
    _json(manifest, _manifest("1" * 64, context_hash=context_expected))
    with pytest.raises(MergeError, match="manifest output SHA-256"):
        export_cds_tsv(gbff, manifest_path=manifest)
    _json(manifest, _manifest(sha256, context_hash=context_expected))
    context.write_bytes(context.read_bytes() + b" ")
    with pytest.raises(MergeError, match="manifest link"):
        export_cds_tsv(gbff, manifest_path=manifest, context_path=context)
    _json(context, _context("2" * 64))
    with pytest.raises(MergeError, match="context report output SHA-256"):
        export_cds_tsv(gbff, context_path=context)
    assert main(["--gbff", str(gbff), "--output", str(gbff)]) == 2
    assert capsys.readouterr().err.startswith("ERROR:")
    assert sha256_bytes(gbff.read_bytes()) == sha256


@pytest.mark.parametrize(
    ("sequence", "policy", "origin", "status"),
    [
        ("MK", "validated-only", "genomically_validated", "passed"),
        ("MM", "import-faa", "imported_faa", "failed"),
    ],
)
def test_exports_core_verified_translation_origin(
    tmp_path: Path,
    sequence: str,
    policy: str,
    origin: str,
    status: str,
) -> None:
    base = write(
        tmp_path / "base.gbff",
        record_bytes("TEST", "T_0001").replace(
            b'                     /translation="MK"\n',
            b'                     /pseudogene="unitary"\n',
        ),
    )
    faa = write(tmp_path / "input.faa", f">T_0001\n{sequence}\n".encode())
    eggnog = write(
        tmp_path / "annotations.tsv",
        eggnog_bytes("T_0001\tseed\t1e-4\t10\t-\t-\t-\t-\t-\t-\thhhhhhhhhhhhh"),
    )
    restored = tmp_path / "restored.gbff"
    manifest = tmp_path / "restored.json"
    restore(
        base, faa, eggnog, restored, manifest_path=manifest, translation_policy=policy
    )
    assert _rows(export_cds_tsv(base))[0]["translation_present"] == "false"
    rows = _rows(export_cds_tsv(restored, manifest_path=manifest))
    assert rows[0]["translation_present"] == "true"
    assert rows[0]["translation_origin"] == origin
    assert rows[0]["translation_genomic_validation_status"] == status


def test_manifest_snapshot_binds_translation_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gbff = write(tmp_path / "input.gbff", record_bytes("R", "L1"))
    manifest_payload = _manifest(sha256_bytes(gbff.read_bytes()))
    manifest_payload["entries"] = [{"entry_type": "translation_restoration"}]
    manifest = _json(tmp_path / "manifest.json", manifest_payload)

    class ReplacementSnapshot(dict):
        manifest_sha256 = "0" * 64

    monkeypatch.setattr(
        cds_export, "load_translation_evidence", lambda *_: ReplacementSnapshot()
    )
    with pytest.raises(MergeError, match="manifest changed"):
        export_cds_tsv(gbff, manifest_path=manifest)


def test_stdout_and_atomic_output_are_identical(
    tmp_path: Path, capsysbinary: pytest.CaptureFixture[bytes]
) -> None:
    gbff = write(tmp_path / "input.gbff", record_bytes("R", "L1"))
    output = tmp_path / "table.tsv"
    assert main(["--gbff", str(gbff)]) == 0
    rendered = capsysbinary.readouterr().out
    assert main(["--gbff", str(gbff), "--output", str(output)]) == 0
    assert rendered == output.read_bytes()
