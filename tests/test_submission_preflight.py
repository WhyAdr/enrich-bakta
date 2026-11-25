from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import pytest
from test_merge_pipeline import record_bytes, write

from enrich_bakta_lib.workflows.submission_preflight import (
    SCHEMA,
    check_submission_gbff,
    main,
)


def test_clean_stdout_and_vocab_snapshot(
    tmp_path: Path, capsysbinary: pytest.CaptureFixture[bytes]
) -> None:
    data = record_bytes(
        "TEST",
        "T_0001",
        db_xrefs=(
            "PDB:12GS",
            "GO:0003677",
            "InterPro:IPR002928",
            "UniProtKB/Swiss-Prot:P12345",
            "ATCC(in host):123456",
            "HGNC:HGNC:2041",
            "BOLD:EPAF263",
        ),
    )
    path = write(tmp_path / "clean.gbff", data)
    assert main(["--gbff", str(path)]) == 0
    captured = capsysbinary.readouterr()
    assert captured.err == b""
    report = json.loads(captured.out)
    assert report["schema"] == SCHEMA
    assert report["input_sha256"] == hashlib.sha256(data).hexdigest()
    assert report["db_xref_count"] == 7
    assert report["finding_count"] == 0
    assert report["vocabulary"]["id"] == "insdc-db-xref-2026-11-25"
    assert report["vocabulary"]["prefix_count"] == 125
    assert report["findings"] == []
    assert path.read_bytes() == data


def test_findings_cover_record_local_features_and_repeated_qualifiers(
    tmp_path: Path,
) -> None:
    data = record_bytes(
        "ONE",
        "T_0001",
        db_xrefs=("pdb:1ABC", "pdb:1ABC", "PDB:", "missing-separator", ":id", "GO: id"),
    ) + record_bytes("TWO", "T_0001", db_xrefs=("Unregistered:X:9",)).replace(
        b'/mol_type="genomic DNA"\n',
        b'/mol_type="genomic DNA"\n                     /db_xref="UnknownSource:1"\n',
    )
    path = write(tmp_path / "sample.gbff", data)
    report = check_submission_gbff(path)
    assert report["record_count"] == 2
    assert report["feature_count"] == 6
    assert report["db_xref_count"] == 8
    assert report["findings_by_code"] == {
        "database_prefix_case": 2,
        "malformed_db_xref": 4,
        "unknown_database_prefix": 2,
    }
    first, duplicate, *rest = report["findings"]
    assert first["record_id"] == "ONE"
    assert first["record_index"] == 0
    assert first["locus_tag"] == "T_0001"
    assert first["feature_type"] == "CDS"
    assert first["feature_location"] == "1..9"
    assert first["feature_uid"] == duplicate["feature_uid"]
    assert first["qualifier"] == "db_xref"
    assert first["qualifier_index"] == 0
    assert duplicate["qualifier_index"] == 1
    assert first["suggested_prefix"] == "PDB"
    assert first["database_prefix"] == "pdb"
    assert first["value"] == duplicate["value"] == "pdb:1ABC"
    assert data[first["feature_byte_offset"] :].startswith(b"     CDS")
    source, last = rest[-2:]
    assert source["record_id"] == last["record_id"] == "TWO"
    assert source["record_index"] == last["record_index"] == 1
    assert source["feature_type"] == "source"
    assert source["locus_tag"] == ""
    assert source["value"] == "UnknownSource:1"
    assert last["value"] == "Unregistered:X:9"
    assert last["database_prefix"] == "Unregistered"
    assert last["feature_uid"] != first["feature_uid"]


def test_tsv_is_deterministic_and_is_written_before_findings_exit(
    tmp_path: Path,
) -> None:
    genome = write(
        tmp_path / "sample.gbff",
        record_bytes("TEST", "T_0001", db_xrefs=("pdb:12GS", "KEGG:K00001")),
    )
    output = tmp_path / "results.tsv"
    flags = ["--gbff", str(genome), "--format", "tsv", "--output", str(output)]
    assert main(flags) == 1
    original = output.read_bytes()
    assert main(flags) == 1
    assert output.read_bytes() == original
    rows = list(csv.DictReader(io.StringIO(original.decode("utf-8")), delimiter="\t"))
    assert [row["code"] for row in rows] == [
        "database_prefix_case",
        "unknown_database_prefix",
    ]
    assert rows[0]["suggested_prefix"] == "PDB"
    assert rows[0]["input_sha256"] == hashlib.sha256(genome.read_bytes()).hexdigest()
    assert rows[0]["vocabulary_id"] == "insdc-db-xref-2026-11-25"


def test_input_error_and_output_collision_preserve_input(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    data = record_bytes("TEST", "T_0001")
    genome = write(tmp_path / "sample.gbff", data)
    assert main(["--gbff", str(genome), "--output", str(genome)]) == 2
    assert "ERROR:" in capsys.readouterr().err
    assert genome.read_bytes() == data
    assert main(["--gbff", str(tmp_path / "missing.gbff")]) == 2
    assert "cannot read submission GBFF" in capsys.readouterr().err
