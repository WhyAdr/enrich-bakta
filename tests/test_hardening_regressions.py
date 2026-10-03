from __future__ import annotations

import io

import pytest
from Bio import SeqIO
from test_merge_pipeline import eggnog_bytes, record_bytes, write

import merge_eggnog_bakta as egg
import merge_engine as core
from enrich_bakta_lib.core.decisions import validate_candidate_ledger

ROW = "T_0001\tseed.1\t1e-20\t100\t-\tabc\t-\t-\t-\t-\thhhhhhhhhhhhh"


def _candidate_fixture():
    data = record_bytes("TEST", "T_0001")
    base = core.parse_genbank_bytes(data)
    feature = next(item for item in base.features if item.feature_type == "CDS")
    insertion = core.qualifier_insertion(
        data, feature, "note", "evidence", "fixture", "evidence", 0
    )
    row = {
        "candidate_id": "candidate:fixture",
        "source_id": "fixture",
        "source_sha256": "a" * 64,
        "target_feature_uids": [core.feature_uid(feature)],
        "field": "note",
        "qualifier": "note",
        "raw_value": "evidence",
        "normalized_value": "evidence",
        "planned_status": "planned",
        "final_status": "emitted",
        "insertion_ids": [core.insertion_uid(insertion)],
        "supporting_candidate_ids": [],
        "confidence_diagnostics": {},
        "reason": "fixture",
    }
    return base, insertion, row


def test_candidate_ledger_rejects_wrong_target_and_value() -> None:
    base, insertion, row = _candidate_fixture()
    wrong_value = dict(row, normalized_value="other")
    with pytest.raises(core.MergeError, match="wrong value"):
        validate_candidate_ledger([wrong_value], [insertion], base=base)

    other_feature = next(
        item for item in base.features if item.locus_tag and item.feature_type != "CDS"
    )
    wrong_target = dict(row, target_feature_uids=[core.feature_uid(other_feature)])
    with pytest.raises(core.MergeError, match="wrong feature"):
        validate_candidate_ledger([wrong_target], [insertion], base=base)


def test_candidate_ledger_rejects_orphaned_final_insertion() -> None:
    base, insertion, _ = _candidate_fixture()
    with pytest.raises(core.MergeError, match="unsupported"):
        validate_candidate_ledger([], [insertion], base=base)


def test_candidate_ledger_requires_functional_support_for_provenance() -> None:
    base, insertion, row = _candidate_fixture()
    provenance = core.qualifier_insertion(
        base.data,
        base.features[2],
        "inference",
        "fixture provenance",
        "fixture provenance",
        "fixture provenance",
        1,
    )
    provenance_row = dict(
        row,
        candidate_id="candidate:provenance",
        qualifier="inference",
        normalized_value="fixture provenance",
        raw_value="fixture provenance",
        insertion_ids=[core.insertion_uid(provenance)],
        confidence_diagnostics={"provenance": True},
    )
    with pytest.raises(core.MergeError, match="functional support"):
        validate_candidate_ledger([provenance_row], [insertion, provenance], base=base)


def test_candidate_ledger_rejects_unknown_support_reference() -> None:
    base, insertion, row = _candidate_fixture()
    row["supporting_candidate_ids"] = ["candidate:missing"]
    with pytest.raises(core.MergeError, match="unknown support"):
        validate_candidate_ledger([row], [insertion], base=base)


def compound(last):
    return record_bytes("TEST", "T_0001").replace(
        b"     CDS             1..9",
        b"     CDS             join(1..3,\n                     " + last + b")",
    )


def test_A01_changed_location_continuation_rejected():
    left, right = compound(b"7..9"), compound(b"4..6")
    core.validate_genbank_semantics(left)
    core.validate_genbank_semantics(right)
    with pytest.raises(core.MergeError):
        core.strict_parity_check(
            core.parse_genbank_bytes(left), core.parse_genbank_bytes(right)
        )


def test_A01_equivalent_location_wrapping_accepted():
    left = compound(b"7..9")
    right = left.replace(b"join(1..3,\n                     7..9)", b"join(1..3,7..9)")
    core.strict_parity_check(
        core.parse_genbank_bytes(left), core.parse_genbank_bytes(right)
    )


def test_A02_free_text_uses_biopython_semantics():
    data = record_bytes("TEST", "T_0001").replace(
        b'/product="test protein"', b'/product="test\n                     protein"'
    )
    bio = list(SeqIO.parse(io.StringIO(data.decode()), "genbank"))[0].features[2]
    raw = core.parse_genbank_bytes(data).features[2]
    assert raw.values("product") == bio.qualifiers["product"]


def test_A09_duplicate_identical_locus_rejected():
    data = record_bytes("TEST", "T_0001").replace(
        b'/protein_id="gnl|Bakta|test"',
        b'/locus_tag="T_0001"\n                     /protein_id="gnl|Bakta|test"',
    )
    with pytest.raises(core.MergeError):
        core.validate_faa_gbff(
            core.parse_genbank_bytes(data),
            {"T_0001": "MK"},
            ["T_0001"],
            source_name="test",
        )


def workbook_bytes(name):
    openpyxl = pytest.importorskip("openpyxl")
    workbook = openpyxl.Workbook()
    worksheet = workbook.active
    worksheet.title = "annotations"
    lines = eggnog_bytes(ROW.replace("\tabc\t", f"\t{name}\t")).decode().splitlines()
    worksheet.append(lines[1].split("\t"))
    worksheet.append(lines[2].split("\t"))
    stream = io.BytesIO()
    workbook.save(stream)
    workbook.close()
    return stream.getvalue()


def test_A10_xlsx_uses_captured_snapshot(tmp_path):
    path = write(tmp_path / "table.xlsx", workbook_bytes("changed"))
    captured = workbook_bytes("original")
    table = egg.parse_eggnog_path(path, expected_version="3.0.0-beta6", data=captured)
    assert table.hits[0].preferred_name == "original"


def test_A11_manifest_directory_failure_preserves_output(tmp_path):
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    output = write(tmp_path / "out.gbff", b"previous output")
    manifest = tmp_path / "manifest.json"
    manifest.mkdir()
    with pytest.raises((core.MergeError, OSError)):
        core.finalize_merge(
            base_path=base,
            base_data=base.read_bytes(),
            output_path=output,
            other_inputs=[],
            insertions=[],
            manifest_path=manifest,
            metadata={},
        )
    assert output.read_bytes() == b"previous output"
