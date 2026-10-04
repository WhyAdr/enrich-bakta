from __future__ import annotations

import json

import pytest
from test_merge_pipeline import eggnog_bytes, record_bytes, write

from enrich_bakta_lib.workflows.enrich import enrich


def run_conflict(tmp_path, policy):
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    fold = write(
        tmp_path / "fold.gbff",
        record_bytes("TEST", "T_0001", gene="xyz", ec_numbers=("1.2.3.4",)),
    )
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    egg = write(
        tmp_path / "egg.tsv",
        eggnog_bytes(
            "T_0001\tseed\t1e-4\t10\t-\tabc\t-\t1.2.3.4\tK00001\t-\thhhhhhhhhhhhh"
        ),
    )
    manifest = tmp_path / "out.json"
    enrich(
        bakta_path=base,
        faa_path=faa,
        eggnog_path=egg,
        baktfold_path=fold,
        output_path=tmp_path / "out.gbff",
        manifest_path=manifest,
        gene_conflict_policy=policy,
    )
    return json.loads(manifest.read_text())


@pytest.mark.parametrize("policy", ["skip", "prefer-eggnog", "prefer-baktfold"])
def test_gene_reconciliation_reason_overrides_planned_insertion(tmp_path, policy):
    payload = run_conflict(tmp_path, policy)
    rejected = [
        row
        for row in payload["decisions"]
        if row["qualifier"] == "gene" and row["final_status"] == "suppressed_conflict"
    ]
    assert len(rejected) == (2 if policy == "skip" else 1)
    expected_reason = (
        "gene_conflict_no_preference"
        if policy == "skip"
        else "gene_conflict_source_preference"
    )
    assert {row["reason_code"] for row in rejected} == {expected_reason}
    by_id = {row["candidate_id"]: row for row in rejected}
    projections = [
        row for row in payload["entries"] if row.get("candidate_id") in by_id
    ]
    assert len(projections) == len(rejected)
    assert all(row["reason_code"] == expected_reason for row in projections)
    assert all(row["planned_reason_code"] == "inserted" for row in projections)


def test_shared_support_source_projection_matches_final_decision(tmp_path):
    payload = run_conflict(tmp_path, "skip")
    shared = {
        row["candidate_id"]: row
        for row in payload["decisions"]
        if row["qualifier"] == "EC_number" and row["final_status"] == "shared_support"
    }
    assert len(shared) == 2
    projections = [
        row
        for row in payload["entries"]
        if row.get("entry_type") in {"baktfold_candidate", "eggnog_candidate"}
        and row.get("field") == "EC"
    ]
    assert len(projections) == 2
    for row in projections:
        decision = shared[row["candidate_id"]]
        assert row["status"] == row["final_status"] == decision["final_status"]
        assert row["reason_code"] == decision["reason_code"] == "shared_support"
        assert row["planned_reason_code"] == "inserted"


@pytest.mark.parametrize("case", ["authoritative", "protein-mismatch", "pair-conflict"])
def test_suppression_total_includes_each_rejection_kind(tmp_path, case):
    base_data = record_bytes(
        "TEST", "T_0001", gene="original" if case == "authoritative" else None
    )
    fold_data = record_bytes(
        "TEST", "T_0001", gene="xyz", ec_numbers=("1.2.3.4",), db_xrefs=("pdb:1ABC",)
    )
    if case == "protein-mismatch":
        fold_data = fold_data.replace(b'/translation="MK"', b'/translation="MM"')
    if case == "pair-conflict":
        fold_data = fold_data.replace(b'/gene="xyz"', b'/gene="other"', 1)
    base = write(tmp_path / "base.gbff", base_data)
    fold = write(tmp_path / "fold.gbff", fold_data)
    manifest = tmp_path / "out.json"
    enrich(
        bakta_path=base,
        baktfold_path=fold,
        output_path=tmp_path / "out.gbff",
        manifest_path=manifest,
    )
    payload = json.loads(manifest.read_text())
    rejected = [
        row
        for row in payload["decisions"]
        if row["final_status"].startswith("suppressed_")
    ]
    assert (
        len(rejected)
        == {
            "authoritative": 1,
            "protein-mismatch": 3,
            "pair-conflict": 2,
        }[case]
    )
    assert payload["metadata"]["suppressed_candidate_count"] == len(rejected)
