from __future__ import annotations

import json
from pathlib import Path

import pytest

from tools.followthrough_acceptance import (
    assert_candidate_manifest,
    assert_expected,
    load_baseline,
)


def test_checked_in_scientific_baseline_matches_runtime_version():
    baseline = load_baseline(Path("docs/data/followthrough-scientific-baseline.json"))
    assert set(baseline["samples"]) == {"C14", "SM"}


def test_expected_value_drift_fails_closed():
    with pytest.raises(AssertionError, match="drift"):
        assert_expected("output hash", "changed", "accepted")


def test_candidate_manifest_semantics_accept_consistent_projection(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "metadata": {"suppressed_candidate_count": 1},
                "decisions": [
                    {
                        "candidate_id": "candidate:1",
                        "final_status": "suppressed_conflict",
                        "reason_code": "gene_conflict_no_preference",
                    }
                ],
                "entries": [
                    {
                        "candidate_id": "candidate:1",
                        "status": "suppressed_conflict",
                        "final_status": "suppressed_conflict",
                        "reason_code": "gene_conflict_no_preference",
                        "planned_reason_code": "inserted",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    assert assert_candidate_manifest(manifest) == {
        "candidate_decisions": 1,
        "source_projections": 1,
        "suppressed_candidates": 1,
    }


@pytest.mark.parametrize(
    "defect", ["projection", "missing-candidate-id", "suppression-total"]
)
def test_candidate_manifest_semantics_reject_residual_defects(tmp_path, defect):
    manifest = tmp_path / "manifest.json"
    payload = {
        "metadata": {"suppressed_candidate_count": 1},
        "decisions": [
            {
                "candidate_id": "candidate:1",
                "final_status": "suppressed_conflict",
                "reason_code": "gene_conflict_no_preference",
            }
        ],
        "entries": [
            {
                "candidate_id": "candidate:1",
                "status": "suppressed_conflict",
                "final_status": "suppressed_conflict",
                "reason_code": "gene_conflict_no_preference",
                "planned_reason_code": "inserted",
            }
        ],
    }
    if defect == "projection":
        payload["entries"][0]["reason_code"] = "inserted"
    elif defect == "missing-candidate-id":
        del payload["entries"][0]["candidate_id"]
    else:
        payload["metadata"]["suppressed_candidate_count"] = 0
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(AssertionError):
        assert_candidate_manifest(manifest)
