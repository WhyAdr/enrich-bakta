from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import pytest
from test_merge_pipeline import record_bytes

from enrich_bakta_lib.core.types import MergeError
from enrich_bakta_lib.workflows.diff import diff_gbff_files, diff_manifest, main


def _gbff_pair(tmp_path: Path) -> tuple[Path, Path]:
    base, enriched = tmp_path / "base.gbff", tmp_path / "enriched.gbff"
    base.write_bytes(record_bytes("TEST", "T_0001", notes=("same",)))
    enriched.write_bytes(
        record_bytes(
            "TEST",
            "T_0001",
            gene="abc",
            notes=("same", "same", "new note"),
            db_xrefs=("GO:0003677",),
            ec_numbers=("1.2.3.4",),
        )
        .replace(b"FEATURES", b"            Enrichment evidence\nFEATURES")
        .replace(
            b'                     /translation="MK"',
            b'                     /inference="similar to AA sequence:test"\n'
            b'                     /translation="MK"',
        )
    )
    return base, enriched


def test_gbff_qualifier_multisets_and_top_level_comment(tmp_path: Path) -> None:
    base, enriched = _gbff_pair(tmp_path)
    result = diff_gbff_files(base, enriched)
    assert result["changed_features"] == 2
    gene, cds = result["deltas"]
    assert gene["added"] == {"gene": ["abc"]}
    assert cds["locus_tag"] == "T_0001"
    assert cds["added"] == {
        "EC_number": ["1.2.3.4"],
        "db_xref": ["GO:0003677"],
        "gene": ["abc"],
        "inference": ["similar to AA sequence:test"],
        "note": ["new note", "same"],
    }
    assert gene["removed"] == cds["removed"] == {}
    assert result["comment_deltas"] == [
        {
            "record_id": "TEST",
            "record_index": 0,
            "added": ["Enrichment evidence"],
            "removed": [],
        }
    ]
    assert diff_gbff_files(enriched, base, qualifiers={"note"})["deltas"][0][
        "removed"
    ] == {"note": ["new note", "same"]}


def test_gbff_fail_on_diff_and_filter(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    base, enriched = _gbff_pair(tmp_path)
    args = ["--base", str(base), "--format", "json", "--fail-on-diff", "--enriched"]
    assert main([*args, str(base)]) == 0
    assert json.loads(capsys.readouterr().out)["changed_features"] == 0
    assert main([*args, str(enriched)]) == 1
    assert main([*args, str(enriched), "--locus", "MISSING"]) == 0
    assert main([*args, str(enriched), "--qualifier", "COMMENT"]) == 1


@pytest.mark.parametrize(
    ("before", "after"),
    [
        (b"atgaaataa", b"atgaaatag"),
        (
            b"FEATURES             Location/Qualifiers\n",
            b"FEATURES             Location/Qualifiers\n     misc_feature    1..9\n",
        ),
        (b"     CDS             1..9", b"     CDS             2..9"),
    ],
)
def test_gbff_rejects_structural_mismatch(
    tmp_path: Path, before: bytes, after: bytes
) -> None:
    base, enriched = _gbff_pair(tmp_path)
    enriched.write_bytes(enriched.read_bytes().replace(before, after))
    with pytest.raises(MergeError, match="structural parity check failed"):
        diff_gbff_files(base, enriched)


def _manifest(tmp_path: Path) -> Path:
    def decision(
        candidate: str,
        source: str,
        qualifier: str,
        status: str,
        uid: str,
        **extra: object,
    ) -> dict[str, object]:
        return {
            "candidate_id": candidate,
            "source_id": source,
            "qualifier": qualifier,
            "final_status": status,
            "normalized_value": "β value",
            "reason_code": "tested_reason",
            "target_feature_uids": [uid],
            **extra,
        }

    path = tmp_path / "manifest.json"
    path.write_text(
        json.dumps(
            {
                "schema": "enrich-bakta.merge-manifest.v2",
                "entries": [
                    {
                        "record": "R1",
                        "locus_tag": "A",
                        "insertion_id": "i",
                    },
                    {"record": "R2", "locus_tag": "B", "target_feature_uids": ["uid2"]},
                ],
                "decisions": [
                    decision("c3", "Baktfold", "gene", "suppressed_conflict", "uid2"),
                    decision(
                        "c1",
                        "KofamScan",
                        "db_xref",
                        "emitted",
                        "uid1",
                        insertion_ids=["i"],
                        supporting_candidate_ids=["c2"],
                    ),
                    decision(
                        "c2",
                        "eggNOG",
                        "db_xref",
                        "shared_support",
                        "uid1",
                        supporting_candidate_ids=["c1"],
                    ),
                    decision(
                        "c4",
                        "Baktfold",
                        "note",
                        "supported_existing",
                        "uid2",
                        confidence_diagnostics={
                            "support_witness": {"value": "already there"}
                        },
                    ),
                ],
            }
        ),
        encoding="utf-8",
    )
    return path


def test_manifest_groups_and_filters(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = _manifest(tmp_path)
    result = diff_manifest(path)
    first, second = result["groups"]
    assert (first["record_id"], first["locus_tag"], first["target_feature_uids"]) == (
        "R1",
        "A",
        ("uid1",),
    )
    assert first["decisions"][0]["supporting_candidate_ids"] == ("c2",)
    assert second["decisions"][0]["reason_code"] == "tested_reason"
    assert second["decisions"][1]["support_witness"] == {"value": "already there"}
    assert result["added_decisions"] == 2

    args = ["--manifest", str(path), "--fail-on-diff", "--format", "json"]
    assert main(args) == 1
    capsys.readouterr()
    filters = [
        "--locus",
        "A",
        "--qualifier",
        "db_xref",
        "--status",
        "shared_support",
        "--source",
    ]
    assert main([*args, *filters, "eggNOG"]) == 1
    result = json.loads(capsys.readouterr().out)
    assert result["added_decisions"] == 1
    assert len(result["groups"][0]["decisions"]) == 1
    assert main([*args, *filters, "Baktfold"]) == 0
    assert json.loads(capsys.readouterr().out)["groups"] == []
    assert main([*args, "--locus", "B"]) == 0
    assert json.loads(capsys.readouterr().out)["groups"][0]["locus_tag"] == "B"


@pytest.mark.parametrize("output_format", ["text", "json", "tsv"])
def test_manifest_rendering_is_deterministic(
    tmp_path: Path, output_format: str
) -> None:
    path = _manifest(tmp_path)
    args = ["--manifest", str(path), "--format", output_format, "--output"]
    first, second = tmp_path / "first", tmp_path / "second"
    assert main([*args, str(first)]) == main([*args, str(second)]) == 0
    assert first.read_bytes() == second.read_bytes()
    rendered = first.read_text(encoding="utf-8")
    assert all(
        value in rendered
        for value in ("supported_existing", "eggNOG", "already there", "β value")
    )
    if output_format == "tsv":
        rows = list(csv.DictReader(io.StringIO(rendered), delimiter="\t"))
        assert len(rows) == 4
        assert rows[0]["candidate_id"] == "c1"


def test_manifest_resolves_suppressed_locus_without_feature_uid(tmp_path: Path) -> None:
    path = _manifest(tmp_path)
    data = json.loads(path.read_bytes())
    data["decisions"][0]["target_feature_uids"] = []
    data["entries"].append({"candidate_id": "c3", "query_id": "C"})
    path.write_text(json.dumps(data), encoding="utf-8")

    result = diff_manifest(path, loci={"C"}, statuses={"suppressed_conflict"})
    assert result["added_decisions"] == 0
    assert result["groups"][0]["locus_tag"] == "C"
    assert result["groups"][0]["decisions"][0]["candidate_id"] == "c3"


def test_cli_rejects_incompatible_modes_and_output_alias(tmp_path: Path) -> None:
    path = _manifest(tmp_path)
    before = path.read_bytes()
    for options in (
        ["--base", str(path)],
        ["--enriched", str(path)],
        ["--output", str(path)],
    ):
        with pytest.raises(SystemExit, match="2"):
            main(["--manifest", str(path), *options])
    assert path.read_bytes() == before
