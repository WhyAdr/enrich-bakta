from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
from test_merge_pipeline import eggnog_bytes, record_bytes, write

from enrich_bakta_lib.core import artifacts
from enrich_bakta_lib.core.types import MergeError, sha256_bytes
from enrich_bakta_lib.workflows.enrich import main


def test_middle_promotion_failure_records_each_artifact(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = write(tmp_path / "first", b"old first")
    second = write(tmp_path / "second", b"old second")
    third = tmp_path / "third"
    journal = tmp_path / "journal.json"
    paths = (first, second, third)
    original_replace = artifacts.os.replace

    def fail_second(src: Path, dst: Path) -> None:
        if dst == first:
            pending = json.loads(journal.read_bytes())
            assert [row["status"] for row in pending["artifacts"]] == ["staged"] * 3
        if dst == second:
            raise OSError("injected promotion failure")
        original_replace(src, dst)

    monkeypatch.setattr(artifacts.os, "replace", fail_second)
    with pytest.raises(MergeError, match="injected promotion failure") as caught:
        artifacts.stage_artifacts(
            [(p, b"new " + p.name.encode()) for p in paths],
            [],
            promotion_journal_path=journal,
        )
    assert str(journal) in str(caught.value)
    payload = json.loads(journal.read_bytes())
    assert payload["schema"] == "enrich-bakta.promotion-journal.v1"
    assert payload["artifacts"] == [
        {
            "destination": str(path),
            "expected_sha256": sha256_bytes(b"new " + path.name.encode()),
            "previous_sha256": sha256_bytes(b"old " + path.name.encode())
            if path != third
            else None,
            "status": status,
        }
        for path, status in zip(paths, ("promoted", "failed", "staged"), strict=True)
    ]
    assert first.read_bytes() == b"new first"
    assert second.read_bytes() == b"old second"
    assert not third.exists()
    assert not list(tmp_path.glob("*.tmp"))


def test_journal_write_failure_prevents_any_promotion(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = write(tmp_path / "target", b"old")
    journal = tmp_path / "journal.json"
    original_replace = artifacts.os.replace

    def fail_journal(src: Path, dst: Path) -> None:
        if dst == journal:
            raise OSError("injected journal write failure")
        original_replace(src, dst)

    monkeypatch.setattr(artifacts.os, "replace", fail_journal)
    with pytest.raises(MergeError, match="injected journal write failure"):
        artifacts.stage_artifacts(
            [(target, b"new")], [], promotion_journal_path=journal
        )
    assert target.read_bytes() == b"old"
    assert not journal.exists()
    assert not list(tmp_path.glob("*.tmp"))


def test_journal_cli_records_hashes_and_rejects_collision(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    eggnog = write(
        tmp_path / "egg.tsv",
        eggnog_bytes(
            "T_0001\tseed\t1e-4\t10\tC\tabc\t-\t1.2.3.4\tK00001\t-\thhhhhhhhhhhhh"
        ),
    )
    output, manifest = tmp_path / "out.gbff", tmp_path / "out.json"
    journal = tmp_path / "journal.json"
    report = tmp_path / "report.md"
    command = [
        "enrich-bakta",
        "--bakta",
        str(base),
        "--faa",
        str(faa),
        "--eggnog",
        str(eggnog),
        "--output",
        str(output),
        "--manifest",
        str(manifest),
        "--summary-report",
        str(report),
        "--promotion-journal",
    ]
    monkeypatch.setattr(sys, "argv", [*command, str(eggnog)])
    assert main() == 2
    assert not output.exists()
    assert not manifest.exists()
    monkeypatch.setattr(sys, "argv", [*command, str(report)])
    assert main() == 2
    assert not report.exists()
    assert not output.exists()
    assert not manifest.exists()
    monkeypatch.setattr(sys, "argv", [*command, str(journal)])
    assert main() == 0
    assert json.loads(journal.read_bytes())["artifacts"] == [
        {
            "destination": str(path),
            "expected_sha256": sha256_bytes(path.read_bytes()),
            "previous_sha256": None,
            "status": "promoted",
        }
        for path in (output, manifest, report)
    ]


def test_journal_update_failure_leaves_the_write_ahead_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = write(tmp_path / "target", b"old")
    journal = tmp_path / "journal.json"
    original_replace = artifacts.os.replace

    def fail_journal_update(src: Path, dst: Path) -> None:
        if dst == journal and journal.exists():
            raise OSError("injected journal update failure")
        original_replace(src, dst)

    monkeypatch.setattr(artifacts.os, "replace", fail_journal_update)
    with pytest.raises(MergeError, match="injected journal update failure"):
        artifacts.stage_artifacts(
            [(target, b"new")], [], promotion_journal_path=journal
        )
    (row,) = json.loads(journal.read_bytes())["artifacts"]
    assert row["status"] == "staged"
    assert row["previous_sha256"] == sha256_bytes(b"old")
    assert row["expected_sha256"] == sha256_bytes(target.read_bytes())
    assert not list(tmp_path.glob("*.tmp"))
