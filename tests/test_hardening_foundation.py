from __future__ import annotations

from pathlib import Path

import pytest
from test_merge_pipeline import record_bytes, write

import merge_engine as core


@pytest.mark.parametrize(
    "location",
    [
        b"join(1..3,7..9)",
        b"complement(join(1..3,7..9))",
        b"join(7..9,1..3)",
        b"<1..>9",
        b"OTHER.1:1..9",
    ],
)
def test_complete_location_roundtrip(location: bytes) -> None:
    base = record_bytes("TEST", "T_0001").replace(
        b"     CDS             1..9", b"     CDS             " + location
    )
    wrapped = base.replace(b"1..3,7..9", b"1..3,\n                     7..9")
    core.strict_parity_check(
        core.parse_genbank_bytes(base), core.parse_genbank_bytes(wrapped)
    )


def test_semantic_roundtrip_rejects_incorrect_payload(tmp_path: Path) -> None:
    data = record_bytes("TEST", "T_0001")
    path = write(tmp_path / "base.gbff", data)
    feature = core.parse_genbank_bytes(data).features[2]
    insertion = core.qualifier_insertion(
        data, feature, "note", "expected", "test", "test", 0
    )
    insertion.payload = insertion.payload.replace(b"expected", b"different")
    with pytest.raises(core.MergeError, match="semantic qualifier round-trip"):
        core.finalize_merge(
            base_path=path,
            base_data=data,
            output_path=tmp_path / "out",
            other_inputs=[],
            insertions=[insertion],
            manifest_path=None,
            metadata={},
        )
    assert not (tmp_path / "out").exists()


def test_all_artifacts_staged_before_replacement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = write(tmp_path / "first", b"old first")
    second = write(tmp_path / "second", b"old second")
    original = core.tempfile.mkstemp
    calls = 0

    def fail_second(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected staging failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(core.tempfile, "mkstemp", fail_second)
    with pytest.raises(core.MergeError, match="injected staging failure"):
        core.stage_artifacts([(first, b"new first"), (second, b"new second")], [])
    assert first.read_bytes() == b"old first"
    assert second.read_bytes() == b"old second"
    assert not list(tmp_path.glob("*.tmp"))


def test_second_promotion_failure_has_documented_boundary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = write(tmp_path / "first", b"old first")
    second = write(tmp_path / "second", b"old second")
    original = core.os.replace
    calls = 0

    def fail_second(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected promotion failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(core.os, "replace", fail_second)
    with pytest.raises(core.MergeError, match="not a cross-file transaction"):
        core.stage_artifacts([(first, b"new first"), (second, b"new second")], [])
    assert first.read_bytes() == b"new first"
    assert second.read_bytes() == b"old second"
    assert not list(tmp_path.glob("*.tmp"))


def test_single_merge_parses_genbank_exactly_twice(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from test_merge_pipeline import eggnog_bytes

    from enrich_bakta_lib.core import genbank as core_genbank
    from enrich_bakta_lib.sources.eggnog import merge as merge_eggnog

    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    egg = write(
        tmp_path / "egg.tsv",
        eggnog_bytes(
            "T_0001\tseed\t1e-4\t10\tC\tabc\t-\t1.2.3.4\tK00001\t-\thhhhhhhhhhhhh"
        ),
    )
    genbank_calls = 0
    original_parse = core_genbank.SeqIO.parse

    def counting_parse(handle, fmt, *args, **kwargs):
        nonlocal genbank_calls
        if fmt == "genbank":
            genbank_calls += 1
        return original_parse(handle, fmt, *args, **kwargs)

    monkeypatch.setattr(core_genbank.SeqIO, "parse", counting_parse)
    merge_eggnog(
        base,
        faa,
        egg,
        tmp_path / "out.gbff",
        manifest_path=tmp_path / "manifest.json",
    )
    assert genbank_calls == 2
