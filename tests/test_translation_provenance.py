from __future__ import annotations

import json
from pathlib import Path

import pytest
from test_merge_pipeline import eggnog_bytes, kofam_bytes, record_bytes, write

from enrich_bakta import enrich
from enrich_bakta_lib.core.merge_engine import (
    MergeError,
    load_translation_evidence,
    parse_genbank_bytes,
)
from enrich_bakta_lib.sources.kofam import merge as merge_kofam
from enrich_bakta_lib.workflows.restore_translations import restore
from graft_baktfold_additions import graft


def _translationless_base(tmp_path: Path) -> Path:
    return write(
        tmp_path / "base.gbff",
        record_bytes("TEST", "T_0001").replace(
            b'                     /translation="MK"\n',
            b'                     /pseudogene="unitary"\n',
        ),
    )


def _eggnog_input(tmp_path: Path) -> Path:
    return write(
        tmp_path / "annotations.tsv",
        eggnog_bytes("T_0001\tseed\t1e-4\t10\t-\t-\t-\t-\t-\t-\thhhhhhhhhhhhh"),
    )


def _kofam_input(tmp_path: Path) -> Path:
    return write(
        tmp_path / "kofam.txt",
        kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha protein"),
    )


def test_restoration_manifest_is_bound_to_exact_output_and_downstream(
    tmp_path: Path,
) -> None:
    base = _translationless_base(tmp_path)
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    eggnog = _eggnog_input(tmp_path)
    restored = tmp_path / "restored.gbff"
    manifest = tmp_path / "restoration.json"

    stats = restore(
        base,
        faa,
        eggnog,
        restored,
        manifest_path=manifest,
        translation_policy="validated-only",
    )

    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert stats["translation_evidence"]["T_0001"]["origin"] == "genomically_validated"
    assert payload["metadata"]["output_sha256"] == stats["output_sha256"]
    assert (
        payload["entries"][0]["translation_evidence"]["bound_output_sha256"]
        == stats["output_sha256"]
    )
    evidence = load_translation_evidence(
        manifest, parse_genbank_bytes(restored.read_bytes(), "restored")
    )
    assert evidence["T_0001"].origin == "genomically_validated"

    next_output = tmp_path / "kofam.gbff"
    merge_kofam(
        restored,
        faa,
        _kofam_input(tmp_path),
        next_output,
        translation_evidence_manifest=manifest,
    )
    assert b"KEGG:K00001" in next_output.read_bytes()


def test_imported_translation_requires_explicit_downstream_opt_in(
    tmp_path: Path,
) -> None:
    base = _translationless_base(tmp_path)
    faa = write(tmp_path / "base.faa", b">T_0001\nMM\n")
    eggnog = _eggnog_input(tmp_path)
    restored = tmp_path / "restored.gbff"
    manifest = tmp_path / "restoration.json"
    restore(
        base,
        faa,
        eggnog,
        restored,
        manifest_path=manifest,
        translation_policy="import-faa",
    )

    with pytest.raises(MergeError, match="allow-imported-translations"):
        merge_kofam(
            restored,
            faa,
            _kofam_input(tmp_path),
            tmp_path / "blocked.gbff",
            translation_evidence_manifest=manifest,
        )
    merge_kofam(
        restored,
        faa,
        _kofam_input(tmp_path),
        tmp_path / "allowed.gbff",
        translation_evidence_manifest=manifest,
        allow_imported_translations=True,
    )


def test_restored_marker_requires_manifest_for_functional_merge(tmp_path: Path) -> None:
    base = _translationless_base(tmp_path)
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    eggnog = _eggnog_input(tmp_path)
    restored = tmp_path / "restored.gbff"
    restore(base, faa, eggnog, restored)

    with pytest.raises(MergeError, match="translation-evidence-manifest"):
        merge_kofam(
            restored,
            faa,
            _kofam_input(tmp_path),
            tmp_path / "blocked.gbff",
        )


def test_baktfold_protein_mismatch_suppresses_protein_derived_transfers(
    tmp_path: Path,
) -> None:
    base = record_bytes("TEST", "T_0001")
    base_path = write(tmp_path / "parity-base.gbff", base)
    source_path = write(
        tmp_path / "baktfold.gbff",
        record_bytes("TEST", "T_0001", gene="abc", db_xrefs=("pdb:1ABC",)).replace(
            b'/translation="MK"', b'/translation="MM"'
        ),
    )
    output = tmp_path / "out.gbff"

    stats = graft(base_path, source_path, output, add_comment_note=False)

    assert stats["parity"]["translation_mismatches"] == ["T_0001"]
    assert stats["translation_mismatch_suppressed_features"] == 2
    assert stats["gene_added"] == 0
    assert stats["xref_added"] == 0
    assert output.read_bytes() == base


def test_conflicting_paired_names_are_suppressed_as_whole_candidates(
    tmp_path: Path,
) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    baktfold = write(
        tmp_path / "baktfold.gbff", record_bytes("TEST", "T_0001", gene="abc")
    )
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    eggnog = write(
        tmp_path / "annotations.tsv",
        eggnog_bytes("T_0001\tseed\t1e-4\t10\t-\txyz\t-\t-\t-\t-\thhhhhhhhhhhhh"),
    )
    output = tmp_path / "out.gbff"
    manifest = tmp_path / "manifest.json"

    enrich(
        bakta_path=base,
        output_path=output,
        faa_path=faa,
        baktfold_path=baktfold,
        eggnog_path=eggnog,
        eggnog_version="3.0.0-beta6",
        manifest_path=manifest,
        add_comment_note=False,
    )

    payload = json.loads(manifest.read_text(encoding="utf-8"))
    decisions = [row for row in payload["decisions"] if row["qualifier"] == "gene"]
    assert decisions
    assert all(row["final_status"] == "suppressed_conflict" for row in decisions)
    assert all(not row["insertion_ids"] for row in decisions)
    assert b'/gene="abc"' not in output.read_bytes()
    assert b'/gene="xyz"' not in output.read_bytes()
