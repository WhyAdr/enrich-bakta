from __future__ import annotations

import json

import pytest
from test_merge_pipeline import eggnog_bytes, kofam_bytes, record_bytes, write

from enrich_bakta import enrich
from enrich_bakta_lib.core.decisions import validate_candidate_ledger
from enrich_bakta_lib.core.merge_engine import (
    MergeError,
    load_translation_evidence,
    parse_genbank_bytes,
)
from enrich_bakta_lib.sources.eggnog import merge as merge_egg
from enrich_bakta_lib.sources.eggnog import parse_eggnog_tsv
from enrich_bakta_lib.sources.kofam import merge as merge_ko
from enrich_bakta_lib.sources.value_rules import validate_preferred_name
from enrich_bakta_lib.workflows.restore_translations import restore

ROW = "T_0001\tseed\t1e-4\t10\t-\tabc\t-\t-\tK00001\t-\thhhhhhhhhhhhh"


def inputs(
    tmp_path, *, sequence=b"atgaaataa", protein=b"MK", extra=b"", location=b"1..9"
):
    data = record_bytes("TEST", "T_0001").replace(
        b'/translation="MK"', b"/pseudo" + extra
    )
    data = data.replace(b"atgaaataa", sequence).replace(
        b"CDS             1..9", b"CDS             " + location
    )
    return (
        write(tmp_path / "base.gbff", data),
        write(tmp_path / "base.faa", b">T_0001\n" + protein + b"\n"),
        write(tmp_path / "egg.tsv", eggnog_bytes(ROW)),
    )


def restored_inputs(tmp_path):
    base, faa, egg = inputs(tmp_path, protein=b"MM")
    out, manifest = tmp_path / "restored.gbff", tmp_path / "restored.json"
    restore(
        base, faa, egg, out, manifest_path=manifest, translation_policy="import-faa"
    )
    return out, faa, egg, manifest


def test_complete_alternative_start_is_validated(tmp_path):
    base, faa, egg = inputs(tmp_path, sequence=b"gtgaaataa")
    restore(base, faa, egg, tmp_path / "out.gbff", translation_policy="validated-only")


@pytest.mark.parametrize(
    "sequence,protein", [(b"atgaaaaaa", b"MKK"), (b"aaaaaataa", b"KK")]
)
def test_invalid_complete_cds_is_not_validated(tmp_path, sequence, protein):
    base, faa, egg = inputs(tmp_path, sequence=sequence, protein=protein)
    with pytest.raises(MergeError):
        restore(
            base, faa, egg, tmp_path / "out.gbff", translation_policy="validated-only"
        )


@pytest.mark.parametrize(
    "extra,location",
    [
        (b"", b"<1..9"),
        (b'\n                     /exception="ribosomal slippage"', b"1..9"),
        (b"", b"1..12"),
    ],
)
def test_unsupported_cds_model_is_not_validated(tmp_path, extra, location):
    base, faa, egg = inputs(tmp_path, extra=extra, location=location)
    with pytest.raises(MergeError):
        restore(
            base, faa, egg, tmp_path / "out.gbff", translation_policy="validated-only"
        )


def test_bad_translation_table_is_contextual_error(tmp_path):
    base, faa, egg = inputs(tmp_path)
    base.write_bytes(
        base.read_bytes().replace(b"/transl_table=11", b"/transl_table=999")
    )
    with pytest.raises(MergeError):
        restore(
            base, faa, egg, tmp_path / "out.gbff", translation_policy="validated-only"
        )


@pytest.mark.parametrize("adapter", ["kofam", "eggnog", "unified"])
def test_restore_merge_ledger_chain_survives(tmp_path, adapter):
    base, faa, egg, parent = restored_inputs(tmp_path)
    out, manifest = tmp_path / "merged.gbff", tmp_path / "merged.json"
    ko = write(tmp_path / "ko.txt", kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha"))
    common = dict(
        manifest_path=manifest,
        translation_evidence_manifest=parent,
        allow_imported_translations=True,
    )
    if adapter == "kofam":
        merge_ko(base, faa, ko, out, **common)
    elif adapter == "eggnog":
        merge_egg(base, faa, egg, out, **common)
    else:
        enrich(
            bakta_path=base, faa_path=faa, eggnog_path=egg, output_path=out, **common
        )
    ledger = load_translation_evidence(manifest, parse_genbank_bytes(out.read_bytes()))
    assert ledger["T_0001"].origin == "imported_faa"


def test_restoration_rerun_requires_parent_ledger(tmp_path):
    base, faa, egg, _ = restored_inputs(tmp_path)
    with pytest.raises(MergeError, match="translation-evidence-manifest"):
        restore(
            base,
            faa,
            egg,
            tmp_path / "rerun.gbff",
            manifest_path=tmp_path / "rerun.json",
            translation_policy="validated-only",
        )


@pytest.mark.parametrize(
    "mutation", ["unknown_origin", "inconsistent_policy", "string_boolean"]
)
def test_invalid_ledger_fields_rejected(tmp_path, mutation):
    base, faa, egg, manifest = restored_inputs(tmp_path)
    payload = json.loads(manifest.read_text())
    evidence = next(
        row["translation_evidence"]
        for row in payload["entries"]
        if row.get("entry_type") == "translation_restoration"
    )
    if mutation == "unknown_origin":
        evidence["origin"] = "imported_faa_typo"
    elif mutation == "inconsistent_policy":
        evidence["origin"] = "genomically_validated"
    else:
        evidence["artifact_sequence_match"] = "false"
    manifest.write_text(json.dumps(payload))
    with pytest.raises(MergeError):
        load_translation_evidence(manifest, parse_genbank_bytes(base.read_bytes()))


def test_kofam_existing_ko_retains_hit_evidence(tmp_path):
    base = write(
        tmp_path / "base.gbff",
        record_bytes("TEST", "T_0001", db_xrefs=("KEGG:K00001",)),
    )
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    ko = write(tmp_path / "ko.txt", kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha"))
    out, manifest = tmp_path / "out.gbff", tmp_path / "out.json"
    merge_ko(base, faa, ko, out, manifest_path=manifest, add_comment_note=False)
    assert b"KofamScan:K00001;threshold=1;score=2;E-value=3e-4" in out.read_bytes()
    payload = json.loads(manifest.read_text())
    ko_decision = next(row for row in payload["decisions"] if row["field"] == "KO")
    assert ko_decision["final_status"] == "supported_existing"


def test_kofam_rerun_decision_is_supported_existing(tmp_path):
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    ko = write(tmp_path / "ko.txt", kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha"))
    first, second, manifest = (
        tmp_path / "first.gbff",
        tmp_path / "second.gbff",
        tmp_path / "second.json",
    )
    merge_ko(base, faa, ko, first)
    merge_ko(first, faa, ko, second, manifest_path=manifest)
    assert first.read_bytes() == second.read_bytes()
    row = next(
        row
        for row in json.loads(manifest.read_text())["decisions"]
        if row["field"] == "KO"
    )
    assert row["final_status"] == "supported_existing"


def test_gene_class_provenance_removed_when_gene_conflicts(tmp_path):
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    fold = write(
        tmp_path / "fold.gbff",
        record_bytes("TEST", "T_0001", gene="xyz", ec_numbers=("1.2.3.4",)),
    )
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    egg = write(tmp_path / "egg.tsv", eggnog_bytes(ROW))
    out = tmp_path / "out.gbff"
    enrich(
        bakta_path=base,
        faa_path=faa,
        eggnog_path=egg,
        baktfold_path=fold,
        output_path=out,
        add_comment_note=False,
    )
    assert b"Baktfold gene-symbol evidence" not in out.read_bytes()
    assert b'/EC_number="1.2.3.4"' in out.read_bytes()


def test_paired_gene_is_one_candidate_with_two_targets(tmp_path):
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    egg = write(tmp_path / "egg.tsv", eggnog_bytes(ROW))
    manifest = tmp_path / "out.json"
    merge_egg(base, faa, egg, tmp_path / "out.gbff", manifest_path=manifest)
    genes = [
        row
        for row in json.loads(manifest.read_text())["decisions"]
        if row["qualifier"] == "gene"
    ]
    assert len(genes) == 1
    assert len(genes[0]["target_feature_uids"]) == 2
    assert len(genes[0]["insertion_ids"]) == 2


def test_two_sources_report_shared_support(tmp_path):
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    fold = write(
        tmp_path / "fold.gbff", record_bytes("TEST", "T_0001", ec_numbers=("1.2.3.4",))
    )
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    egg = write(
        tmp_path / "egg.tsv",
        eggnog_bytes(ROW.replace("\t-\t-\tK00001\t", "\t-\t1.2.3.4\tK00001\t")),
    )
    manifest = tmp_path / "out.json"
    enrich(
        bakta_path=base,
        faa_path=faa,
        eggnog_path=egg,
        baktfold_path=fold,
        output_path=tmp_path / "out.gbff",
        manifest_path=manifest,
    )
    ecs = [
        row
        for row in json.loads(manifest.read_text())["decisions"]
        if row["qualifier"] == "EC_number"
    ]
    assert len(ecs) == 2
    assert all(row["final_status"] == "shared_support" for row in ecs)


def test_different_authoritative_gene_is_not_supported_existing(tmp_path):
    base = write(
        tmp_path / "base.gbff", record_bytes("TEST", "T_0001", gene="original")
    )
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    egg = write(tmp_path / "egg.tsv", eggnog_bytes(ROW))
    manifest = tmp_path / "out.json"
    merge_egg(base, faa, egg, tmp_path / "out.gbff", manifest_path=manifest)
    row = next(
        row
        for row in json.loads(manifest.read_text())["decisions"]
        if row["field"] == "Preferred_name"
    )
    assert row["final_status"] != "supported_existing"


def test_reserved_weekday_query_row_is_not_silently_skipped():
    data = eggnog_bytes(ROW, ROW.replace("T_0001", "##Monday_bad"))
    with pytest.raises(MergeError):
        parse_eggnog_tsv(data)


def test_missing_sentinel_rechecked_after_suffix_cleaning():
    result = validate_preferred_name("NA_123", clean_suffix=True)
    assert result.status != "valid"


def test_cog_candidate_uses_note_qualifier(tmp_path):
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    egg = write(
        tmp_path / "egg.tsv",
        eggnog_bytes(ROW.replace("\t-\tabc\t", "\tCOG0001\tabc\t")),
    )
    manifest = tmp_path / "out.json"
    merge_egg(base, faa, egg, tmp_path / "out.gbff", manifest_path=manifest)
    row = next(
        row
        for row in json.loads(manifest.read_text())["decisions"]
        if row["field"] == "COG_category"
    )
    assert row["qualifier"] == "note"
    assert row["final_status"] == "emitted"


def test_restoration_noop_carries_imported_origin(tmp_path):
    base, faa, egg, parent = restored_inputs(tmp_path)
    out, manifest = tmp_path / "rerun.gbff", tmp_path / "rerun.json"
    stats = restore(
        base,
        faa,
        egg,
        out,
        manifest_path=manifest,
        translation_policy="validated-only",
        translation_evidence_manifest=parent,
    )
    assert stats["restored_translation_count"] == 0
    assert stats["insertions"] == 0
    assert out.read_bytes() == base.read_bytes()
    ledger = load_translation_evidence(manifest, parse_genbank_bytes(out.read_bytes()))
    assert ledger["T_0001"].origin == "imported_faa"


def test_merge_rerun_keeps_imported_origin_and_requires_optin(tmp_path):
    base, faa, egg, parent = restored_inputs(tmp_path)
    first, first_manifest = tmp_path / "first.gbff", tmp_path / "first.json"
    second, second_manifest = tmp_path / "second.gbff", tmp_path / "second.json"
    merge_egg(
        base,
        faa,
        egg,
        first,
        manifest_path=first_manifest,
        translation_evidence_manifest=parent,
        allow_imported_translations=True,
    )
    with pytest.raises(MergeError, match="allow-imported-translations"):
        merge_egg(first, faa, egg, second, translation_evidence_manifest=first_manifest)
    merge_egg(
        first,
        faa,
        egg,
        second,
        manifest_path=second_manifest,
        translation_evidence_manifest=first_manifest,
        allow_imported_translations=True,
    )
    assert first.read_bytes() == second.read_bytes()
    assert (
        load_translation_evidence(
            second_manifest, parse_genbank_bytes(second.read_bytes())
        )["T_0001"].origin
        == "imported_faa"
    )


def test_xlsx_query_alias_and_duplicate_alias_collision(tmp_path):
    import io

    openpyxl = pytest.importorskip("openpyxl")
    from enrich_bakta_lib.sources.eggnog import parse_eggnog_path

    workbook = openpyxl.Workbook()
    worksheet = workbook.active
    worksheet.title = "annotations"
    lines = eggnog_bytes(ROW).decode().splitlines()
    worksheet.append(lines[1].replace("#query", "query").split("\t"))
    worksheet.append(lines[2].split("\t"))
    stream = io.BytesIO()
    workbook.save(stream)
    path = write(tmp_path / "egg.xlsx", stream.getvalue())
    assert (
        parse_eggnog_path(path, expected_version="3.0.0-beta6").hits[0].query_id
        == "T_0001"
    )
    worksheet.cell(1, 2, "#query")
    stream = io.BytesIO()
    workbook.save(stream)
    workbook.close()
    path.write_bytes(stream.getvalue())
    with pytest.raises(MergeError, match="duplicate columns"):
        parse_eggnog_path(path, expected_version="3.0.0-beta6")


def test_paired_candidate_cannot_claim_only_one_insertion(tmp_path):
    from enrich_bakta_lib.core.merge_engine import (
        feature_uid,
        insertion_uid,
        qualifier_insertion,
    )

    base = parse_genbank_bytes(record_bytes("TEST", "T_0001"))
    gene, cds = base.features[1:]
    insertion = qualifier_insertion(base.data, cds, "gene", "abc", "eggNOG", "abc", 0)
    row = {
        "candidate_id": "candidate:test",
        "source_id": "eggNOG",
        "qualifier": "gene",
        "normalized_value": "abc",
        "final_status": "emitted",
        "target_feature_uids": [feature_uid(gene), feature_uid(cds)],
        "insertion_ids": [insertion_uid(insertion)],
        "supporting_candidate_ids": [],
    }
    with pytest.raises(MergeError, match="every target"):
        validate_candidate_ledger([row], [insertion], base=base)


def test_baktfold_retains_translation_lineage(tmp_path):
    from enrich_bakta_lib.sources.baktfold import graft

    base, _, _, parent = restored_inputs(tmp_path)
    fold = write(
        tmp_path / "fold.gbff",
        base.read_bytes().replace(
            b'/product="test protein"',
            b'/EC_number="1.2.3.4"\n                     /product="test protein"',
        ),
    )
    out, manifest = tmp_path / "fold-out.gbff", tmp_path / "fold-out.json"
    with pytest.raises(MergeError, match="allow-imported-translations"):
        graft(base, fold, out, translation_evidence_manifest=parent)
    graft(
        base,
        fold,
        out,
        manifest_path=manifest,
        translation_evidence_manifest=parent,
        allow_imported_translations=True,
    )
    assert (
        load_translation_evidence(manifest, parse_genbank_bytes(out.read_bytes()))[
            "T_0001"
        ].origin
        == "imported_faa"
    )


def test_parent_manifest_hash_uses_the_loaded_snapshot(tmp_path, monkeypatch):
    import enrich_bakta_lib.sources.eggnog as adapter
    from enrich_bakta_lib.core.merge_engine import sha256_bytes

    base, faa, egg, parent = restored_inputs(tmp_path)
    expected = sha256_bytes(parent.read_bytes())
    original = adapter.load_translation_evidence

    def replace_after_loading(path, document):
        ledger = original(path, document)
        path.write_bytes(b"a later disk version")
        return ledger

    monkeypatch.setattr(adapter, "load_translation_evidence", replace_after_loading)
    manifest = tmp_path / "out.json"
    adapter.merge(
        base,
        faa,
        egg,
        tmp_path / "out.gbff",
        manifest_path=manifest,
        translation_evidence_manifest=parent,
        allow_imported_translations=True,
    )
    assert (
        json.loads(manifest.read_text())["metadata"][
            "translation_evidence_parent_sha256"
        ]
        == expected
    )
