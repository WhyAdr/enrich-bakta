from __future__ import annotations

import io
from pathlib import Path

import pytest
from Bio import SeqIO

from enrich_bakta import enrich
from graft_baktfold_additions import graft
from merge_eggnog_bakta import merge as merge_eggnog
from merge_eggnog_bakta import parse_eggnog_path, parse_eggnog_tsv
from merge_engine import (
    Insertion,
    MergeError,
    apply_insertions,
    format_qualifier,
    parse_genbank_bytes,
    reconcile_insertions,
    strict_parity_check,
)
from merge_kofamscan_bakta import merge, parse_kofam_table
from restore_bakta_translations import restore


def record_bytes(
    record_id: str,
    locus_tag: str,
    *,
    gene: str | None = None,
    notes: tuple[str, ...] = (),
    db_xrefs: tuple[str, ...] = (),
    ec_numbers: tuple[str, ...] = (),
    newline: bytes = b"\n",
) -> bytes:
    qualifiers = [f'                     /locus_tag="{locus_tag}"']
    if gene is not None:
        qualifiers.append(f'                     /gene="{gene}"')
    cds_qualifiers = [
        *(f'                     /note="{value}"' for value in notes),
        *(f'                     /db_xref="{value}"' for value in db_xrefs),
        *(f'                     /EC_number="{value}"' for value in ec_numbers),
        '                     /product="test protein"',
        f'                     /locus_tag="{locus_tag}"',
        '                     /protein_id="gnl|Bakta|test"',
        '                     /translation="MK"',
        "                     /codon_start=1",
        "                     /transl_table=11",
    ]
    if gene is not None:
        cds_qualifiers.append(f'                     /gene="{gene}"')
    lines = [
        f"LOCUS       {record_id:<16}{9:>12} bp    DNA     linear   BCT 01-JAN-2000",
        "DEFINITION  synthetic test record.",
        f"ACCESSION   {record_id}",
        f"VERSION     {record_id}",
        "KEYWORDS    .",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "            other sequences; artificial sequences.",
        "COMMENT     Software: v0.1.0",
        "FEATURES             Location/Qualifiers",
        "     source          1..9",
        '                     /organism="synthetic construct"',
        '                     /mol_type="genomic DNA"',
        "     gene            1..9",
        *qualifiers,
        "     CDS             1..9",
        *cds_qualifiers,
        "ORIGIN",
        "        1 atgaaataa",
        "//",
    ]
    return newline.join(line.encode("ascii") for line in lines) + newline


def write(path: Path, data: bytes) -> Path:
    path.write_bytes(data)
    return path


def kofam_bytes(*rows: str) -> bytes:
    header = (
        "# gene name           KO     thrshld  score   E-value KO definition\n"
        "#-------------------- ------ ------- ------ --------- ---------------------\n"
    )
    return (header + "\n".join(rows) + "\n").encode("ascii")


def eggnog_bytes(*rows: str, version: str = "3.0.0-beta6") -> bytes:
    header = (
        "#query\tseed_ortholog\tevalue\tscore\tCOG_category\tPreferred_name"
        "\tGOs\tEC\tKEGG_ko\tCAZy\tannotation_confidence\n"
    )
    return (f"## emapper-{version}\n" + header + "\n".join(rows) + "\n").encode("utf-8")


def test_baktfold_exact_prefixes_blank_gene_and_separate_provenance(
    tmp_path: Path,
) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001", gene=""))
    source = write(
        tmp_path / "source.gbff",
        record_bytes(
            "TEST",
            "T_0001",
            gene="abc",
            db_xrefs=(
                "pdb:1ABC",
                "pdbx:false",
                "cathode:false",
                "afdb_v6_extra:false",
                "EC:1.2.3.-",
            ),
        ),
    )
    output = tmp_path / "out.gbff"
    stats = graft(base, source, output, add_comment_note=False)
    text = output.read_text()
    assert stats["gene_added"] == 2  # paired gene and CDS features
    assert stats["ec_added"] == 1
    assert stats["xref_added"] == 1
    assert '/db_xref="pdb:1ABC"' in text
    assert "pdbx:false" not in text
    assert "cathode:false" not in text
    assert "afdb_v6_extra:false" not in text
    assert "Baktfold gene-symbol evidence" in text
    assert "Baktfold functional-annotation evidence" in text
    assert "protein structure similarity:Baktfold Foldseek" in text
    assert stats["self_check"] is True


def test_baktfold_rejects_ambiguous_gene_and_writes_nothing(tmp_path: Path) -> None:
    base = record_bytes("TEST", "T_0001")
    source = record_bytes("TEST", "T_0001", gene="abc").replace(
        b'/gene="abc"', b'/gene="abc"\n                     /gene="xyz"', 1
    )
    base_path = write(tmp_path / "base.gbff", base)
    source_path = write(tmp_path / "source.gbff", source)
    output = tmp_path / "out.gbff"
    with pytest.raises(MergeError, match="ambiguous gene symbols"):
        graft(base_path, source_path, output, add_comment_note=False)
    assert not output.exists()


@pytest.mark.parametrize(
    "mutation, message",
    [
        (
            lambda data: data.replace(
                b"     CDS             1..9", b"     CDS             2..9"
            ),
            "feature 2 differs",
        ),
        (
            lambda data: data.replace(
                b"     CDS             1..9",
                b'     CDS             1..9\n                     /locus_tag="T_0001"',
            ),
            "ambiguous locus_tag",
        ),
    ],
)
def test_parity_failures_prevent_output(tmp_path: Path, mutation, message: str) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    source = write(tmp_path / "source.gbff", mutation(record_bytes("TEST", "T_0001")))
    output = tmp_path / "out.gbff"
    with pytest.raises(MergeError, match=message):
        graft(base, source, output, add_comment_note=False)
    assert not output.exists()


def test_record_local_parity_catches_reordered_multirecord_input() -> None:
    first = record_bytes("ONE", "A_0001")
    second = record_bytes("TWO", "B_0001")
    base = parse_genbank_bytes(first + second)
    source = parse_genbank_bytes(second + first)
    with pytest.raises(MergeError, match="record 0 ID differs"):
        strict_parity_check(base, source)


def test_kofam_multiple_hits_existing_ko_manifest_and_no_ec_inference(
    tmp_path: Path,
) -> None:
    base = write(
        tmp_path / "base.gbff",
        record_bytes("TEST", "T_0001", notes=("KEGG:K00001",)),
    )
    faa = write(tmp_path / "base.faa", b">T_0001 test\nMK\n")
    kofam = write(
        tmp_path / "kofam.txt",
        kofam_bytes(
            "* T_0001 K00001 10.0 20.0 1e-5 existing enzyme [EC:1.2.3.4]",
            "* T_0001 K00002 11.0 21.0 0 second assignment [EC:9.9.9.9]",
        ),
    )
    output = tmp_path / "out.gbff"
    manifest = tmp_path / "manifest.json"
    stats = merge(
        base,
        faa,
        kofam,
        output,
        manifest_path=manifest,
        kofamscan_version="1.3.0",
        add_comment_note=False,
    )
    text = output.read_text()
    assert stats["kofam"]["hits"] == 2
    assert stats["kofam"]["hit_cds"] == 1
    assert stats["kofam"]["new_gene_ko_pairs"] == 1
    assert stats["kofam"]["already_represented_pairs"] == 1
    assert text.count('/db_xref="KEGG:K00002"') == 1
    assert '/db_xref="KEGG:K00001"' not in text
    assert text.count('/inference="profile:KofamScan:1.3.0"') == 1
    assert "/EC_number" not in text
    parsed_manifest = __import__("json").loads(manifest.read_text())
    hit_rows = [
        row for row in parsed_manifest["entries"] if row["entry_type"] == "kofam_hit"
    ]
    assert len(hit_rows) == 2
    assert "[EC:1.2.3.4]" in hit_rows[0]["definition"]


def test_combined_merge_is_idempotent_and_preserves_crlf(tmp_path: Path) -> None:
    base_data = record_bytes("TEST", "T_0001", newline=b"\r\n")
    base = write(tmp_path / "base.gbff", base_data)
    source = write(
        tmp_path / "source.gbff",
        record_bytes(
            "TEST", "T_0001", gene="abc", db_xrefs=("pdb:1ABC",), newline=b"\r\n"
        ),
    )
    faa = write(tmp_path / "base.faa", b">T_0001\r\nMK\r\n")
    kofam = write(
        tmp_path / "kofam.txt", kofam_bytes("* T_0001 K00002 1 2 3e-4 alpha protein")
    )
    first = tmp_path / "first.gbff"
    second = tmp_path / "second.gbff"
    merge(
        base,
        faa,
        kofam,
        first,
        baktfold_path=source,
        kofamscan_version="1",
    )
    rerun = merge(
        first,
        faa,
        kofam,
        second,
        baktfold_path=source,
        kofamscan_version="1",
    )
    assert first.read_bytes() == second.read_bytes()
    assert b"\n" not in first.read_bytes().replace(b"\r\n", b"")
    assert first.read_bytes().endswith(b"\r\n")
    assert rerun["insertions"] == 0


def test_no_additions_preserve_legacy_bytes_and_final_newline_style(
    tmp_path: Path,
) -> None:
    data = record_bytes("TEST", "T_0001").replace(
        b"COMMENT     Software", b"COMMENT     legacy=\xff; Software"
    )
    data = data.removesuffix(b"\n")
    base = write(tmp_path / "base.gbff", data)
    source = write(tmp_path / "source.gbff", data)
    output = tmp_path / "out.gbff"
    stats = graft(base, source, output, add_comment_note=False)
    assert stats["insertions"] == 0
    assert output.read_bytes() == data
    assert not output.read_bytes().endswith(b"\n")


def test_duplicate_feature_key_is_rejected() -> None:
    data = record_bytes("TEST", "T_0001")
    duplicate = b'     CDS             1..9\n                     /locus_tag="T_0001"\n'
    source_data = data.replace(b"ORIGIN\n", duplicate + b"ORIGIN\n")
    base = parse_genbank_bytes(data)
    source = parse_genbank_bytes(source_data)
    with pytest.raises(MergeError, match="duplicate feature key"):
        strict_parity_check(base, source)


def test_output_path_collision_and_missing_kofam_id_fail_before_write(
    tmp_path: Path,
) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    source = write(tmp_path / "source.gbff", record_bytes("TEST", "T_0001", gene="abc"))
    with pytest.raises(MergeError, match="output path must differ"):
        graft(base, source, base, add_comment_note=False)

    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    kofam = write(
        tmp_path / "kofam.txt", kofam_bytes("* WRONG K00001 1 2 3e-4 protein name")
    )
    output = tmp_path / "missing.gbff"
    with pytest.raises(MergeError, match="do not map exactly"):
        merge(base, faa, kofam, output, add_comment_note=False)
    assert not output.exists()


def test_parser_numeric_validation_and_insertion_allowlist() -> None:
    with pytest.raises(MergeError, match="invalid KO identifier"):
        parse_kofam_table(kofam_bytes("* T_0001 K1234 1 2 3e-4 protein name"))
    with pytest.raises(MergeError, match="invalid threshold"):
        parse_kofam_table(kofam_bytes("* T_0001 K00001 nope 2 3e-4 protein name"))
    forbidden = Insertion(0, b"x", "R", "CDS", "T", "product", "x", "test", "x")
    with pytest.raises(MergeError, match="allowlist rejected"):
        apply_insertions(b"base", [forbidden])


def test_quoted_formatter_wraps_and_biopython_parses() -> None:
    value = 'quoted "value" ' + "x" * 100
    rendered = format_qualifier("note", value, b"\n")
    assert all(len(line) <= 80 for line in rendered.splitlines())
    data = record_bytes("TEST", "T_0001").replace(b"ORIGIN\n", rendered + b"ORIGIN\n")
    records = list(SeqIO.parse(io.StringIO(data.decode()), "genbank"))
    assert len(records) == 1


def test_eggnog_header_driven_parser_and_merge_provenance(tmp_path: Path) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    eggnog = write(
        tmp_path / "annotations.tsv",
        eggnog_bytes(
            "T_0001\tseed.1\t1e-20\t50\tCOG0001\tname_12\tGO:0000001"
            "\tEC:1.2.3.4\tko:K00001\tGT2|Glycosyltransferase Family 2.\thhhhhhhhhhhhh"
        ),
    )
    output = tmp_path / "out.gbff"
    manifest = tmp_path / "manifest.json"
    stats = merge_eggnog(
        base,
        faa,
        eggnog,
        output,
        manifest_path=manifest,
        eggnog_version="3.0.0-beta6",
        clean_gene_suffix=True,
        add_comment_note=False,
    )
    text = output.read_text()
    assert stats["eggnog"]["emitted_values"] == 7  # paired gene plus five CDS values
    assert text.count('/gene="name"') == 2
    assert '/db_xref="GO:0000001"' in text
    assert '/EC_number="1.2.3.4"' in text
    assert '/db_xref="KEGG:K00001"' in text
    assert '/note="COG:COG0001"' in text
    assert '/db_xref="CAZy:GT2"' in text
    parsed = list(SeqIO.parse(io.StringIO(text), "genbank"))
    cds = next(feature for feature in parsed[0].features if feature.type == "CDS")
    assert cds.qualifiers["inference"] == [
        "DESCRIPTION:similar to AA sequence:eggNOG:seed.1"
    ]
    entries = __import__("json").loads(manifest.read_text())["entries"]
    gene = next(row for row in entries if row.get("field") == "Preferred_name")
    assert gene["raw_value"] == "name_12"
    assert gene["normalized_value"] == "name"


def test_eggnog_rejects_bad_fields_and_version() -> None:
    row = "T_0001\tseed\t1e-4\t10\tZ\t-\tGO:1\t-\t-\tGT2bad|wrong\thhhhhhhhhhhhh"
    with pytest.raises(MergeError, match="invalid GO"):
        parse_eggnog_tsv(eggnog_bytes(row))
    valid = "T_0001\tseed\t1e-4\t10\tZ\t-\t-\t-\t-\t-\thhhhhhhhhhhhh"
    with pytest.raises(MergeError, match="version mismatch"):
        parse_eggnog_tsv(eggnog_bytes(valid), expected_version="4")


def test_eggnog_xlsx_requires_explicit_version_before_optional_import(
    tmp_path: Path,
) -> None:
    with pytest.raises(MergeError, match="--eggnog-version is required"):
        parse_eggnog_path(tmp_path / "annotations.xlsx")


def test_eggnog_confidence_and_existing_values_are_not_reinserted(
    tmp_path: Path,
) -> None:
    base = write(
        tmp_path / "base.gbff",
        record_bytes("TEST", "T_0001", db_xrefs=("GO:0000001", "KEGG:K00001")),
    )
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    eggnog = write(
        tmp_path / "annotations.tsv",
        eggnog_bytes(
            "T_0001\tseed\t1e-4\t10\tCOG0001\t-\tGO:0000001\t1.2.3.4"
            "\tK00001\t-\thhlhhhhhhhhhh"
        ),
    )
    output = tmp_path / "out.gbff"
    stats = merge_eggnog(
        base,
        faa,
        eggnog,
        output,
        min_confidence="high",
        add_comment_note=False,
    )
    text = output.read_text()
    assert text.count('/db_xref="GO:0000001"') == 1
    assert text.count('/db_xref="KEGG:K00001"') == 1
    assert "/EC_number" not in text
    assert '/note="COG:COG0001"' in text  # COG deliberately has no confidence code
    assert stats["eggnog"]["confidence_filtered"] == 1


def test_reconciliation_collapses_duplicates_and_handles_gene_conflicts() -> None:
    shared = Insertion(
        10,
        b"",
        "TEST",
        "CDS",
        "T_0001",
        "db_xref",
        "KEGG:K00001",
        "KofamScan",
        "K00001",
        1,
    )
    duplicate = Insertion(
        10,
        b"",
        "TEST",
        "CDS",
        "T_0001",
        "db_xref",
        "KEGG:K00001",
        "eggNOG",
        "K00001",
        2,
    )
    baktfold_gene = Insertion(
        10, b"", "TEST", "CDS", "T_0001", "gene", "bakt", "Baktfold", "bakt", 3
    )
    eggnog_gene = Insertion(
        10, b"", "TEST", "CDS", "T_0001", "gene", "egg", "eggNOG", "egg", 4
    )
    reconciled, rows, stats = reconcile_insertions(
        [shared, duplicate, baktfold_gene, eggnog_gene]
    )
    assert len(reconciled) == 1
    assert reconciled[0].qualifier == "db_xref"
    assert stats["gene_conflicts"] == 1
    assert any(row["status"] == "exact_duplicate_collapsed" for row in rows)
    preferred, _, _ = reconcile_insertions(
        [baktfold_gene, eggnog_gene], gene_conflict_policy="prefer-eggnog"
    )
    assert [item.value for item in preferred] == ["egg"]


def test_unified_merge_reconciles_eggnog_and_kofam(tmp_path: Path) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    kofam = write(tmp_path / "kofam.txt", kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha"))
    eggnog = write(
        tmp_path / "annotations.tsv",
        eggnog_bytes("T_0001\tseed\t1e-4\t10\t-\t-\t-\t-\tK00001\t-\thhhhhhhhhhhhh"),
    )
    output = tmp_path / "out.gbff"
    manifest = tmp_path / "manifest.json"
    enrich(
        bakta_path=base,
        output_path=output,
        faa_path=faa,
        kofamscan_path=kofam,
        eggnog_path=eggnog,
        manifest_path=manifest,
        add_comment_note=False,
    )
    assert output.read_text().count('/db_xref="KEGG:K00001"') == 1
    entries = __import__("json").loads(manifest.read_text())["entries"]
    assert any(row.get("status") == "exact_duplicate_collapsed" for row in entries)


def test_translation_restoration_is_limited_to_eggnog_pseudogenes(
    tmp_path: Path,
) -> None:
    base_data = record_bytes("TEST", "T_0001").replace(
        b'                     /translation="MK"\n',
        b'                     /pseudogene="unitary"\n',
    )
    base = write(tmp_path / "base.gbff", base_data)
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    eggnog = write(
        tmp_path / "annotations.tsv",
        eggnog_bytes("T_0001\tseed\t1e-4\t10\t-\t-\t-\t-\t-\t-\thhhhhhhhhhhhh"),
    )
    output = tmp_path / "restored.gbff"
    stats = restore(base, faa, eggnog, output, eggnog_version="3.0.0-beta6")
    assert stats["restored_translation_count"] == 1
    assert output.read_bytes().count(b'/translation="MK"') == 1
