"""Tests for InterProScan streaming TSV parser and identity adapter (Phase B)."""

from __future__ import annotations

from pathlib import Path

import pytest

from enrich_bakta_lib.core.merge_engine import (
    MergeError,
    parse_genbank_bytes,
)
from enrich_bakta_lib.sources.interproscan import (
    MAX_RAW_LINE_BYTES,
    build_parser,
    main,
    merge,
    plan_interproscan,
    stream_interproscan_tsv,
)

FIXTURES_DIR = Path(__file__).parent / "fixtures" / "interproscan"

DUMMY_SEQ = (
    "MKGEHHYQVAVTWQGNLGAGTEDYRAYGRDHLIGAVGKADIAGSADPAFRGDATRWNPED"
    "LLVASLSACHKLWYLHLCATSGISVLAYQDNAVGVMREDAARGGFFTSVTLRPQVTVRDG"
    "DDLQLAQQLHEQAHHLCFIANSVNFPVACEPQAEYATR"
)
DUMMY_MD5 = "d73bd3e8f74ce6911a99f09cf9ad0077"
DUMMY_LEN = 158
DUMMY_QUERY = "OJMJKD_006675"


def record_bytes(
    record_id: str,
    locus_tag: str,
    translation: str | None = DUMMY_SEQ,
    *,
    comment_header: str = "COMMENT     Software: v0.1.0",
) -> bytes:
    """Construct a minimal valid GenBank flatfile byte sequence with a CDS."""
    cds_qualifiers = [
        '                     /product="test protein"',
        f'                     /locus_tag="{locus_tag}"',
        '                     /protein_id="gnl|Bakta|test"',
    ]
    if translation is not None:
        cds_qualifiers.append(f'                     /translation="{translation}"')
    cds_qualifiers.extend(
        [
            "                     /codon_start=1",
            "                     /transl_table=11",
        ]
    )

    lines = [
        f"LOCUS       {record_id:<16}{474:>12} bp    DNA     linear   BCT 01-JAN-2000",
        "DEFINITION  synthetic test record.",
        f"ACCESSION   {record_id}",
        f"VERSION     {record_id}",
        "KEYWORDS    .",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "            other sequences; artificial sequences.",
        comment_header,
        "FEATURES             Location/Qualifiers",
        "     source          1..474",
        '                     /organism="synthetic construct"',
        '                     /mol_type="genomic DNA"',
        "     gene            1..474",
        f'                     /locus_tag="{locus_tag}"',
        "     CDS             1..474",
        *cds_qualifiers,
        "ORIGIN",
        "        1 " + "atg" * 158,
        "//",
    ]
    return b"\n".join(line.encode("ascii") for line in lines) + b"\n"


def faa_bytes(query_id: str = DUMMY_QUERY, seq: str = DUMMY_SEQ) -> bytes:
    """Construct a minimal FAA FASTA byte sequence."""
    return f">{query_id} test protein\n{seq}\n".encode("ascii")


def test_layout_parsing_all_profiles() -> None:
    """Verify stream_interproscan_tsv correctly parses valid fixtures across all 4 layouts."""
    # 1. ipr-go-pathways
    sha, hits, diag = stream_interproscan_tsv(
        FIXTURES_DIR / "valid_ipr_go_pathways.tsv", "ipr-go-pathways"
    )
    assert len(sha) == 64
    assert len(hits) == 3
    assert diag["row_count"] == 3
    assert diag["unintegrated_count"] == 1
    assert diag["integrated_count"] == 2
    assert diag["pathway_counters"]["total_occurrences"] == 3
    assert "MetaCyc" in diag["pathway_counters"]["by_database"]
    assert "Reactome" in diag["pathway_counters"]["by_database"]
    assert hits[0].interpro_accession is None
    assert hits[1].interpro_accession == "IPR003718"
    assert hits[2].go_terms == ("GO:0006950", "GO:0009405")

    # 2. ipr-go
    sha_go, hits_go, diag_go = stream_interproscan_tsv(
        FIXTURES_DIR / "valid_ipr_go.tsv", "ipr-go"
    )
    assert len(hits_go) == 3
    assert hits_go[2].go_terms == ("GO:0006950", "GO:0009405")
    assert hits_go[2].pathways == ()

    # 3. ipr-pathways
    sha_pw, hits_pw, diag_pw = stream_interproscan_tsv(
        FIXTURES_DIR / "valid_ipr_pathways.tsv", "ipr-pathways"
    )
    assert len(hits_pw) == 3
    assert hits_pw[1].pathways == ("MetaCyc: PWY-5292",)
    assert hits_pw[1].go_terms == ()

    # 4. ipr-only
    sha_only, hits_only, diag_only = stream_interproscan_tsv(
        FIXTURES_DIR / "valid_ipr_only.tsv", "ipr-only"
    )
    assert len(hits_only) == 3
    assert hits_only[1].interpro_accession == "IPR003718"
    assert hits_only[1].go_terms == ()
    assert hits_only[1].pathways == ()


def test_1mib_line_ceiling_failure(tmp_path: Path) -> None:
    """Verify that lines exceeding 1 MiB without newline fail immediately."""
    oversized = tmp_path / "oversized.tsv"
    # Write a chunk exceeding MAX_RAW_LINE_BYTES with no newline
    with oversized.open("wb") as handle:
        handle.write(b"X" * (MAX_RAW_LINE_BYTES + 100))

    with pytest.raises(MergeError, match=r"exceeds maximum allowed size of 1 MiB"):
        stream_interproscan_tsv(oversized, "ipr-go-pathways")


def test_embedded_nul_byte_failure(tmp_path: Path) -> None:
    """Verify that lines containing NUL characters trigger MergeError."""
    nul_tsv = tmp_path / "nul.tsv"
    nul_tsv.write_bytes(
        b"OJMJKD_006675\td73bd3e8f74ce6911a99f09cf9ad0077\t158\tPAN\x00THER\tPTHR42830\tDESC\t4\t154\t1.5E-25\tT\t25-08-2026\t-\t-\n"
    )
    with pytest.raises(MergeError, match=r"contains embedded NUL character"):
        stream_interproscan_tsv(nul_tsv, "ipr-go-pathways")


def test_invalid_utf8_failure(tmp_path: Path) -> None:
    """Verify that invalid UTF-8 bytes trigger MergeError."""
    bad_utf8 = tmp_path / "bad_utf8.tsv"
    bad_utf8.write_bytes(
        b"OJMJKD_006675\td73bd3e8f74ce6911a99f09cf9ad0077\t158\tPANTHER\tPTHR42830\tDESC\xff\xfe\t4\t154\t1.5E-25\tT\t25-08-2026\t-\t-\n"
    )
    with pytest.raises(MergeError, match=r"contains invalid UTF-8"):
        stream_interproscan_tsv(bad_utf8, "ipr-go-pathways")


def test_conflicting_repeated_query_md5_length_failure(tmp_path: Path) -> None:
    """Verify conflicting MD5 or length on repeated queries raises MergeError."""
    # Conflicting MD5
    conflict_md5 = tmp_path / "conflict_md5.tsv"
    conflict_md5.write_text(
        f"{DUMMY_QUERY}\t{DUMMY_MD5}\t{DUMMY_LEN}\tPANTHER\tP1\tD1\t1\t10\t1.0\tT\t25-08-2026\t-\t-\n"
        f"{DUMMY_QUERY}\t00000000000000000000000000000000\t{DUMMY_LEN}\tPfam\tPF00001\tD2\t1\t10\t1.0\tT\t25-08-2026\t-\t-\n",
        encoding="utf-8",
    )
    with pytest.raises(MergeError, match=r"\[query_conflict\].*conflicting MD5"):
        stream_interproscan_tsv(conflict_md5, "ipr-go-pathways")

    # Conflicting length
    conflict_len = tmp_path / "conflict_len.tsv"
    conflict_len.write_text(
        f"{DUMMY_QUERY}\t{DUMMY_MD5}\t{DUMMY_LEN}\tPANTHER\tP1\tD1\t1\t10\t1.0\tT\t25-08-2026\t-\t-\n"
        f"{DUMMY_QUERY}\t{DUMMY_MD5}\t{DUMMY_LEN + 10}\tPfam\tPF00001\tD2\t1\t10\t1.0\tT\t25-08-2026\t-\t-\n",
        encoding="utf-8",
    )
    with pytest.raises(MergeError, match=r"\[query_conflict\].*conflicting.*length"):
        stream_interproscan_tsv(conflict_len, "ipr-go-pathways")


def test_antifam_qc_tracking(tmp_path: Path) -> None:
    """Verify AntiFam hits are tracked in diagnostics qc_flags."""
    antifam_tsv = tmp_path / "antifam.tsv"
    antifam_tsv.write_text(
        f"{DUMMY_QUERY}\t{DUMMY_MD5}\t{DUMMY_LEN}\tAntiFam\tANF00001\tShadow ORFs\t1\t50\t0.001\tT\t25-08-2026\t-\t-\n",
        encoding="utf-8",
    )
    _, hits, diag = stream_interproscan_tsv(antifam_tsv, "ipr-go-pathways")
    assert len(diag["anti_fam_hits"]) == 1
    assert diag["anti_fam_hits"][0]["signature_accession"] == "ANF00001"
    assert diag["anti_fam_hits"][0]["query_id"] == DUMMY_QUERY


def test_missing_queries_in_faa_vs_tsv_vs_gbff(tmp_path: Path) -> None:
    """Verify that query IDs present in TSV must match FAA and GenBank CDS."""
    base_data = record_bytes("REC1", DUMMY_QUERY)
    base = parse_genbank_bytes(base_data, "base")

    faa_file = tmp_path / "test.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    # 1. TSV query missing from FAA
    tsv_extra_query = tmp_path / "tsv_extra.tsv"
    tsv_extra_query.write_text(
        f"UNKNOWN_QUERY\t{DUMMY_MD5}\t{DUMMY_LEN}\tPANTHER\tP1\tD1\t1\t10\t1.0\tT\t25-08-2026\t-\t-\n",
        encoding="utf-8",
    )
    with pytest.raises(
        MergeError,
        match=r"(?s)InterProScan query IDs do not map exactly.*missing_from_faa",
    ):
        plan_interproscan(base, faa_file, tsv_extra_query)

    # 2. TSV query present in FAA but missing from GBFF
    faa_different = tmp_path / "different.faa"
    faa_different.write_bytes(faa_bytes("DIFFERENT_QUERY"))
    tsv_different = tmp_path / "tsv_diff.tsv"
    tsv_different.write_text(
        f"DIFFERENT_QUERY\t{DUMMY_MD5}\t{DUMMY_LEN}\tPANTHER\tP1\tD1\t1\t10\t1.0\tT\t25-08-2026\t-\t-\n",
        encoding="utf-8",
    )
    with pytest.raises(
        MergeError,
        match=r"(?s)InterProScan query IDs do not map exactly.*missing_from_gbff",
    ):
        plan_interproscan(base, faa_different, tsv_different)


def test_sequence_md5_length_mismatch_between_faa_and_tsv(tmp_path: Path) -> None:
    """Verify that sequence length or MD5 mismatch between FAA and TSV is rejected."""
    base = parse_genbank_bytes(record_bytes("REC1", DUMMY_QUERY), "base")
    # FAA has modified sequence (first char changed)
    bad_seq = "A" + DUMMY_SEQ[1:]
    bad_faa = tmp_path / "bad.faa"
    bad_faa.write_bytes(faa_bytes(DUMMY_QUERY, bad_seq))

    valid_tsv = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"
    with pytest.raises(
        MergeError,
        match=r"(?s)InterProScan FAA protein validation failed.*sequence_mismatches",
    ):
        plan_interproscan(base, bad_faa, valid_tsv)


def test_missing_translation_failure(tmp_path: Path) -> None:
    """Verify that a queried CDS lacking /translation in GenBank flatfile is rejected."""
    # Construct base GenBank record where CDS has NO translation qualifier
    no_trans_base_data = record_bytes("REC1", DUMMY_QUERY, translation=None)
    base = parse_genbank_bytes(no_trans_base_data, "base")

    faa_file = tmp_path / "test.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    valid_tsv = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"
    with pytest.raises(
        MergeError,
        match=r"(?s)InterProScan FAA/GBFF protein validation failed:.*missing_or_ambiguous_translation.*"
        + DUMMY_QUERY,
    ):
        plan_interproscan(base, faa_file, valid_tsv)


def test_known_legacy_restoration_preflight_failure(tmp_path: Path) -> None:
    """Verify base GenBank with legacy normalize_baktfold comment without evidence marker is rejected."""
    legacy_header = "COMMENT     normalize_baktfold.py: restore Bakta provenance"
    legacy_base_data = record_bytes("REC1", DUMMY_QUERY, comment_header=legacy_header)
    base = parse_genbank_bytes(legacy_base_data, "base")

    faa_file = tmp_path / "test.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    valid_tsv = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"
    with pytest.raises(
        MergeError, match=r"known legacy restoration provenance detected"
    ):
        plan_interproscan(base, faa_file, valid_tsv)


def test_standalone_merge_and_cli(tmp_path: Path) -> None:
    """Verify standalone merge() function and CLI parser."""
    # Test parser construction
    parser = build_parser()
    assert parser.prog == "merge_interproscan_bakta"

    # Minimal standalone merge execution
    base_file = tmp_path / "base.gbff"
    base_file.write_bytes(record_bytes("REC1", DUMMY_QUERY))

    faa_file = tmp_path / "sample.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    tsv_file = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"
    output_gbff = tmp_path / "output.gbff"
    manifest_json = tmp_path / "manifest.json"
    context_json = tmp_path / "context.json"

    res = merge(
        bakta_path=base_file,
        faa_path=faa_file,
        interproscan_path=tsv_file,
        output_path=output_gbff,
        manifest_path=manifest_json,
        context_report_path=context_json,
        merge_timestamp="2026-10-04T12:00:00Z",
    )

    assert res["self_check"] is True
    assert output_gbff.is_file()
    assert manifest_json.is_file()
    assert context_json.is_file()

    # CLI main() invocation with bad args returns 1 or raises SystemExit
    with pytest.raises(SystemExit):
        main(["--help"])

    rc = main(
        [
            "--bakta",
            str(base_file),
            "--faa",
            str(faa_file),
            "--interproscan",
            str(tsv_file),
            "--output",
            str(base_file),  # colliding output path should fail cleanly
        ]
    )
    assert rc == 1
