"""Tests for unified workflow integration with InterProScan (Phase D)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from enrich_bakta_lib.core.merge_engine import (
    MergeError,
    sha256_file,
)
from enrich_bakta_lib.workflows.enrich import enrich, main

FIXTURES_DIR = Path(__file__).parent / "fixtures" / "interproscan"

DUMMY_SEQ = (
    "MKGEHHYQVAVTWQGNLGAGTEDYRAYGRDHLIGAVGKADIAGSADPAFRGDATRWNPED"
    "LLVASLSACHKLWYLHLCATSGISVLAYQDNAVGVMREDAARGGFFTSVTLRPQVTVRDG"
    "DDLQLAQQLHEQAHHLCFIANSVNFPVACEPQAEYATR"
)
DUMMY_QUERY = "OJMJKD_006675"


def record_bytes(
    record_id: str,
    locus_tag: str,
    translation: str = DUMMY_SEQ,
) -> bytes:
    """Construct a minimal valid GenBank flatfile byte sequence with a CDS."""
    cds_qualifiers = [
        '                     /product="test protein"',
        f'                     /locus_tag="{locus_tag}"',
        '                     /protein_id="gnl|Bakta|test"',
        f'                     /translation="{translation}"',
        "                     /codon_start=1",
        "                     /transl_table=11",
    ]

    lines = [
        f"LOCUS       {record_id:<16}{474:>12} bp    DNA     linear   BCT 01-JAN-2000",
        "DEFINITION  synthetic test record.",
        f"ACCESSION   {record_id}",
        f"VERSION     {record_id}",
        "KEYWORDS    .",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "            other sequences; artificial sequences.",
        "COMMENT     Software: v0.1.0",
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


def test_unified_workflow_interproscan_alone(tmp_path: Path) -> None:
    """Verify unified enrich() works with --interproscan alone and produces interproscan context."""
    base_file = tmp_path / "base.gbff"
    base_file.write_bytes(record_bytes("REC1", DUMMY_QUERY))

    faa_file = tmp_path / "sample.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    tsv_file = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"
    out_gbff = tmp_path / "out.gbff"
    manifest = tmp_path / "manifest.json"
    context = tmp_path / "context.json"

    res = enrich(
        bakta_path=base_file,
        faa_path=faa_file,
        interproscan_path=tsv_file,
        output_path=out_gbff,
        manifest_path=manifest,
        context_report_path=context,
    )

    assert res["self_check"] is True
    assert out_gbff.is_file()
    assert manifest.is_file()
    assert context.is_file()

    ctx_data = json.loads(context.read_text(encoding="utf-8"))
    assert ctx_data["schema"] == "enrich-bakta.interproscan-context.v1"
    assert ctx_data["metadata"]["operation"] == "unified-enrichment"

    man_data = json.loads(manifest.read_text(encoding="utf-8"))
    assert man_data["metadata"]["context_report_sha256"] == sha256_file(context)


def test_unified_workflow_context_report_guard(tmp_path: Path) -> None:
    """Verify --context-report without eggnog or interproscan is rejected."""
    base_file = tmp_path / "base.gbff"
    base_file.write_bytes(record_bytes("REC1", DUMMY_QUERY))

    with pytest.raises(
        MergeError, match=r"--context-report requires --eggnog or --interproscan"
    ):
        enrich(
            bakta_path=base_file,
            output_path=tmp_path / "out.gbff",
            baktfold_path=base_file,
            context_report_path=tmp_path / "context.json",
        )


def test_unified_cli_interproscan(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify unified CLI main() parses --interproscan arguments properly."""
    base_file = tmp_path / "base.gbff"
    base_file.write_bytes(record_bytes("REC1", DUMMY_QUERY))

    faa_file = tmp_path / "sample.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    tsv_file = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"
    out_gbff = tmp_path / "out.gbff"

    monkeypatch.setattr(
        "sys.argv",
        [
            "enrich-bakta",
            "--bakta",
            str(base_file),
            "--faa",
            str(faa_file),
            "--interproscan",
            str(tsv_file),
            "--output",
            str(out_gbff),
        ],
    )
    rc = main()
    assert rc == 0
    assert out_gbff.is_file()
