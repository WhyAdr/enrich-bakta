"""Tests for InterProScan ledger reconciliation and merge engine integration (Phase C)."""

from __future__ import annotations

import json
from pathlib import Path

from enrich_bakta_lib.core.decisions import (
    build_candidate_ledger,
)
from enrich_bakta_lib.core.merge_engine import (
    cds_by_locus,
    insertion_uid,
    parse_genbank_bytes,
    qualifier_insertion,
    reconcile_insertions,
    sha256_file,
)
from enrich_bakta_lib.sources.interproscan import (
    merge as merge_interproscan,
)
from enrich_bakta_lib.sources.interproscan import (
    plan_interproscan,
)

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
    *,
    db_xrefs: tuple[str, ...] = (),
    notes: tuple[str, ...] = (),
) -> bytes:
    """Construct a minimal valid GenBank flatfile byte sequence with a CDS."""
    cds_qualifiers = [
        *(f'                     /db_xref="{val}"' for val in db_xrefs),
        *(f'                     /note="{val}"' for val in notes),
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


def test_shared_support_reconciliation(tmp_path: Path) -> None:
    """Verify exact GO proposal between eggNOG and InterProScan results in shared_support."""
    base_data = record_bytes("REC1", DUMMY_QUERY)
    base = parse_genbank_bytes(base_data, "base")
    cds = cds_by_locus(base)[DUMMY_QUERY]

    # Create two insertions proposing the exact same GO term: GO:0006950
    # One from eggNOG, one from InterProScan
    ins_eggnog = qualifier_insertion(
        base.data,
        cds,
        "db_xref",
        "GO:0006950",
        "eggNOG",
        "GO:0006950",
        order=0,
        candidate_role="functional_proposal",
        evidence_class="GO",
    )
    ins_interpro = qualifier_insertion(
        base.data,
        cds,
        "db_xref",
        "GO:0006950",
        "InterProScan",
        "GO:0006950",
        order=1,
        candidate_role="functional_proposal",
        evidence_class="GO",
    )

    all_insertions = [ins_eggnog, ins_interpro]
    reconciled, _recon_rows, _recon_stats = reconcile_insertions(all_insertions)

    # Exactly one insertion survives deduplication
    assert len(reconciled) == 1
    surviving = reconciled[0]

    candidate_rows, ledger_counts = build_candidate_ledger(
        base,
        all_insertions,
        reconciled,
        [],
        source_hashes={"eggNOG": "a" * 64, "InterProScan": "b" * 64},
    )

    # Both candidates must be represented in decisions
    decisions = [
        r for r in candidate_rows if r.get("entry_type") == "candidate_decision"
    ]
    assert len(decisions) == 2

    c_eggnog = next(d for d in decisions if d["source_id"] == "eggNOG")
    c_interpro = next(d for d in decisions if d["source_id"] == "InterProScan")

    # Both must have shared_support
    assert c_eggnog["final_status"] == "shared_support"
    assert c_eggnog["reason_code"] == "shared_support"
    assert c_interpro["final_status"] == "shared_support"
    assert c_interpro["reason_code"] == "shared_support"

    # Both must cross-reference each other
    assert c_eggnog["supporting_candidate_ids"] == [c_interpro["candidate_id"]]
    assert c_interpro["supporting_candidate_ids"] == [c_eggnog["candidate_id"]]

    # Both must point to the single surviving insertion
    surviving_uid = insertion_uid(surviving)
    assert c_eggnog["insertion_ids"] == [surviving_uid]
    assert c_interpro["insertion_ids"] == [surviving_uid]


def test_provenance_pruning_when_all_proposals_supported_existing(
    tmp_path: Path,
) -> None:
    """Verify InterProScan producer provenance is suppressed when all proposals are already in base."""
    # Base already has the exact InterPro accession and GO terms from valid_ipr_go_pathways.tsv
    # Hit 1: unintegrated (no proposals)
    # Hit 2: IPR003718, Pfam PF02566
    # Hit 3: IPR036102, GO:0006950, GO:0009405
    base_data = record_bytes(
        "REC1",
        DUMMY_QUERY,
        db_xrefs=(
            "InterPro:IPR003718",
            "InterPro:IPR036102",
            "GO:0006950",
            "GO:0009405",
        ),
        notes=(
            "PFAM:PF02566.20",
        ),  # Versioned Pfam note matches canonical bare PF02566
    )
    base = parse_genbank_bytes(base_data, "base")

    faa_file = tmp_path / "test.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    tsv_file = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"

    plan = plan_interproscan(base, faa_file, tsv_file, add_comment_note=False)

    # In the plan, since all proposals are already supported_existing,
    # new_qualifier_features was empty, so no producer provenance insertion was planned!
    reconciled, _, _ = reconcile_insertions(plan.insertions)
    candidate_rows, counts = build_candidate_ledger(
        base,
        plan.insertions,
        reconciled,
        plan.evidence_rows,
        source_hashes={"InterProScan": plan.source_sha256},
    )

    # Check evidence rows: all candidates must be supported_existing
    for row in plan.evidence_rows:
        assert row["status"] == "supported_existing"

    # There should be 0 insertions planned and 0 candidate insertions
    assert len(reconciled) == 0
    assert counts["candidate_insertion_count"] == 0
    assert counts["accepted_candidate_count"] == 5  # 2 IPR, 2 GO, 1 Pfam
    assert counts["emitted_candidate_count"] == 0


def test_provenance_emitted_when_new_annotation_survives(tmp_path: Path) -> None:
    """Verify InterProScan producer provenance is emitted when at least one new annotation is inserted."""
    # Base has NO annotations
    base_data = record_bytes("REC1", DUMMY_QUERY)
    base = parse_genbank_bytes(base_data, "base")

    faa_file = tmp_path / "test.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    tsv_file = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"

    plan = plan_interproscan(base, faa_file, tsv_file)
    reconciled, _, _ = reconcile_insertions(plan.insertions)
    candidate_rows, counts = build_candidate_ledger(
        base,
        plan.insertions,
        reconciled,
        plan.evidence_rows,
        source_hashes={"InterProScan": plan.source_sha256},
    )

    decisions = [
        r for r in candidate_rows if r.get("entry_type") == "candidate_decision"
    ]
    prov_decision = next(
        d for d in decisions if d["candidate_role"] == "producer_provenance"
    )

    assert prov_decision["final_status"] == "emitted"
    assert len(prov_decision["supporting_candidate_ids"]) > 0


def test_suppressed_candidate_count_aggregation(tmp_path: Path) -> None:
    """Verify suppressed_candidate_count in ledger aggregates all suppressed_* final statuses."""
    base_data = record_bytes("REC1", DUMMY_QUERY)
    base = parse_genbank_bytes(base_data, "base")
    cds = cds_by_locus(base)[DUMMY_QUERY]

    # Two conflicting gene symbols from different sources with default skip policy
    ins1 = qualifier_insertion(
        base.data,
        cds,
        "gene",
        "geneA",
        "eggNOG",
        "geneA",
        order=0,
        candidate_role="functional_proposal",
    )
    ins2 = qualifier_insertion(
        base.data,
        cds,
        "gene",
        "geneB",
        "Baktfold",
        "geneB",
        order=1,
        candidate_role="functional_proposal",
    )

    reconciled, _, _ = reconcile_insertions([ins1, ins2], gene_conflict_policy="skip")
    assert len(reconciled) == 0

    candidate_rows, counts = build_candidate_ledger(
        base,
        [ins1, ins2],
        reconciled,
        [],
        source_hashes={"eggNOG": "a" * 64, "Baktfold": "b" * 64},
        gene_conflict_policy="skip",
    )

    decisions = [
        r for r in candidate_rows if r.get("entry_type") == "candidate_decision"
    ]
    for d in decisions:
        assert d["final_status"] == "suppressed_conflict"

    assert counts["suppressed_candidate_count"] == 2


def test_context_report_sha256_in_manifest(tmp_path: Path) -> None:
    """Verify finalize_merge embeds context_report_sha256 into manifest metadata."""
    base_file = tmp_path / "base.gbff"
    base_file.write_bytes(record_bytes("REC1", DUMMY_QUERY))

    faa_file = tmp_path / "sample.faa"
    faa_file.write_bytes(faa_bytes(DUMMY_QUERY))

    tsv_file = FIXTURES_DIR / "valid_ipr_go_pathways.tsv"
    output_gbff = tmp_path / "output.gbff"
    manifest_json = tmp_path / "manifest.json"
    context_json = tmp_path / "context.json"

    merge_interproscan(
        bakta_path=base_file,
        faa_path=faa_file,
        interproscan_path=tsv_file,
        output_path=output_gbff,
        manifest_path=manifest_json,
        context_report_path=context_json,
    )

    assert manifest_json.is_file()
    assert context_json.is_file()

    manifest_data = json.loads(manifest_json.read_text(encoding="utf-8"))
    expected_context_sha = sha256_file(context_json)

    assert "context_report_sha256" in manifest_data["metadata"]
    assert manifest_data["metadata"]["context_report_sha256"] == expected_context_sha
