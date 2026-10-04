"""Phase A: Tests for InterProScan shared semantics, member-note witnesses, and layouts."""

from __future__ import annotations

from pathlib import Path

import pytest
from test_merge_pipeline import record_bytes

from enrich_bakta_lib.core.decisions import (
    CandidateDecision,
    validate_candidate_ledger,
)
from enrich_bakta_lib.core.merge_engine import (
    MergeError,
    feature_uid,
    parse_genbank_bytes,
)
from enrich_bakta_lib.sources.interproscan import (
    parse_member_dbs,
    validate_tsv_row_layout,
)
from enrich_bakta_lib.sources.value_rules import (
    find_member_note_witness,
    parse_standalone_member_note,
    validate_go_term,
    validate_interpro_accession,
    validate_pfam_accession,
    validate_tigrfam_accession,
)


def test_parse_standalone_member_note_positive() -> None:
    # Bare Pfam
    match = parse_standalone_member_note("PFAM:PF00001")
    assert match is not None
    assert match.db == "Pfam"
    assert match.accession == "PF00001"
    assert match.version is None
    assert match.canonical_note == "PFAM:PF00001"
    assert match.raw_note == "PFAM:PF00001"

    # Versioned Pfam
    match = parse_standalone_member_note("PFAM:PF12345.33")
    assert match is not None
    assert match.db == "Pfam"
    assert match.accession == "PF12345"
    assert match.version == "33"
    assert match.canonical_note == "PFAM:PF12345"
    assert match.raw_note == "PFAM:PF12345.33"

    # Bare TIGRFAM
    match = parse_standalone_member_note("TIGRFAM:TIGR00001")
    assert match is not None
    assert match.db == "TIGRFAM"
    assert match.accession == "TIGR00001"
    assert match.version is None
    assert match.canonical_note == "TIGRFAM:TIGR00001"

    # Versioned TIGRFAM
    match = parse_standalone_member_note("TIGRFAM:TIGR00002.5")
    assert match is not None
    assert match.db == "TIGRFAM"
    assert match.accession == "TIGR00002"
    assert match.version == "5"
    assert match.canonical_note == "TIGRFAM:TIGR00002"


@pytest.mark.parametrize(
    "invalid_note",
    [
        "NOT_PFAM:PF00001",
        "prefix_PFAM:PF00001",
        "MY_TIGRFAM:TIGR00001",
        "PFAM:PF00001.0",  # 0 is not a positive version
        "PFAM:PF00001.",
        "PFAM:PF00001.abc",
        "PFAM:PF00001.1.2",
        "PFAM:PF00001-v1",
        "similar to PFAM:PF00001",
        "no PFAM:PF00001 found",
        "PFAM:PF00001 (score: 5)",
        "PFAMs: PF00001",  # composite eggNOG note
        "PFAM:PF00001, PFAM:PF00002",
        "PFAM:PF001",  # too short
        "PFAM:PF123456",  # too long
        "TIGRFAM:TIGR1",
        "TIGRFAM:TIGR123456",
        "",
        "-",
    ],
)
def test_parse_standalone_member_note_negative(invalid_note: str) -> None:
    assert parse_standalone_member_note(invalid_note) is None


def test_find_member_note_witness() -> None:
    notes = [
        "some product note",
        "PFAM:PF00005.33",
        "TIGRFAM:TIGR00042.1",
        "PFAMs: PF00196, PF04545",
    ]
    # Positive witness lookup returns exact raw note
    assert find_member_note_witness(notes, "PFAM:PF00005") == "PFAM:PF00005.33"
    assert find_member_note_witness(notes, "TIGRFAM:TIGR00042") == "TIGRFAM:TIGR00042.1"

    # Absent or negative returns None
    assert find_member_note_witness(notes, "PFAM:PF00006") is None
    # Composite eggNOG note is NOT a witness
    assert find_member_note_witness(notes, "PFAM:PF00196") is None


def test_accession_validators() -> None:
    assert validate_interpro_accession("IPR003718") == "IPR003718"
    assert validate_interpro_accession("IPR0037189") is None
    assert validate_interpro_accession("IPR00371") is None
    assert validate_interpro_accession("PF00001") is None

    assert validate_go_term("GO:0005524") == "GO:0005524"
    assert validate_go_term("GO:000552") is None
    assert validate_go_term("GO:00055240") is None
    assert validate_go_term("GO:0005524(IEA)") is None

    assert validate_pfam_accession("PF00005") == ("PF00005", None)
    assert validate_pfam_accession("PF00005.33") == ("PF00005", "33")
    assert validate_pfam_accession("PF00005.0") is None
    assert validate_pfam_accession("PF005") is None

    assert validate_tigrfam_accession("TIGR00001") == ("TIGR00001", None)
    assert validate_tigrfam_accession("TIGR00001.2") == ("TIGR00001", "2")
    assert validate_tigrfam_accession("TIGR001") is None


def test_pfam_version_witness_survives_ledger_validation(tmp_path: Path) -> None:
    # Base feature has versioned note: PFAM:PF00005.33
    base_data = record_bytes("TEST", "T_0001", notes=("PFAM:PF00005.33",))
    base = parse_genbank_bytes(base_data, "base")
    feature = next(f for f in base.features if f.feature_type == "CDS")
    f_uid = feature_uid(feature)

    # Candidate decision claims supported_existing for canonical PFAM:PF00005
    candidate = CandidateDecision(
        candidate_id="candidate:0123456789abcdef0123456789abcdef",
        source_id="InterProScan",
        source_sha256="0" * 64,
        target_feature_uids=(f_uid,),
        field="Pfam",
        qualifier="note",
        raw_value="PF00005",
        normalized_value="PFAM:PF00005",
        planned_status="existing",
        candidate_role="substantive_evidence",
        evidence_class="member:Pfam",
        reason_code="value_already_present",
        final_status="supported_existing",
        confidence_diagnostics={"support_witness": "PFAM:PF00005.33"},
        reason="Pfam match supported by existing base note",
    )

    # validate_candidate_ledger should accept the versioned witness against bare canonical
    validate_candidate_ledger(
        [candidate.as_dict()],
        final_insertions=(),
        base=base,
        output=base,
    )


def test_pfam_version_witness_rejects_unsupported_or_malformed(tmp_path: Path) -> None:
    # Base feature has free-text note: similar to PFAM:PF00005
    base_data = record_bytes("TEST", "T_0001", notes=("similar to PFAM:PF00005",))
    base = parse_genbank_bytes(base_data, "base")
    feature = next(f for f in base.features if f.feature_type == "CDS")
    f_uid = feature_uid(feature)

    candidate = CandidateDecision(
        candidate_id="candidate:0123456789abcdef0123456789abcdef",
        source_id="InterProScan",
        source_sha256="0" * 64,
        target_feature_uids=(f_uid,),
        field="Pfam",
        qualifier="note",
        raw_value="PF00005",
        normalized_value="PFAM:PF00005",
        planned_status="existing",
        candidate_role="substantive_evidence",
        evidence_class="member:Pfam",
        reason_code="value_already_present",
        final_status="supported_existing",
        confidence_diagnostics={},
        reason="",
    )

    with pytest.raises(MergeError, match="absent from the final output"):
        validate_candidate_ledger(
            [candidate.as_dict()],
            final_insertions=(),
            base=base,
            output=base,
        )


def test_parse_member_dbs() -> None:
    assert parse_member_dbs("Pfam,TIGRFAM") == ("Pfam", "TIGRFAM")
    assert parse_member_dbs("Pfam, TIGRFAM") == ("Pfam", "TIGRFAM")
    assert parse_member_dbs("Pfam") == ("Pfam",)
    assert parse_member_dbs("") == ()
    assert parse_member_dbs(["Pfam", "TIGRFAM"]) == ("Pfam", "TIGRFAM")

    with pytest.raises(
        MergeError, match="unsupported InterProScan member database 'CDD'"
    ):
        parse_member_dbs("Pfam,CDD")

    with pytest.raises(
        MergeError, match="unsupported InterProScan member database 'SUPERFAMILY'"
    ):
        parse_member_dbs("SUPERFAMILY")


def test_validate_tsv_row_layout_all_profiles() -> None:
    fixtures_dir = Path(__file__).parent / "fixtures" / "interproscan"

    # 1. ipr-go-pathways (13 unintegrated, 15 integrated)
    lines = (
        (fixtures_dir / "valid_ipr_go_pathways.tsv")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    assert len(lines) == 3
    # Row 1: unintegrated (13 cols)
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[0].split("\t"), "ipr-go-pathways", 1
    )
    assert not is_int and ipr is None and desc is None and not gos and not paths
    # Row 2: integrated (15 cols, no GO, has pathways)
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[1].split("\t"), "ipr-go-pathways", 2
    )
    assert (
        is_int
        and ipr == "IPR003718"
        and desc == "OsmC/Ohr family"
        and not gos
        and paths == ("MetaCyc: PWY-5292",)
    )
    # Row 3: integrated (15 cols, has GO and pathways)
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[2].split("\t"), "ipr-go-pathways", 3
    )
    assert is_int and ipr == "IPR036102" and gos == ("GO:0006950", "GO:0009405")
    assert paths == ("MetaCyc: PWY-5292", "Reactome: R-BTA-1234")

    # 2. ipr-go (13 unintegrated, 14 integrated)
    lines = (
        (fixtures_dir / "valid_ipr_go.tsv")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[0].split("\t"), "ipr-go", 1
    )
    assert not is_int
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[2].split("\t"), "ipr-go", 3
    )
    assert (
        is_int
        and ipr == "IPR036102"
        and gos == ("GO:0006950", "GO:0009405")
        and not paths
    )

    # 3. ipr-pathways (13 unintegrated, 14 integrated)
    lines = (
        (fixtures_dir / "valid_ipr_pathways.tsv")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[0].split("\t"), "ipr-pathways", 1
    )
    assert not is_int
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[2].split("\t"), "ipr-pathways", 3
    )
    assert is_int and ipr == "IPR036102" and not gos and len(paths) == 2

    # 4. ipr-only (13 unintegrated, 13 integrated)
    lines = (
        (fixtures_dir / "valid_ipr_only.tsv")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[0].split("\t"), "ipr-only", 1
    )
    assert not is_int
    is_int, ipr, desc, gos, paths = validate_tsv_row_layout(
        lines[1].split("\t"), "ipr-only", 2
    )
    assert is_int and ipr == "IPR003718" and not gos and not paths


def test_validate_tsv_row_layout_negative_fixtures() -> None:
    fixtures_dir = Path(__file__).parent / "fixtures" / "interproscan"

    # Native 11-column export (lookup disabled) -> unsupported_tsv_layout
    lines_11 = (
        (fixtures_dir / "invalid_11_cols_lookup_off.tsv")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    with pytest.raises(
        MergeError,
        match="11 columns \\(lookup disabled\\); this native layout is unsupported",
    ):
        validate_tsv_row_layout(lines_11[0].split("\t"), "ipr-go-pathways", 1)

    # Incoherent placeholder (IPR accession '-' but description present)
    lines_incoherent = (
        (fixtures_dir / "invalid_placeholder_incoherent.tsv")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    with pytest.raises(MergeError, match="placeholder_incoherence"):
        validate_tsv_row_layout(lines_incoherent[0].split("\t"), "ipr-go-pathways", 1)
    # IPR accession present but description is '-'
    with pytest.raises(MergeError, match="placeholder_incoherence"):
        validate_tsv_row_layout(lines_incoherent[1].split("\t"), "ipr-go-pathways", 2)

    # Ambiguous 14 cols when ipr-go-pathways expects 15
    lines_14 = (
        (fixtures_dir / "invalid_ambiguous_14_cols.tsv")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    with pytest.raises(MergeError, match="expected 15 for layout 'ipr-go-pathways'"):
        validate_tsv_row_layout(lines_14[0].split("\t"), "ipr-go-pathways", 1)

    # Invalid layout name
    with pytest.raises(
        MergeError, match="unsupported InterProScan TSV layout 'bad-layout'"
    ):
        validate_tsv_row_layout(lines_14[0].split("\t"), "bad-layout", 1)
