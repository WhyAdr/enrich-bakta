from __future__ import annotations

from pathlib import Path

import pytest
from test_merge_pipeline import eggnog_bytes

from enrich_bakta_lib.core.merge_engine import MergeError
from enrich_bakta_lib.sources.eggnog import parse_eggnog_tsv
from enrich_bakta_lib.sources.value_rules import (
    structured_note_tokens,
    validate_ec_value,
    validate_preferred_name,
    validate_structural_xref,
)

ROW = "T_0001\tseed\t1e-4\t10\t-\t-\t-\t-\t-\t-\thhhhhhhhhhhhh"


def test_unknown_eggnog_version_does_not_inherit_v3_contract() -> None:
    with pytest.raises(MergeError, match="unsupported eggNOG producer version"):
        parse_eggnog_tsv(eggnog_bytes(ROW, version="99.0.0"))


def test_comment_free_table_requires_explicit_schema() -> None:
    data = eggnog_bytes(ROW).split(b"\n", 1)[1]
    with pytest.raises(MergeError, match="metadata is absent"):
        parse_eggnog_tsv(data)
    table = parse_eggnog_tsv(data, schema_id="eggnog-v3.0.0-beta6")
    assert table.schema_id == "eggnog-v3.0.0-beta6"
    assert table.schema_selection_source == "explicit_schema"


def test_explicit_schema_can_audit_unknown_producer() -> None:
    table = parse_eggnog_tsv(
        eggnog_bytes(ROW, version="99.0.0"),
        schema_id="eggnog-v3.0.0-beta6",
    )
    assert table.schema_selection_source == "explicit_schema"
    assert table.producer_version_raw == "99.0.0"


def test_empty_and_repeated_schema_metadata_fail_closed() -> None:
    empty = eggnog_bytes(ROW).replace(
        b"## emapper-3.0.0-beta6\n",
        b"## emapper-3.0.0-beta6\n## confidence field order:\n",
    )
    with pytest.raises(MergeError, match="legend is empty"):
        parse_eggnog_tsv(empty)
    repeated = eggnog_bytes(ROW).replace(
        b"## emapper-3.0.0-beta6\n",
        b"## emapper-3.0.0-beta6\n## emapper-3.0.0-beta6\n",
    )
    with pytest.raises(MergeError, match="duplicate producer"):
        parse_eggnog_tsv(repeated)


def test_tabbed_reserved_hash_row_after_header_is_not_silently_skipped() -> None:
    data = (
        eggnog_bytes(ROW)
        .replace(
            b"#query\tseed_ortholog",
            b"#query\tseed_ortholog",
        )
        .replace(
            b"T_0001\tseed",
            b"##T_0001\tseed",
        )
    )
    with pytest.raises(MergeError, match="reserved '#'"):
        parse_eggnog_tsv(data)


def test_actual_emapper_command_metadata_with_tabs_is_accepted() -> None:
    data = eggnog_bytes(ROW).replace(
        b"## emapper-3.0.0-beta6\n",
        b"## emapper-3.0.0-beta6\n## /usr/bin/emapper.py\t--cpu\t4\n",
    )
    assert parse_eggnog_tsv(data).version == "3.0.0-beta6"


def test_pinned_official_v3_fixture_uses_explicit_auditable_schema() -> None:
    fixture = (
        Path(__file__).parent
        / "fixtures"
        / "official"
        / "eggnog-v3"
        / "test_diamond.emapper.annotations"
    )
    table = parse_eggnog_tsv(
        fixture.read_bytes(),
        schema_id="eggnog-v3.0.0-beta6",
    )
    assert table.schema_id == "eggnog-v3.0.0-beta6"
    assert table.schema_selection_source == "explicit_schema"
    assert table.producer_version_raw == "v3.0.0-beta5-2-ga3fdb5b"


def test_owner_beta6_fixture_uses_automatic_registry_selection() -> None:
    fixture = (
        Path(__file__).parent
        / "fixtures"
        / "official"
        / "eggnog-v3"
        / "beta6-owner-sample.annotations"
    )
    table = parse_eggnog_tsv(fixture.read_bytes())
    assert table.schema_id == "eggnog-v3.0.0-beta6"
    assert table.schema_selection_source == "version_registry"


def test_value_validation_retains_raw_normalized_status_and_reason() -> None:
    ec = validate_ec_value(" EC:1.2.3.- ", source="fixture")
    assert ec.raw == " EC:1.2.3.- "
    assert ec.normalized == "1.2.3.-"
    assert ec.status == "valid"
    assert ec.reason == "partial EC grammar"

    invalid = validate_ec_value("3.5.1.n3", source="Baktfold")
    assert invalid.normalized == "3.5.1.n3"
    assert invalid.status == "invalid"
    assert "invalid EC value" in invalid.reason

    missing = validate_preferred_name("-", clean_suffix=True)
    assert missing.status == "missing"
    assert missing.normalized is None

    structural = validate_structural_xref("pdb: 1ABC", source="Baktfold")
    assert structural.normalized == "pdb:1ABC"
    assert structural.status == "valid"
    assert validate_structural_xref("KEGG:K00001", source="Baktfold").status == (
        "not_applicable"
    )


@pytest.mark.parametrize(
    ("note", "expected"),
    [
        ("GO:0000001", ("GO:0000001",)),
        ("GO:0000001, KEGG:K00001", ("GO:0000001", "KEGG:K00001")),
        ("notGO:0000001x", ()),
        ("KEGG:K000010", ()),
        ("CAZy:GH10_2", ("CAZy:GH10_2",)),
        ("CAZy:GH10", ("CAZy:GH10",)),
    ],
)
def test_structured_note_tokens_require_complete_affirmative_tokens(note, expected):
    assert structured_note_tokens(note) == expected
