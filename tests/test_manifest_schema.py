from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator
from jsonschema.exceptions import ValidationError
from test_merge_pipeline import eggnog_bytes, kofam_bytes, record_bytes, write

from enrich_bakta_lib.sources.baktfold import graft
from enrich_bakta_lib.sources.eggnog import merge as merge_eggnog
from enrich_bakta_lib.sources.interproscan import merge as merge_interproscan
from enrich_bakta_lib.sources.kofam import merge as merge_kofam
from enrich_bakta_lib.workflows.enrich import enrich
from enrich_bakta_lib.workflows.restore_translations import restore

SCHEMA_PATH = Path(__file__).parents[1] / "schemas" / "merge-manifest.v2.schema.json"
IPS_CONTEXT_SCHEMA_PATH = (
    Path(__file__).parents[1] / "schemas" / "interproscan-context.v1.schema.json"
)
CONTEXT_SCHEMA_PATH = Path(__file__).parents[1] / "schemas" / "context.v1.schema.json"


def _validator() -> Draft202012Validator:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


def _ips_context_validator() -> Draft202012Validator:
    schema = json.loads(IPS_CONTEXT_SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


def _context_validator() -> Draft202012Validator:
    schema = json.loads(CONTEXT_SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


def _assert_valid(path: Path) -> None:
    _validator().validate(json.loads(path.read_text(encoding="utf-8")))


def test_checked_in_merge_manifest_schema_is_valid() -> None:
    _validator()
    _ips_context_validator()
    _context_validator()


def test_all_manifest_writers_conform_to_checked_in_schema(tmp_path: Path) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    fold = write(
        tmp_path / "fold.gbff",
        record_bytes(
            "TEST",
            "T_0001",
            gene="abc",
            ec_numbers=("1.2.3.4",),
            db_xrefs=("pdb:1ABC",),
        ),
    )
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    ko = write(tmp_path / "ko.txt", kofam_bytes("* T_0001 K00001 1 2 3e-4 alpha"))
    egg = write(
        tmp_path / "egg.tsv",
        eggnog_bytes(
            "T_0001\tseed\t1e-4\t10\tC\tabc\t-\t1.2.3.4\tK00001\t-\thhhhhhhhhhhhh"
        ),
    )
    mk_md5 = hashlib.md5(b"MK").hexdigest()
    ips = write(
        tmp_path / "ips.tsv",
        f"T_0001\t{mk_md5}\t2\tPfam\tPF00001\tdesc\t1\t2\t1.0\tT\t25-08-2026\t-\t-\n".encode(
            "utf-8"
        ),
    )

    manifests: list[Path] = []
    cases = [
        (
            "baktfold",
            lambda output, manifest: graft(base, fold, output, manifest_path=manifest),
        ),
        (
            "kofam",
            lambda output, manifest: merge_kofam(
                base, faa, ko, output, manifest_path=manifest
            ),
        ),
        (
            "eggnog",
            lambda output, manifest: merge_eggnog(
                base, faa, egg, output, manifest_path=manifest
            ),
        ),
        (
            "interproscan",
            lambda output, manifest: merge_interproscan(
                bakta_path=base,
                faa_path=faa,
                interproscan_path=ips,
                output_path=output,
                manifest_path=manifest,
                interproscan_version="5.59-91.0",
            ),
        ),
        (
            "unified",
            lambda output, manifest: enrich(
                bakta_path=base,
                output_path=output,
                faa_path=faa,
                baktfold_path=fold,
                kofamscan_path=ko,
                eggnog_path=egg,
                interproscan_path=ips,
                interproscan_version="5.59-91.0",
                manifest_path=manifest,
            ),
        ),
    ]
    for label, writer in cases:
        manifest = tmp_path / f"{label}.json"
        writer(tmp_path / f"{label}.gbff", manifest)
        manifests.append(manifest)

    translationless = write(
        tmp_path / "translationless.gbff",
        record_bytes("RESTORE", "T_0001").replace(b'/translation="MK"', b"/pseudo"),
    )
    restoration = tmp_path / "restoration.json"
    restore(
        translationless,
        faa,
        egg,
        tmp_path / "restored.gbff",
        manifest_path=restoration,
    )
    manifests.append(restoration)

    for manifest in manifests:
        _assert_valid(manifest)


def test_interproscan_context_schema_validation(tmp_path: Path) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    mk_md5 = hashlib.md5(b"MK").hexdigest()
    ips = write(
        tmp_path / "ips.tsv",
        f"T_0001\t{mk_md5}\t2\tPfam\tPF00001\tdesc\t1\t2\t1.0\tT\t25-08-2026\t-\t-\n".encode(
            "utf-8"
        ),
    )
    out_gbff = tmp_path / "ips_out.gbff"
    out_manifest = tmp_path / "ips_manifest.json"
    out_context = tmp_path / "ips_context.json"

    merge_interproscan(
        bakta_path=base,
        faa_path=faa,
        interproscan_path=ips,
        output_path=out_gbff,
        manifest_path=out_manifest,
        interproscan_version="5.59-91.0",
        context_report_path=out_context,
    )

    _assert_valid(out_manifest)
    ctx_data = json.loads(out_context.read_text(encoding="utf-8"))
    _ips_context_validator().validate(ctx_data)


def test_schema_negative_rejects_corrupted_interproscan_entries(tmp_path: Path) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    mk_md5 = hashlib.md5(b"MK").hexdigest()
    ips = write(
        tmp_path / "ips.tsv",
        f"T_0001\t{mk_md5}\t2\tPfam\tPF00001\tdesc\t1\t2\t1.0\tT\t25-08-2026\t-\t-\n".encode(
            "utf-8"
        ),
    )
    out_manifest = tmp_path / "manifest.json"
    merge_interproscan(
        bakta_path=base,
        faa_path=faa,
        interproscan_path=ips,
        output_path=tmp_path / "out.gbff",
        manifest_path=out_manifest,
        interproscan_version="5.59-91.0",
    )

    validator = _validator()
    valid_data = json.loads(out_manifest.read_text(encoding="utf-8"))
    validator.validate(valid_data)

    # 1. Corrupt normalized_value to integer
    corrupted_data = json.loads(json.dumps(valid_data))
    corrupted_data["entries"][0]["normalized_value"] = 12345
    with pytest.raises(ValidationError):
        validator.validate(corrupted_data)

    # 2. Corrupt target_feature_uids to empty list
    corrupted_data2 = json.loads(json.dumps(valid_data))
    corrupted_data2["entries"][0]["target_feature_uids"] = []
    with pytest.raises(ValidationError):
        validator.validate(corrupted_data2)

    # 3. Corrupt final_status to invented status
    corrupted_data3 = json.loads(json.dumps(valid_data))
    corrupted_data3["entries"][0]["final_status"] = "invented_status"
    with pytest.raises(ValidationError):
        validator.validate(corrupted_data3)

    # 4. Corrupt metadata interproscan_version to number
    corrupted_data4 = json.loads(json.dumps(valid_data))
    corrupted_data4["metadata"]["interproscan_version"] = 5.59
    with pytest.raises(ValidationError):
        validator.validate(corrupted_data4)


def test_schema_negative_rejects_corrupted_interproscan_context(tmp_path: Path) -> None:
    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
    faa = write(tmp_path / "base.faa", b">T_0001\nMK\n")
    mk_md5 = hashlib.md5(b"MK").hexdigest()
    ips = write(
        tmp_path / "ips.tsv",
        f"T_0001\t{mk_md5}\t2\tGene3D\tG3D001\tdesc\t1\t2\t1.0\tT\t25-08-2026\t-\t-\n".encode(
            "utf-8"
        ),
    )
    out_context = tmp_path / "context.json"
    merge_interproscan(
        bakta_path=base,
        faa_path=faa,
        interproscan_path=ips,
        output_path=tmp_path / "out.gbff",
        interproscan_version="5.59-91.0",
        context_report_path=out_context,
    )

    validator = _ips_context_validator()
    valid_data = json.loads(out_context.read_text(encoding="utf-8"))
    validator.validate(valid_data)

    # 1. Negative pathway occurrences
    corrupted1 = json.loads(json.dumps(valid_data))
    corrupted1["entries"][0]["pathway_occurrences"] = -1
    with pytest.raises(ValidationError):
        validator.validate(corrupted1)

    # 2. More than 50 context signatures
    corrupted2 = json.loads(json.dumps(valid_data))
    corrupted2["entries"][0]["context_signatures"] = ["sig"] * 51
    with pytest.raises(ValidationError):
        validator.validate(corrupted2)

    # 3. Invalid SHA256 in metadata
    corrupted3 = json.loads(json.dumps(valid_data))
    corrupted3["metadata"]["interproscan_sha256"] = "invalid_hash"
    with pytest.raises(ValidationError):
        validator.validate(corrupted3)
