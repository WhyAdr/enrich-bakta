from __future__ import annotations

import json
from pathlib import Path

from jsonschema import Draft202012Validator
from test_merge_pipeline import eggnog_bytes, kofam_bytes, record_bytes, write

from enrich_bakta_lib.sources.baktfold import graft
from enrich_bakta_lib.sources.eggnog import merge as merge_eggnog
from enrich_bakta_lib.sources.kofam import merge as merge_kofam
from enrich_bakta_lib.workflows.enrich import enrich
from enrich_bakta_lib.workflows.restore_translations import restore

SCHEMA_PATH = Path(__file__).parents[1] / "schemas" / "merge-manifest.v2.schema.json"


def _validator() -> Draft202012Validator:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


def _assert_valid(path: Path) -> None:
    _validator().validate(json.loads(path.read_text(encoding="utf-8")))


def test_checked_in_merge_manifest_schema_is_valid() -> None:
    _validator()


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
            "unified",
            lambda output, manifest: enrich(
                bakta_path=base,
                output_path=output,
                faa_path=faa,
                baktfold_path=fold,
                kofamscan_path=ko,
                eggnog_path=egg,
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
