#!/usr/bin/env python3
"""Run real InterProScan C14, SM, and BK71A acceptance checks."""

from __future__ import annotations

import argparse
import json
import time
import traceback
from collections.abc import Callable
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator

from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    MergeError,
    load_translation_evidence,
    parse_genbank_bytes,
    sha256_file,
)
from enrich_bakta_lib.sources.interproscan import merge as merge_interproscan
from enrich_bakta_lib.workflows.enrich import enrich

DEFAULT_BASELINE_PATH = Path("docs/data/interproscan-scientific-baseline.json")

EXPECTED_C14_IMPORTED_ORIGINS = (
    "OJMJKD_002410",
    "OJMJKD_010260",
    "OJMJKD_012060",
    "OJMJKD_020450",
)


def _validate_json_schema(schema_path: Path, payload: dict[str, Any]) -> None:
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    Draft202012Validator(schema).validate(payload)


EXPECTED_C14 = {
    "output_sha256": "472b0bd60d3990dda6307c89de7153d2ab3fbfd9f9324180e8889cb1b47bc9e1",
    "insertions": 34080,
    "candidate_count": 37942,
    "emitted_candidate_count": 34079,
    "accepted_candidate_count": 37942,
    "context_report_sha256": "59cad1cb6b7f97adb5467745c8d085916941cebe1d49c07c210b28bd81da0843",
    "four_source_output_sha256": "fbf4bb5eae6a35b3405469c443ee36570f04254bf128f2b925c9020b7bbadb0e",
    "four_source_insertions": 74910,
}

EXPECTED_SM = {
    "output_sha256": "0df5dce4a2aa615e313eaf858d948e8ff8b8698cc11ffc4bca444a2860d4e471",
    "insertions": 33542,
    "candidate_count": 37412,
    "emitted_candidate_count": 33540,
    "accepted_candidate_count": 37412,
    "context_report_sha256": "17eb69dc2710adc9808ba7052b4734af0ff63664e43970f97252f30252abcb9c",
    "four_source_output_sha256": "69388307bedc377ce9acdf6ccce838e9af9eca38dff51f13ed9b76263c8d1035",
    "four_source_insertions": 73889,
}

EXPECTED_SM_MISSING_TRANSLATIONS = (
    "NNLMDO_002705",
    "NNLMDO_007000",
    "NNLMDO_007005",
    "NNLMDO_013725",
    "NNLMDO_023610",
)

EXPECTED_SM_IMPORTED_ORIGINS = (
    "NNLMDO_002705",
    "NNLMDO_007000",
    "NNLMDO_007005",
    "NNLMDO_011850",
    "NNLMDO_013725",
    "NNLMDO_023610",
)


def _check(
    label: str,
    action: Callable[[], dict[str, Any]],
    *,
    results: list[dict[str, Any]],
    report: Path,
    tracebacks: Path,
) -> None:
    started = time.monotonic()
    try:
        row: dict[str, Any] = {
            "label": label,
            "status": "passed",
            "result": action(),
        }
    except Exception as exc:
        row = {
            "label": label,
            "status": "failed",
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
        (tracebacks / f"{label}.traceback.txt").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
    row["seconds"] = round(time.monotonic() - started, 2)
    results.append(row)
    report.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(row), flush=True)


def run_checks(
    root: Path,
    output: Path,
    baseline_path: Path | None = None,
) -> int:
    output.mkdir(parents=True, exist_ok=True)
    report = output / "interproscan-real-checks.json"
    results: list[dict[str, Any]] = []

    c14_expected = dict(EXPECTED_C14)
    sm_expected = dict(EXPECTED_SM)

    effective_baseline = baseline_path or (root / DEFAULT_BASELINE_PATH)
    if not effective_baseline.is_file():
        raise FileNotFoundError(f"Scientific baseline not found: {effective_baseline}")
    base_json = json.loads(effective_baseline.read_text(encoding="utf-8"))
    if base_json.get("tool") != {"name": "enrich-bakta", "version": TOOL_VERSION}:
        raise AssertionError(
            f"scientific baseline tool identity does not match enrich-bakta {TOOL_VERSION}"
        )
    c14_exp = base_json.get("samples", {}).get("C14", {}).get("expected", {})
    if "interproscan_only" in c14_exp:
        c14_expected.update(c14_exp["interproscan_only"])
    if "four_source" in c14_exp:
        c14_expected["four_source_output_sha256"] = c14_exp["four_source"][
            "output_sha256"
        ]
        c14_expected["four_source_insertions"] = c14_exp["four_source"]["insertions"]
    sm_exp = base_json.get("samples", {}).get("SM", {}).get("expected", {})
    if "restored_interproscan_only" in sm_exp:
        sm_expected.update(sm_exp["restored_interproscan_only"])
    if "four_source" in sm_exp:
        sm_expected["four_source_output_sha256"] = sm_exp["four_source"][
            "output_sha256"
        ]
        sm_expected["four_source_insertions"] = sm_exp["four_source"]["insertions"]

    manifest_schema = root / "schemas" / "merge-manifest.v2.schema.json"
    ips_context_schema = root / "schemas" / "interproscan-context.v1.schema.json"
    combined_context_schema = root / "schemas" / "context.v1.schema.json"

    # 1. C14 standalone enrichment
    def c14_check() -> dict[str, Any]:
        c14_gbff = root / "data" / "C14" / "bakta" / "C14-NMZ.gbff"
        c14_faa = root / "data" / "C14" / "bakta" / "C14-NMZ.faa"
        c14_tsv = (
            root
            / "data"
            / "C14"
            / "evidence"
            / "interproscan"
            / "C14-NMZ.interproscan.tsv"
        )
        out_gbff = output / "C14-interpro-enriched.gbff"
        manifest_path = output / "C14-interpro-enriched.manifest.json"
        context_path = output / "C14-interpro-context.json"

        stats = enrich(
            bakta_path=c14_gbff,
            faa_path=c14_faa,
            interproscan_path=c14_tsv,
            interproscan_version="5.59-91.0",
            output_path=out_gbff,
            manifest_path=manifest_path,
            context_report_path=context_path,
        )

        assert stats["self_check"] is True
        assert stats["output_sha256"] == c14_expected["output_sha256"]
        assert stats["insertions"] == c14_expected["insertions"]

        man_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        _validate_json_schema(manifest_schema, man_data)
        ctx_data = json.loads(context_path.read_text(encoding="utf-8"))
        _validate_json_schema(ips_context_schema, ctx_data)
        for key in (
            "candidate_count",
            "emitted_candidate_count",
            "accepted_candidate_count",
            "context_report_sha256",
        ):
            if man_data["metadata"][key] != c14_expected[key]:
                raise AssertionError(
                    f"C14 {key} mismatch: {man_data['metadata'][key]} != {c14_expected[key]}"
                )

        ctx_sha = sha256_file(context_path)
        assert ctx_sha == c14_expected["context_report_sha256"]
        return {
            "output_sha256": stats["output_sha256"],
            "insertions": stats["insertions"],
            "context_report_sha256": ctx_sha,
        }

    _check(
        "c14-interpro-acceptance",
        c14_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    # 2. C14 No-op: Output as input
    def c14_noop_check() -> dict[str, Any]:
        enriched_gbff = output / "C14-interpro-enriched.gbff"
        c14_faa = root / "data" / "C14" / "bakta" / "C14-NMZ.faa"
        c14_tsv = (
            root
            / "data"
            / "C14"
            / "evidence"
            / "interproscan"
            / "C14-NMZ.interproscan.tsv"
        )
        noop_gbff = output / "C14-interpro-noop.gbff"

        stats = enrich(
            bakta_path=enriched_gbff,
            faa_path=c14_faa,
            interproscan_path=c14_tsv,
            interproscan_version="5.59-91.0",
            output_path=noop_gbff,
        )
        assert stats["self_check"] is True
        assert stats["insertions"] == 0, (
            f"Expected 0 insertions, got {stats['insertions']}"
        )
        assert stats["output_sha256"] == c14_expected["output_sha256"], (
            f"No-op output SHA changed: {stats['output_sha256']} != {c14_expected['output_sha256']}"
        )
        return {"insertions": 0, "output_sha256": stats["output_sha256"]}

    _check(
        "c14-interpro-noop",
        c14_noop_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    # 3. C14 Determinism: Fresh run produces identical file digests
    def c14_determinism_check() -> dict[str, Any]:
        c14_gbff = root / "data" / "C14" / "bakta" / "C14-NMZ.gbff"
        c14_faa = root / "data" / "C14" / "bakta" / "C14-NMZ.faa"
        c14_tsv = (
            root
            / "data"
            / "C14"
            / "evidence"
            / "interproscan"
            / "C14-NMZ.interproscan.tsv"
        )
        fresh_gbff = output / "C14-fresh-determinism.gbff"
        fresh_manifest = output / "C14-fresh-determinism.manifest.json"
        fresh_context = output / "C14-fresh-determinism.context.json"

        stats = enrich(
            bakta_path=c14_gbff,
            faa_path=c14_faa,
            interproscan_path=c14_tsv,
            interproscan_version="5.59-91.0",
            output_path=fresh_gbff,
            manifest_path=fresh_manifest,
            context_report_path=fresh_context,
        )
        assert stats["output_sha256"] == c14_expected["output_sha256"]
        assert sha256_file(fresh_context) == c14_expected["context_report_sha256"]
        return {
            "output_sha256": stats["output_sha256"],
            "context_sha256": sha256_file(fresh_context),
        }

    _check(
        "c14-interpro-determinism",
        c14_determinism_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    # 4. SM pristine rejection
    def sm_pristine_check() -> dict[str, Any]:
        sm_gbff = root / "data" / "SM" / "bakta" / "SM-NMZ.gbff"
        sm_faa = root / "data" / "SM" / "bakta" / "SM-NMZ.faa"
        sm_tsv = (
            root
            / "data"
            / "SM"
            / "evidence"
            / "interproscan"
            / "SM-NMZ.interproscan.tsv"
        )
        out_gbff = output / "SM-pristine-fail.gbff"

        try:
            enrich(
                bakta_path=sm_gbff,
                faa_path=sm_faa,
                interproscan_path=sm_tsv,
                interproscan_version="5.59-91.0",
                output_path=out_gbff,
            )
        except MergeError as exc:
            msg = str(exc)
            for query in EXPECTED_SM_MISSING_TRANSLATIONS:
                if query not in msg:
                    raise AssertionError(
                        f"Expected missing query {query} in error message: {msg}"
                    )
            return {"rejected_as_expected": True, "error_snippet": msg[:200]}
        raise AssertionError(
            "SM pristine run unexpectedly succeeded without translation evidence!"
        )

    _check(
        "sm-pristine-rejection",
        sm_pristine_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    # 5. SM restored acceptance
    def sm_restored_check() -> dict[str, Any]:
        sm_restored_gbff = (
            root
            / "data"
            / "SM"
            / "derived"
            / "restored"
            / "SM-NMZ-eggnog-restored.gbff"
        )
        sm_restored_manifest = (
            root
            / "data"
            / "SM"
            / "derived"
            / "restored"
            / "SM-NMZ-eggnog-restored.manifest.json"
        )
        sm_faa = root / "data" / "SM" / "bakta" / "SM-NMZ.faa"
        sm_tsv = (
            root
            / "data"
            / "SM"
            / "evidence"
            / "interproscan"
            / "SM-NMZ.interproscan.tsv"
        )
        out_gbff = output / "SM-interpro-enriched.gbff"
        manifest_path = output / "SM-interpro-enriched.manifest.json"
        context_path = output / "SM-interpro-context.json"

        stats = enrich(
            bakta_path=sm_restored_gbff,
            faa_path=sm_faa,
            interproscan_path=sm_tsv,
            interproscan_version="5.59-91.0",
            translation_evidence_manifest=sm_restored_manifest,
            allow_imported_translations=True,
            output_path=out_gbff,
            manifest_path=manifest_path,
            context_report_path=context_path,
        )

        assert stats["self_check"] is True
        assert stats["output_sha256"] == sm_expected["output_sha256"]
        assert stats["insertions"] == sm_expected["insertions"]

        man_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        _validate_json_schema(manifest_schema, man_data)
        ctx_data = json.loads(context_path.read_text(encoding="utf-8"))
        _validate_json_schema(ips_context_schema, ctx_data)
        for key in (
            "candidate_count",
            "emitted_candidate_count",
            "accepted_candidate_count",
            "context_report_sha256",
        ):
            if man_data["metadata"][key] != sm_expected[key]:
                raise AssertionError(
                    f"SM {key} mismatch: {man_data['metadata'][key]} != {sm_expected[key]}"
                )

        te = man_data["metadata"].get("translation_evidence", {})
        assert set(te.keys()) == set(EXPECTED_SM_IMPORTED_ORIGINS), (
            f"Exact translation evidence mismatch: {set(te.keys())} != {set(EXPECTED_SM_IMPORTED_ORIGINS)}"
        )
        for query in EXPECTED_SM_IMPORTED_ORIGINS:
            if te[query].get("origin") != "imported_faa":
                raise AssertionError(f"Expected origin imported_faa for {query}")
        reloaded = load_translation_evidence(
            manifest_path, parse_genbank_bytes(out_gbff.read_bytes(), "SM output")
        )
        assert set(reloaded.keys()) == set(EXPECTED_SM_IMPORTED_ORIGINS)

        return {
            "output_sha256": stats["output_sha256"],
            "insertions": stats["insertions"],
            "imported_translations_carried_forward": len(EXPECTED_SM_IMPORTED_ORIGINS),
        }

    _check(
        "sm-restored-acceptance",
        sm_restored_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    # 6. SM standalone restored run
    def sm_standalone_check() -> dict[str, Any]:
        sm_restored_gbff = (
            root
            / "data"
            / "SM"
            / "derived"
            / "restored"
            / "SM-NMZ-eggnog-restored.gbff"
        )
        sm_restored_manifest = (
            root
            / "data"
            / "SM"
            / "derived"
            / "restored"
            / "SM-NMZ-eggnog-restored.manifest.json"
        )
        sm_faa = root / "data" / "SM" / "bakta" / "SM-NMZ.faa"
        sm_tsv = (
            root
            / "data"
            / "SM"
            / "evidence"
            / "interproscan"
            / "SM-NMZ.interproscan.tsv"
        )
        out_gbff = output / "SM-standalone-restored.gbff"
        manifest_path = output / "SM-standalone-restored.manifest.json"
        context_path = output / "SM-standalone-restored.context.json"

        stats = merge_interproscan(
            bakta_path=sm_restored_gbff,
            faa_path=sm_faa,
            interproscan_path=sm_tsv,
            interproscan_version="5.59-91.0",
            translation_evidence_manifest=sm_restored_manifest,
            allow_imported_translations=True,
            output_path=out_gbff,
            manifest_path=manifest_path,
            context_report_path=context_path,
        )

        assert stats["self_check"] is True
        assert stats["output_sha256"] == sm_expected["output_sha256"]
        assert stats["insertions"] == sm_expected["insertions"]
        man_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        _validate_json_schema(manifest_schema, man_data)
        ctx_data = json.loads(context_path.read_text(encoding="utf-8"))
        _validate_json_schema(ips_context_schema, ctx_data)
        reloaded = load_translation_evidence(
            manifest_path,
            parse_genbank_bytes(out_gbff.read_bytes(), "SM standalone output"),
        )
        assert set(reloaded.keys()) == set(EXPECTED_SM_IMPORTED_ORIGINS)
        return {
            "output_sha256": stats["output_sha256"],
            "insertions": stats["insertions"],
        }

    _check(
        "sm-standalone-restored",
        sm_standalone_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    # 7. BK71A lineage rejection
    def bk71a_check() -> dict[str, Any]:
        bk_restored = (
            root / "data" / "BK71A" / "derived" / "restored" / "BK71A-restored.gbff"
        )
        bk_tsv = (
            root
            / "data"
            / "BK71A"
            / "evidence"
            / "interproscan"
            / "BK71A.interproscan.tsv"
        )
        c14_faa = root / "data" / "C14" / "bakta" / "C14-NMZ.faa"
        out_gbff = output / "BK71A-fail.gbff"

        try:
            enrich(
                bakta_path=bk_restored,
                faa_path=c14_faa,
                interproscan_path=bk_tsv,
                interproscan_version="5.59-91.0",
                output_path=out_gbff,
            )
        except MergeError as exc:
            msg = str(exc)
            if "known legacy restoration provenance detected" not in msg:
                raise AssertionError(
                    f"Expected legacy provenance preflight rejection: {msg}"
                )
            return {"rejected_as_expected": True, "error_snippet": msg[:200]}
        raise AssertionError("BK71A run unexpectedly succeeded!")

    _check(
        "bk71a-lineage-rejection",
        bk71a_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    # 8. C14 four-source enrichment
    def c14_four_source_check() -> dict[str, Any]:
        c14_restored_gbff = (
            root
            / "data"
            / "C14"
            / "derived"
            / "restored"
            / "C14-NMZ-eggnog-restored.gbff"
        )
        c14_restored_manifest = (
            root
            / "data"
            / "C14"
            / "derived"
            / "restored"
            / "C14-NMZ-eggnog-restored.manifest.json"
        )
        c14_faa = root / "data" / "C14" / "bakta" / "C14-NMZ.faa"

        stats = enrich(
            bakta_path=c14_restored_gbff,
            faa_path=c14_faa,
            baktfold_path=root / "data" / "C14" / "baktfold" / "C14_NMZ.gbff",
            baktfold_invalid_ec_policy="skip",
            kofamscan_path=root / "data" / "C14" / "evidence" / "C14-KofamKOALA.txt",
            eggnog_path=root
            / "data"
            / "C14"
            / "evidence"
            / "C14-NMZ-query.emapper.annotations",
            interproscan_path=root
            / "data"
            / "C14"
            / "evidence"
            / "interproscan"
            / "C14-NMZ.interproscan.tsv",
            interproscan_version="5.59-91.0",
            translation_evidence_manifest=c14_restored_manifest,
            allow_imported_translations=True,
            output_path=output / "C14-four-source.gbff",
            manifest_path=output / "C14-four-source.manifest.json",
            context_report_path=output / "C14-four-source.context.json",
        )
        assert stats["self_check"] is True
        assert stats["insertions"] == c14_expected["four_source_insertions"]
        assert stats["output_sha256"] == c14_expected["four_source_output_sha256"]
        man_path = output / "C14-four-source.manifest.json"
        ctx_path = output / "C14-four-source.context.json"
        out_path = output / "C14-four-source.gbff"
        man_data = json.loads(man_path.read_text(encoding="utf-8"))
        _validate_json_schema(manifest_schema, man_data)
        ctx_data = json.loads(ctx_path.read_text(encoding="utf-8"))
        _validate_json_schema(combined_context_schema, ctx_data)
        _validate_json_schema(ips_context_schema, ctx_data["interproscan"])
        reloaded = load_translation_evidence(
            man_path, parse_genbank_bytes(out_path.read_bytes(), "C14 four-source")
        )
        assert set(reloaded.keys()) == set(EXPECTED_C14_IMPORTED_ORIGINS)
        return {
            "output_sha256": stats["output_sha256"],
            "insertions": stats["insertions"],
        }

    _check(
        "c14-four-source-acceptance",
        c14_four_source_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    # 9. SM four-source enrichment
    def sm_four_source_check() -> dict[str, Any]:
        sm_restored_gbff = (
            root
            / "data"
            / "SM"
            / "derived"
            / "restored"
            / "SM-NMZ-eggnog-restored.gbff"
        )
        sm_restored_manifest = (
            root
            / "data"
            / "SM"
            / "derived"
            / "restored"
            / "SM-NMZ-eggnog-restored.manifest.json"
        )
        sm_faa = root / "data" / "SM" / "bakta" / "SM-NMZ.faa"

        stats = enrich(
            bakta_path=sm_restored_gbff,
            faa_path=sm_faa,
            baktfold_path=root / "data" / "SM" / "baktfold" / "SM_NMZ.gbff",
            baktfold_invalid_ec_policy="skip",
            kofamscan_path=root / "data" / "SM" / "evidence" / "SM-KofamKOALA.txt",
            eggnog_path=root
            / "data"
            / "SM"
            / "evidence"
            / "SM-NMZ-query.emapper.annotations",
            interproscan_path=root
            / "data"
            / "SM"
            / "evidence"
            / "interproscan"
            / "SM-NMZ.interproscan.tsv",
            interproscan_version="5.59-91.0",
            translation_evidence_manifest=sm_restored_manifest,
            allow_imported_translations=True,
            output_path=output / "SM-four-source.gbff",
            manifest_path=output / "SM-four-source.manifest.json",
            context_report_path=output / "SM-four-source.context.json",
        )
        assert stats["self_check"] is True
        assert stats["insertions"] == sm_expected["four_source_insertions"]
        assert stats["output_sha256"] == sm_expected["four_source_output_sha256"]
        man_path = output / "SM-four-source.manifest.json"
        ctx_path = output / "SM-four-source.context.json"
        out_path = output / "SM-four-source.gbff"
        man_data = json.loads(man_path.read_text(encoding="utf-8"))
        _validate_json_schema(manifest_schema, man_data)
        ctx_data = json.loads(ctx_path.read_text(encoding="utf-8"))
        _validate_json_schema(combined_context_schema, ctx_data)
        _validate_json_schema(ips_context_schema, ctx_data["interproscan"])
        reloaded = load_translation_evidence(
            man_path, parse_genbank_bytes(out_path.read_bytes(), "SM four-source")
        )
        assert set(reloaded.keys()) == set(EXPECTED_SM_IMPORTED_ORIGINS)
        return {
            "output_sha256": stats["output_sha256"],
            "insertions": stats["insertions"],
        }

    _check(
        "sm-four-source-acceptance",
        sm_four_source_check,
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path, default=Path(".test-output/interproscan-gates")
    )
    parser.add_argument(
        "--baseline", type=Path, default=None, help="Path to versioned baseline JSON"
    )
    args = parser.parse_args()
    return run_checks(args.root.resolve(), args.output.resolve(), args.baseline)


if __name__ == "__main__":
    raise SystemExit(main())
