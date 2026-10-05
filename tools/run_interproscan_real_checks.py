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

from enrich_bakta_lib.core.merge_engine import (
    MergeError,
    sha256_file,
)
from enrich_bakta_lib.workflows.enrich import enrich

EXPECTED_C14 = {
    "output_sha256": "472b0bd60d3990dda6307c89de7153d2ab3fbfd9f9324180e8889cb1b47bc9e1",
    "insertions": 34080,
    "candidate_count": 37942,
    "emitted_candidate_count": 34079,
    "accepted_candidate_count": 37942,
    "context_report_sha256": "2113c50f83e35df548f6cbc0087750d925b30153edfb7aa04b64b093e4a9adc9",
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


def run_checks(root: Path, output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    report = output / "interproscan-real-checks.json"
    results: list[dict[str, Any]] = []

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
        assert stats["output_sha256"] == EXPECTED_C14["output_sha256"]
        assert stats["insertions"] == EXPECTED_C14["insertions"]

        man_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        for key in (
            "candidate_count",
            "emitted_candidate_count",
            "accepted_candidate_count",
            "context_report_sha256",
        ):
            if man_data["metadata"][key] != EXPECTED_C14[key]:
                raise AssertionError(
                    f"C14 {key} mismatch: {man_data['metadata'][key]} != {EXPECTED_C14[key]}"
                )

        ctx_sha = sha256_file(context_path)
        assert ctx_sha == EXPECTED_C14["context_report_sha256"]
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

    # 2. SM pristine rejection
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

    # 3. SM restored acceptance
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
        man_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        te = man_data["metadata"].get("translation_evidence", {})
        for query in EXPECTED_SM_IMPORTED_ORIGINS:
            if query not in te:
                raise AssertionError(
                    f"Expected imported translation {query} in manifest metadata"
                )
            if te[query].get("origin") != "imported_faa":
                raise AssertionError(f"Expected origin imported_faa for {query}")

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

    # 4. BK71A lineage rejection
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
        c14_faa = (
            root / "data" / "C14" / "bakta" / "C14-NMZ.faa"
        )  # any FAA to reach base preflight
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

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path, default=Path(".test-output/interproscan-gates")
    )
    args = parser.parse_args()
    return run_checks(args.root.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
