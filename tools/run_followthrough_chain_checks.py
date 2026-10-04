#!/usr/bin/env python3
"""Run published C14/SM restoration and unified output-as-input chains."""

from __future__ import annotations

import argparse
import json
import time
import traceback
from collections.abc import Callable
from pathlib import Path
from typing import Any

from followthrough_acceptance import (
    assert_candidate_manifest,
    assert_expected,
    assert_imported_evidence,
    load_baseline,
    validate_baseline_inputs,
)

from enrich_bakta_lib.core.merge_engine import (
    MergeError,
)
from enrich_bakta_lib.sources.eggnog import merge as merge_eggnog
from enrich_bakta_lib.workflows.enrich import enrich
from enrich_bakta_lib.workflows.restore_translations import restore


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
    except Exception as exc:  # pragma: no cover - exercised only by acceptance failures
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


def run(root: Path, output: Path, baseline_path: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    report = output / "chain-checks.json"
    results: list[dict[str, Any]] = []
    baseline = load_baseline(baseline_path)
    _check(
        "baseline-inputs",
        lambda: validate_baseline_inputs(root, baseline),
        results=results,
        report=report,
        tracebacks=output,
    )
    if results[-1]["status"] != "passed":
        return 1

    for sample in ("C14", "SM"):
        expected = baseline["samples"][sample]["expected"]
        folder = root / "data" / sample
        faa = folder / "bakta" / f"{sample}-NMZ.faa"
        egg = folder / "evidence" / f"{sample}-NMZ-query.emapper.annotations"
        fold = folder / "baktfold" / f"{sample}_NMZ.gbff"
        ko = folder / "evidence" / f"{sample}-KofamKOALA.txt"
        restored = output / f"{sample}-restored.gbff"
        restoration_manifest = output / f"{sample}-restored.json"

        def restoration_noop(
            faa: Path = faa,
            egg: Path = egg,
            restored: Path = restored,
            restoration_manifest: Path = restoration_manifest,
            sample: str = sample,
        ) -> dict[str, Any]:
            rerun = output / f"{sample}-restored-noop.gbff"
            ledger_path = output / f"{sample}-restored-noop.json"
            stats = restore(
                restored,
                faa,
                egg,
                rerun,
                manifest_path=ledger_path,
                translation_policy="validated-only",
                translation_evidence_manifest=restoration_manifest,
            )
            if stats["insertions"] != 0 or rerun.read_bytes() != restored.read_bytes():
                raise AssertionError("restoration output-as-input rerun changed bytes")
            evidence = assert_imported_evidence(
                ledger_path, rerun, expected["imported_query_ids"]
            )
            return {
                "insertions": 0,
                "ledger_entries": len(evidence["query_ids"]),
                **evidence,
            }

        _check(
            f"{sample}-restoration-noop",
            restoration_noop,
            results=results,
            report=report,
            tracebacks=output,
        )
        egg_output = output / f"{sample}-egg.gbff"
        egg_manifest = output / f"{sample}-egg.json"

        def eggnog_noop(
            faa: Path = faa,
            egg: Path = egg,
            egg_output: Path = egg_output,
            egg_manifest: Path = egg_manifest,
            sample: str = sample,
        ) -> dict[str, Any]:
            rerun = output / f"{sample}-egg-noop.gbff"
            ledger_path = output / f"{sample}-egg-noop.json"
            stats = merge_eggnog(
                egg_output,
                faa,
                egg,
                rerun,
                manifest_path=ledger_path,
                translation_evidence_manifest=egg_manifest,
                allow_imported_translations=True,
            )
            if (
                stats["insertions"] != 0
                or rerun.read_bytes() != egg_output.read_bytes()
            ):
                raise AssertionError("eggNOG output-as-input rerun changed bytes")
            evidence = assert_imported_evidence(
                ledger_path, rerun, expected["imported_query_ids"]
            )
            return {
                "insertions": 0,
                "ledger_entries": len(evidence["query_ids"]),
                "manifest_semantics": assert_candidate_manifest(ledger_path),
                **evidence,
            }

        _check(
            f"{sample}-egg-noop",
            eggnog_noop,
            results=results,
            report=report,
            tracebacks=output,
        )
        unified = output / f"{sample}-unified.gbff"
        unified_manifest = output / f"{sample}-unified.json"

        def unified_check(
            faa: Path = faa,
            egg: Path = egg,
            fold: Path = fold,
            ko: Path = ko,
            restored: Path = restored,
            restoration_manifest: Path = restoration_manifest,
            unified: Path = unified,
            unified_manifest: Path = unified_manifest,
            sample: str = sample,
        ) -> dict[str, Any]:
            stats = enrich(
                bakta_path=restored,
                faa_path=faa,
                eggnog_path=egg,
                baktfold_path=fold,
                kofamscan_path=ko,
                kofamscan_version="1.3.0",
                output_path=unified,
                manifest_path=unified_manifest,
                translation_evidence_manifest=restoration_manifest,
                allow_imported_translations=True,
                baktfold_invalid_ec_policy="skip",
            )
            evidence = assert_imported_evidence(
                unified_manifest, unified, expected["imported_query_ids"]
            )
            try:
                enrich(
                    bakta_path=unified,
                    faa_path=faa,
                    eggnog_path=egg,
                    output_path=output / f"{sample}-blocked.gbff",
                    translation_evidence_manifest=unified_manifest,
                )
            except MergeError as exc:
                if "allow-imported-translations" not in str(exc):
                    raise
            else:
                raise AssertionError("unified imported-protein opt-in was not enforced")
            result = {
                "output_sha256": stats["output_sha256"],
                "insertions": stats["insertions"],
                "ledger_entries": len(evidence["query_ids"]),
            }
            assert_expected(
                f"{sample} unified",
                {key: result[key] for key in ("output_sha256", "insertions")},
                expected["unified"],
            )
            result["manifest_semantics"] = assert_candidate_manifest(unified_manifest)
            result["translation_evidence"] = evidence
            return result

        _check(
            f"{sample}-unified",
            unified_check,
            results=results,
            report=report,
            tracebacks=output,
        )

        def unified_noop(
            faa: Path = faa,
            egg: Path = egg,
            fold: Path = fold,
            ko: Path = ko,
            unified: Path = unified,
            unified_manifest: Path = unified_manifest,
            sample: str = sample,
        ) -> dict[str, Any]:
            rerun = output / f"{sample}-unified-noop.gbff"
            ledger_path = output / f"{sample}-unified-noop.json"
            stats = enrich(
                bakta_path=unified,
                faa_path=faa,
                eggnog_path=egg,
                baktfold_path=fold,
                kofamscan_path=ko,
                kofamscan_version="1.3.0",
                output_path=rerun,
                manifest_path=ledger_path,
                translation_evidence_manifest=unified_manifest,
                allow_imported_translations=True,
                baktfold_invalid_ec_policy="skip",
            )
            if stats["insertions"] != 0 or rerun.read_bytes() != unified.read_bytes():
                raise AssertionError("unified output-as-input rerun changed bytes")
            assert_expected(
                f"{sample} unified rerun hash",
                stats["output_sha256"],
                expected["unified"]["output_sha256"],
            )
            evidence = assert_imported_evidence(
                ledger_path, rerun, expected["imported_query_ids"]
            )
            return {
                "insertions": 0,
                "ledger_entries": len(evidence["query_ids"]),
                "manifest_semantics": assert_candidate_manifest(ledger_path),
                **evidence,
            }

        _check(
            f"{sample}-unified-noop",
            unified_noop,
            results=results,
            report=report,
            tracebacks=output,
        )
    return 1 if any(row["status"] != "passed" for row in results) else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output", type=Path, default=Path(".test-output/followthrough")
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=Path("docs/data/followthrough-scientific-baseline.json"),
    )
    args = parser.parse_args()
    return run(args.root.resolve(), args.output.resolve(), args.baseline.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
