#!/usr/bin/env python3
"""Run published C14/SM equivalence and first-pass enrichment checks."""

from __future__ import annotations

import argparse
import json
import time
import traceback
from collections.abc import Callable
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.merge_engine import (
    load_translation_evidence,
    parse_genbank_bytes,
)
from enrich_bakta_lib.sources.eggnog import merge as merge_eggnog
from enrich_bakta_lib.sources.eggnog import parse_eggnog_path
from enrich_bakta_lib.sources.kofam import merge as merge_kofam
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


def run(root: Path, output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    report = output / "real-checks.json"
    results: list[dict[str, Any]] = []
    normalized_fields = (
        "query_id",
        "seed_ortholog",
        "evalue",
        "score",
        "cog_category",
        "preferred_name",
        "gos",
        "ec",
        "ec_partial",
        "kegg_ko",
        "cazy",
        "pfams",
        "eggnog_ogs",
        "confidence",
    )

    for sample in ("C14", "SM"):
        folder = root / "data" / sample
        prefix = f"{sample}-NMZ"
        base = folder / "bakta" / f"{prefix}.gbff"
        faa = folder / "bakta" / f"{prefix}.faa"
        fold = folder / "baktfold" / f"{sample}_NMZ.gbff"
        ko = folder / "evidence" / f"{sample}-KofamKOALA.txt"
        egg = folder / "evidence" / f"{prefix}-query.emapper.annotations"
        xlsx = folder / "evidence" / f"{prefix}-emapper_annotations.xlsx"

        def table_check(egg: Path = egg, xlsx: Path = xlsx) -> dict[str, Any]:
            tsv = parse_eggnog_path(egg)
            excel = parse_eggnog_path(xlsx, expected_version="3.0.0-beta6")
            left = {
                hit.query_id: tuple(getattr(hit, name) for name in normalized_fields)
                for hit in tsv.hits
            }
            right = {
                hit.query_id: tuple(getattr(hit, name) for name in normalized_fields)
                for hit in excel.hits
            }
            if left != right:
                raise AssertionError("normalized TSV and XLSX hit fields differ")
            return {
                "rows": len(left),
                "schema": tsv.schema_id,
                "tsv_xlsx_equal": True,
            }

        _check(
            f"{sample}-tables",
            table_check,
            results=results,
            report=report,
            tracebacks=output,
        )
        combined = output / f"{sample}-bk.gbff"

        def combined_check(
            base: Path = base,
            faa: Path = faa,
            ko: Path = ko,
            fold: Path = fold,
            combined: Path = combined,
            sample: str = sample,
        ) -> dict[str, Any]:
            stats = merge_kofam(
                base,
                faa,
                ko,
                combined,
                baktfold_path=fold,
                manifest_path=output / f"{sample}-bk.json",
                baktfold_invalid_ec_policy="skip",
            )
            return {
                "output_sha256": stats["output_sha256"],
                "insertions": stats["insertions"],
                "hits": stats["kofam"]["hits"],
            }

        _check(
            f"{sample}-bk",
            combined_check,
            results=results,
            report=report,
            tracebacks=output,
        )

        def combined_rerun(
            faa: Path = faa,
            ko: Path = ko,
            fold: Path = fold,
            combined: Path = combined,
            sample: str = sample,
        ) -> dict[str, Any]:
            rerun = output / f"{sample}-bk-rerun.gbff"
            stats = merge_kofam(
                combined,
                faa,
                ko,
                rerun,
                baktfold_path=fold,
                manifest_path=output / f"{sample}-bk-rerun.json",
                baktfold_invalid_ec_policy="skip",
            )
            if stats["insertions"] != 0 or rerun.read_bytes() != combined.read_bytes():
                raise AssertionError(
                    "Baktfold/Kofam output-as-input rerun changed bytes"
                )
            return {"output_sha256": stats["output_sha256"], "insertions": 0}

        _check(
            f"{sample}-bk-rerun",
            combined_rerun,
            results=results,
            report=report,
            tracebacks=output,
        )
        restored = output / f"{sample}-restored.gbff"
        restoration_manifest = output / f"{sample}-restored.json"

        def restoration_check(
            base: Path = base,
            faa: Path = faa,
            egg: Path = egg,
            restored: Path = restored,
            restoration_manifest: Path = restoration_manifest,
        ) -> dict[str, Any]:
            stats = restore(
                base,
                faa,
                egg,
                restored,
                manifest_path=restoration_manifest,
                translation_policy="import-faa",
            )
            return {
                "count": stats["restored_translation_count"],
                "output_sha256": stats["output_sha256"],
            }

        _check(
            f"{sample}-restore",
            restoration_check,
            results=results,
            report=report,
            tracebacks=output,
        )
        enriched = output / f"{sample}-egg.gbff"
        egg_manifest = output / f"{sample}-egg.json"

        def eggnog_check(
            faa: Path = faa,
            egg: Path = egg,
            restored: Path = restored,
            restoration_manifest: Path = restoration_manifest,
            enriched: Path = enriched,
            egg_manifest: Path = egg_manifest,
        ) -> dict[str, Any]:
            stats = merge_eggnog(
                restored,
                faa,
                egg,
                enriched,
                manifest_path=egg_manifest,
                translation_evidence_manifest=restoration_manifest,
                allow_imported_translations=True,
            )
            return {
                "output_sha256": stats["output_sha256"],
                "insertions": stats["insertions"],
                "emitted_values": stats["eggnog"]["emitted_values"],
            }

        _check(
            f"{sample}-egg",
            eggnog_check,
            results=results,
            report=report,
            tracebacks=output,
        )

        def chain_check(
            enriched: Path = enriched, egg_manifest: Path = egg_manifest
        ) -> dict[str, Any]:
            ledger = load_translation_evidence(
                egg_manifest, parse_genbank_bytes(enriched.read_bytes())
            )
            return {"entries": len(ledger)}

        _check(
            f"{sample}-egg-chain",
            chain_check,
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
    args = parser.parse_args()
    return run(args.root.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
