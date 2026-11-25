"""Artifact collision checks, manifest/sidecar serialization, atomic staging, and finalize_merge."""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.genbank import (
    parse_genbank_bytes,
    validate_genbank_semantics,
)
from enrich_bakta_lib.core.lineage import (
    _bind_metadata_translation_evidence,
    _bind_translation_evidence,
    has_translation_evidence_marker,
)
from enrich_bakta_lib.core.splicer import (
    apply_insertions,
    insertion_rows,
    verify_insertion_semantics,
)
from enrich_bakta_lib.core.types import (
    TOOL_VERSION,
    Insertion,
    MergeError,
    sha256_bytes,
    sha256_file,
)


def paths_collide(output: Path, inputs: Iterable[Path]) -> bool:
    output_resolved = output.resolve()
    for input_path in inputs:
        if output_resolved == input_path.resolve():
            return True
        if output.exists() and input_path.exists():
            try:
                if os.path.samefile(output, input_path):
                    return True
            except OSError:
                pass
    return False


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except OSError:
        try:
            temporary.unlink(missing_ok=True)
        finally:
            raise


def manifest_bytes(
    path: Path, metadata: dict[str, Any], rows: list[dict[str, Any]]
) -> bytes:
    if path.suffix.lower() == ".json":
        decisions = [
            row for row in rows if row.get("entry_type") == "candidate_decision"
        ]
        entries = [row for row in rows if row.get("entry_type") != "candidate_decision"]
        payload = (
            json.dumps(
                {
                    "schema": "enrich-bakta.merge-manifest.v2",
                    "metadata": metadata,
                    "entries": entries,
                    "decisions": decisions,
                },
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
            ).encode("utf-8")
            + b"\n"
        )
        return payload
    preferred_columns = [
        "entry_type",
        "query_id",
        "record",
        "feature_type",
        "locus_tag",
        "ko",
        "threshold",
        "score",
        "e_value",
        "definition",
        "ko_already_present",
        "emitted_qualifiers",
        "qualifier",
        "value",
        "source",
        "source_value",
        "base_offset",
        "output_offset",
        "emitted_sha256",
        "order",
    ]
    present = {key for row in rows for key in row}
    columns = [column for column in preferred_columns if column in present]
    columns.extend(sorted(present - set(columns)))
    output = ["# schema\tenrich-bakta.merge-manifest.v2"]
    output.extend(
        f"# {key}\t{_manifest_cell(metadata[key])}" for key in sorted(metadata)
    )
    output.append("\t".join(columns))
    for row in rows:
        output.append(
            "\t".join(_manifest_cell(row.get(column, "")) for column in columns)
        )
    return ("\n".join(output) + "\n").encode("utf-8")


def _manifest_cell(value: object) -> str:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(
            value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        )
    return str(value).replace("\t", " ").replace("\n", " ")


def write_manifest(
    path: Path, metadata: dict[str, Any], rows: list[dict[str, Any]]
) -> None:
    atomic_write(path, manifest_bytes(path, metadata, rows))


def json_sidecar_bytes(payload: dict[str, Any]) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False).encode(
            "utf-8"
        )
        + b"\n"
    )


def write_json_sidecar(path: Path, payload: dict[str, Any]) -> None:
    atomic_write(path, json_sidecar_bytes(payload))


def preflight_artifact_paths(outputs: Iterable[Path], inputs: Iterable[Path]) -> None:
    paths, input_paths = list(outputs), list(inputs)
    for i, path in enumerate(paths):
        if paths_collide(path, [*input_paths, *paths[:i]]):
            raise MergeError(f"artifact destination aliases an input/output: {path}")
        if path.exists() and not path.is_file():
            raise MergeError(f"artifact destination is not a regular file: {path}")
        for parent in path.parents:
            if parent.exists() and not parent.is_dir():
                raise MergeError(f"artifact parent is not a directory: {parent}")
        resolved = path.resolve()
        for other in [*input_paths, *paths[:i]]:
            other_resolved = other.resolve()
            if resolved in other_resolved.parents or other_resolved in resolved.parents:
                raise MergeError(
                    f"ancestor/descendant artifact path relationship: {path} / {other}"
                )


def stage_artifacts(
    artifacts: list[tuple[Path, bytes]], inputs: Iterable[Path]
) -> None:
    """Stage every artifact before the first replacement; promotions are sequential."""
    input_paths = list(inputs)
    paths = [path for path, _ in artifacts]
    preflight_artifact_paths(paths, input_paths)
    staged: list[tuple[Path, Path]] = []
    try:
        for path, data in artifacts:
            path.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(
                prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
            )
            temporary = Path(name)
            staged.append((temporary, path))
            with os.fdopen(fd, "wb") as handle:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
        preflight_artifact_paths(paths, input_paths)
        for temporary, path in staged:
            os.replace(temporary, path)
    except OSError as exc:
        raise MergeError(
            "artifact staging/promotion failed; promotions are not a cross-file transaction: "
            f"{exc}"
        ) from exc
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)


def finalize_merge(
    *,
    base_path: Path,
    base_data: bytes,
    output_path: Path,
    other_inputs: Iterable[Path],
    insertions: list[Insertion],
    manifest_path: Path | None,
    metadata: dict[str, Any],
    evidence_rows: list[dict[str, Any]] | None = None,
    extra_allowed_qualifiers: Iterable[str] = (),
    sidecar_path: Path | None = None,
    sidecar_payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if has_translation_evidence_marker(base_data) and not metadata.get(
        "translation_evidence"
    ):
        raise MergeError(
            "restored input requires --translation-evidence-manifest and preserved translation evidence"
        )
    inputs = [base_path, *other_inputs]
    if paths_collide(output_path, inputs):
        raise MergeError("output path must differ from every input path")
    if manifest_path is not None and paths_collide(
        manifest_path, [*inputs, output_path, *([sidecar_path] if sidecar_path else [])]
    ):
        raise MergeError("manifest path must differ from inputs, output, and sidecar")
    if sidecar_path is None and sidecar_payload is not None:
        raise MergeError("sidecar payload requires a sidecar path")
    if sidecar_path is not None and sidecar_payload is None:
        raise MergeError("sidecar path requires a sidecar payload")
    if sidecar_path is not None and paths_collide(
        sidecar_path,
        [*inputs, output_path, *([manifest_path] if manifest_path else [])],
    ):
        raise MergeError("sidecar path must differ from inputs, output, and manifest")
    if sha256_file(base_path) != sha256_bytes(base_data):
        raise MergeError(
            "base file changed on disk since planning; re-run the merge before writing"
        )
    merged, applied = apply_insertions(
        base_data, insertions, extra_allowed_qualifiers=extra_allowed_qualifiers
    )
    validate_genbank_semantics(merged, "merged output")
    parsed_output = parse_genbank_bytes(merged, "merged output")
    verify_insertion_semantics(
        parse_genbank_bytes(base_data, "Bakta input"), parsed_output, applied
    )
    if evidence_rows is not None:
        candidate_rows = [
            row
            for row in evidence_rows
            if row.get("entry_type") == "candidate_decision"
        ]
        if candidate_rows:
            from enrich_bakta_lib.core.decisions import validate_candidate_ledger

            validate_candidate_ledger(
                candidate_rows,
                insertions,
                base=parse_genbank_bytes(base_data, "Bakta input"),
                output=parsed_output,
            )
    output_sha256 = sha256_bytes(merged)
    metadata = _bind_metadata_translation_evidence(
        {
            **metadata,
            "tool_version": TOOL_VERSION,
            "base_sha256": sha256_bytes(base_data),
            "output_sha256": output_sha256,
            "insertions": len(applied),
            "records": len(parsed_output.records),
        },
        output_sha256,
    )
    sidecar_bytes: bytes | None = None
    if sidecar_payload is not None:
        sidecar_metadata = sidecar_payload.get("metadata")
        if isinstance(sidecar_metadata, dict):
            sidecar_metadata.update(
                {
                    "base_sha256": metadata["base_sha256"],
                    "output_sha256": metadata["output_sha256"],
                    "insertions": metadata["insertions"],
                    "records": metadata["records"],
                    "tool_version": TOOL_VERSION,
                }
            )
        sidecar_bytes = json_sidecar_bytes(sidecar_payload)
        metadata["context_report_sha256"] = sha256_bytes(sidecar_bytes)
    artifacts = [(output_path, merged)]
    if manifest_path is not None:
        rows = insertion_rows(applied)
        if evidence_rows is not None:
            rows = [*_bind_translation_evidence(evidence_rows, output_sha256), *rows]
        represented = {
            row.get("query_id")
            for row in rows
            if row.get("entry_type") == "translation_restoration"
        }
        for query_id, evidence in sorted(
            metadata.get("translation_evidence", {}).items()
        ):
            if query_id not in represented:
                rows.append(
                    {
                        "entry_type": "translation_restoration",
                        "query_id": query_id,
                        "status": "carried_forward",
                        "translation_evidence": evidence,
                    }
                )
        artifacts.append((manifest_path, manifest_bytes(manifest_path, metadata, rows)))
    if sidecar_path is not None and sidecar_bytes is not None:
        artifacts.append((sidecar_path, sidecar_bytes))
    stage_artifacts(artifacts, inputs)
    return {
        "output_sha256": metadata["output_sha256"],
        "base_sha256": metadata["base_sha256"],
        "insertions": len(applied),
        "self_check": True,
        "manifest": str(manifest_path) if manifest_path else None,
        "sidecar": str(sidecar_path) if sidecar_path else None,
    }
