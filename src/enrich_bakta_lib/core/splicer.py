"""Qualifier and COMMENT byte formatting, splicing, reversibility, and semantic verify."""

from __future__ import annotations

import collections
import os
import re
import textwrap
from collections.abc import Iterable
from dataclasses import asdict
from typing import Any

from enrich_bakta_lib.core.types import (
    QUALIFIER_INDENT,
    AppliedInsertion,
    Insertion,
    MergeError,
    RawDocument,
    RawFeature,
    RawRecord,
    insertion_evidence_class,
    insertion_role,
    insertion_uid,
    sha256_bytes,
)


def newline_for_offset(data: bytes, offset: int) -> bytes:
    previous_start = data.rfind(b"\n", 0, max(0, offset - 1)) + 1
    previous = data[previous_start:offset]
    if previous.endswith(b"\r\n"):
        return b"\r\n"
    if previous.endswith(b"\n"):
        return b"\n"
    first = re.search(rb"\r\n|\n|\r", data)
    return first.group(0) if first else os.linesep.encode("ascii")


def format_qualifier(key: str, value: str, newline: bytes) -> bytes:
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
        raise MergeError(f"invalid qualifier name: {key!r}")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise MergeError(f"qualifier /{key} value contains a control character")
    escaped = value.replace('"', '""')
    payload = f'/{key}="{escaped}"'
    width = 80 - len(QUALIFIER_INDENT)
    chunks = [payload]
    if key != "inference":
        candidate_chunks = textwrap.wrap(
            payload,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
            replace_whitespace=False,
            drop_whitespace=True,
        ) or [payload]
        # A continuation line beginning with '/' is parsed as a new qualifier.
        # Keeping this uncommon value on one long line preserves its semantics.
        if not any(chunk.startswith("/") for chunk in candidate_chunks[1:]):
            chunks = candidate_chunks
    return b"".join(
        QUALIFIER_INDENT + chunk.encode("utf-8") + newline for chunk in chunks
    )


def qualifier_insertion(
    base_data: bytes,
    feature: RawFeature,
    key: str,
    value: str,
    source: str,
    source_value: str,
    order: int,
    *,
    candidate_role: str = "",
    evidence_class: str = "",
) -> Insertion:
    newline = newline_for_offset(base_data, feature.end)
    return Insertion(
        offset=feature.end,
        payload=format_qualifier(key, value, newline),
        record=feature.record_id,
        feature_type=feature.feature_type,
        locus_tag=feature.locus_tag or "",
        qualifier=key,
        value=value,
        source=source,
        source_value=source_value,
        order=order,
        candidate_role=candidate_role,
        evidence_class=evidence_class,
    )


def comment_insertion(
    base: RawDocument,
    record: RawRecord,
    marker: str,
    lines: list[str],
    source: str,
    order: int,
) -> Insertion | None:
    if record.features_offset is None:
        raise MergeError(f"record {record.record_id!r} has no FEATURES offset")
    record_prefix = base.data[record.start : record.features_offset]
    marker_bytes = marker.encode("utf-8")
    if marker_bytes in record_prefix:
        return None
    newline = newline_for_offset(base.data, record.features_offset)
    has_comment = re.search(rb"(?m)^COMMENT\b", record_prefix) is not None
    wrapped_lines = [
        piece
        for line in [marker, *lines]
        for piece in (
            textwrap.wrap(
                line,
                width=68,
                break_long_words=False,
                break_on_hyphens=False,
            )
            or [""]
        )
    ]
    rendered: list[bytes] = []
    for index, line in enumerate(wrapped_lines):
        if index == 0 and not has_comment:
            rendered.append(b"COMMENT     " + line.encode("utf-8") + newline)
        else:
            rendered.append(b"            " + line.encode("utf-8") + newline)
    return Insertion(
        offset=record.features_offset,
        payload=b"".join(rendered),
        record=record.record_id,
        feature_type="",
        locus_tag="",
        qualifier="COMMENT",
        value=marker,
        source=source,
        source_value=" | ".join(lines),
        order=order,
        candidate_role="context",
        evidence_class="context",
    )


def apply_insertions(
    base_data: bytes,
    insertions: list[Insertion],
    *,
    extra_allowed_qualifiers: Iterable[str] = (),
) -> tuple[bytes, list[AppliedInsertion]]:
    """Apply a fully computed allowlist and prove exact reversibility."""
    allowed_qualifiers = {
        "gene",
        "EC_number",
        "db_xref",
        "note",
        "inference",
        "COMMENT",
    }
    allowed_qualifiers.update(extra_allowed_qualifiers)
    ordered = sorted(insertions, key=lambda item: (item.offset, item.order))
    cursor = 0
    output = bytearray()
    applied: list[AppliedInsertion] = []
    for insertion in ordered:
        if insertion.qualifier not in allowed_qualifiers:
            raise MergeError(f"insertion allowlist rejected /{insertion.qualifier}")
        if not 0 <= insertion.offset <= len(base_data) or insertion.offset < cursor:
            raise MergeError(
                f"invalid or unsorted insertion offset: {insertion.offset}"
            )
        output.extend(base_data[cursor : insertion.offset])
        output_offset = len(output)
        output.extend(insertion.payload)
        applied.append(
            AppliedInsertion(insertion=insertion, output_offset=output_offset)
        )
        cursor = insertion.offset
    output.extend(base_data[cursor:])
    result = bytes(output)

    reconstructed = bytearray()
    reconstructed_cursor = 0
    for item in applied:
        start = item.output_offset
        end = start + len(item.insertion.payload)
        if result[start:end] != item.insertion.payload:
            raise MergeError(
                "byte-preservation audit could not locate a recorded insertion"
            )
        reconstructed.extend(result[reconstructed_cursor:start])
        reconstructed_cursor = end
    reconstructed.extend(result[reconstructed_cursor:])
    if bytes(reconstructed) != base_data or sha256_bytes(
        bytes(reconstructed)
    ) != sha256_bytes(base_data):
        raise MergeError(
            "byte-preservation audit failed after removing recorded insertions"
        )
    return result, applied


def verify_insertion_semantics(
    base: RawDocument, output: RawDocument, applied: list[AppliedInsertion]
) -> None:
    """Prove semantic values and targets as well as exact byte reversibility."""
    if len(base.records) != len(output.records):
        raise MergeError("semantic round-trip changed record count")
    targets = {feature.end: feature for feature in base.features}
    additions: dict[int, list[Insertion]] = collections.defaultdict(list)
    for item in applied:
        insertion = item.insertion
        if insertion.qualifier == "COMMENT":
            continue
        target = targets.get(insertion.offset)
        if target is None or (
            target.record_id,
            target.feature_type,
            target.locus_tag or "",
        ) != (insertion.record, insertion.feature_type, insertion.locus_tag):
            raise MergeError("insertion does not identify its exact target feature")
        additions[target.start].append(insertion)
    for before_record, after_record in zip(base.records, output.records, strict=True):
        if before_record.sequence != after_record.sequence or len(
            before_record.features
        ) != len(after_record.features):
            raise MergeError("semantic round-trip changed sequence or feature count")
        for before, after in zip(
            before_record.features, after_record.features, strict=True
        ):
            if (before.feature_type, before.location_key) != (
                after.feature_type,
                after.location_key,
            ):
                raise MergeError("semantic round-trip changed feature location/type")
            expected = {key: list(values) for key, values in before.qualifiers.items()}
            for insertion in additions.get(before.start, []):
                expected.setdefault(insertion.qualifier, []).append(insertion.value)
            if expected != after.qualifiers:
                raise MergeError(
                    f"semantic qualifier round-trip failed for {before.record_id}/{before.locus_tag}"
                )


def insertion_rows(applied: list[AppliedInsertion]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in applied:
        insertion = item.insertion
        row = asdict(insertion)
        row.pop("payload")
        row["candidate_role"] = insertion_role(insertion)
        row["evidence_class"] = insertion_evidence_class(insertion)
        row["insertion_id"] = insertion_uid(insertion)
        row["entry_type"] = "insertion"
        row["base_offset"] = row.pop("offset")
        row["output_offset"] = item.output_offset
        row["emitted_sha256"] = sha256_bytes(insertion.payload)
        rows.append(row)
    return rows
