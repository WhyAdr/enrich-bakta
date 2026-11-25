"""Raw and semantic GenBank/FAA parsing, feature indexing, and structural parity."""

from __future__ import annotations

import collections
import io
import json
import re
import warnings
from collections.abc import Iterable
from typing import Any

from Bio import BiopythonParserWarning, SeqIO
from Bio.SeqRecord import SeqRecord

from enrich_bakta_lib.core.types import (
    QUALIFIER_INDENT,
    MergeError,
    RawDocument,
    RawFeature,
    RawRecord,
    sha256_bytes,
)

_FEATURE_RE = re.compile(rb"^ {5}(\S+)\s+(.+)$")
_QUALIFIER_RE = re.compile(rb"^ {21}/([^=\s]+)(?:=(.*))?$")
_LOCUS_LENGTH_RE = re.compile(rb"^LOCUS\s+\S+\s+(\d+)\s+bp\b")


def _split_raw_lines(data: bytes) -> list[tuple[int, int, bytes, bytes]]:
    result: list[tuple[int, int, bytes, bytes]] = []
    offset = 0
    for line in data.splitlines(keepends=True):
        body = line.rstrip(b"\r\n")
        newline = line[len(body) :]
        result.append((offset, offset + len(line), body, newline))
        offset += len(line)
    if offset < len(data) or not result:
        result.append((offset, len(data), data[offset:], b""))
    return result


def _decode(value: bytes, label: str) -> str:
    try:
        return value.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise MergeError(f"{label} is not valid UTF-8 at byte {exc.start}") from exc


def parse_genbank_bytes(data: bytes, label: str = "GenBank input") -> RawDocument:
    """Locate records/features without altering the original bytes."""
    lines = _split_raw_lines(data)
    records: list[RawRecord] = []
    current_record: dict[str, Any] | None = None
    current_feature: dict[str, Any] | None = None
    in_origin = False

    def flush_feature(end: int) -> None:
        nonlocal current_feature
        if current_feature is None or current_record is None:
            return
        feature = RawFeature(
            record_index=current_record["index"],
            record_id=current_record["record_id"],
            feature_type=current_feature["feature_type"],
            location=current_feature["location"],
            start=current_feature["start"],
            end=end,
        )
        current_record["features"].append(feature)
        current_feature = None

    def flush_record(end: int) -> None:
        nonlocal current_record
        if current_record is None:
            return
        flush_feature(end)
        records.append(
            RawRecord(
                index=current_record["index"],
                record_id=current_record["record_id"],
                accessions=tuple(current_record["accessions"]),
                declared_length=current_record["declared_length"],
                topology=current_record["topology"],
                start=current_record["start"],
                end=end,
                features_offset=current_record["features_offset"],
                sequence=b"".join(current_record["sequence"]).upper(),
                features=current_record["features"],
            )
        )
        current_record = None

    for start, end, body, _newline in lines:
        if body.startswith(b"LOCUS"):
            flush_record(start)
            fields = body.split()
            if len(fields) < 2:
                raise MergeError(f"{label}: malformed LOCUS line at byte {start}")
            record_id = _decode(fields[1], "LOCUS identifier")
            length_match = _LOCUS_LENGTH_RE.match(body)
            topology = next(
                (
                    _decode(token, "topology")
                    for token in fields
                    if token in {b"linear", b"circular"}
                ),
                None,
            )
            current_record = {
                "index": len(records),
                "record_id": record_id,
                "accessions": [],
                "declared_length": int(length_match.group(1)) if length_match else None,
                "topology": topology,
                "start": start,
                "features_offset": None,
                "features": [],
                "sequence": [],
            }
            current_feature = None
            in_origin = False
            continue
        if current_record is None:
            continue
        if body.startswith(b"ACCESSION"):
            current_record["accessions"].extend(
                _decode(token, "ACCESSION") for token in body.split()[1:]
            )
            continue
        if body.startswith(b"FEATURES"):
            flush_feature(start)
            current_record["features_offset"] = start
            continue
        if body.startswith(b"ORIGIN"):
            flush_feature(start)
            in_origin = True
            continue
        if body.startswith(b"//"):
            flush_record(end)
            in_origin = False
            continue
        if in_origin:
            current_record["sequence"].append(re.sub(rb"[^A-Za-z]", b"", body))
            continue
        feature_match = _FEATURE_RE.match(body)
        if feature_match:
            flush_feature(start)
            current_feature = {
                "feature_type": _decode(feature_match.group(1), "feature type"),
                "location": _decode(feature_match.group(2).strip(), "feature location"),
                "start": start,
                "has_qualifier": False,
            }
            continue
        if current_feature is not None and body.startswith(QUALIFIER_INDENT):
            continuation = body[len(QUALIFIER_INDENT) :].strip()
            if continuation.startswith(b"/"):
                current_feature["has_qualifier"] = True
            elif not current_feature["has_qualifier"]:
                current_feature["location"] += _decode(
                    continuation, "feature location continuation"
                )

    flush_record(len(data))
    if not records:
        raise MergeError(f"{label}: no GenBank LOCUS records found")
    for record in records:
        if record.features_offset is None:
            raise MergeError(
                f"{label}: record {record.record_id!r} has no FEATURES table"
            )
        if (
            record.declared_length is not None
            and len(record.sequence) != record.declared_length
        ):
            raise MergeError(
                f"{label}: record {record.record_id!r} declares {record.declared_length} bp "
                f"but ORIGIN contains {len(record.sequence)} bp"
            )
    semantic_records = _semantic_records(data, label)
    if len(records) != len(semantic_records):
        raise MergeError(f"{label}: raw/semantic record counts differ")
    for raw_record, semantic_record in zip(records, semantic_records, strict=True):
        if len(raw_record.features) != len(semantic_record.features):
            raise MergeError(f"{label}: raw/semantic feature counts differ")
        for raw, semantic in zip(
            raw_record.features, semantic_record.features, strict=True
        ):
            if raw.feature_type != semantic.type or semantic.location is None:
                raise MergeError(
                    f"{label}: raw/semantic feature alignment failed at byte {raw.start}"
                )
            raw.semantic_location = semantic.location
            raw.qualifiers = {
                key: list(values) for key, values in semantic.qualifiers.items()
            }
    return RawDocument(data=data, records=records)


def _semantic_records(data: bytes, label: str) -> list[SeqRecord]:
    """Parse one authoritative semantic view, keeping raw bytes separately."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        # The splice engine never rewrites source bytes. Latin-1 provides a
        # lossless one-codepoint-per-byte semantic view for legacy comments.
        text = data.decode("latin-1")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            records = list(SeqIO.parse(io.StringIO(text), "genbank"))
    except BiopythonParserWarning as exc:
        raise MergeError(f"{label}: Biopython GenBank warning: {exc}") from exc
    except Exception as exc:  # Biopython exposes multiple parser exception types
        raise MergeError(f"{label}: Biopython GenBank parsing failed: {exc}") from exc
    if not records:
        raise MergeError(f"{label}: Biopython found no GenBank records")
    return records


def validate_genbank_semantics(data: bytes, label: str = "GenBank input") -> int:
    return len(_semantic_records(data, label))


def feature_key(feature: RawFeature) -> tuple[str, str, str] | None:
    values = feature.values("locus_tag")
    if not values:
        return None
    if len(values) != 1 or not values[0].strip():
        raise MergeError(
            f"feature {feature.record_id}/{feature.feature_type} at byte {feature.start} "
            f"has ambiguous locus_tag values: {feature.locus_tags!r}"
        )
    return feature.record_id, feature.feature_type, values[0]


def feature_uid(feature: RawFeature) -> str:
    """Return a stable identity for a semantically located feature."""
    payload = json.dumps(
        {
            "record_id": feature.record_id,
            "feature_type": feature.feature_type,
            "locus_tag": feature.locus_tag or "",
            "location": feature.location_key,
        },
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return f"feature:{sha256_bytes(payload)[:32]}"


def index_unique_features(
    features: Iterable[RawFeature], *, feature_type: str | None = None, label: str
) -> dict[tuple[str, str, str], RawFeature]:
    index: dict[tuple[str, str, str], RawFeature] = {}
    for feature in features:
        if feature_type is not None and feature.feature_type != feature_type:
            continue
        key = feature_key(feature)
        if key is None:
            continue
        if key in index:
            raise MergeError(f"{label}: duplicate feature key {key!r}")
        index[key] = feature
    return index


def strict_parity_check(base: RawDocument, source: RawDocument) -> dict[str, Any]:
    """Validate record-local identity before source annotations are consulted."""
    problems: list[str] = []
    diagnostics: dict[str, Any] = {"translation_mismatches": []}
    if len(base.records) != len(source.records):
        problems.append(
            f"record count differs: base={len(base.records)} source={len(source.records)}"
        )

    # Duplicate keyed features are unsafe even when both files contain them.
    try:
        index_unique_features(base.features, label="base")
        index_unique_features(source.features, label="source")
    except MergeError as exc:
        problems.append(str(exc))

    for record_index, (left, right) in enumerate(zip(base.records, source.records)):
        prefix = f"record {record_index}"
        for name, left_value, right_value in (
            ("ID", left.record_id, right.record_id),
            ("accessions", left.accessions, right.accessions),
            ("declared length", left.declared_length, right.declared_length),
            ("topology", left.topology, right.topology),
            ("sequence SHA-256", left.sequence_sha256, right.sequence_sha256),
        ):
            if left_value != right_value:
                problems.append(
                    f"{prefix} {name} differs: base={left_value!r} source={right_value!r}"
                )
        left_counts = collections.Counter(
            feature.feature_type for feature in left.features
        )
        right_counts = collections.Counter(
            feature.feature_type for feature in right.features
        )
        if left_counts != right_counts:
            problems.append(
                f"{prefix} feature-type counts differ: base={dict(left_counts)} source={dict(right_counts)}"
            )
        if len(left.features) != len(right.features):
            problems.append(
                f"{prefix} feature count differs: base={len(left.features)} source={len(right.features)}"
            )
        for feature_index, (base_feature, source_feature) in enumerate(
            zip(left.features, right.features)
        ):
            base_tuple = (
                base_feature.feature_type,
                tuple(base_feature.locus_tags),
                base_feature.location_key,
            )
            source_tuple = (
                source_feature.feature_type,
                tuple(source_feature.locus_tags),
                source_feature.location_key,
            )
            if base_tuple != source_tuple:
                problems.append(
                    f"{prefix} feature {feature_index} differs: "
                    f"base={base_tuple!r} source={source_tuple!r}"
                )
            if base_feature.feature_type == "CDS":
                left_translation = normalize_protein(
                    "".join(base_feature.values("translation"))
                )
                right_translation = normalize_protein(
                    "".join(source_feature.values("translation"))
                )
                if (
                    left_translation
                    and right_translation
                    and left_translation != right_translation
                ):
                    diagnostics["translation_mismatches"].append(
                        base_feature.locus_tag
                        or f"{base_feature.record_id}:{feature_index}"
                    )
            if len(problems) >= 25:
                break
        if len(problems) >= 25:
            break
    if problems:
        report = {"status": "failed", "problems": problems, **diagnostics}
        raise MergeError(
            "structural parity check failed:\n" + json.dumps(report, indent=2)
        )
    return {
        "status": "passed",
        "records": len(base.records),
        "features": len(base.features),
        **diagnostics,
    }


def normalize_protein(sequence: str) -> str:
    value = re.sub(r"\s+", "", sequence).upper()
    return value.removesuffix("*")


def parse_faa(data: bytes) -> dict[str, str]:
    """Parse a UTF-8 protein FAA without silently normalizing identifiers."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise MergeError(f"FAA is not UTF-8 at byte {exc.start}") from exc
    proteins: dict[str, str] = {}
    try:
        for record in SeqIO.parse(io.StringIO(text), "fasta"):
            if record.id in proteins:
                raise MergeError(f"FAA contains duplicate record ID {record.id!r}")
            sequence = normalize_protein(str(record.seq))
            if not sequence:
                raise MergeError(
                    f"FAA record {record.id!r} has an empty protein sequence"
                )
            proteins[record.id] = sequence
    except MergeError:
        raise
    except Exception as exc:
        raise MergeError(f"FAA parsing failed: {exc}") from exc
    if not proteins:
        raise MergeError("FAA contains no protein records")
    return proteins


def cds_by_locus(base: RawDocument) -> dict[str, RawFeature]:
    """Build a global CDS locus index, rejecting ambiguous record membership."""
    result: dict[str, RawFeature] = {}
    for feature in base.features:
        if feature.feature_type != "CDS":
            continue
        key = feature_key(feature)
        if key is None:
            continue
        locus_tag = key[2]
        if locus_tag in result:
            other = result[locus_tag]
            raise MergeError(
                f"GBFF duplicate CDS locus_tag {locus_tag!r} in "
                f"{other.record_id!r} and {feature.record_id!r}"
            )
        result[locus_tag] = feature
    return result


def protein_sha256(sequence: str) -> str:
    normalized = normalize_protein(sequence)
    try:
        encoded = normalized.encode("ascii")
    except UnicodeEncodeError as exc:
        raise MergeError(
            f"protein sequence contains non-ASCII characters: {normalized[:32]!r}"
        ) from exc
    return sha256_bytes(encoded)
