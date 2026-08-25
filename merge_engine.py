#!/usr/bin/env python3
"""Shared byte-preserving GenBank annotation merge primitives."""

from __future__ import annotations

import collections
import hashlib
import io
import json
import os
import re
import tempfile
import textwrap
import warnings
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from Bio import BiopythonParserWarning, SeqIO

QUALIFIER_INDENT = b" " * 21
TOOL_VERSION = "0.2.0"
_FEATURE_RE = re.compile(rb"^ {5}(\S+)\s+(.+)$")
_QUALIFIER_RE = re.compile(rb"^ {21}/([^=\s]+)(?:=(.*))?$")
_LOCUS_LENGTH_RE = re.compile(rb"^LOCUS\s+\S+\s+(\d+)\s+bp\b")


class MergeError(ValueError):
    """Raised when an input cannot be merged safely."""


@dataclass
class RawFeature:
    record_index: int
    record_id: str
    feature_type: str
    location: str
    start: int
    end: int
    qualifiers: dict[str, list[str]] = field(default_factory=dict)

    def values(self, key: str) -> list[str]:
        return list(self.qualifiers.get(key, []))

    @property
    def locus_tags(self) -> list[str]:
        return [value for value in self.values("locus_tag") if value]

    @property
    def locus_tag(self) -> str | None:
        values = list(dict.fromkeys(self.locus_tags))
        return values[0] if len(values) == 1 else None


@dataclass
class RawRecord:
    index: int
    record_id: str
    accessions: tuple[str, ...]
    declared_length: int | None
    topology: str | None
    start: int
    end: int
    features_offset: int | None
    sequence: bytes
    features: list[RawFeature] = field(default_factory=list)

    @property
    def sequence_sha256(self) -> str:
        return sha256_bytes(self.sequence)


@dataclass
class RawDocument:
    data: bytes
    records: list[RawRecord]

    @property
    def features(self) -> list[RawFeature]:
        return [feature for record in self.records for feature in record.features]


@dataclass
class Insertion:
    offset: int
    payload: bytes
    record: str
    feature_type: str
    locus_tag: str
    qualifier: str
    value: str
    source: str
    source_value: str
    order: int = 0


@dataclass
class AppliedInsertion:
    insertion: Insertion
    output_offset: int


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _parse_qualifiers(block: bytes) -> dict[str, list[str]]:
    qualifiers: dict[str, list[str]] = collections.defaultdict(list)
    current_key: str | None = None
    pieces: list[str] = []

    def flush() -> None:
        nonlocal current_key, pieces
        if current_key is None:
            return
        raw = "".join(piece.strip() for piece in pieces)
        if len(raw) >= 2 and raw.startswith('"') and raw.endswith('"'):
            raw = raw[1:-1].replace('""', '"')
        qualifiers[current_key].append(raw)
        current_key = None
        pieces = []

    for raw_line in block.splitlines()[1:]:
        match = _QUALIFIER_RE.match(raw_line)
        if match:
            flush()
            current_key = _decode(match.group(1), "qualifier name")
            pieces = [_decode(match.group(2) or b"", f"/{current_key} value")]
        elif current_key is not None and raw_line.startswith(QUALIFIER_INDENT):
            pieces.append(
                _decode(
                    raw_line[len(QUALIFIER_INDENT) :], f"/{current_key} continuation"
                )
            )
    flush()
    return dict(qualifiers)


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
        block = data[current_feature["start"] : end]
        feature = RawFeature(
            record_index=current_record["index"],
            record_id=current_record["record_id"],
            feature_type=current_feature["feature_type"],
            location=current_feature["location"],
            start=current_feature["start"],
            end=end,
            qualifiers=_parse_qualifiers(block),
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
            }

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
    return RawDocument(data=data, records=records)


def validate_genbank_semantics(data: bytes, label: str = "GenBank input") -> int:
    """Require canonical Biopython parsing and return the record count."""
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
    return len(records)


def feature_key(feature: RawFeature) -> tuple[str, str, str] | None:
    if not feature.locus_tags:
        return None
    unique = list(dict.fromkeys(feature.locus_tags))
    if len(unique) != 1 or len(feature.locus_tags) != 1:
        raise MergeError(
            f"feature {feature.record_id}/{feature.feature_type} at byte {feature.start} "
            f"has ambiguous locus_tag values: {feature.locus_tags!r}"
        )
    return feature.record_id, feature.feature_type, unique[0]


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
                base_feature.location,
            )
            source_tuple = (
                source_feature.feature_type,
                tuple(source_feature.locus_tags),
                source_feature.location,
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
        if feature.feature_type != "CDS" or not feature.locus_tag:
            continue
        if feature.locus_tag in result:
            other = result[feature.locus_tag]
            raise MergeError(
                f"GBFF duplicate CDS locus_tag {feature.locus_tag!r} in "
                f"{other.record_id!r} and {feature.record_id!r}"
            )
        result[feature.locus_tag] = feature
    return result


def validate_faa_gbff(
    base: RawDocument,
    proteins: dict[str, str],
    query_ids: Iterable[str],
    *,
    source_name: str,
) -> tuple[dict[str, RawFeature], dict[str, Any]]:
    """Prove that every evidence query maps exactly to an identical GBFF CDS."""
    cds = cds_by_locus(base)
    unique_ids = set(query_ids)
    missing_faa = sorted(unique_ids - proteins.keys())
    missing_gbff = sorted(unique_ids - cds.keys())
    if missing_faa or missing_gbff:
        raise MergeError(
            f"{source_name} query IDs do not map exactly: "
            + json.dumps(
                {
                    "missing_from_faa": missing_faa[:25],
                    "missing_from_gbff": missing_gbff[:25],
                },
                indent=2,
            )
        )
    mismatches: list[str] = []
    missing_translations: list[str] = []
    for query_id in sorted(unique_ids):
        translations = cds[query_id].values("translation")
        if len(translations) != 1 or not normalize_protein(translations[0]):
            missing_translations.append(query_id)
        elif normalize_protein(translations[0]) != proteins[query_id]:
            mismatches.append(query_id)
    if missing_translations or mismatches:
        raise MergeError(
            f"{source_name} FAA/GBFF protein validation failed: "
            + json.dumps(
                {
                    "missing_or_ambiguous_translation": missing_translations[:25],
                    "sequence_mismatches": mismatches[:25],
                },
                indent=2,
            )
        )
    return cds, {
        "faa_records": len(proteins),
        "faa_gbff_sequence_matches": len(unique_ids),
        "missing_faa_query_ids": 0,
        "missing_gbff_query_ids": 0,
    }


def protein_sha256(sequence: str) -> str:
    return sha256_bytes(normalize_protein(sequence).encode("ascii"))


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
    if any(character in value for character in "\r\n\t"):
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


def insertion_rows(applied: list[AppliedInsertion]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in applied:
        insertion = item.insertion
        row = asdict(insertion)
        row.pop("payload")
        row["entry_type"] = "insertion"
        row["base_offset"] = row.pop("offset")
        row["output_offset"] = item.output_offset
        row["emitted_sha256"] = sha256_bytes(insertion.payload)
        rows.append(row)
    return rows


def reconcile_insertions(
    insertions: Iterable[Insertion], *, gene_conflict_policy: str = "skip"
) -> tuple[list[Insertion], list[dict[str, Any]], dict[str, Any]]:
    """Collapse exact overlaps and resolve paired gene-symbol disagreements."""
    if gene_conflict_policy not in {"skip", "prefer-eggnog", "prefer-baktfold"}:
        raise MergeError(f"invalid gene conflict policy {gene_conflict_policy!r}")
    ordered = sorted(
        insertions,
        key=lambda item: (
            item.record,
            item.locus_tag,
            item.feature_type,
            item.qualifier,
            item.value,
            item.order,
            item.source,
        ),
    )
    rows: list[dict[str, Any]] = []
    gene_values: dict[tuple[str, str], set[str]] = collections.defaultdict(set)
    for item in ordered:
        if item.qualifier == "gene":
            gene_values[(item.record, item.locus_tag)].add(item.value)
    conflicts = {key: values for key, values in gene_values.items() if len(values) > 1}
    preferred_source = {
        "prefer-eggnog": "eggnog",
        "prefer-baktfold": "baktfold",
    }.get(gene_conflict_policy)
    filtered: list[Insertion] = []
    for item in ordered:
        conflict_values = conflicts.get((item.record, item.locus_tag))
        if (
            item.qualifier == "gene"
            and conflict_values
            and (preferred_source is None or item.source.lower() != preferred_source)
        ):
            continue
        filtered.append(item)
    unique: dict[tuple[str, str, str, str, str], Insertion] = {}
    sources: dict[tuple[str, str, str, str, str], set[str]] = collections.defaultdict(
        set
    )
    for item in filtered:
        key = (
            item.record,
            item.feature_type,
            item.locus_tag,
            item.qualifier,
            item.value,
        )
        sources[key].add(item.source)
        unique.setdefault(key, item)
    for key, supported_by in sorted(sources.items()):
        if len(supported_by) > 1:
            rows.append(
                {
                    "entry_type": "reconciliation",
                    "record": key[0],
                    "feature_type": key[1],
                    "locus_tag": key[2],
                    "qualifier": key[3],
                    "value": key[4],
                    "status": "exact_duplicate_collapsed",
                    "supporting_sources": " | ".join(sorted(supported_by)),
                }
            )
    for (record, locus_tag), values in sorted(conflicts.items()):
        rows.append(
            {
                "entry_type": "reconciliation",
                "record": record,
                "locus_tag": locus_tag,
                "qualifier": "gene",
                "value": " | ".join(sorted(values)),
                "status": "gene_conflict_skipped"
                if preferred_source is None
                else gene_conflict_policy,
                "supporting_sources": " | ".join(
                    sorted(
                        {
                            item.source
                            for item in ordered
                            if item.record == record
                            and item.locus_tag == locus_tag
                            and item.qualifier == "gene"
                        }
                    )
                ),
            }
        )
    reconciled = sorted(unique.values(), key=lambda item: (item.offset, item.order))
    return (
        reconciled,
        rows,
        {
            "input_insertions": len(ordered),
            "reconciled_insertions": len(reconciled),
            "exact_duplicates_collapsed": len(filtered) - len(reconciled),
            "gene_conflicts": len(conflicts),
            "gene_conflict_policy": gene_conflict_policy,
        },
    )


def write_manifest(
    path: Path, metadata: dict[str, Any], rows: list[dict[str, Any]]
) -> None:
    if path.suffix.lower() == ".json":
        payload = (
            json.dumps(
                {
                    "schema": "enrich-bakta.merge-manifest.v1",
                    "metadata": metadata,
                    "entries": rows,
                },
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
            ).encode("utf-8")
            + b"\n"
        )
        atomic_write(path, payload)
        return
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
    output = [f"# {key}\t{metadata[key]}" for key in sorted(metadata)]
    output.append("\t".join(columns))
    for row in rows:
        output.append(
            "\t".join(
                str(row.get(column, "")).replace("\t", " ").replace("\n", " ")
                for column in columns
            )
        )
    atomic_write(path, ("\n".join(output) + "\n").encode("utf-8"))


def write_json_sidecar(path: Path, payload: dict[str, Any]) -> None:
    data = (
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False).encode(
            "utf-8"
        )
        + b"\n"
    )
    atomic_write(path, data)


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
    metadata = {
        **metadata,
        "tool_version": TOOL_VERSION,
        "base_sha256": sha256_bytes(base_data),
        "output_sha256": sha256_bytes(merged),
        "insertions": len(applied),
        "records": len(parsed_output.records),
    }
    if sidecar_payload is not None:
        sidecar_metadata = sidecar_payload.get("metadata")
        if isinstance(sidecar_metadata, dict):
            sidecar_metadata.update(
                {
                    "base_sha256": metadata["base_sha256"],
                    "output_sha256": metadata["output_sha256"],
                    "insertions": metadata["insertions"],
                    "records": metadata["records"],
                }
            )
    atomic_write(output_path, merged)
    if manifest_path is not None:
        rows = insertion_rows(applied)
        if evidence_rows is not None:
            rows = [*evidence_rows, *rows]
        write_manifest(manifest_path, metadata, rows)
    if sidecar_path is not None and sidecar_payload is not None:
        write_json_sidecar(sidecar_path, sidecar_payload)
    return {
        "output_sha256": metadata["output_sha256"],
        "base_sha256": metadata["base_sha256"],
        "insertions": len(applied),
        "self_check": True,
        "manifest": str(manifest_path) if manifest_path else None,
        "sidecar": str(sidecar_path) if sidecar_path else None,
    }
