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
from typing import Any, Mapping

from Bio import BiopythonParserWarning, SeqIO
from Bio.SeqFeature import CompoundLocation, SimpleLocation
from Bio.SeqRecord import SeqRecord

QUALIFIER_INDENT = b" " * 21
TOOL_VERSION = "0.3.0"
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
    semantic_location: SimpleLocation | CompoundLocation | None = None

    @property
    def location_key(self) -> str:
        return (
            str(self.semantic_location)
            if self.semantic_location is not None
            else self.location
        )

    def values(self, key: str) -> list[str]:
        return list(self.qualifiers.get(key, []))

    @property
    def locus_tags(self) -> list[str]:
        return [value for value in self.values("locus_tag") if value]

    @property
    def locus_tag(self) -> str | None:
        values = self.values("locus_tag")
        return values[0] if len(values) == 1 and values[0] else None


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


@dataclass(frozen=True)
class TranslationEvidence:
    feature_uid: str
    protein_sha256: str
    origin: str
    artifact_sequence_match: bool
    genomic_validation_status: str
    validation_reason: str
    original_base_sha256: str
    faa_sha256: str
    restoration_policy: str | None
    producer_lineage: dict[str, object] = field(default_factory=dict)
    bound_output_sha256: str | None = None

    def as_dict(self) -> dict[str, object]:
        return {
            "feature_uid": self.feature_uid,
            "protein_sha256": self.protein_sha256,
            "origin": self.origin,
            "artifact_sequence_match": self.artifact_sequence_match,
            "genomic_validation_status": self.genomic_validation_status,
            "validation_reason": self.validation_reason,
            "original_base_sha256": self.original_base_sha256,
            "faa_sha256": self.faa_sha256,
            "restoration_policy": self.restoration_policy,
            "producer_lineage": self.producer_lineage,
            "bound_output_sha256": self.bound_output_sha256,
        }


def insertion_uid(insertion: Insertion) -> str:
    """Return a deterministic identity for one planned insertion."""
    payload = json.dumps(
        {
            "record": insertion.record,
            "feature_type": insertion.feature_type,
            "locus_tag": insertion.locus_tag,
            "qualifier": insertion.qualifier,
            "value": insertion.value,
            "source": insertion.source,
            "source_value": insertion.source_value,
            "order": insertion.order,
            "payload_sha256": sha256_bytes(insertion.payload),
        },
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return f"insertion:{sha256_bytes(payload)[:32]}"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_input_bytes(path: Path, label: str) -> bytes:
    """Read an input artifact with a stable public-workflow error type."""
    try:
        return path.read_bytes()
    except OSError as exc:
        raise MergeError(f"cannot read {label} input {path}: {exc}") from exc


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


def validate_faa_gbff(
    base: RawDocument,
    proteins: dict[str, str],
    query_ids: Iterable[str],
    *,
    source_name: str,
    translation_evidence: Mapping[str, TranslationEvidence] | None = None,
    allow_imported_translations: bool = False,
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
    if translation_evidence is not None:
        base_hash = sha256_bytes(base.data)
        required_query_ids: Iterable[str] = getattr(
            translation_evidence, "required_query_ids", frozenset()
        )
        for query_id in sorted(unique_ids):
            evidence = translation_evidence.get(query_id)
            if evidence is None:
                # Restoration manifests are intentionally partial: they bind
                # imported/restored targets, while untouched translated CDSs
                # remain covered by the FAA/GBFF equality check above.
                if query_id in required_query_ids:
                    raise MergeError(
                        f"{source_name} translation evidence ledger is missing {query_id!r}"
                    )
                continue
            expected_hash = (
                evidence.bound_output_sha256 or evidence.original_base_sha256
            )
            if expected_hash != base_hash:
                raise MergeError(
                    f"{source_name} translation evidence ledger is bound to a different base"
                )
            if evidence.feature_uid != feature_uid(cds[query_id]):
                raise MergeError(
                    f"{source_name} translation evidence ledger disagrees for {query_id!r}"
                )
            if evidence.protein_sha256 != protein_sha256(proteins[query_id]):
                raise MergeError(
                    f"{source_name} translation evidence ledger disagrees for {query_id!r}"
                )
            if evidence.origin == "imported_faa" and not allow_imported_translations:
                raise MergeError(
                    f"{source_name} uses imported translation evidence for {query_id!r}; "
                    "pass --allow-imported-translations explicitly"
                )
    return cds, {
        "faa_records": len(proteins),
        "faa_gbff_sequence_matches": len(unique_ids),
        "missing_faa_query_ids": 0,
        "missing_gbff_query_ids": 0,
    }


TRANSLATION_EVIDENCE_MARKER = "##enrich-bakta:translation-evidence:v1##"


class TranslationEvidenceLedger(dict[str, TranslationEvidence]):
    """Structured ledger with the restored query IDs it is expected to cover."""

    def __init__(
        self,
        values: Mapping[str, TranslationEvidence],
        *,
        required_query_ids: Iterable[str] = (),
    ) -> None:
        super().__init__(values)
        self.required_query_ids = frozenset(required_query_ids)


def has_translation_evidence_marker(data: bytes) -> bool:
    return TRANSLATION_EVIDENCE_MARKER.encode("ascii") in data


def load_translation_evidence(
    path: Path, base: RawDocument
) -> dict[str, TranslationEvidence]:
    """Load and bind a JSON restoration manifest to the current GBFF bytes."""
    if path.suffix.lower() != ".json":
        raise MergeError(
            "translation evidence must use a JSON merge manifest; TSV manifests "
            "cannot preserve the structured evidence ledger"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise MergeError(
            f"cannot read translation evidence manifest {path}: {exc}"
        ) from exc
    if not isinstance(payload, dict) or payload.get("schema") not in {
        "enrich-bakta.merge-manifest.v1",
        "enrich-bakta.merge-manifest.v2",
    }:
        raise MergeError(f"unsupported translation evidence manifest schema: {path}")
    metadata = payload.get("metadata")
    entries = payload.get("entries")
    if not isinstance(metadata, dict) or not isinstance(entries, list):
        raise MergeError(
            "translation evidence manifest has invalid metadata or entries"
        )
    current_hash = sha256_bytes(base.data)
    if metadata.get("output_sha256") != current_hash:
        raise MergeError(
            "translation evidence manifest does not describe the current GBFF bytes"
        )
    cds = cds_by_locus(base)
    result: dict[str, TranslationEvidence] = {}
    for entry in entries:
        if (
            not isinstance(entry, dict)
            or entry.get("entry_type") != "translation_restoration"
        ):
            continue
        query_id = entry.get("query_id")
        raw = entry.get("translation_evidence")
        if not isinstance(query_id, str) or not isinstance(raw, dict):
            raise MergeError(
                "translation restoration manifest entry lacks structured evidence"
            )
        if query_id in result:
            raise MergeError(f"translation evidence manifest repeats {query_id!r}")
        if query_id not in cds:
            raise MergeError(
                f"translation evidence manifest references unknown CDS {query_id!r}"
            )
        required = (
            "feature_uid",
            "protein_sha256",
            "origin",
            "artifact_sequence_match",
            "genomic_validation_status",
            "validation_reason",
            "original_base_sha256",
            "faa_sha256",
            "restoration_policy",
            "bound_output_sha256",
        )
        if any(key not in raw for key in required):
            raise MergeError(f"translation evidence for {query_id!r} is incomplete")
        try:
            evidence = TranslationEvidence(
                feature_uid=str(raw["feature_uid"]),
                protein_sha256=str(raw["protein_sha256"]),
                origin=str(raw["origin"]),
                artifact_sequence_match=bool(raw["artifact_sequence_match"]),
                genomic_validation_status=str(raw["genomic_validation_status"]),
                validation_reason=str(raw["validation_reason"]),
                original_base_sha256=str(raw["original_base_sha256"]),
                faa_sha256=str(raw["faa_sha256"]),
                restoration_policy=(
                    None
                    if raw["restoration_policy"] is None
                    else str(raw["restoration_policy"])
                ),
                producer_lineage=(
                    raw["producer_lineage"]
                    if isinstance(raw["producer_lineage"], dict)
                    else {}
                ),
                bound_output_sha256=str(raw["bound_output_sha256"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise MergeError(
                f"translation evidence for {query_id!r} has invalid fields"
            ) from exc
        if evidence.bound_output_sha256 != current_hash:
            raise MergeError(
                f"translation evidence for {query_id!r} is bound to another output"
            )
        if evidence.feature_uid != feature_uid(cds[query_id]):
            raise MergeError(
                f"translation evidence for {query_id!r} has the wrong feature identity"
            )
        result[query_id] = evidence
    if not result:
        raise MergeError(
            "translation evidence manifest contains no restoration entries"
        )
    metadata_evidence = metadata.get("translation_evidence")
    required_query_ids: set[str] = set()
    if metadata_evidence is not None:
        if not isinstance(metadata_evidence, dict) or any(
            not isinstance(query_id, str) or not isinstance(value, dict)
            for query_id, value in metadata_evidence.items()
        ):
            raise MergeError(
                "translation evidence manifest has invalid metadata ledger"
            )
        required_query_ids = set(metadata_evidence)
        if required_query_ids != set(result):
            raise MergeError(
                "translation evidence manifest metadata and entries disagree"
            )
    return TranslationEvidenceLedger(result, required_query_ids=required_query_ids)


def protein_sha256(sequence: str) -> str:
    normalized = normalize_protein(sequence)
    try:
        encoded = normalized.encode("ascii")
    except UnicodeEncodeError as exc:
        raise MergeError(
            f"protein sequence contains non-ASCII characters: {normalized[:32]!r}"
        ) from exc
    return sha256_bytes(encoded)


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
        row["insertion_id"] = insertion_uid(insertion)
        row["entry_type"] = "insertion"
        row["base_offset"] = row.pop("offset")
        row["output_offset"] = item.output_offset
        row["emitted_sha256"] = sha256_bytes(insertion.payload)
        rows.append(row)
    return rows


def _is_provenance_source(source: str) -> bool:
    return "provenance" in source.lower()


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

    def source_family(source: str) -> str:
        return source.lower().removesuffix(" provenance")

    surviving_functional = {
        (
            item.record,
            item.locus_tag,
            source_family(item.source),
        )
        for item in filtered
        if item.qualifier != "COMMENT" and not _is_provenance_source(item.source)
    }
    provenance_filtered: list[Insertion] = []
    for item in filtered:
        if (
            item.qualifier != "COMMENT"
            and _is_provenance_source(item.source)
            and (
                item.record,
                item.locus_tag,
                source_family(item.source),
            )
            not in surviving_functional
        ):
            rows.append(
                {
                    "entry_type": "reconciliation",
                    "record": item.record,
                    "feature_type": item.feature_type,
                    "locus_tag": item.locus_tag,
                    "qualifier": item.qualifier,
                    "value": item.value,
                    "source": item.source,
                    "status": "provenance_without_surviving_annotation_suppressed",
                    "supporting_sources": item.source,
                }
            )
            continue
        provenance_filtered.append(item)
    filtered = provenance_filtered
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


def _bind_translation_evidence(
    rows: Iterable[dict[str, Any]], output_sha256: str
) -> list[dict[str, Any]]:
    bound_rows: list[dict[str, Any]] = []
    for row in rows:
        copied = dict(row)
        evidence = copied.get("translation_evidence")
        if isinstance(evidence, dict):
            evidence = dict(evidence)
            evidence["bound_output_sha256"] = output_sha256
            copied["translation_evidence"] = evidence
        bound_rows.append(copied)
    return bound_rows


def _bind_metadata_translation_evidence(
    metadata: dict[str, Any], output_sha256: str
) -> dict[str, Any]:
    copied = dict(metadata)
    evidence = copied.get("translation_evidence")
    if isinstance(evidence, dict):
        copied["translation_evidence"] = {
            query_id: (
                {
                    **value,
                    "bound_output_sha256": output_sha256,
                }
                if isinstance(value, dict)
                else value
            )
            for query_id, value in evidence.items()
        }
    return copied


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
    artifacts = [(output_path, merged)]
    if manifest_path is not None:
        rows = insertion_rows(applied)
        if evidence_rows is not None:
            rows = [*_bind_translation_evidence(evidence_rows, output_sha256), *rows]
        artifacts.append((manifest_path, manifest_bytes(manifest_path, metadata, rows)))
    if sidecar_path is not None and sidecar_payload is not None:
        artifacts.append((sidecar_path, json_sidecar_bytes(sidecar_payload)))
    stage_artifacts(artifacts, inputs)
    return {
        "output_sha256": metadata["output_sha256"],
        "base_sha256": metadata["base_sha256"],
        "insertions": len(applied),
        "self_check": True,
        "manifest": str(manifest_path) if manifest_path else None,
        "sidecar": str(sidecar_path) if sidecar_path else None,
    }
