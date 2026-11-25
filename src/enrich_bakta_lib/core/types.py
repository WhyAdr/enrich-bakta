"""Core data models, constants, hashing, and insertion identity primitives."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from Bio.SeqFeature import CompoundLocation, SimpleLocation

QUALIFIER_INDENT = b" " * 21
TOOL_VERSION = "0.5.0"


class MergeError(ValueError):
    """Raised when an input cannot be merged safely."""


CANDIDATE_ROLES = frozenset(
    {
        "functional_proposal",
        "substantive_evidence",
        "producer_provenance",
        "context",
    }
)


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
    text_encoding: str = "utf-8"

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
    candidate_role: str = ""
    evidence_class: str = ""


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


class TranslationEvidenceLedger(dict[str, TranslationEvidence]):
    """Structured ledger with the restored query IDs it is expected to cover."""

    def __init__(
        self,
        values: Mapping[str, TranslationEvidence],
        *,
        required_query_ids: Iterable[str] = (),
        manifest_sha256: str = "",
    ) -> None:
        super().__init__(values)
        self.required_query_ids = frozenset(required_query_ids)
        self.manifest_sha256 = manifest_sha256


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
            "candidate_role": insertion_role(insertion),
            "evidence_class": insertion_evidence_class(insertion),
            "payload_sha256": sha256_bytes(insertion.payload),
        },
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return f"insertion:{sha256_bytes(payload)[:32]}"


def insertion_role(insertion: Insertion) -> str:
    """Return the explicit semantic role of an insertion."""
    if insertion.candidate_role:
        if insertion.candidate_role not in CANDIDATE_ROLES:
            raise MergeError(
                f"insertion has unsupported candidate role {insertion.candidate_role!r}"
            )
        return insertion.candidate_role
    if insertion.qualifier == "COMMENT":
        return "context"
    if insertion.source in {
        "Baktfold provenance",
        "KofamScan provenance",
        "eggNOG provenance",
    }:
        return "producer_provenance"
    if insertion.source == "KofamScan hit evidence":
        return "substantive_evidence"
    return "functional_proposal"


def insertion_evidence_class(insertion: Insertion) -> str:
    """Return the evidence class used for support-link reconciliation."""
    if insertion.evidence_class:
        return insertion.evidence_class
    if insertion.qualifier == "COMMENT":
        return "context"
    if insertion.source.startswith("KofamScan"):
        return "ko"
    if insertion.source.startswith("eggNOG"):
        return "functional_annotation"
    if insertion.source.startswith("Baktfold"):
        if insertion.source == "Baktfold provenance" and insertion.source_value in {
            "gene",
            "ec",
            "structural",
        }:
            return insertion.source_value
        return {
            "gene": "gene",
            "EC_number": "ec",
            "db_xref": "structural",
        }.get(insertion.qualifier, "functional_annotation")
    return "functional_annotation"
