"""Uniform merge context, normalized source plan, and adapter protocol."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

from enrich_bakta_lib.core.types import Insertion, RawDocument, TranslationEvidence


@dataclass(frozen=True)
class MergeContext:
    """Single-snapshot input bundle shared across every adapter in a merge."""

    bakta_path: Path
    base_data: bytes
    base_doc: RawDocument
    base_sha256: str
    faa_path: Path | None = None
    faa_data: bytes = b""
    faa_sha256: str = ""
    proteins: Mapping[str, str] | None = None
    translation_evidence_manifest: Path | None = None
    translation_evidence: Mapping[str, TranslationEvidence] | None = None
    allow_imported_translations: bool = False
    add_comment_note: bool = True
    add_feature_provenance: bool = True
    merge_timestamp: str | None = None
    clean_gene_suffix: bool = False


@dataclass(frozen=True)
class SourcePlan:
    """Normalized output returned by every source adapter."""

    source_id: str
    source_sha256: str
    insertions: list[Insertion]
    evidence_rows: list[dict[str, Any]]
    stats: dict[str, Any]
    manifest_metadata: dict[str, Any]
    context_report: dict[str, Any] | None = None
    context_key: str | None = None
    input_paths: tuple[Path, ...] = ()


class SourceAdapter(Protocol):
    """Protocol implemented by every evidence source adapter."""

    source_id: str

    def plan(self, ctx: MergeContext, *, starting_order: int = 0) -> SourcePlan:
        """Validate source inputs against ctx and return a deterministic SourcePlan."""
        ...
