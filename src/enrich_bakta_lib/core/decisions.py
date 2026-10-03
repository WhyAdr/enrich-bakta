"""Typed, deterministic candidate decisions for merge manifests."""

from __future__ import annotations

import json
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from typing import Any

from enrich_bakta_lib.core.merge_engine import (
    Insertion,
    MergeError,
    RawDocument,
    cds_by_locus,
    feature_uid,
    insertion_uid,
    sha256_bytes,
)


@dataclass
class CandidateDecision:
    candidate_id: str
    source_id: str
    source_sha256: str
    target_feature_uids: tuple[str, ...]
    field: str
    qualifier: str
    raw_value: str
    normalized_value: str
    planned_status: str
    final_status: str | None = None
    insertion_ids: tuple[str, ...] = ()
    supporting_candidate_ids: tuple[str, ...] = ()
    confidence_diagnostics: dict[str, Any] = dataclass_field(default_factory=dict)
    reason: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "entry_type": "candidate_decision",
            "candidate_id": self.candidate_id,
            "source_id": self.source_id,
            "source_sha256": self.source_sha256,
            "target_feature_uids": list(self.target_feature_uids),
            "field": self.field,
            "qualifier": self.qualifier,
            "raw_value": self.raw_value,
            "normalized_value": self.normalized_value,
            "planned_status": self.planned_status,
            "final_status": self.final_status,
            "insertion_ids": list(self.insertion_ids),
            "supporting_candidate_ids": list(self.supporting_candidate_ids),
            "confidence_diagnostics": self.confidence_diagnostics,
            "reason": self.reason,
        }


def stable_id(prefix: str, payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return f"{prefix}:{sha256_bytes(encoded)[:32]}"


def _source_id(value: str) -> str:
    lowered = value.lower()
    if lowered.startswith("baktfold"):
        return "Baktfold"
    if lowered.startswith("kofamscan") or lowered == "kofam":
        return "KofamScan"
    if lowered.startswith("eggnog"):
        return "eggNOG"
    return value


def _source_hash(source_id: str, source_hashes: Mapping[str, str]) -> str:
    return source_hashes.get(source_id) or source_hashes.get(source_id.lower(), "")


def _feature_index(base: RawDocument) -> dict[tuple[str, str, str], str]:
    return {
        (feature.record_id, feature.feature_type, feature.locus_tag or ""): feature_uid(
            feature
        )
        for feature in base.features
        if feature.locus_tag
    }


def _insertion_feature_uid(
    insertion: Insertion, index: Mapping[tuple[str, str, str], str]
) -> str | None:
    return index.get((insertion.record, insertion.feature_type, insertion.locus_tag))


def _semantic_key(insertion: Insertion) -> tuple[str, str, str, str, str]:
    return (
        insertion.record,
        insertion.feature_type,
        insertion.locus_tag,
        insertion.qualifier,
        insertion.value,
    )


def _is_provenance_insertion(insertion: Insertion) -> bool:
    return "provenance" in insertion.source.lower()


def _row_source(row: Mapping[str, Any]) -> str | None:
    entry_type = str(row.get("entry_type", ""))
    if entry_type.startswith("eggnog"):
        return "eggNOG"
    if entry_type == "kofam_hit":
        return "KofamScan"
    source = row.get("source")
    return _source_id(str(source)) if source else None


def _row_field_and_value(row: Mapping[str, Any]) -> tuple[str, str, str]:
    if row.get("entry_type") == "kofam_hit":
        ko = str(row.get("ko", ""))
        return "KO", "db_xref", f"KEGG:{ko}"
    if row.get("entry_type") == "baktfold_invalid":
        return "EC", "EC_number", str(row.get("normalized_value", ""))
    field = str(row.get("field", row.get("qualifier", "")))
    qualifier = {
        "Preferred_name": "gene",
        "GOs": "db_xref",
        "EC": "EC_number",
        "KEGG_ko": "db_xref",
        "CAZy": "db_xref",
        "PFAMs": "note",
        "eggNOG_OGs": "note",
    }.get(field, str(row.get("qualifier", field)))
    value = str(row.get("normalized_value", row.get("value", "")))
    return field, qualifier, value


def build_candidate_ledger(
    base: RawDocument,
    planned_insertions: Iterable[Insertion],
    final_insertions: Iterable[Insertion],
    evidence_rows: Iterable[dict[str, Any]],
    *,
    source_hashes: Mapping[str, str],
    gene_conflict_policy: str = "skip",
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Reconstruct typed planned/final decisions from all adapter plans."""
    planned = [item for item in planned_insertions if item.qualifier != "COMMENT"]
    final = [item for item in final_insertions if item.qualifier != "COMMENT"]
    index = _feature_index(base)
    final_by_semantic: dict[tuple[str, str, str, str, str], list[str]] = defaultdict(
        list
    )
    for item in final:
        final_by_semantic[_semantic_key(item)].append(insertion_uid(item))
    cds_index = cds_by_locus(base)
    valid_uids = set(index.values())
    planned_by_value: dict[tuple[str, str, str], list[Insertion]] = defaultdict(list)
    for item in planned:
        planned_by_value[(_source_id(item.source), item.qualifier, item.value)].append(
            item
        )
    covered: set[int] = set()
    candidates: list[CandidateDecision] = []

    def matching_insertions(
        source_id: str,
        qualifier: str,
        value: str,
        raw_value: str,
        target_uids: set[str],
    ) -> list[Insertion]:
        matches = []
        for item in planned_by_value.get((source_id, qualifier, value), []):
            item_uid = _insertion_feature_uid(item, index)
            if item_uid not in target_uids and target_uids:
                continue
            if raw_value and item.source_value != raw_value:
                continue
            matches.append(item)
        return matches

    for row in evidence_rows:
        entry_type = str(row.get("entry_type", ""))
        if entry_type not in {
            "baktfold_invalid",
            "eggnog_candidate",
            "eggnog_pair_conflict",
            "kofam_hit",
        }:
            continue
        source_id = _row_source(row)
        if source_id is None:
            continue
        source_sha = _source_hash(source_id, source_hashes)
        field_name, qualifier, normalized = _row_field_and_value(row)
        raw_value = str(row.get("raw_value", row.get("ko", "")))
        target_uids: set[str] = set()
        query_id = str(row.get("query_id", ""))
        if query_id in cds_index:
            target_uid = feature_uid(cds_index[query_id])
            target_uids.add(target_uid)
        target_uids.update(
            value
            for value in row.get("target_feature_uids", [])
            if isinstance(value, str) and value in valid_uids
        )
        if not target_uids and row.get("record") is not None:
            record = str(row.get("record", ""))
            feature_type = str(row.get("feature_type", ""))
            locus_tag = str(row.get("locus_tag", ""))
            target_location = str(row.get("target_location", ""))
            for feature in base.features:
                if (
                    feature.record_id == record
                    and feature.feature_type == feature_type
                    and (feature.locus_tag or "") == locus_tag
                    and (not target_location or feature.location_key == target_location)
                ):
                    target_uids.add(feature_uid(feature))
        if not target_uids:
            target_uids.update(
                item_uid
                for item in planned_by_value.get((source_id, qualifier, normalized), [])
                if (item_uid := _insertion_feature_uid(item, index)) is not None
            )
        matches = (
            []
            if entry_type in {"baktfold_invalid", "eggnog_pair_conflict"}
            else matching_insertions(
                source_id, qualifier, normalized, raw_value, target_uids
            )
        )
        covered.update(id(item) for item in matches)
        if not target_uids:
            continue
        row_status = str(row.get("status", "planned"))
        emitted = bool(row.get("emitted_qualifiers"))
        planned_status = row_status or "planned"
        final_ids = sorted(
            {
                insertion_id
                for item in matches
                for insertion_id in final_by_semantic.get(_semantic_key(item), [])
            }
        )
        if emitted:
            final_status = "emitted" if final_ids else "suppressed_conflict"
        elif row_status in {"existing", "existing_gene"}:
            final_status = "supported_existing"
        else:
            final_status = row_status or "not_applicable"
        payload = {
            "entry_type": entry_type,
            "source_sha256": source_sha,
            "row_number": row.get("row_number", row.get("kofam_row", 0)),
            "source_role": row.get("source_role", ""),
            "field": field_name,
            "raw_value": raw_value,
            "normalized_value": normalized,
            "target_feature_uids": sorted(target_uids),
        }
        candidates.append(
            CandidateDecision(
                candidate_id=stable_id("candidate", payload),
                source_id=source_id,
                source_sha256=source_sha,
                target_feature_uids=tuple(sorted(target_uids)),
                field=field_name,
                qualifier=qualifier,
                raw_value=raw_value,
                normalized_value=normalized,
                planned_status=planned_status,
                final_status=final_status,
                insertion_ids=tuple(final_ids),
                confidence_diagnostics={
                    key: row[key]
                    for key in (
                        "confidence_code",
                        "confidence_status",
                        "score",
                        "e_value",
                    )
                    if key in row
                },
                reason=str(row.get("reason", "")),
            )
        )

    grouped: dict[tuple[str, str, str, str, str, str], list[Insertion]] = defaultdict(
        list
    )
    for item in planned:
        if id(item) in covered:
            continue
        source_id = _source_id(item.source)
        grouped[
            (
                source_id,
                _source_hash(source_id, source_hashes),
                item.record,
                item.qualifier,
                item.value,
                item.source_value,
            )
        ].append(item)
    for (
        source_id,
        source_sha,
        record,
        qualifier,
        value,
        raw_value,
    ), items in sorted(grouped.items()):
        group_target_uids = sorted(
            item_uid
            for item in items
            if (item_uid := _insertion_feature_uid(item, index)) is not None
        )
        if not group_target_uids:
            continue
        final_ids = sorted(
            {
                insertion_id
                for item in items
                for insertion_id in final_by_semantic.get(_semantic_key(item), [])
            }
        )
        payload = {
            "source_sha256": source_sha,
            "record": record,
            "field": qualifier,
            "raw_value": raw_value,
            "normalized_value": value,
            "target_feature_uids": group_target_uids,
        }
        candidates.append(
            CandidateDecision(
                candidate_id=stable_id("candidate", payload),
                source_id=source_id,
                source_sha256=source_sha,
                target_feature_uids=tuple(group_target_uids),
                field=qualifier,
                qualifier=qualifier,
                raw_value=raw_value,
                normalized_value=value,
                planned_status="planned",
                final_status=(
                    "emitted"
                    if final_ids
                    else (
                        "suppressed_conflict"
                        if qualifier == "gene" and gene_conflict_policy == "skip"
                        else "not_applicable"
                    )
                ),
                insertion_ids=tuple(final_ids),
                confidence_diagnostics=(
                    {"provenance": True}
                    if any(_is_provenance_insertion(item) for item in items)
                    else {}
                ),
            )
        )

    by_semantic: dict[tuple[str, str, str], list[CandidateDecision]] = defaultdict(list)
    for candidate in candidates:
        if candidate.final_status not in {"emitted", "shared_support"}:
            continue
        for target_uid in candidate.target_feature_uids:
            by_semantic[
                (target_uid, candidate.qualifier, candidate.normalized_value)
            ].append(candidate)
    for candidate in candidates:
        supporting = {
            other.candidate_id
            for target_uid in candidate.target_feature_uids
            for other in by_semantic.get(
                (target_uid, candidate.qualifier, candidate.normalized_value), []
            )
            if other.candidate_id != candidate.candidate_id
        }
        candidate.supporting_candidate_ids = tuple(sorted(supporting))
        if len(supporting) > 1 and candidate.final_status == "emitted":
            candidate.final_status = "shared_support"

    functional_support: dict[tuple[str, str], set[str]] = defaultdict(set)
    for candidate in candidates:
        if candidate.final_status in {
            "emitted",
            "shared_support",
        } and not candidate.confidence_diagnostics.get("provenance"):
            for target_uid in candidate.target_feature_uids:
                functional_support[(candidate.source_id, target_uid)].add(
                    candidate.candidate_id
                )
    for candidate in candidates:
        if candidate.confidence_diagnostics.get("provenance"):
            supporting = {
                candidate_id
                for target_uid in candidate.target_feature_uids
                for candidate_id in functional_support.get(
                    (candidate.source_id, target_uid), set()
                )
            }
            candidate.supporting_candidate_ids = tuple(sorted(supporting))

    candidates.sort(key=lambda item: item.candidate_id)
    rows = [candidate.as_dict() for candidate in candidates]
    counts = {
        "candidate_count": len(candidates),
        "accepted_candidate_count": sum(
            candidate.final_status
            in {"emitted", "shared_support", "supported_existing"}
            for candidate in candidates
        ),
        "emitted_candidate_count": sum(
            candidate.final_status in {"emitted", "shared_support"}
            for candidate in candidates
        ),
        "suppressed_candidate_count": sum(
            candidate.final_status == "suppressed_conflict" for candidate in candidates
        ),
        "candidate_insertion_count": len(final),
        "candidate_emitting_feature_count": len(
            {
                target_uid
                for candidate in candidates
                if candidate.final_status in {"emitted", "shared_support"}
                for target_uid in candidate.target_feature_uids
            }
        ),
        "candidate_ledger_schema": "enrich-bakta.candidate-ledger.v2",
    }
    return rows, counts


def validate_candidate_ledger(
    rows: Iterable[Mapping[str, Any]],
    final_insertions: Iterable[Insertion],
    *,
    base: RawDocument | None = None,
) -> None:
    """Check candidate-to-insertion references before artifacts are promoted."""
    final_items = [item for item in final_insertions if item.qualifier != "COMMENT"]
    insertion_by_id = {insertion_uid(item): item for item in final_items}
    if len(insertion_by_id) != len(final_items):
        raise MergeError("candidate ledger validation found duplicate insertion IDs")
    feature_index = _feature_index(base) if base is not None else {}
    valid_feature_uids = set(feature_index.values()) if base is not None else None
    ledger_rows = list(rows)
    seen: set[str] = set()
    referenced_final_ids: set[str] = set()
    for row in ledger_rows:
        candidate_id = row.get("candidate_id")
        if not isinstance(candidate_id, str) or not candidate_id:
            raise MergeError("candidate ledger contains a missing candidate_id")
        if candidate_id in seen:
            raise MergeError(f"candidate ledger repeats {candidate_id!r}")
        seen.add(candidate_id)
        status = row.get("final_status")
        qualifier = row.get("qualifier")
        normalized_value = row.get("normalized_value")
        if not isinstance(qualifier, str) or not qualifier:
            raise MergeError(f"candidate {candidate_id!r} has no qualifier")
        if not isinstance(normalized_value, str):
            raise MergeError(
                f"candidate {candidate_id!r} has an invalid normalized value"
            )
        targets = row.get("target_feature_uids", [])
        if (
            not isinstance(targets, list)
            or any(not isinstance(value, str) or not value for value in targets)
            or len(set(targets)) != len(targets)
        ):
            raise MergeError(
                f"candidate {candidate_id!r} has invalid target feature UIDs"
            )
        if valid_feature_uids is not None and not set(targets).issubset(
            valid_feature_uids
        ):
            raise MergeError(
                f"candidate {candidate_id!r} references a feature outside the base"
            )
        referenced = row.get("insertion_ids", [])
        if (
            not isinstance(referenced, list)
            or any(not isinstance(value, str) for value in referenced)
            or len(set(referenced)) != len(referenced)
        ):
            raise MergeError(f"candidate {candidate_id!r} has invalid insertion IDs")
        supporting = row.get("supporting_candidate_ids", [])
        if (
            not isinstance(supporting, list)
            or any(not isinstance(value, str) or not value for value in supporting)
            or len(set(supporting)) != len(supporting)
        ):
            raise MergeError(
                f"candidate {candidate_id!r} has invalid supporting candidate IDs"
            )
        if status in {"emitted", "shared_support"}:
            if not referenced:
                raise MergeError(
                    f"accepted candidate {candidate_id!r} has no insertion support"
                )
            for reference in referenced:
                insertion = insertion_by_id.get(reference)
                if insertion is None:
                    raise MergeError(
                        f"accepted candidate {candidate_id!r} references a missing insertion"
                    )
                insertion_target = _insertion_feature_uid(insertion, feature_index)
                if base is not None and insertion_target not in set(targets):
                    raise MergeError(
                        f"candidate {candidate_id!r} references an insertion for the wrong feature"
                    )
                if insertion.qualifier != qualifier:
                    raise MergeError(
                        f"candidate {candidate_id!r} references the wrong qualifier"
                    )
                if insertion.value != normalized_value:
                    raise MergeError(
                        f"candidate {candidate_id!r} references the wrong value"
                    )
                referenced_final_ids.add(reference)
        if status == "suppressed_conflict" and referenced:
            raise MergeError(
                f"suppressed candidate {candidate_id!r} still references insertions"
            )

    by_id = {str(row["candidate_id"]): row for row in ledger_rows}
    for row in ledger_rows:
        candidate_id = str(row["candidate_id"])
        for supporting_id in row.get("supporting_candidate_ids", []):
            if supporting_id not in by_id:
                raise MergeError(
                    f"candidate {candidate_id!r} references unknown support {supporting_id!r}"
                )
            if supporting_id == candidate_id:
                raise MergeError(f"candidate {candidate_id!r} cannot support itself")

    for row in ledger_rows:
        if row.get("final_status") not in {"emitted", "shared_support"}:
            continue
        provenance_refs = [
            insertion_by_id[reference]
            for reference in row.get("insertion_ids", [])
            if _is_provenance_insertion(insertion_by_id[reference])
        ]
        if not provenance_refs:
            continue
        supporting_rows = [
            by_id[supporting_id]
            for supporting_id in row.get("supporting_candidate_ids", [])
            if by_id[supporting_id].get("final_status") in {"emitted", "shared_support"}
            and not by_id[supporting_id]
            .get("confidence_diagnostics", {})
            .get("provenance")
            and by_id[supporting_id].get("source_id") == row.get("source_id")
            and set(by_id[supporting_id].get("target_feature_uids", []))
            & set(row.get("target_feature_uids", []))
        ]
        if not supporting_rows:
            raise MergeError(
                f"provenance candidate {row['candidate_id']!r} has no functional support"
            )

    orphaned = set(insertion_by_id) - referenced_final_ids
    if orphaned:
        raise MergeError(
            "candidate ledger leaves final insertions unsupported: "
            + ", ".join(sorted(orphaned))
        )
