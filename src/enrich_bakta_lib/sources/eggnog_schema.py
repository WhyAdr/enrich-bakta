"""Versioned eggNOG-mapper confidence contracts."""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass

from enrich_bakta_lib.core.merge_engine import MergeError

CONFIDENCE_SCORED_FIELDS = (
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
)
CONFIDENCE_FIELD_INDEX = {
    field: index for index, field in enumerate(CONFIDENCE_SCORED_FIELDS)
}
CONFIDENCE_ORDER_PREFIX = "## confidence field order:"


@dataclass(frozen=True)
class EggnogSchemaContract:
    schema_id: str
    producer_version_key: str
    confidence_field_order: tuple[str, ...]
    allow_missing_legend: bool = True


SCHEMA_V3_BETA6 = "eggnog-v3.0.0-beta6"
SCHEMA_REGISTRY: dict[str, EggnogSchemaContract] = {
    SCHEMA_V3_BETA6: EggnogSchemaContract(
        schema_id=SCHEMA_V3_BETA6,
        producer_version_key="3.0.0-beta6",
        confidence_field_order=CONFIDENCE_SCORED_FIELDS,
    )
}

_VERSION_KEY_RE = re.compile(r"^v?(?P<key>\d+\.\d+\.\d+(?:-[A-Za-z0-9]+)?)(?:[-+].*)?$")


def producer_version_key(value: str | None) -> str | None:
    if value is None:
        return None
    match = _VERSION_KEY_RE.fullmatch(value.strip())
    return match.group("key") if match else None


def resolve_eggnog_contract(
    *,
    producer_version: str | None,
    requested_version: str | None,
    requested_schema: str | None,
    declared_order: Iterable[str] | None,
    columns: Iterable[str],
) -> tuple[EggnogSchemaContract, str]:
    """Resolve a known schema without silently inheriting a positional contract."""
    raw_producer = producer_version.strip() if producer_version else None
    raw_requested = requested_version.strip() if requested_version else None
    producer_key = producer_version_key(raw_producer)
    requested_key = producer_version_key(raw_requested)
    if raw_producer and producer_key is None:
        if requested_schema is None:
            raise MergeError(
                f"eggNOG producer version is not parseable: {raw_producer!r}"
            )
    if raw_requested and requested_key is None:
        if requested_schema is None:
            if raw_producer is not None:
                raise MergeError(
                    f"eggNOG version mismatch: declared={raw_producer!r}, requested={raw_requested!r}"
                )
            raise MergeError(
                f"requested eggNOG version is not parseable: {raw_requested!r}"
            )
    if producer_key and requested_key and producer_key != requested_key:
        raise MergeError(
            f"eggNOG version mismatch: declared={raw_producer!r}, requested={raw_requested!r}"
        )

    effective_key = producer_key or requested_key
    selection_source = "version_registry"
    if requested_schema is not None:
        try:
            contract = SCHEMA_REGISTRY[requested_schema]
        except KeyError as exc:
            raise MergeError(f"unsupported eggNOG schema {requested_schema!r}") from exc
        known_keys = {
            candidate.producer_version_key for candidate in SCHEMA_REGISTRY.values()
        }
        if (
            effective_key in known_keys
            and effective_key != contract.producer_version_key
        ):
            raise MergeError(
                "explicit eggNOG schema conflicts with the declared/requested "
                f"producer version {raw_producer or raw_requested!r}"
            )
        selection_source = "explicit_schema"
    else:
        if effective_key is None:
            raise MergeError(
                "eggNOG producer metadata is absent; provide --eggnog-schema"
            )
        matching = [
            candidate
            for candidate in SCHEMA_REGISTRY.values()
            if candidate.producer_version_key == effective_key
        ]
        if not matching:
            raise MergeError(
                "unsupported eggNOG producer version; provide an explicit compatible "
                "--eggnog-schema if policy permits it"
            )
        contract = matching[0]

    if declared_order is not None:
        order = tuple(declared_order)
        if not order:
            raise MergeError("eggNOG confidence field-order legend is empty")
        if order != contract.confidence_field_order:
            raise MergeError(
                "eggNOG confidence field order does not match the selected schema: "
                f"{order!r}"
            )
    elif not contract.allow_missing_legend:
        raise MergeError(
            "selected eggNOG schema requires a confidence field-order legend"
        )

    header_order = tuple(
        column for column in columns if column in CONFIDENCE_FIELD_INDEX
    )
    expected_header_order = tuple(
        column for column in contract.confidence_field_order if column in header_order
    )
    if header_order != expected_header_order:
        raise MergeError(
            "eggNOG scored columns conflict with the selected confidence field-order contract"
        )
    return contract, selection_source
