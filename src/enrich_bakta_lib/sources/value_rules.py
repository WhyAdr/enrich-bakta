"""Field-specific validation and affirmative structured-token rules."""

from __future__ import annotations

import re
from dataclasses import dataclass

from enrich_bakta_lib.core.merge_engine import MergeError

EC_FULL_RE = re.compile(r"\d+\.\d+\.\d+\.\d+\Z")
EC_PARTIAL_RE = re.compile(r"(?:\d+\.\d+\.\d+\.-|\d+\.\d+\.-\.-|\d+\.-\.-\.-)\Z")
STRUCTURED_TOKEN_RE = re.compile(
    r"(?:GO:\d{7}|KEGG:K\d{5}|CAZy:(?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?)\Z"
)
STRUCTURED_TOKEN_SEARCH_RE = re.compile(
    r"(?:GO:\d{7}|KEGG:K\d{5}|CAZy:(?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?)"
)
KOFAM_NOTE_RE = re.compile(
    r"KofamScan:(K\d{5});threshold=[^;]+;score=[^;]+;E-value=[^;]+\Z"
)
MISSING_VALUES = {"", "-", "NA"}


@dataclass(frozen=True)
class ValueValidationResult:
    """Auditable result for one raw annotation value."""

    raw: str
    normalized: str | None
    status: str
    reason: str


def validate_ec_value(
    raw: str, *, source: str, allow_partial: bool = True
) -> ValueValidationResult:
    normalized = raw.strip()
    if not normalized:
        return ValueValidationResult(raw, None, "missing", "empty EC value")
    if normalized.lower().startswith("ec:"):
        normalized = normalized[3:]
    if EC_FULL_RE.fullmatch(normalized):
        return ValueValidationResult(raw, normalized, "valid", "full EC grammar")
    if allow_partial and EC_PARTIAL_RE.fullmatch(normalized):
        return ValueValidationResult(raw, normalized, "valid", "partial EC grammar")
    return ValueValidationResult(
        raw,
        normalized,
        "invalid",
        f"{source}: invalid EC value {raw!r}",
    )


def validate_structural_xref(value: str, *, source: str) -> ValueValidationResult:
    prefix, separator, identifier = value.partition(":")
    if prefix not in {"afdb_v6", "cath", "pdb"} or not separator:
        return ValueValidationResult(
            value,
            None,
            "not_applicable",
            "xref prefix is not a supported structural database",
        )
    normalized_identifier = identifier.strip()
    if not normalized_identifier or any(
        ord(char) < 32 or ord(char) == 127 for char in normalized_identifier
    ):
        return ValueValidationResult(
            value,
            None,
            "invalid",
            f"{source}: structural xref has an empty/invalid identifier {value!r}",
        )
    return ValueValidationResult(
        value,
        f"{prefix}:{normalized_identifier}",
        "valid",
        "supported structural xref",
    )


def validate_preferred_name(value: str, *, clean_suffix: bool) -> ValueValidationResult:
    candidate = value.strip()
    if candidate in MISSING_VALUES:
        return ValueValidationResult(
            value, None, "missing", "preferred name is missing"
        )
    if clean_suffix:
        candidate = re.sub(r"_\d+\Z", "", candidate)
    if not candidate:
        return ValueValidationResult(
            value, None, "invalid", "preferred name is empty after normalization"
        )
    return ValueValidationResult(value, candidate, "valid", "preferred name accepted")


def _parse_structured_note(value: str) -> tuple[str, ...]:
    kofam = KOFAM_NOTE_RE.fullmatch(value)
    if kofam:
        return (f"KEGG:{kofam.group(1)}",)
    parts = [part.strip() for part in re.split(r"[,;]", value)]
    if not parts or any(not STRUCTURED_TOKEN_RE.fullmatch(part) for part in parts):
        return ()
    return tuple(dict.fromkeys(parts))


def validate_structured_note(note: str) -> ValueValidationResult:
    value = note.strip()
    if not value:
        return ValueValidationResult(note, None, "missing", "note is empty")
    tokens = _parse_structured_note(value)
    if not tokens:
        return ValueValidationResult(
            note,
            None,
            "invalid",
            "note does not contain only recognized structured tokens",
        )
    return ValueValidationResult(
        note,
        ", ".join(tokens),
        "valid",
        "recognized structured tokens",
    )


def validated_ec_values(
    values: list[str], *, source: str, allow_partial: bool = True
) -> set[str]:
    result: set[str] = set()
    for raw in values:
        validation = validate_ec_value(raw, source=source, allow_partial=allow_partial)
        if validation.status == "missing":
            continue
        if validation.status != "valid" or validation.normalized is None:
            raise MergeError(validation.reason)
        result.add(validation.normalized)
    return result


def structural_xref(value: str, *, source: str) -> str | None:
    validation = validate_structural_xref(value, source=source)
    if validation.status == "not_applicable":
        return None
    if validation.status != "valid" or validation.normalized is None:
        raise MergeError(validation.reason)
    return validation.normalized


def normalized_preferred_name(value: str, *, clean_suffix: bool) -> str | None:
    validation = validate_preferred_name(value, clean_suffix=clean_suffix)
    return validation.normalized if validation.status == "valid" else None


def structured_note_tokens(note: str) -> tuple[str, ...]:
    """Return tokens only from an affirmative, recognized note grammar."""
    value = note.strip()
    return _parse_structured_note(value) if value else ()


def structured_tokens_in_text(value: str) -> tuple[str, ...]:
    """Extract only complete tokens from an already approved text field."""
    return tuple(
        dict.fromkeys(
            token
            for token in STRUCTURED_TOKEN_SEARCH_RE.findall(value)
            if STRUCTURED_TOKEN_RE.fullmatch(token)
        )
    )
