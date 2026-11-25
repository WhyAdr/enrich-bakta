"""Compare InterProScan export occurrences without feeding them into enrichment."""

from __future__ import annotations

import argparse
import codecs
import hashlib
import json
import re
import sys
from collections import Counter
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any, BinaryIO
from urllib.parse import unquote

from enrich_bakta_lib.core.artifacts import (
    json_sidecar_bytes,
    preflight_artifact_paths,
    stage_artifacts,
)
from enrich_bakta_lib.core.types import TOOL_VERSION, MergeError
from enrich_bakta_lib.sources.interproscan import (
    INTERPROSCAN_PROFILES,
    MAX_RAW_LINE_BYTES,
    TSV_LAYOUT_PROFILES,
    validate_tsv_row_layout,
)

_SCHEMA = "enrich-bakta.interproscan-audit.v1"
_SAMPLE_LIMIT = 50
_CHUNK_SIZE = 1024 * 1024
_ANALYSIS_NORM_RE = re.compile(r"[^a-z0-9]+")
_WHITESPACE = re.compile(r"[ \t\r\n]*")
MatchKey = tuple[str, str, str, int, int]


def normalize_analysis_name(name: str) -> str:
    return _ANALYSIS_NORM_RE.sub("", name.lower())


class _ChunkedJson:
    """Decode UTF-8 and track exact input bytes as the JSON cursor advances."""

    def __init__(self, handle: BinaryIO, chunk_size: int) -> None:
        self.handle = handle
        self.chunk_size = chunk_size
        self.text = ""
        self.pos = 0
        self.eof = False
        self.hasher = hashlib.sha256()
        self.byte_count = 0
        self.decoder = codecs.getincrementaldecoder("utf-8")()
        self.json_decoder = json.JSONDecoder(parse_constant=self._invalid_constant)

    @staticmethod
    def _invalid_constant(token: str) -> Any:
        raise MergeError(f"invalid JSON constant {token!r}")

    def refill(self) -> None:
        raw = self.handle.read(self.chunk_size)
        self.eof = not raw
        self.hasher.update(raw)
        self.byte_count += len(raw)
        self.text = self.text[self.pos :] + self.decoder.decode(raw, final=self.eof)
        self.pos = 0

    def peek(self) -> str:
        while True:
            whitespace = _WHITESPACE.match(self.text, self.pos)
            assert whitespace is not None
            self.pos = whitespace.end()
            if self.pos < len(self.text):
                return self.text[self.pos]
            if self.eof:
                return ""
            self.refill()

    def consume(self, token: str) -> None:
        if self.peek() != token:
            raise MergeError(f"expected JSON token {token!r}")
        self.pos += 1

    def value(self) -> Any:
        self.peek()
        while True:
            try:
                obj, end = self.json_decoder.raw_decode(self.text, self.pos)
                if end < len(self.text) or self.eof:
                    self.pos = end
                    return obj
            except json.JSONDecodeError as exc:
                if self.eof:
                    raise MergeError(f"malformed JSON: {exc.msg}") from exc
            self.refill()


class InterProScanJsonReader:
    """Iterate the top-level results array, buffering one protein at a time."""

    def __init__(self, path: Path, *, chunk_size: int = _CHUNK_SIZE) -> None:
        if chunk_size < 1:
            raise ValueError("chunk_size must be positive")
        self.path = path
        self.chunk_size = chunk_size
        self.producer_version: str | None = None
        self.sha256 = ""
        self.byte_count = 0

    def __iter__(self) -> Iterator[dict[str, Any]]:
        self.producer_version = None
        self.sha256 = ""
        self.byte_count = 0
        with _audit_input(self.path, "JSON") as handle:
            reader = _ChunkedJson(handle, self.chunk_size)
            reader.consume("{")
            fields: set[str] = set()
            while reader.peek() != "}":
                if fields:
                    reader.consume(",")
                field = reader.value()
                if not isinstance(field, str) or field in fields:
                    raise MergeError("expected a distinct top-level JSON field name")
                fields.add(field)
                reader.consume(":")
                if field == "results":
                    reader.consume("[")
                    first = True
                    while reader.peek() != "]":
                        if not first:
                            reader.consume(",")
                        protein = reader.value()
                        if not isinstance(protein, dict):
                            raise MergeError(
                                "results array must contain protein objects"
                            )
                        first = False
                        yield protein
                    reader.consume("]")
                else:
                    value = reader.value()
                    if field == "interproscan-version":
                        if not isinstance(value, str) or not value:
                            raise MergeError("invalid interproscan-version value")
                        self.producer_version = value
            reader.consume("}")
            if "results" not in fields:
                raise MergeError("missing top-level results array")
            if reader.peek():
                raise MergeError("unexpected data after top-level JSON object")
            self.sha256 = reader.hasher.hexdigest()
            self.byte_count = reader.byte_count


@contextmanager
def _audit_input(path: Path, kind: str) -> Iterator[BinaryIO]:
    try:
        with path.open("rb") as handle:
            yield handle
    except (MergeError, OSError, UnicodeError, RecursionError) as exc:
        raise MergeError(f"InterProScan {kind} {path}: {exc}") from exc


class _LineInput:
    def __init__(self, handle: BinaryIO) -> None:
        self.handle = handle
        self.hasher = hashlib.sha256()
        self.byte_count = 0

    def __iter__(self) -> Iterator[tuple[int, str]]:
        for number, raw in enumerate(
            iter(lambda: self.handle.readline(MAX_RAW_LINE_BYTES + 2), b""), 1
        ):
            self._update(raw)
            content = raw.removesuffix(b"\n").removesuffix(b"\r")
            if len(content) > MAX_RAW_LINE_BYTES:
                raise MergeError(f"line {number} exceeds maximum allowed size of 1 MiB")
            yield number, content.decode("utf-8")

    def _update(self, raw: bytes) -> None:
        self.hasher.update(raw)
        self.byte_count += len(raw)

    def drain(self) -> None:
        for raw in iter(lambda: self.handle.read(_CHUNK_SIZE), b""):
            self._update(raw)

    def metadata(self, path: Path) -> dict[str, Any]:
        return _file_metadata(path, self.byte_count, self.hasher.hexdigest())


def _file_metadata(path: Path, byte_count: int, sha256: str) -> dict[str, Any]:
    return {"path": str(path), "size_bytes": byte_count, "sha256": sha256}


def _match_key(
    query: str, analysis: str, signature: str, start: Any, end: Any
) -> MatchKey:
    normalized = normalize_analysis_name(analysis)
    if not query or not normalized or not signature:
        raise MergeError("empty query, analysis, or signature accession")
    for value in (start, end):
        if type(value) is not int and (
            not isinstance(value, str) or not re.fullmatch(r"[0-9]+", value)
        ):
            raise MergeError(f"invalid coordinate {value!r}")
    try:
        start, end = int(start), int(end)
    except ValueError as exc:
        raise MergeError(f"invalid coordinates {start!r}..{end!r}") from exc
    if not 1 <= start <= end:
        raise MergeError(f"invalid coordinates {start}..{end}")
    return query, normalized, signature, start, end


def _read_tsv(path: Path, layout: str) -> tuple[Counter[MatchKey], dict[str, Any]]:
    counts: Counter[MatchKey] = Counter()
    with _audit_input(path, "TSV") as handle:
        stream = _LineInput(handle)
        for number, line in stream:
            if not line:
                continue
            fields = line.split("\t")
            validate_tsv_row_layout(fields, layout, number)
            try:
                counts[
                    _match_key(
                        fields[0].strip(),
                        fields[3].strip(),
                        fields[4].strip(),
                        fields[6].strip(),
                        fields[7].strip(),
                    )
                ] += 1
            except MergeError as exc:
                raise MergeError(f"line {number}: {exc}") from exc
        return counts, stream.metadata(path)


class _Occurrences:
    def __init__(self, expected: Counter[MatchKey]) -> None:
        self.remaining = expected.copy()
        self.extras: Counter[str] = Counter()
        self.panther_pairs: Counter[tuple[MatchKey, MatchKey]] = Counter()

    def count(self, key: MatchKey, alternate_signature: str | None = None) -> None:
        query, analysis, signature, start, end = key
        if (
            analysis == "panther"
            and alternate_signature
            and alternate_signature.startswith(f"{signature}:")
        ):
            alternate = (query, analysis, alternate_signature, start, end)
            if alternate in self.remaining:
                self.panther_pairs[key, alternate] += 1
                return
        if self.remaining[key]:
            self.remaining[key] -= 1
        else:
            self.extras[analysis] += 1

    def summary(self) -> dict[str, Any]:
        # Count exact rows first. A family/subfamily occurrence can then satisfy
        # its remaining subfamily first, or the shared family as a fallback.
        for (family, subfamily), count in sorted(self.panther_pairs.items()):
            subfamily_matches = min(count, self.remaining[subfamily])
            self.remaining[subfamily] -= subfamily_matches
            count -= subfamily_matches
            family_matches = min(count, self.remaining[family])
            if family_matches:
                self.remaining[family] -= family_matches
            count -= family_matches
            if count:
                self.extras["panther"] += count
        samples = []
        for (query, analysis, signature, start, end), occurrences in sorted(
            (+self.remaining).items()
        )[:_SAMPLE_LIMIT]:
            samples.append(
                {
                    "query_id": query,
                    "analysis": analysis,
                    "signature_accession": signature,
                    "start": start,
                    "end": end,
                    "occurrences": occurrences,
                }
            )
        return {
            "missing_tsv_occurrences": {
                "count": self.remaining.total(),
                "samples": samples,
            },
            "extra_occurrences_by_analysis": dict(self.extras),
        }


def _required_field(obj: Any, key: str, expected_type: type[Any]) -> Any:
    value = obj.get(key) if isinstance(obj, dict) else None
    if not isinstance(value, expected_type):
        raise MergeError(f"field {key!r} must be a {expected_type.__name__}")
    return value


def _audit_json(
    reader: InterProScanJsonReader, expected: Counter[MatchKey], tsv_queries: set[str]
) -> dict[str, Any]:
    occurrences = _Occurrences(expected)
    protein_objects = total_xrefs = multi_xrefs = md5_mismatches = 0
    multi_samples: list[dict[str, Any]] = []
    unmatched_queries: set[str] = set()
    libraries: dict[str, str] = {}
    for protein_objects, protein in enumerate(reader, 1):
        try:
            sequence = _required_field(protein, "sequence", str)
            md5 = _required_field(protein, "md5", str).lower()
            actual_md5 = hashlib.md5(sequence.encode("ascii"), usedforsecurity=False)
            md5_mismatches += actual_md5.hexdigest() != md5
            xrefs = _required_field(protein, "xref", list)
            query_ids = sorted(_required_field(xref, "id", str) for xref in xrefs)
            if not query_ids or not all(query_ids) or len(set(query_ids)) != len(xrefs):
                raise MergeError("xref must contain distinct, nonempty query IDs")
            total_xrefs += len(query_ids)
            unmatched_queries.update(q for q in query_ids if q not in tsv_queries)
            if len(query_ids) > 1:
                multi_xrefs += 1
                if len(multi_samples) < _SAMPLE_LIMIT:
                    multi_samples.append({"md5": md5, "query_ids": query_ids})

            for match in _required_field(protein, "matches", list):
                signature = _required_field(match, "signature", dict)
                accession = _required_field(signature, "accession", str)
                release = _required_field(signature, "signatureLibraryRelease", dict)
                library = _required_field(release, "library", str)
                library_version = _required_field(release, "version", str)
                if libraries.setdefault(library, library_version) != library_version:
                    raise MergeError(f"inconsistent library version for {library!r}")
                alternate = match.get("model-ac")
                for location in _required_field(match, "locations", list):
                    start = _required_field(location, "start", int)
                    end = _required_field(location, "end", int)
                    for query in query_ids:
                        key = _match_key(query, library, accession, start, end)
                        occurrences.count(
                            key, alternate if isinstance(alternate, str) else None
                        )
        except (MergeError, UnicodeError) as exc:
            raise MergeError(
                f"InterProScan JSON {reader.path}, protein {protein_objects}: {exc}"
            ) from exc
    return {
        "protein_objects": protein_objects,
        "total_xrefs": total_xrefs,
        "multi_xref_proteins": {"count": multi_xrefs, "samples": multi_samples},
        "md5_mismatches": md5_mismatches,
        "unmatched_queries_count": len(unmatched_queries),
        "library_versions": libraries,
        **occurrences.summary(),
    }


def _audit_gff3(
    path: Path, expected: Counter[MatchKey]
) -> tuple[str | None, dict[str, Any], dict[str, Any]]:
    occurrences = _Occurrences(expected)
    polypeptides = match_rows = 0
    version = None
    with _audit_input(path, "GFF3") as handle:
        stream = _LineInput(handle)
        for number, line in stream:
            if line == "##FASTA":
                stream.drain()
                break
            if line.startswith("##interproscan-version "):
                declared = line.partition(" ")[2].strip()
                if not declared:
                    raise MergeError(f"line {number}: empty producer version")
                if version is not None and version != declared:
                    raise MergeError(f"line {number}: conflicting producer version")
                version = declared
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise MergeError(f"line {number}: expected 9 GFF3 columns")
            if fields[2] == "polypeptide":
                polypeptides += 1
            if fields[2] == "protein_match":
                match_rows += 1
                _, marker, attributes = (";" + fields[8]).partition(";Name=")
                if not marker:
                    raise MergeError(f"line {number}: missing Name attribute")
                signature = unquote(attributes.split(";", 1)[0], errors="strict")
                try:
                    occurrences.count(
                        _match_key(
                            unquote(fields[0], errors="strict"),
                            fields[1],
                            signature,
                            fields[3],
                            fields[4],
                        )
                    )
                except MergeError as exc:
                    raise MergeError(f"line {number}: {exc}") from exc
        data = stream.metadata(path)
    return (
        version,
        data,
        {
            "polypeptide_count": polypeptides,
            "match_rows": match_rows,
            **occurrences.summary(),
        },
    )


def audit_interproscan(
    tsv_path: Path,
    version: str,
    *,
    layout: str = "ipr-go-pathways",
    json_path: Path | None = None,
    gff3_path: Path | None = None,
) -> dict[str, Any]:
    if json_path is None and gff3_path is None:
        raise MergeError("provide at least one of --json and --gff3")
    profile = INTERPROSCAN_PROFILES.get(version)
    if profile is None:
        raise MergeError(f"unsupported InterProScan producer version {version!r}")
    if layout not in profile.supported_layouts:
        raise MergeError(f"unsupported InterProScan TSV layout {layout!r}")

    expected, tsv_meta = _read_tsv(tsv_path, layout)
    queries = {query for query, *_ in expected}
    files = {"tsv": tsv_meta}
    json_version = gff3_version = None
    report: dict[str, Any] = {}
    missing_count = bad_md5s = 0
    version_mismatch = False
    if json_path is not None:
        reader = InterProScanJsonReader(json_path)
        json_report = _audit_json(reader, expected, queries)
        report["json_audit"] = json_report
        json_version = reader.producer_version
        version_mismatch |= json_version != version
        missing_count += json_report["missing_tsv_occurrences"]["count"]
        bad_md5s = json_report["md5_mismatches"]
        files["json"] = _file_metadata(json_path, reader.byte_count, reader.sha256)
    if gff3_path is not None:
        gff3_version, files["gff3"], gff3_report = _audit_gff3(gff3_path, expected)
        report["gff3_audit"] = gff3_report
        version_mismatch |= gff3_version != version
        missing_count += gff3_report["missing_tsv_occurrences"]["count"]

    return {
        "schema": _SCHEMA,
        "inputs": {
            **files,
            "asserted_version": version,
            "tsv_layout": layout,
        },
        "version_check": {
            "asserted_version": version,
            "json_version": json_version,
            "gff3_version": gff3_version,
            "version_mismatch": version_mismatch,
        },
        "tsv_summary": {
            "rows": expected.total(),
            "unique_queries": len(queries),
            "unique_match_keys": len(expected),
            "occurrences_total": expected.total(),
        },
        **report,
        "status": "failed"
        if version_mismatch or missing_count or bad_md5s
        else "passed",
    }


def render_markdown(report: dict[str, Any]) -> str:
    version = report["version_check"]
    lines = [
        "# InterProScan cross-format audit",
        "",
        f"**Status** {report['status']}",
        f"**Asserted producer version** {version['asserted_version']}",
        f"**Producer version mismatch** {version['version_mismatch']}",
        f"**TSV occurrences** {report['tsv_summary']['occurrences_total']:,}",
    ]
    for key, label in (("json_audit", "JSON"), ("gff3_audit", "GFF3")):
        data = report.get(key)
        if data is None:
            continue
        lines.extend(["", f"## {label}", ""])
        if key == "json_audit":
            lines.extend(
                [
                    f"- Protein objects: {data['protein_objects']:,}",
                    f"- Multi-xref proteins: {data['multi_xref_proteins']['count']:,}",
                    f"- MD5 mismatches: {data['md5_mismatches']:,}",
                ]
            )
        else:
            lines.append(f"- Protein match rows: {data['match_rows']:,}")
        lines.append(
            f"- Missing TSV occurrences: {data['missing_tsv_occurrences']['count']:,}"
        )
        lines.extend(
            f"- Extra {analysis} occurrences: {count:,}"
            for analysis, count in sorted(data["extra_occurrences_by_analysis"].items())
        )
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument("--tsv", type=Path, required=True)
    parser.add_argument("--interproscan-version", required=True)
    parser.add_argument(
        "--interproscan-tsv-layout",
        default="ipr-go-pathways",
        choices=TSV_LAYOUT_PROFILES,
    )
    parser.add_argument("--json", type=Path)
    parser.add_argument("--gff3", type=Path)
    parser.add_argument(
        "--output", type=Path, help="JSON destination (default: stdout)"
    )
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--strict", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args(argv)

    inputs = [path for path in (args.tsv, args.json, args.gff3) if path is not None]
    outputs = [path for path in (args.output, args.markdown) if path is not None]
    try:
        preflight_artifact_paths(outputs, inputs)
        report = audit_interproscan(
            args.tsv,
            args.interproscan_version,
            layout=args.interproscan_tsv_layout,
            json_path=args.json,
            gff3_path=args.gff3,
        )
        data = json_sidecar_bytes(report)
        artifacts = []
        if args.output is not None:
            artifacts.append((args.output, data))
        if args.markdown is not None:
            artifacts.append((args.markdown, render_markdown(report).encode("utf-8")))
        stage_artifacts(artifacts, inputs)
        if args.output is None:
            sys.stdout.write(data.decode("utf-8"))
    except (MergeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return int(args.strict and report["status"] == "failed")


if __name__ == "__main__":
    raise SystemExit(main())
