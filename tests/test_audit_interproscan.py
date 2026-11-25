from __future__ import annotations

import hashlib
import json
from pathlib import Path
from urllib.parse import quote

import pytest
from jsonschema import Draft202012Validator

from enrich_bakta_lib.core.types import MergeError
from enrich_bakta_lib.workflows.audit_interproscan import (
    InterProScanJsonReader,
    main,
)

_VERSION = "5.59-91.0"
_SEQUENCE = "MNGIGK"
_MD5 = hashlib.md5(_SEQUENCE.encode(), usedforsecurity=False).hexdigest()
_SCHEMA = Path("schemas/interproscan-audit.v1.schema.json")


def _match(library, accession, start=1, end=3, **kw):
    return {
        "signature": {
            "accession": accession,
            "signatureLibraryRelease": {"library": library, "version": "1.0"},
        },
        "locations": [{"start": start, "end": end}],
        **kw,
    }


@pytest.fixture
def exports(tmp_path):
    tsv, js, gff = (tmp_path / f"sample.{ext}" for ext in ("tsv", "json", "gff3"))
    proteins = [
        {
            "sequence": _SEQUENCE,
            "md5": _MD5,
            "xref": [{"id": "a", "label": "ÂµÎ²"}, {"id": "b"}],
            "matches": [
                _match("PFAM", "PF001"),
                _match("PROSITE_PATTERNS", "PS001", 2, 4),
                _match("MOBIDB_LITE", "mobi;1", 2, 5),
                _match("PANTHER", "PTHR001", 3, 5, **{"model-ac": "PTHR001:SF2"}),
                _match("PANTHER", "PTHR002", 4, 5, **{"model-ac": "PTHR002:SF5"}),
                _match("PIRSR", "PIRSR001", 4, 6),
            ],
        },
        {"sequence": _SEQUENCE, "md5": _MD5, "xref": [{"id": "c"}], "matches": []},
    ]
    rows = [
        ("Pfam", "PF001", 1, 3),
        ("ProSitePatterns", "PS001", 2, 4),
        ("MobiDBLite", "mobi;1", 2, 5),
        ("PANTHER", "PTHR001", 3, 5),
        ("PANTHER", "PTHR002:SF5", 4, 5),
    ]
    tsv.write_text(
        "".join(
            f"{q}\t{_MD5}\t6\t{db}\t{acc}\t-\t{start}\t{end}\t-\tT\t01-01-2025\t-\t-\n"
            for q in ("a", "b")
            for db, acc, start, end in rows
        ),
        encoding="utf-8",
    )
    js.write_text(
        json.dumps(
            {"results": proteins, "interproscan-version": _VERSION}, ensure_ascii=False
        ),
        encoding="utf-8",
    )
    gff.write_bytes(
        (
            f"##gff-version 3\n##interproscan-version {_VERSION}\n"
            + "".join(
                f"{q}\t.\tpolypeptide\t1\t6\t.\t+\t.\tID={q}\n" for q in ("a", "b")
            )
            + "".join(
                f"{q}\t{db}\tprotein_match\t{start}\t{end}\t.\t+\t.\t"
                f"oldName=ignored;Name={quote(acc)};other=x\n"
                for q in ("a", "b")
                for db, acc, start, end in [*rows, ("PIRSR", "PIRSR001", 4, 6)]
            )
            + "##FASTA\n>a\nMNGIGK\n>b\nMNGIGK\n"
        ).encode()
        + b"x" * (1024 * 1024 + 1)  # No feature-row line limit applies to FASTA.
    )
    return tsv, js, gff


def _run(exports, tmp_path, *options):
    tsv, js, gff = exports
    output = tmp_path / "audit.json"
    command = [
        "--tsv",
        str(tsv),
        "--interproscan-version",
        _VERSION,
        "--output",
        str(output),
    ]
    for flag, path in (("--json", js), ("--gff3", gff)):
        if path is not None:
            command.extend([flag, str(path)])
    code = main([*command, *options])
    payload = json.loads(output.read_bytes())
    return code, payload


def _validate_schema(payload):
    schema = json.loads(_SCHEMA.read_bytes())
    Draft202012Validator.check_schema(schema)
    Draft202012Validator(schema).validate(payload)


def test_audit_all_exports(tmp_path, exports):
    markdown = tmp_path / "report.md"
    code, payload = _run(exports, tmp_path, "--markdown", str(markdown))
    assert code == 0
    assert payload["status"] == "passed"
    assert payload["tsv_summary"] == {
        "rows": 10,
        "unique_queries": 2,
        "unique_match_keys": 10,
        "occurrences_total": 10,
    }
    for kind, source in zip(("tsv", "json", "gff3"), exports, strict=True):
        metadata = payload["inputs"][kind]
        assert metadata["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
        assert metadata["size_bytes"] == source.stat().st_size
    for kind in ("json_audit", "gff3_audit"):
        assert payload[kind]["extra_occurrences_by_analysis"] == {"pirsr": 2}
        assert payload[kind]["missing_tsv_occurrences"] == {"count": 0, "samples": []}
    summary = payload["json_audit"]
    assert summary["multi_xref_proteins"] == {
        "count": 1,
        "samples": [{"md5": _MD5, "query_ids": ["a", "b"]}],
    }
    assert summary["total_xrefs"] == 3
    assert summary["unmatched_queries_count"] == 1
    assert summary["md5_mismatches"] == 0
    assert "Extra pirsr occurrences: 2" in markdown.read_text("utf-8")
    original = (tmp_path / "audit.json").read_bytes()
    assert _run(exports, tmp_path)[0] == 0
    assert (tmp_path / "audit.json").read_bytes() == original
    _validate_schema(payload)


@pytest.mark.parametrize("kind", ["json", "gff3"])
def test_version_mismatch_and_nonstrict_exit(exports, tmp_path, kind):
    tsv, js, gff = exports
    selected = js if kind == "json" else gff
    selected.write_bytes(selected.read_bytes().replace(_VERSION.encode(), b"5.60-92.0"))
    selected_exports = (
        tsv,
        js if kind == "json" else None,
        gff if kind == "gff3" else None,
    )
    code, report = _run(selected_exports, tmp_path)
    assert code == 1
    assert report["status"] == "failed"
    assert report["version_check"]["version_mismatch"]
    assert _run(selected_exports, tmp_path, "--no-strict") == (0, report)
    _validate_schema(report)


def test_panther_family_and_subfamily_match_without_order_dependence(exports, tmp_path):
    tsv, js, _ = exports
    source = json.loads(js.read_bytes())
    for query in ("a", "b"):
        with tsv.open("ab") as handle:
            handle.write(
                f"{query}\t{_MD5}\t6\tPANTHER\tPTHR001:SF2\t-\t3\t5\t-\tT\t01-01-2025\t-\t-\n".encode()
            )
    matches = source["results"][0]["matches"]
    matches.append(_match("PANTHER", "PTHR001", 3, 5, **{"model-ac": "PTHR001:SF3"}))
    js.write_text(json.dumps(source), encoding="utf-8")
    assert _run((tsv, js, None), tmp_path)[0] == 0
    matches.reverse()
    js.write_text(json.dumps(source), encoding="utf-8")
    assert _run((tsv, js, None), tmp_path)[0] == 0


def test_missing_occurrences_respect_tsv_multiplicity(exports, tmp_path):
    tsv, js, gff = exports
    with tsv.open("ab") as handle:
        handle.write(tsv.read_bytes().splitlines()[0] + b"\n")
    code, report = _run((tsv, js, gff), tmp_path)
    assert code == 1
    for kind in ("json_audit", "gff3_audit"):
        assert report[kind]["missing_tsv_occurrences"] == {
            "count": 1,
            "samples": [
                {
                    "query_id": "a",
                    "analysis": "pfam",
                    "signature_accession": "PF001",
                    "start": 1,
                    "end": 3,
                    "occurrences": 1,
                }
            ],
        }


def test_md5_mismatch(exports, tmp_path):
    tsv, js, _ = exports
    data = json.loads(js.read_bytes())
    data["results"][0]["sequence"] += "W"
    js.write_text(json.dumps(data), encoding="utf-8")
    code, payload = _run((tsv, js, None), tmp_path)
    assert code == 1
    assert payload["json_audit"]["md5_mismatches"] == 1


def test_json_reader_split_utf8_and_trailing_hash(exports):
    _, path, _ = exports
    payload = json.loads(path.read_bytes())
    path.write_bytes(path.read_bytes() + b" \r\n " * 40)
    for chunk_size in (1, 64):
        reader = InterProScanJsonReader(path, chunk_size=chunk_size)
        assert list(reader) == payload["results"]
        assert reader.producer_version == _VERSION
        assert reader.sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
        assert reader.byte_count == path.stat().st_size


@pytest.mark.parametrize(
    "content",
    [
        b'{"results": [ {} {} ]}',
        b'{"results": [ {}, ]}',
        b'{"results": [ {} ]',
        b'{"results": [ {} ]} trailing',
        b'{"other": {"results": []}}',
        b'{"results": [], "results": []}',
    ],
)
def test_json_reader_rejects_malformed_or_missing_results(tmp_path, content):
    path = tmp_path / "sample.json"
    path.write_bytes(content)
    with pytest.raises(MergeError):
        list(InterProScanJsonReader(path, chunk_size=4))


def test_cli_rejects_missing_comparison_and_collisions(exports):
    tsv, js, _ = exports
    args = ["--tsv", str(tsv), "--interproscan-version", _VERSION]
    assert main(args) == 1
    original = js.read_bytes()
    assert main([*args, "--json", str(js), "--output", str(js)]) == 1
    assert js.read_bytes() == original
