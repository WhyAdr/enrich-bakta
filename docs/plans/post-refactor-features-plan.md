# Post-Refactor Features Implementation Plan (Codex Handoff)

- **Date:** 2026-11-25
- **Baseline:** `main` at `5ebee37` (`0.4.1` + audit-hardening + 4-phase architectural refactor; 222 unit tests, 23 typed source modules, 42 published LFS artifacts)
- **Target release:** `0.5.0` (bump version only in the final release phase after all gates pass)
- **Execution recipient:** Codex agents
- **Source backlog:** [deferred-features-notes.md](deferred-features-notes.md)

---

## 1. Executive Summary & Scope Split

Now that `core/` is decomposed (`types.py`, `genbank.py`, `lineage.py`, `splicer.py`, `reconcile.py`, `artifacts.py`), GenBank parsing is down to `2` parses per merge, and all four adapters share `MergeContext` / `SourcePlan` / `workflows/pipeline.py`, this plan implements the full post-refactor feature suite cleanly in **five phases**:

| Phase | Deliverable | New / Updated Modules |
|---|---|---|
| **Phase 0** | Small pipeline polish (`PipelineResult` dataclass replacing mutable `candidate_counts_out`) | `workflows/pipeline.py`, `sources/{baktfold,kofam,eggnog,interproscan}.py`, `workflows/enrich.py` |
| **Phase 1** | **Enrichment Summary Report** (`--summary-report` inline flag + standalone `enrich-bakta-report` CLI; deterministic Markdown & self-contained HTML) | `workflows/summary_report.py`, `workflows/pipeline.py`, `workflows/enrich.py`, `core/artifacts.py`, `pyproject.toml`, `tests/test_summary_report.py` |
| **Phase 2** | **InterProScan Cross-Format Audit Tool** (`enrich-bakta-audit-interproscan` streaming TSV vs. JSON/GFF3 parity & MD5/multi-xref checker + JSON schema) | `workflows/audit_interproscan.py`, `schemas/interproscan-audit.v1.schema.json`, `pyproject.toml`, `tests/test_audit_interproscan.py` |
| **Phase 3** | **Locus-by-Locus Diff Tool** (`enrich-bakta-diff` supporting manifest mode and pairwise GBFF mode in `text`/`json`/`tsv`) | `workflows/diff.py`, `pyproject.toml`, `tests/test_diff.py` |
| **Phase 4** | **Producer-Profile Registry & Cross-File Promotion Journal** (typed `InterProScanProfile` registry ready for future `5.6x`/`5.7x` verification; opt-in write-ahead promotion recovery journal in `core/artifacts.py`; adapter onboarding template doc) | `sources/interproscan.py`, `core/artifacts.py`, `docs/architecture/adapter-onboarding.md`, `tests/test_promotion_journal.py` |
| **Phase 5** | **Release & Verification Gate** (bump to `0.5.0`, align baselines/docs/CI wheel smoke checks, run full unit + real-data C14/SM/BK71A verification) | `pyproject.toml`, `core/types.py`, `README.md`, `CHANGELOG.md`, baselines, `.github/workflows/ci.yml` |

### Hard Scientific & Engineering Constraints (Do Not Violate)

1. **Zero GBFF byte drift on existing enrichment runs:** Adding `--summary-report`, `audit-interproscan`, `enrich-bakta-diff`, `PipelineResult`, or the `InterProScanProfile` registry MUST NOT alter any existing C14 or SM GBFF `output_sha256` or insertion count. *(Only the final `0.5.0` version bump in Phase 5 may touch `tool_version` in `followthrough-scientific-baseline.json` / `interproscan-scientific-baseline.json`, and even then all GBFF outputs and annotation sets must be verified identical.)*
2. **Zero new mandatory runtime dependencies:** Use Python stdlib + `biopython>=1.83,<2` only (`html.escape`, `json.JSONDecoder`, `hashlib`, `argparse`, `dataclasses`, `collections`). Do **not** add `ijson`, `jinja2`, or `pandas` to `[project.dependencies]`.
3. **Never fabricate unverified bioinformatics profiles or external adapters:** In Phase 4, introduce the `InterProScanProfile` registry with **`5.59-91.0` as the sole verified runtime profile** (plus documented schema/test hooks for future `5.6x`/`5.7x` or `eggNOG 2.1.x`). Do **not** register unverified producer versions or speculative `dbCAN`/`AMRFinderPlus`/`DefenseFinder`/`antiSMASH` parsers without real validated fixtures; instead document their `SourceAdapter` onboarding contract in `docs/architecture/adapter-onboarding.md`.
4. **Keep JSON/GFF3 strictly out of runtime enrichment:** `audit-interproscan` is a standalone audit/QC command. Runtime InterProScan enrichment remains TSV-only.
5. **Preserve determinism & scientific caveats:** No wall-clock timestamps unless explicitly supplied; always state explicitly that cross-source GO/KO overlaps can share upstream databases and are not independent experimental proof, and that `import-faa` translations remain imported evidence rather than genomic validation.

---

## 2. Phase 0 — Pipeline Cleanup (`PipelineResult`)

Before adding `--summary-report` to `workflows/pipeline.py`, replace the mutable out-parameter `candidate_counts_out: dict[str, Any] | None = None` with a frozen `PipelineResult` dataclass so callers receive `final`, `plans`, and `candidate_counts` cleanly.

### Rough Pseudocode — `workflows/pipeline.py`

```python
@dataclass(frozen=True)
class PipelineResult:
    """Structured return value from execute_merge_pipeline()."""

    final: dict[str, Any]
    plans: tuple[SourcePlan, ...]
    candidate_counts: dict[str, Any]


def execute_merge_pipeline(
    ctx: MergeContext,
    adapters: Sequence[SourceAdapter],
    *,
    operation: str,
    output_path: Path,
    manifest_path: Path | None = None,
    context_report_path: Path | None = None,
    summary_report_path: Path | None = None,  # wired in Phase 1
    gene_conflict_policy: str = "skip",
    reconcile: bool = True,
    include_clean_gene_suffix_policy: bool = False,
    extra_metadata: Mapping[str, Any] | None = None,
    metadata_override: (
        Callable[[dict[str, Any], list[SourcePlan], dict[str, Any]], dict[str, Any]]
        | None
    ) = None,
) -> PipelineResult:
    ...
    return PipelineResult(
        final=final,
        plans=tuple(plans),
        candidate_counts=candidate_counts,
    )
```

Update the five callers (`baktfold.graft`, `kofam.merge`, `eggnog.merge`, `interproscan.merge`, `workflows/enrich.py`) to unpack `result.final`, `result.plans`, and `result.candidate_counts`. Keep every standalone return dict identical!

---

## 3. Phase 1 — Enrichment Summary Report (`--summary-report` & `enrich-bakta-report`)

### Contract & CLI Surface

1. **Inline flag on unified CLI & API (`enrich_bakta_lib.workflows.enrich`):**
   - `--summary-report PATH` (`summary_report_path: Path | None = None`)
   - Format inferred from suffix: `.md` / `.markdown` → Markdown; `.html` / `.htm` → self-contained HTML (inline CSS, zero external scripts/fonts/CDNs, `html.escape()` on every dynamic string). Reject unsupported suffixes loudly with `MergeError`.
   - Must participate in **collision preflight and `stage_artifacts()`** alongside GBFF, manifest, and context sidecar!
   - Works even if `--manifest` was not written to disk (builds the in-memory `enrich-bakta.merge-manifest.v2` dictionary from `metadata`, `rows`, and `sidecar_payload` right inside `finalize_merge()` or immediately before staging, so the report always reflects finalized hashes, insertion counts, and bound translation evidence!).
   - Does **not** mutate GBFF bytes or existing manifest/context SHA-256 digests (`summary_report` is an output view, not hashed back into `manifest["metadata"]` unless explicitly desired — **keep existing manifest/context hashes untouched**!).
2. **Standalone offline CLI (`enrich-bakta-report` → `enrich_bakta_lib.workflows.summary_report:main`):**
   ```text
   enrich-bakta-report --manifest ENRICHED.manifest.json [--context-report ENRICHED.context.json] --output REPORT.html [--format auto|md|html]
   ```
   - Validates input schema (`enrich-bakta.merge-manifest.v2`) and optional context schema (`enrich-bakta.context.v1`, `enrich-bakta.interproscan-context.v1`, or `enrich-bakta.eggnog-context.v1`).
   - If both `--manifest` and `--context-report` are passed and `metadata.context_report_sha256` is present in the manifest, verify the context file's SHA-256 matches `metadata.context_report_sha256` (fail closed on mismatch!).
   - Atomic write via `stage_artifacts()` with path collision checks against all inputs.

### Report Sections & Exact Metrics

Compute everything in one deterministic `SummaryMetrics` dataclass (`build_summary_metrics(manifest_payload, context_payload=None)`):

1. **Run Overview & Provenance Bindings:** `tool_version`, `operation`, `base_sha256`, `output_sha256`, `records`, total `insertions` (functional vs. provenance vs. comment), `merge_timestamp` (or `"omitted (deterministic)"`), source versions & hashes (`baktfold`, `kofamscan`, `eggnog`, `interproscan`).
2. **Per-Source Annotation Gain Table:** Rows = each `source_id` present in `decisions` (sorted); columns = `emitted`, `shared_support`, `supported_existing`, `suppressed_*`, `filtered_confidence`, `skipped_partial_ec`, `invalid_value`, `total`.
3. **Per-Qualifier Breakdown Table:** Rows = `(qualifier, source_id)` across `gene`, `EC_number`, `db_xref`, `note`, `inference`; columns = `emitted`, `shared_support`, `supported_existing`, `suppressed`.
4. **Cross-Source Overlap Summary:**
   - Group `decisions` where `final_status == "shared_support"` by `(target_feature_uids, qualifier, normalized_value)` so each overlapping annotation is counted **once** with the sorted pair/set of corroborating `source_id`s (e.g., `eggNOG + InterProScan` on `db_xref=GO:...`, `KofamScan + eggNOG` on `db_xref=KEGG:...`).
   - Include explicit mandatory caveat text: *"Shared support across bioinformatic tools may stem from common upstream reference databases and does not constitute independent experimental confirmation."*
5. **QC, Conflicts & Suppressions:**
   - **AntiFam QC alerts:** Pull from InterProScan context `metadata.qc_flags.anti_fam` (directly or under `context["interproscan"]` in combined envelope), listing `query_id`, `signature_accession`, `signature_description`, `score`, `start..end`.
   - **Gene conflicts & authoritative base-name collisions:** Counts + bounded sample table (first 50 sorted by `(record, locus_tag)`) of `suppressed_conflict` (`qualifier == "gene"`) and `suppressed_authoritative_gene`, plus `gene_conflict_policy`.
   - **Skipped partial ECs & invalid EC tokens:** Count of `skipped_partial_ec` and `invalid_value` decisions (including `baktfold_invalid_ec_policy`).
6. **Translation Lineage Audit Table:**
   - Pull from `metadata.translation_evidence` (plus `translation_evidence_parent_sha256` and `allow_imported_translations`).
   - List exact restored/carried-forward loci, `origin` (`imported_faa` vs. `genomically_validated`), `restoration_policy`, `genomic_validation_status`, and `bound_output_sha256`.
   - Include explicit caveat if any `imported_faa` entry exists: *"Imported FAA translations reflect explicit protein-artifact restoration, not independent genomic CDS validation."*
7. **Context Sidecar Highlights (when context is available):**
   - eggNOG context hits/values counts.
   - InterProScan pathway occurrence totals by database (with caveat that Reactome taxon labels are database hits, not bacterial pathway completeness calls), `support_digest`, `match_evidence_digest`, `pathway_digest`, and count of queries with `context_signatures_truncated == true`.

### Rough Pseudocode — `workflows/summary_report.py`

```python
# src/enrich_bakta_lib/workflows/summary_report.py
from __future__ import annotations

import collections
import html
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.artifacts import stage_artifacts
from enrich_bakta_lib.core.types import MergeError, read_input_bytes, sha256_bytes


@dataclass(frozen=True)
class SummaryMetrics:
    operation: str
    tool_version: str
    base_sha256: str
    output_sha256: str
    records: int
    insertions: int
    functional_insertions: int
    provenance_insertions: int
    by_source: dict[str, dict[str, int]]
    by_qualifier_source: dict[tuple[str, str], dict[str, int]]
    cross_source_overlaps: dict[tuple[str, str], int]  # (qualifier, "InterProScan + eggNOG") -> unique count
    suppressed_samples: list[dict[str, Any]]
    suppressed_total: int
    anti_fam_hits: list[dict[str, Any]]
    translation_evidence: dict[str, dict[str, Any]]
    has_imported_translations: bool
    interproscan_context_summary: dict[str, Any] | None
    eggnog_context_summary: dict[str, Any] | None


def resolve_report_format(path: Path, explicit: str = "auto") -> str:
    if explicit in {"md", "html"}:
        return explicit
    suffix = path.suffix.lower()
    if suffix in {".md", ".markdown"}:
        return "md"
    if suffix in {".html", ".htm"}:
        return "html"
    raise MergeError(
        f"unsupported summary report extension {path.suffix!r}; use .md or .html"
    )


def build_summary_metrics(
    manifest: dict[str, Any],
    context: dict[str, Any] | None = None,
) -> SummaryMetrics:
    if manifest.get("schema") != "enrich-bakta.merge-manifest.v2":
        raise MergeError("summary report requires enrich-bakta.merge-manifest.v2 JSON")
    metadata = manifest.get("metadata", {})
    decisions = manifest.get("decisions", [])

    by_source: dict[str, dict[str, int]] = collections.defaultdict(
        lambda: collections.defaultdict(int)
    )
    overlap_groups: dict[tuple[tuple[str, ...], str, str], set[str]] = (
        collections.defaultdict(set)
    )
    for dec in decisions:
        src = str(dec["source_id"])
        status = str(dec["final_status"])
        by_source[src][status] += 1
        by_source[src]["total"] += 1
        if status == "shared_support":
            key = (
                tuple(dec["target_feature_uids"]),
                str(dec["qualifier"]),
                str(dec["normalized_value"]),
            )
            overlap_groups[key].add(src)

    cross_source_overlaps: dict[tuple[str, str], int] = collections.defaultdict(int)
    for (_uids, qualifier, _val), srcs in overlap_groups.items():
        if len(srcs) > 1:
            cross_source_overlaps[(qualifier, " + ".join(sorted(srcs)))] += 1

    # Extract InterProScan / eggNOG context regardless of envelope vs standalone
    ips_ctx = None
    egg_ctx = None
    if context is not None:
        schema = context.get("schema")
        if schema == "enrich-bakta.context.v1":
            ips_ctx = context.get("interproscan")
            egg_ctx = context.get("eggnog")
        elif schema == "enrich-bakta.interproscan-context.v1":
            ips_ctx = context
        elif schema == "enrich-bakta.eggnog-context.v1":
            egg_ctx = context
        else:
            raise MergeError(f"unsupported context report schema: {schema!r}")
    ...


def render_markdown(metrics: SummaryMetrics) -> str: ...
def render_html(metrics: SummaryMetrics) -> str: ...  # escape ALL dynamic fields with html.escape()
```

### Required Tests (`tests/test_summary_report.py`)
- Inline `--summary-report` generation in `.md` and `.html` on multi-source merge (verify byte-identical GBFF and unchanged `manifest` + `context` SHA-256 compared to a run without `--summary-report`).
- Standalone `enrich-bakta-report` offline generation from saved manifest + combined context envelope, standalone InterProScan context, and standalone eggNOG context.
- Context SHA-256 mismatch failure (fail closed when `--context-report` bytes disagree with `manifest["metadata"]["context_report_sha256"]`).
- HTML escaping test (inject `<script>alert(1)</script>` into a synthetic AntiFam description or reason and confirm zero unescaped `<script>` tags in HTML output).
- Deterministic byte equality across two runs with identical inputs.
- Collision test (`--summary-report` aliasing an input, output, manifest, or context path fails before writing).

---

## 4. Phase 2 — InterProScan Cross-Format Audit Tool (`enrich-bakta-audit-interproscan`)

### Motivation & Exact Corpus Baselines

The published corpus contains **~950 MB JSON** and **~170 MB GFF3** exports alongside the runtime **~150 MB TSV** for C14, SM, and BK71A. From the audited 5.59-91.0 corpus baseline (`docs/plans/enrich-bakta-interproscan-luna-handoff-2026-10-04.md`):

| Sample | TSV rows | JSON protein objects | JSON multi-xref objects | JSON total query xrefs | Extra GFF3/JSON match occurrences beyond TSV |
|---|---:|---:|---:|---:|---:|
| C14 | 40,234 | 4,634 | 2 | 4,636 | 1,944 |
| SM | 39,684 | 4,573 | 3 | 4,576 | 1,850 |
| BK71A | 32,791 | 4,122 | 7 | 4,131 | 1,968 |

### Format Quirks Confirmed on Published 5.59-91.0 Exports

1. **Analysis name casing/punctuation differs between TSV/GFF3 and JSON:**
   - TSV & GFF3 use `Pfam`, `Gene3D`, `ProSitePatterns`, `MobiDBLite`, `SUPERFAMILY`, `PANTHER`, `PRINTS`.
   - JSON `signature.signatureLibraryRelease.library` uses `PFAM`, `GENE3D`, `PROSITE_PATTERNS`, `MOBIDB_LITE`, `SUPERFAMILY`, `PANTHER`, `PRINTS`.
   - Normalize both sides with: `re.sub(r"[^a-z0-9]+", "", name.lower())` (so `Pfam` ↔ `PFAM`, `ProSitePatterns` ↔ `PROSITE_PATTERNS`, `MobiDBLite` ↔ `MOBIDB_LITE`).
2. **PANTHER signature accession:**
   - In JSON, `match["signature"]["accession"]` is `"PTHR42830"` while `match["accession"]` / `match["model-ac"]` may be `"PTHR42830:SF2"`. Match TSV/GFF3 against `signature["accession"]` (and also accept `model-ac` if a TSV row records a subfamily accession).
3. **JSON structure (`results` array):**
   ```json
   {
     "interproscan-version": "5.59-91.0",
     "results": [
       {
         "sequence": "MKGE...",
         "md5": "d73bd3e8...",
         "matches": [ { "signature": { "accession": "...", "signatureLibraryRelease": { "library": "PFAM", "version": "35.0" } }, "locations": [ { "start": 51, "end": 149 } ] } ],
         "xref": [ { "name": "...", "id": "OJMJKD_006675" } ]
       }
     ]
   }
   ```
   **Never call `json.load()` on the whole file!** Stream the `"results"` array one protein dict at a time using a pure-stdlib chunked `json.JSONDecoder().raw_decode()` iterator so peak Python memory stays under a few tens of MB even on the 951 MB C14 JSON.
4. **GFF3 structure:**
   - Header directives: `##gff-version 3`, `##interproscan-version 5.59-91.0`, `##sequence-region ...`.
   - Stop scanning immediately at `##FASTA` (do not read trailing FASTA sequences as GFF rows).
   - Feature type `protein_match`: `col[0]` = query ID, `col[1]` = analysis, `col[3]` = start (int), `col[4]` = end (int), `col[8]` attributes → parse `Name=<signature_accession>` (URL-unescaped with `urllib.parse.unquote`). Enforce `MAX_RAW_LINE_BYTES = 1 MiB` just like TSV streaming.

### CLI Surface & Output Contract

```text
enrich-bakta-audit-interproscan \
  --tsv SAMPLE.interproscan.tsv \
  --interproscan-version 5.59-91.0 \
  [--interproscan-tsv-layout ipr-go-pathways] \
  [--json SAMPLE.interproscan.json] \
  [--gff3 SAMPLE.interproscan.gff3] \
  [--output AUDIT.json] \
  [--markdown AUDIT.md] \
  [--strict]
```

- Require at least one of `--json` or `--gff3`.
- Emit deterministic JSON report adhering to new checked-in schema `schemas/interproscan-audit.v1.schema.json` (`schema: "enrich-bakta.interproscan-audit.v1"`):
  - `inputs`: paths, sizes, SHA-256 hashes (hashed during the same streaming pass!), asserted version & TSV layout.
  - `version_check`: `asserted_version`, `json_version`, `gff3_version`, `version_mismatch: bool`.
  - `tsv_summary`: `rows`, `unique_queries`, `unique_match_keys`, `occurrences_total`.
  - `json_audit` (when `--json` supplied): `protein_objects`, `total_xrefs`, `multi_xref_proteins` (count + bounded sample list up to 50 with `md5` and sorted `query_ids`), `md5_mismatches` (recalculated `md5(sequence)` vs declared `md5`), `unmatched_queries_count` (queries in JSON absent from TSV), `missing_tsv_occurrences` (count + bounded sample up to 50), `extra_occurrences_by_analysis` (e.g., `{"pirsr": ...}`), `library_versions` (`{library: version}`).
  - `gff3_audit` (when `--gff3` supplied): `polypeptide_count`, `match_rows`, `missing_tsv_occurrences` (count + bounded sample up to 50), `extra_occurrences_by_analysis`.
  - `status`: `"passed"` if `version_mismatch == False`, `md5_mismatches == 0`, and `missing_tsv_occurrences == 0` across all provided formats; otherwise `"failed"`. With `--strict` (default on), exit nonzero (`1`) when `status == "failed"`.

### Rough Pseudocode — Streaming JSON `results` Iterator (`workflows/audit_interproscan.py`)

```python
# src/enrich_bakta_lib/workflows/audit_interproscan.py
import hashlib
import json
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

_ANALYSIS_NORM_RE = re.compile(r"[^a-z0-9]+")


def normalize_analysis_name(name: str) -> str:
    return _ANALYSIS_NORM_RE.sub("", name.strip().lower())


def stream_interproscan_json(
    json_path: Path,
    *,
    chunk_size: int = 1024 * 1024,
) -> tuple[str, str | None, Iterator[dict[str, Any]]]:
    """Stream InterProScan JSON protein result objects with bounded memory.

    Returns a state object or generator yielding one result dict at a time
    while updating a SHA-256 hasher over the exact bytes read.
    """
    # Implementation note for Codex:
    # Return a small context-managed reader class (InterProScanJsonReader) so the
    # caller can iterate `for protein in reader:` and read `reader.sha256` and
    # `reader.producer_version` after EOF without loading the 1 GB file.
    ...


class InterProScanJsonReader:
    """Pure-stdlib incremental reader for InterProScan JSON `results` array."""

    def __init__(self, path: Path, *, chunk_size: int = 1024 * 1024) -> None:
        self.path = path
        self.chunk_size = chunk_size
        self.producer_version: str | None = None
        self.sha256: str = ""

    def __iter__(self) -> Iterator[dict[str, Any]]:
        hasher = hashlib.sha256()
        decoder = json.JSONDecoder()
        with self.path.open("rb") as handle:
            buf = ""
            eof = False

            def refill() -> bool:
                nonlocal buf, eof
                raw = handle.read(self.chunk_size)
                if not raw:
                    eof = True
                    return False
                hasher.update(raw)
                buf += raw.decode("utf-8")
                return True

            # 1. Read prefix until '"results"' array '[' opens, extracting "interproscan-version"
            while '"results"' not in buf and refill():
                pass
            ver_match = re.search(r'"interproscan-version"\s*:\s*"([^"]+)"', buf)
            if ver_match:
                self.producer_version = ver_match.group(1)
            results_idx = buf.find('"results"')
            if results_idx < 0:
                raise MergeError(f"InterProScan JSON {self.path} lacks 'results' array")
            bracket_idx = buf.find("[", results_idx)
            while bracket_idx < 0 and refill():
                bracket_idx = buf.find("[", results_idx)
            if bracket_idx < 0:
                raise MergeError(f"InterProScan JSON {self.path} has malformed 'results' array")
            buf = buf[bracket_idx + 1 :]

            # 2. Incrementally decode each protein object inside [...]
            pos = 0
            while True:
                while True:
                    # Skip whitespace and commas
                    while pos < len(buf) and buf[pos] in " \t\r\n,":
                        pos += 1
                    if pos < len(buf):
                        break
                    buf = ""
                    pos = 0
                    if not refill():
                        raise MergeError(f"InterProScan JSON {self.path} ended unexpectedly")
                if buf[pos] == "]":
                    # Drain any remaining bytes into hasher so SHA-256 covers the full file!
                    while refill():
                        buf = ""
                    self.sha256 = hasher.hexdigest()
                    return
                while True:
                    try:
                        obj, end = decoder.raw_decode(buf, pos)
                        pos = end
                        # Periodically compact buffer so memory stays O(1 protein)
                        if pos > self.chunk_size:
                            buf = buf[pos:]
                            pos = 0
                        yield obj
                        break
                    except json.JSONDecodeError:
                        if eof:
                            raise MergeError(
                                f"InterProScan JSON {self.path} contains malformed JSON"
                            )
                        refill()
```

> **Memory note:** Match occurrence keys `(query, norm_analysis, sig_acc, start, end)` can be tracked as `collections.Counter` of TSV occurrences (~40,000 keys → < 5 MB RAM) and decremented/checked while streaming JSON and GFF3! Do **not** store all JSON/GFF3 matches in a giant set — check each streamed JSON/GFF3 key against the small TSV counter and only count extras per normalized analysis!

### Required Tests (`tests/test_audit_interproscan.py`)
- Synthetic TSV + JSON + GFF3 happy path (0 missing occurrences, version agreement, 0 MD5 mismatches, multi-xref detection, extra `PIRSR` occurrence counted under `extra_occurrences_by_analysis`).
- Version mismatch between `--interproscan-version` and JSON/GFF3 header fails closed (`status == "failed"`, exit code `1`).
- Missing TSV occurrence in JSON or GFF3 fails closed with exact missing key in report.
- Corrupted JSON sequence MD5 (`md5 != md5(sequence)`) fails closed.
- Analysis name normalization (`Pfam` ↔ `PFAM`, `ProSitePatterns` ↔ `PROSITE_PATTERNS`, `MobiDBLite` ↔ `MOBIDB_LITE`) and PANTHER `signature.accession` match.
- Small `chunk_size=64` test proving `InterProScanJsonReader` handles arbitrary mid-token buffer boundaries and hashes the exact full file bytes (`reader.sha256 == sha256_file(path)`).
- Schema validation of output against `schemas/interproscan-audit.v1.schema.json`.
- Real-data smoke verification on C14 (and optional SM/BK71A) confirming the exact baseline numbers (`C14`: 4,634 JSON protein objects, 2 multi-xref, 4,636 xrefs, 0 missing TSV occurrences, 0 MD5 mismatches, 1,944 extra occurrences!).

---

## 5. Phase 3 — Locus-by-Locus Diff Tool (`enrich-bakta-diff`)

### Contract & CLI Surface

```text
# Mode A: Manifest mode (decisions + provenance + suppressions)
enrich-bakta-diff --manifest ENRICHED.manifest.json \
  [--locus LOCUS_TAG ...] [--source SOURCE ...] [--qualifier QUAL ...] \
  [--status STATUS ...] [--format text|json|tsv] [--output PATH] [--fail-on-diff]

# Mode B: Pairwise GBFF mode (direct raw/semantic qualifier comparison)
enrich-bakta-diff --base BAKTA.gbff --enriched ENRICHED.gbff \
  [--locus LOCUS_TAG ...] [--qualifier QUAL ...] \
  [--format text|json|tsv] [--output PATH] [--fail-on-diff]
```

- Enforce mutually exclusive mode selection: either `--manifest` OR (`--base` and `--enriched`).
- Reuse `core/genbank.py:parse_genbank_bytes()` and `strict_parity_check()` in pairwise GBFF mode so sequence/record/coordinate/feature-order mismatches fail loudly rather than producing a bogus diff!
- In pairwise GBFF mode, compute per-feature multiset differences (`collections.Counter(after.values(q)) - collections.Counter(before.values(q))` for added, reverse for removed) across every feature (grouped by `(record_id, locus_tag, feature_type, location_key)`), plus top-level record `COMMENT` additions.
- In manifest mode, group `decisions` by `(record, locus_tag, target_feature_uids)` (mapping UIDs via `entries` when available) and display:
  - Added (`emitted`, `shared_support`) with `source_id`, `qualifier`, `normalized_value`, and `supporting_candidate_ids`.
  - Corroborated existing (`supported_existing`) with `support_witness`.
  - Suppressed/skipped (`suppressed_*`, `filtered_confidence`, `skipped_partial_ec`, `invalid_value`) with `reason_code`.
- `--fail-on-diff`: exits `1` when any matching locus has additions/removals/emitted decisions (useful in CI no-op checks); exits `0` when zero diffs exist.

### Rough Pseudocode — `workflows/diff.py`

```python
# src/enrich_bakta_lib/workflows/diff.py
from __future__ import annotations

import collections
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.genbank import parse_genbank_bytes, strict_parity_check
from enrich_bakta_lib.core.types import MergeError, read_input_bytes, sha256_bytes


@dataclass(frozen=True)
class FeatureQualifierDelta:
    record_id: str
    locus_tag: str
    feature_type: str
    location: str
    added: dict[str, list[str]]
    removed: dict[str, list[str]]


def diff_gbff_files(
    base_path: Path,
    enriched_path: Path,
    *,
    loci: set[str] | None = None,
    qualifiers: set[str] | None = None,
) -> dict[str, Any]:
    base_bytes = read_input_bytes(base_path, "base GBFF")
    enriched_bytes = read_input_bytes(enriched_path, "enriched GBFF")
    base_doc = parse_genbank_bytes(base_bytes, "base GBFF")
    enriched_doc = parse_genbank_bytes(enriched_bytes, "enriched GBFF")
    strict_parity_check(base_doc, enriched_doc)  # guarantees identical record/feature alignment

    deltas: list[FeatureQualifierDelta] = []
    for before_rec, after_rec in zip(base_doc.records, enriched_doc.records, strict=True):
        for before, after in zip(before_rec.features, after_rec.features, strict=True):
            locus = before.locus_tag or ""
            if loci and locus not in loci:
                continue
            keys = sorted(set(before.qualifiers) | set(after.qualifiers))
            added: dict[str, list[str]] = {}
            removed: dict[str, list[str]] = {}
            for key in keys:
                if qualifiers and key not in qualifiers:
                    continue
                b_cnt = collections.Counter(before.values(key))
                a_cnt = collections.Counter(after.values(key))
                inc = sorted((a_cnt - b_cnt).elements())
                dec = sorted((b_cnt - a_cnt).elements())
                if inc:
                    added[key] = inc
                if dec:
                    removed[key] = dec
            if added or removed:
                deltas.append(
                    FeatureQualifierDelta(
                        record_id=before.record_id,
                        locus_tag=locus,
                        feature_type=before.feature_type,
                        location=before.location_key,
                        added=added,
                        removed=removed,
                    )
                )
    return {
        "schema": "enrich-bakta.diff.v1",
        "mode": "gbff",
        "base_sha256": sha256_bytes(base_bytes),
        "enriched_sha256": sha256_bytes(enriched_bytes),
        "changed_features": len(deltas),
        "deltas": [d.__dict__ for d in deltas],
    }


def diff_manifest(
    manifest_path: Path,
    *,
    loci: set[str] | None = None,
    sources: set[str] | None = None,
    qualifiers: set[str] | None = None,
    statuses: set[str] | None = None,
) -> dict[str, Any]:
    ...
```

### Required Tests (`tests/test_diff.py`)
- Pairwise GBFF diff on base vs. enriched GBFF shows exact added qualifiers (`db_xref`, `note`, `gene`, `EC_number`, `inference`) and zero removed qualifiers.
- Pairwise GBFF diff on identical files returns `changed_features == 0` and exits `0` even with `--fail-on-diff`.
- Pairwise GBFF diff fails closed when sequence or feature count differs (`strict_parity_check`).
- Manifest mode diff groups by locus, filters cleanly by `--locus`, `--source`, `--qualifier`, and `--status`, and renders deterministic `text`, `json`, and `tsv`.
- `--fail-on-diff` exits `1` when diffs exist and `0` when filtered diff is empty.

---

## 6. Phase 4 — Producer-Profile Registry, Adapter Onboarding Guide & Optional Promotion Journal

### 6.1 InterProScan Producer-Profile Registry (`sources/interproscan.py`)

Replace hardcoded `SUPPORTED_INTERPROSCAN_VERSIONS = frozenset({"5.59-91.0"})` with a typed profile registry so future verified releases (e.g., `5.6x`/`5.7x` with source-decorated GO syntax like `GO:0005524(InterPro)`) can be added cleanly once real exports are audited, without scattering version `if/else` checks:

```python
@dataclass(frozen=True)
class InterProScanProfile:
    version: str
    supported_layouts: tuple[str, ...]
    allow_decorated_go: bool = False  # False for 5.59-91.0 (bare GO:ddddddd only)
    verified_on_corpus: bool = True


INTERPROSCAN_PROFILES: dict[str, InterProScanProfile] = {
    "5.59-91.0": InterProScanProfile(
        version="5.59-91.0",
        supported_layouts=TSV_LAYOUT_PROFILES,
        allow_decorated_go=False,
        verified_on_corpus=True,
    ),
}
SUPPORTED_INTERPROSCAN_VERSIONS = frozenset(INTERPROSCAN_PROFILES)
```

Keep `SUPPORTED_INTERPROSCAN_VERSIONS` exported for backwards compatibility; continue rejecting any version outside `INTERPROSCAN_PROFILES`.

### 6.2 Opt-In Cross-File Promotion Journal (`core/artifacts.py`)

Currently `stage_artifacts()` stages all `.tmp` files before the first `os.replace()` and documents clearly that multi-file replacement is sequential, not a cross-file transaction. Add an optional, opt-in `promotion_journal_path: Path | None = None` parameter to `stage_artifacts()` (and `--promotion-journal` or automatic crash-safe journal recording when requested):

- Writes a deterministic JSON journal (`enrich-bakta.promotion-journal.v1`) recording planned destinations, expected SHA-256 hashes, pre-existing target SHA-256 (or `null`), and per-file promotion status (`staged` → `promoted` / `failed`).
- If a mid-sequence `os.replace()` fails (e.g., second artifact fails after first was replaced), the journal preserves exact state (`promoted` vs. `unpromoted` files and their hashes) on disk and surfaces the journal path inside `MergeError`, enabling deterministic audit or operator recovery without claiming magical OS-level multi-file atomicity.
- Default `promotion_journal_path=None` preserves 100% existing behavior and zero extra files by default.

### 6.3 Adapter Onboarding Guide (`docs/architecture/adapter-onboarding.md`)

Write a concise (~100 lines) developer guide showing how future adapters (`dbCAN`, `AMRFinderPlus`, `DefenseFinder`, `antiSMASH`, or future eggNOG/InterProScan profiles) plug into:
- `sources/base.py` (`MergeContext`, `SourcePlan`, `SourceAdapter`),
- `sources/value_rules.py` (strict ASCII regex grammars),
- `core/decisions.py` (`candidate_role`, `evidence_class`, exact target UIDs, producer provenance mapping),
- `schemas/merge-manifest.v2.schema.json` (conditional `allOf` entry typing), and
- `workflows/pipeline.py` + `workflows/enrich.py`.

---

## 7. Phase 5 — Packaging, Docs, and `0.5.0` Release Gate

1. Register new console entry points in `pyproject.toml`:
   ```toml
   [project.scripts]
   enrich-bakta = "enrich_bakta_lib.workflows.enrich:main"
   enrich-bakta-baktfold = "enrich_bakta_lib.sources.baktfold:main"
   enrich-bakta-kofam = "enrich_bakta_lib.sources.kofam:main"
   enrich-bakta-eggnog = "enrich_bakta_lib.sources.eggnog:main"
   enrich-bakta-interproscan = "enrich_bakta_lib.sources.interproscan:main"
   enrich-bakta-restore-translations = "enrich_bakta_lib.workflows.restore_translations:main"
   enrich-bakta-report = "enrich_bakta_lib.workflows.summary_report:main"
   enrich-bakta-audit-interproscan = "enrich_bakta_lib.workflows.audit_interproscan:main"
   enrich-bakta-diff = "enrich_bakta_lib.workflows.diff:main"
   ```
   All three new CLIs must support `--version` (`%(prog)s {TOOL_VERSION}`) and `--help`, and be verified inside the isolated wheel test in `.github/workflows/ci.yml`!
2. Bump version to `0.5.0` in `pyproject.toml`, `src/enrich_bakta_lib/core/types.py`, `README.md`, `CHANGELOG.md`, `tests/test_compatibility.py`, `tests/test_merge_pipeline.py`, and review/update `tool.version` (and any `tool_version`-affected `context_report_sha256`) in `docs/data/followthrough-scientific-baseline.json`, `docs/data/interproscan-scientific-baseline.json`, and `tools/run_interproscan_real_checks.py` — **explicitly confirming every C14/SM GBFF `output_sha256` and insertion count stayed 100% identical**!
3. Update `README.md` with executable CLI examples for `--summary-report`, `enrich-bakta-report`, `enrich-bakta-audit-interproscan`, and `enrich-bakta-diff`.
4. Move `docs/plans/post-refactor-features-plan.md` to completed in `docs/plans/README.md` upon completion.

---

## 8. Verification Commands

Run after each phase (fast gate) and in full at Phase 5 closeout (on Windows sandboxed shells, pass `--basetemp .test-output/pytest-tmp` to pytest):

```bash
python -m ruff format --check .
python -m ruff check .
python -m mypy
git diff --check
python -m pytest -q -p no:cacheprovider --basetemp .test-output/pytest-tmp
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
python tools/run_followthrough_real_checks.py
python tools/run_followthrough_chain_checks.py
python tools/run_interproscan_real_checks.py
python -m build
```

Real-data smoke checks for the new tools at Phase 5 (write outputs under `.test-output/v050-smoke/`):

```bash
# 1. Cross-format audit on published C14 InterProScan TSV + JSON + GFF3
enrich-bakta-audit-interproscan \
  --tsv data/C14/evidence/interproscan/C14-NMZ.interproscan.tsv \
  --json data/C14/evidence/interproscan/C14-NMZ.interproscan.json \
  --gff3 data/C14/evidence/interproscan/C14-NMZ.interproscan.gff3 \
  --interproscan-version 5.59-91.0 \
  --output .test-output/v050-smoke/C14-audit.json \
  --markdown .test-output/v050-smoke/C14-audit.md

# 2. Offline summary report from C14 InterProScan manifest + context
enrich-bakta-report \
  --manifest .test-output/interproscan-gates/C14-interpro-enriched.manifest.json \
  --context-report .test-output/interproscan-gates/C14-interpro-context.json \
  --output .test-output/v050-smoke/C14-summary.html

# 3. Pairwise GBFF diff and manifest diff on C14
enrich-bakta-diff \
  --base data/C14/bakta/C14-NMZ.gbff \
  --enriched .test-output/interproscan-gates/C14-interpro-enriched.gbff \
  --format json --output .test-output/v050-smoke/C14-gbff-diff.json
```

---

## 9. Codex Execution Prompt

> Implement the feature plan in `docs/plans/post-refactor-features-plan.md` phase by phase (Phases 0–5), making one focused, well-tested commit per phase. Keep zero new mandatory runtime dependencies (Python stdlib + Biopython only), preserve 100% byte-identical GBFF outputs on C14 and SM across all existing enrichment workflows, and do not fabricate unverified bioinformatics producer profiles or external tool adapters without real fixtures. Verify every phase with `ruff format --check .`, `ruff check .`, `mypy`, `git diff --check`, and `pytest`; at Phase 5 run the full published dataset manifest check, real-data C14/SM/BK71A scientific gates, real C14 `enrich-bakta-audit-interproscan` / `enrich-bakta-report` / `enrich-bakta-diff` smoke runs, and isolated wheel installation checks (`--help`, `--version` for all console scripts). Report exact commit SHAs, test counts, real-data C14 audit metrics (protein objects, multi-xref count, missing TSV occurrences, extra PIRSR occurrences, peak memory/runtime), and any baseline hash changes reviewed for the `0.5.0` bump. Do not push or tag unless explicitly instructed.
