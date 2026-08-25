# Architectural and Implementation Review: `enrich-bakta`

**Repository:** `WhyAdr/enrich-bakta` (`d:\W\baktfold`)
**Review Date:** 2026-08-25
**Reviewed Commits:** `7f47838` -> `3ec54ea` -> `9ea1058` -> `48edcd4` (HEAD)
**Branch:** `main`, clean working tree, synchronized with `origin/main`

---

## Executive Summary

`enrich-bakta` is a byte-preserving annotation enrichment toolkit for bacterial
genomes annotated by [Bakta](https://github.com/oschwengers/bakta). It grafts
complementary evidence from three upstream sources onto Bakta GenBank flat files
(`.gbff`) without re-serializing the file through Biopython or any other
round-trip parser:

- **Baktfold** — Foldseek-based structural alignments contributing gene symbols,
  EC numbers, and structural cross-references (`afdb_v6:`, `cath:`, `pdb:`).
- **KofamScan / KofamKOALA** — KEGG profile HMM assignments contributing KO
  identifiers and per-hit score provenance.
- **eggNOG-mapper** — Orthology-based functional transfers contributing GO terms,
  EC numbers, KEGG KO, COG IDs, CAZy families, PFAM domains,
  eggNOG ortholog groups, and higher-order metabolic context.

The central design principle is that the original Bakta `.gbff` is treated as an
**immutable raw byte stream**. The merge engine splices new qualifier lines at
computed byte offsets, then proves reversibility by stripping the recorded
insertions and comparing both the raw bytes and SHA-256 digest against the
original. Biopython parses the spliced output as an independent syntax check, but
never participates in serialization.

### Core Design Properties

1. **No re-serialization.** Insertions are offset-based byte splices. Existing
   bytes, line endings (LF or CRLF), comments, and qualifier ordering are never
   touched.
2. **Cryptographic reversibility audit.** After splicing, the engine strips the
   recorded insertions and verifies
   `reconstructed_bytes == base_data and sha256(reconstructed) == sha256(base_data)`.
3. **Pre-write identity verification.** Every evidence query must map to a unique
   FAA protein and a unique GBFF CDS `locus_tag` with an identical normalized
   amino acid sequence. Mismatches abort before any output is created.
4. **Non-destructive annotation policy.** Existing Bakta `/gene`, `/product`,
   `/translation`, `/protein_id`, `/EC_number`, `/db_xref`, `/note`, and
   `/inference` qualifiers are never replaced or removed.
5. **Deterministic output.** No timestamp is generated unless `--merge-timestamp`
   is explicitly supplied, so identical inputs always produce identical output.
6. **Multi-source reconciliation.** The unified CLI (`enrich_bakta.py`) collapses
   exact duplicate qualifiers across sources and resolves gene-name conflicts via
   an explicit policy. Standalone commands bypass reconciliation and operate
   directly on a single evidence source.
7. **Separation of feature-level and higher-order evidence.** Direct CDS
   annotations go into the GBFF; metabolic pathways, BRITE hierarchies, and
   taxonomic context are routed to a deterministic JSON sidecar.

---

## Component Architecture

### Data Flow Overview

```
                             +---------------------+
                             |  Bakta Base (.gbff)  |
                             |  (immutable bytes)   |
                             +----------+----------+
                                        |
       +------------------+-------------+-------------------+
       |                  |                                  |
       v                  v                                  v
+-------------+   +-----------------+              +-------------------+
|  Baktfold   |   | KofamScan/KOALA |              |   eggNOG-mapper   |
|   (.gbff)   |   | (FAA + table)   |              |  (FAA + TSV/XLSX) |
+------+------+   +--------+--------+              +---------+---------+
       |                   |                                  |
       v                   v                                  v
 strict_parity_check  validate_faa_gbff              validate_faa_gbff
       |                   |                                  |
       v                   v                                  v
 plan_baktfold_       plan_kofam_                    plan_eggnog_
 additions            additions                      additions
       |                   |                                  |
       +-------------------+----------------------------------+
                           |
              +------------+-------------+
              |                          |
    (unified CLI only)           (standalone CLIs)
              |                          |
              v                          |
   reconcile_insertions                  |
              |                          |
              +------------+-------------+
                           |
                           v
                  apply_insertions
                  (byte splice + reversibility audit)
                           |
                           v
                  validate_genbank_semantics
                  (Biopython parse check)
                           |
                           v
                     atomic_write
                  (tempfile + fsync + rename)
                           |
              +------------+------------+
              |                         |
              v                         v
    Enriched Output GBFF      Manifest + Sidecar
    (byte-preserved)          (JSON or TSV)
```

> **Important distinction:** `reconcile_insertions` is called **only** by the
> unified multi-source CLI (`enrich_bakta.py`). The three standalone commands
> (`graft_baktfold_additions.py`, `merge_kofamscan_bakta.py`,
> `merge_eggnog_bakta.py`) pass their insertions directly to `apply_insertions`
> via `finalize_merge` without cross-source reconciliation.

---

### Module Inventory

| Module | Purpose | Key Exports |
|---|---|---|
| [`merge_engine.py`](file:///d:/W/baktfold/merge_engine.py) | Foundation: raw-byte GenBank parser, offset-based splicer, reversibility auditor, FAA/GBFF validator, cross-source reconciler, manifest writer, and atomic file I/O. | `RawDocument`, `RawRecord`, `RawFeature`, `Insertion`, `AppliedInsertion`, `parse_genbank_bytes`, `validate_genbank_semantics`, `apply_insertions`, `strict_parity_check`, `validate_faa_gbff`, `parse_faa`, `cds_by_locus`, `reconcile_insertions`, `finalize_merge`, `atomic_write`, `write_manifest`, `write_json_sidecar` |
| [`graft_baktfold_additions.py`](file:///d:/W/baktfold/graft_baktfold_additions.py) | Grafts Baktfold structural xrefs, gene symbols, and EC numbers onto Bakta. | `plan_baktfold_additions`, `graft`, `detect_baktfold_version` |
| [`merge_kofamscan_bakta.py`](file:///d:/W/baktfold/merge_kofamscan_bakta.py) | Merges KofamScan/KOALA profile hits. Supports optional `--baktfold` for one-pass combined merge. | `KofamHit`, `parse_kofam_table`, `plan_kofam_additions`, `merge` |
| [`merge_eggnog_bakta.py`](file:///d:/W/baktfold/merge_eggnog_bakta.py) | Merges eggNOG-mapper evidence (TSV or optional XLSX). Supports `--context-report` sidecar. | `EggnogHit`, `EggnogTable`, `parse_eggnog_tsv`, `parse_eggnog_xlsx`, `parse_eggnog_path`, `plan_eggnog_additions`, `merge` |
| [`enrich_bakta.py`](file:///d:/W/baktfold/enrich_bakta.py) | Unified multi-source CLI: reconciles all three sources in one pass. | `enrich`, `main` |
| [`restore_bakta_translations.py`](file:///d:/W/baktfold/restore_bakta_translations.py) | Restores `/translation` for translationless pseudogene CDSs from a matched Bakta FAA. | `plan_translation_restoration`, `restore`, `prepare_restoration` |
| [`normalize_baktfold.py`](file:///d:/W/baktfold/normalize_baktfold.py) | Legacy normalization utility. Attempted to restore Bakta qualifier formatting (e.g., `/transl_table`, `/protein_id` namespace, `/inference` strings, COMMENT blocks, LOCUS dates) from Baktfold-processed files. **Not imported by any other module; superseded by the byte-preserving insertion approach.** | `normalize`, `normalize_with_overrides`, `count_features` |

---

## Merge Engine Internals (`merge_engine.py`)

### Raw-Byte Document Model

The engine parses GenBank bytes into offset-tracked structures without
normalizing content:

- **`RawFeature`**: Captures `record_index`, `record_id`, `feature_type`,
  `location` text, byte `start`/`end` offsets, and parsed `qualifiers` dict.
- **`RawRecord`**: Captures `record_id`, `accessions`, `declared_length`,
  `topology`, byte `start`/`end`, `features_offset` (byte offset of the
  `FEATURES` line), and uppercase-normalized `sequence` bytes (extracted from
  `ORIGIN`).
- **`RawDocument`**: Wraps the complete original `data: bytes` alongside the
  parsed `RawRecord` list.

The parser validates declared sequence length against actual ORIGIN content and
rejects files with no LOCUS records or missing FEATURES tables.

### Insertion Data Model

```python
@dataclass
class Insertion:
    offset: int          # Byte offset in original file (insertion point)
    payload: bytes       # Formatted qualifier lines (indented, with matching newlines)
    record: str          # Target record ID
    feature_type: str    # Feature type (CDS, gene, etc.)
    locus_tag: str       # Locus tag
    qualifier: str       # Qualifier key (e.g., "db_xref", "note")
    value: str           # Qualifier value
    source: str          # Provenance label (e.g., "eggNOG", "KofamScan")
    source_value: str    # Raw upstream value
    order: int = 0       # Deterministic ordering tie-breaker
```

Multiple insertions at the same byte offset are supported: after splicing,
`cursor` advances to `insertion.offset` (not past the payload), so the next
insertion at the same offset copies zero original bytes and appends its payload
immediately after.

### Qualifier Formatting

[`format_qualifier`](file:///d:/W/baktfold/merge_engine.py#L543-L564) wraps
qualifier lines to 80 columns at the standard 21-character GenBank indent, with
two notable behaviors:

- **`/inference` is never wrapped.** The INSDC inference format contains
  structured colons and version strings that would be corrupted by mid-line
  breaks.
- **Newline detection** uses
  [`newline_for_offset`](file:///d:/W/baktfold/merge_engine.py#L532-L540),
  which inspects the bytes immediately preceding the insertion point to determine
  whether to emit `\r\n` or `\n`. This is how CRLF and LF styles are preserved
  per-location within a file.

### Qualifier Allowlist

[`apply_insertions`](file:///d:/W/baktfold/merge_engine.py#L640-L695) enforces
a hard-coded allowlist of qualifier keys that may be inserted:
`gene`, `EC_number`, `db_xref`, `note`, `inference`, and `COMMENT`. An
`extra_allowed_qualifiers` parameter extends this set (used by
`restore_bakta_translations.py` to add `translation`). Any insertion with a
qualifier outside this set raises `MergeError`.

### Reversibility Audit

After splicing, the engine performs an in-memory reverse-reconstruction:

1. Iterate over applied insertions in order.
2. Verify that `output[output_offset : output_offset + len(payload)] == payload`.
3. Copy the inter-insertion segments into a reconstruction buffer.
4. Assert `bytes(reconstructed) == base_data` **and**
   `sha256(reconstructed) == sha256(base_data)`.

If either check fails, the merge aborts with `MergeError` before any file is
written.

### Biopython Semantic Validation

[`validate_genbank_semantics`](file:///d:/W/baktfold/merge_engine.py#L292-L306)
parses the spliced output through `Bio.SeqIO.parse(..., "genbank")` as a syntax
check. It first attempts UTF-8 decoding; if that fails, it falls back to
**Latin-1** to provide a lossless one-codepoint-per-byte view for files
containing legacy non-ASCII bytes in comments (the test suite exercises this path
with embedded `\xff` bytes).

### Atomic Output

[`atomic_write`](file:///d:/W/baktfold/merge_engine.py#L712-L728) creates a
temporary file in the same directory as the target, writes + `fsync`s, then uses
`os.replace` for an atomic rename. Path collision detection
([`paths_collide`](file:///d:/W/baktfold/merge_engine.py#L698-L709)) prevents
the output from overwriting any input, manifest, or sidecar path.

---

## Annotation Source Details

### Baktfold (`graft_baktfold_additions.py`)

**Validation:** [`strict_parity_check`](file:///d:/W/baktfold/merge_engine.py#L337-L428)
requires exact match of record IDs, accessions, declared lengths, topologies,
source sequence SHA-256s, ordered feature-type counts, per-feature
`(type, locus_tags, location)` tuples, and unique feature keys. CDS translation
differences are recorded as diagnostics but do not block the merge; only the
structural checks are hard gates.

**Evidence categories and provenance:**

| Category | Source qualifier | Emitted qualifier | Provenance |
|---|---|---|---|
| Structural xrefs | `/db_xref` with exact prefix in `{afdb_v6, cath, pdb}` | `/db_xref="<prefix>:<id>"` | `/inference="protein structure similarity:Baktfold Foldseek:<version>"` |
| Gene symbols | `/gene` (only when Bakta has blank `/gene` on both paired `gene` and `CDS` features) | `/gene="<symbol>"` on both features | `/note="Baktfold gene-symbol evidence:v<version>"` |
| EC numbers | `/EC_number` and `/db_xref="EC:..."` (both representations) | `/EC_number="<value>"` | `/note="Baktfold functional-annotation evidence:v<version>"` |

Key policies:
- A non-empty Bakta `/gene` is **never** overwritten.
- Ambiguous Baktfold gene symbols (multiple distinct values on one feature)
  raise `MergeError` and abort the merge.
- EC values are collected from both `/EC_number` and `/db_xref="EC:..."` for
  deduplication. Partial/wildcard ECs (e.g., `1.2.3.-`) are accepted — the
  Baktfold graft does not require fully-specified four-part numeric EC numbers.
- The `--force` flag is accepted by the CLI but immediately errors with
  `"--force is no longer supported"`. Failed identity checks cannot be bypassed.

### KofamScan / KofamKOALA (`merge_kofamscan_bakta.py`)

**Parsing:** The table parser is header-driven (requires the standard
`# gene name  KO  thrshld  score  E-value  KO definition` header line). Only
rows with a leading `*` (above-threshold assignments) are retained. Numeric
fields (threshold, score, E-value) are validated as finite `Decimal` values.

**Protein integrity:** Every hit query ID must exist exactly once in the FAA and
map to a single GBFF CDS `locus_tag`. The normalized FAA sequence must equal the
CDS `/translation`. Missing hit IDs in either the FAA or GBFF are fatal.
Non-hit FAA proteins (those with no Kofam row) are expected and valid.

**Emitted qualifiers per hit:**

| Qualifier | Value | Condition |
|---|---|---|
| `/db_xref` | `"KEGG:Kxxxxx"` | KO not already present in Bakta `/db_xref` or `/note` |
| `/note` | `"KofamScan:Kxxxxx;threshold=...;score=...;E-value=..."` | Always (per-hit provenance) |
| `/inference` | `"profile:KofamScan[:VERSION]"` | Once per CDS, not per hit |

Key policies:
- Multiple KO assignments per CDS are all retained.
- `[EC:...]` text in KO definitions is **never** promoted to `/EC_number`.
- KO definitions are **never** promoted to `/product`.
- Existing Bakta qualifiers are never replaced.
- The `--baktfold` flag enables a one-pass combined Baktfold + KofamScan merge
  through a single `finalize_merge` call.

### eggNOG-mapper (`merge_eggnog_bakta.py`)

**Parsing:** Two input formats are supported:

- **TSV** (dependency-free): Header-driven via the `#query` column header.
  Extracts `## emapper-VERSION` metadata. Validates unique query IDs, column
  width consistency, required columns, UTF-8 encoding, and non-empty data.
- **XLSX** (requires optional `openpyxl`): Must contain exactly one sheet named
  `annotations` with no formulas. Requires an explicit `--eggnog-version` because
  the `## emapper-*` metadata is not preserved in XLSX. If `openpyxl` is not
  installed, a descriptive error is raised before attempting import.

**Confidence filtering:** The 13-character `annotation_confidence` vector is
validated to contain only `l`/`m`/`h`/`-` codes. Each feature-level field maps
to a specific confidence position:

| Index | Field | Index | Field |
|---|---|---|---|
| 0 | `Preferred_name` | 7 | `KEGG_rclass` |
| 1 | `GOs` | 8 | `BRITE` |
| 2 | `EC` | 9 | `KEGG_TC` |
| 3 | `KEGG_ko` | 10 | `CAZy` |
| 4 | `KEGG_Pathway` | 11 | `BiGG_Reaction` |
| 5 | `KEGG_Module` | 12 | `PFAMs` |
| 6 | `KEGG_Reaction` | | |

Fields without a confidence position (`COG_category`, `eggNOG_OGs`) are never
confidence-filtered. The `--min-eggnog-confidence` threshold accepts `low`
(default), `medium`, or `high`.

**Feature-level qualifier mapping:**

| eggNOG field | Qualifier | Confidence index | Notes |
|---|---|---|---|
| `Preferred_name` | `/gene` | 0 | Only on blank, location-matched gene/CDS pairs. `--clean-gene-suffix` strips one terminal `_digits` from the planned value while preserving the raw value in the manifest. |
| `GOs` | `/db_xref="GO:ddddddd"` | 1 | Comma-split; validated as `GO:` + seven digits. |
| `EC` | `/EC_number="X.Y.Z.W"` | 2 | Comma-split; optional `ec:` prefix stripped; requires fully-specified four-part numeric EC (unlike Baktfold, partial ECs are rejected). |
| `KEGG_ko` | `/db_xref="KEGG:Kddddd"` | 3 | Comma-split; optional `ko:` prefix stripped; deduplicated against existing Bakta notes/xrefs. |
| `COG_category` | `/note="COG:COGdddd"` | None | Only validated `COG\d{4}` identifiers are emitted. Alphabetic functional-category letters (e.g., `J`, `KL`) are recognized as valid COG_category values but are **not** converted to COG IDs. |
| `CAZy` | `/db_xref="CAZy:FAMILY"` | 10 | Each token is split at `\|` to strip description text; the family code (`GH`, `GT`, `PL`, `CE`, `AA`, `CBM` + digits, optional `_subfamily`) is then fully validated. |
| `PFAMs` | `/note="PFAM:TOKEN"` | 12 | Comma-split; validated as non-empty whitespace-free tokens. |
| `eggNOG_OGs` | `/note="eggNOG_OG:TOKEN"` | None | Comma-split; no confidence filtering. |

**Inference:** One `/inference="DESCRIPTION:similar to AA sequence:eggNOG:<seed_ortholog>"`
per CDS with at least one emitted value.

**Higher-order context sidecar:** The following fields are **never** injected as
feature qualifiers. When `--context-report` is supplied, they are written to a
deterministic JSON sidecar (schema `enrich-bakta.eggnog-context.v1`) grouped by
query/CDS. They are also recorded as `sidecar_only` entries in the manifest:

- `Description` (eggNOG functional description)
- `KEGG_Pathway`, `KEGG_Module`, `KEGG_Reaction`, `KEGG_rclass`
- `BRITE`, `KEGG_TC`, `BiGG_Reaction`
- `tax_ceiling`, `max_annot_lvl`
- `farthest_donor_taxid`, `farthest_donor_lineage`

---

## Multi-Source Reconciliation (`enrich_bakta.py`)

The unified CLI assigns non-overlapping order tiers to each source's insertions:

| Source | Starting order | Typical range |
|---|---|---|
| Baktfold | `0` | 0 .. 999,999 |
| KofamScan | `1,000,000` | 1,000,000 .. 1,999,999 |
| eggNOG | `2,000,000` | 2,000,000 .. 2,999,999 |

[`reconcile_insertions`](file:///d:/W/baktfold/merge_engine.py#L745-L846)
then applies three rules:

1. **Exact duplicate collapse.** If multiple sources propose identical
   `(record, feature_type, locus_tag, qualifier, value)` tuples, the qualifier
   is emitted once. All supporting sources are recorded in a
   `exact_duplicate_collapsed` manifest entry.
2. **Multi-value coexistence.** Different EC numbers or KO assignments from
   different sources for the same CDS coexist as separate qualifiers.
3. **Gene symbol conflict resolution.** If Baktfold and eggNOG propose different
   gene names for the same blank locus:
   - `--gene-conflict-policy skip` (default): Neither is emitted;
     `gene_conflict_skipped` is recorded.
   - `--gene-conflict-policy prefer-eggnog`: Emits the eggNOG preferred name.
   - `--gene-conflict-policy prefer-baktfold`: Emits the Baktfold gene symbol.

The unified CLI also requires:
- At least one evidence source.
- `--faa` whenever KofamScan or eggNOG is supplied.
- Version flags only with their corresponding source.
- `--context-report` only with `--eggnog`.

---

## Translationless Pseudogene Repair (`restore_bakta_translations.py`)

Bakta's FAA can retain predicted protein sequences for pseudogene CDSs whose
GBFF features carry `/pseudogene="..."` instead of `/translation`. Because the
eggNOG merge mandates exact FAA/GBFF translation identity, these pseudogenes
would fail validation without prior restoration.

`restore_bakta_translations.py` provides a controlled, two-step workflow:

1. **Validate** (`--dry-run`): Confirms every eggNOG query maps to the FAA and
   GBFF, that targetted CDSs are marked `/pseudogene` and lack `/translation`,
   and that CDSs which already have `/translation` are not touched.
2. **Restore**: Splices `/translation="<FAA_sequence>"` into an intermediate
   `RESTORED.gbff` via the standard `finalize_merge` path (using
   `extra_allowed_qualifiers=("translation",)` to extend the allowlist). The
   original `BAKTA.gbff` is left untouched. Each restoration is recorded with
   the protein SHA-256 in the manifest.

---

## Git History & Evolution

The repository was built in a single day (2026-08-18) across four commits:

```
48edcd4  Add PFAM and eggNOG context evidence
9ea1058  Add eggNOG enrichment pipeline
3ec54ea  Implement byte-preserving KofamScan and Baktfold merge engine with comprehensive test suite
7f47838  Initial commit: enrich-bakta architecture, grafting pipeline, and merge roadmap
```

### Commit 1: `7f47838` — Initial Architecture

- Created `graft_baktfold_additions.py` (586 lines) as the first annotation
  grafting script.
- Created `normalize_baktfold.py` (441 lines) — an attempt to restore Bakta
  qualifier formatting from Baktfold-processed files by line-level text
  transformations. This approach was fragile (required exact knowledge of
  Bakta's COMMENT layout, qualifier spelling, and protein_id namespace) and was
  superseded by the byte-preserving insertion approach in the next commit.
- Created `merge-koala-baktfold-plan.md` (283 lines) documenting the audit
  findings, the structural parity requirements, and the planned KofamScan merger.

### Commit 2: `3ec54ea` — Byte-Preserving Core Engine

- **Created `merge_engine.py`** (739 lines) with raw-byte parsing,
  offset-based qualifier splicing, SHA-256 reversibility audit, FAA/GBFF
  sequence validation, atomic file I/O, and TSV/JSON manifest writing.
- **Rewrote `graft_baktfold_additions.py`** to delegate all splicing and
  validation to `merge_engine.py`. Removed the earlier line-level approach.
- **Created `merge_kofamscan_bakta.py`** (509 lines) with header-driven Kofam
  table parsing, per-hit evidence tracking, and optional `--baktfold` one-pass
  integration.
- **Created `tests/test_merge_pipeline.py`** with synthetic parity, preservation,
  qualifier, KO, failure, CRLF, and idempotence tests.

### Commit 3: `9ea1058` — eggNOG Pipeline & Unified Orchestrator

- **Created `merge_eggnog_bakta.py`** (663 lines) with TSV and optional XLSX
  parsing, 13-character confidence vector enforcement, paired gene/CDS emission,
  and standalone CLI.
- **Created `enrich_bakta.py`** (206 lines) as the canonical multi-source CLI.
- **Created `restore_bakta_translations.py`** (216 lines) for pseudogene
  translation restoration.
- **Extended `merge_engine.py`** (+233 lines) with `validate_faa_gbff`,
  `cds_by_locus`, `reconcile_insertions`, `write_json_sidecar`, and refactored
  shared protein validation out of `merge_kofamscan_bakta.py`.
- **Refactored `merge_kofamscan_bakta.py`** (-97 lines) to use the shared
  protein validation helper.
- Extended test suite with eggNOG parser, provenance, confidence, existing-value,
  reconciliation, unified merge, and pseudogene restoration tests.

### Commit 4: `48edcd4` — PFAM, eggNOG OGs & Context Sidecar

- **Extended `merge_eggnog_bakta.py`** (+207 lines) with `PFAMs` mapped to
  `/note="PFAM:..."` (confidence index 12), `eggNOG_OGs` mapped to
  `/note="eggNOG_OG:..."` (no confidence position), and higher-order context
  collection (`_collect_context`) producing per-query grouped JSON sidecar
  entries.
- **Extended `merge_engine.py`** (+39 lines) with `write_json_sidecar` and
  sidecar path collision checks in `finalize_merge`.
- **Extended `enrich_bakta.py`** (+13 lines) to pass `--context-report` through
  to the eggNOG planner and tag the sidecar with `unified-enrichment` operation.
- Added PFAM/OG note tests, confidence-position-12 test, and context sidecar
  tests to the test suite.

---

## Workspace Data Inventory

### Tracked Files (committed to `main`)

| Category | Files |
|---|---|
| Core pipeline | `merge_engine.py`, `graft_baktfold_additions.py`, `merge_kofamscan_bakta.py`, `merge_eggnog_bakta.py`, `enrich_bakta.py`, `restore_bakta_translations.py` |
| Legacy | `normalize_baktfold.py` (not imported anywhere) |
| Tests | `tests/test_merge_pipeline.py` |
| Documentation | `README.md`, `merge-koala-baktfold-plan.md`, `integrate-eggnog-plan.md` |
| Configuration | `.gitignore` |

### Untracked Data (present in workspace, excluded by `.gitignore`)

Three bacterial isolate datasets are present locally:

| Dataset | Organism / Description | Records | Bakta CDSs | Files |
|---|---|---|---|---|
| **C14** | Single circular chromosome | 1 | 4,636 | `C14-NMZ_bakta/` (Bakta outputs), `C14-baktfold/` (Baktfold outputs), `C14-KofamKOALA.txt`, `C14-NMZ-query.emapper.annotations`, `C14-NMZ-emapper_annotations.xlsx` |
| **SM** | Two circular replicons | 2 | 4,576 | `SM-NMZ_bakta/` (Bakta outputs), `SM-baktfold/` (Baktfold outputs), `SM-KofamKOALA.txt`, `SM-NMZ-query.emapper.annotations`, `SM-NMZ-emapper_annotations.xlsx` |
| **BK71A** | *Bacillus halotolerans* chromosome (4,133,448 bp) | 1 | — | `BK71A-restored.gbff` (a single restored GBFF; appears to be output from `normalize_baktfold.py`) |

Additionally, enriched eggNOG outputs already exist in the Bakta directories:
- `C14-NMZ_bakta/C14-NMZ-eggnog-enriched.gbff` (~14 MB) with manifest (~47 MB)
- `C14-NMZ_bakta/C14-NMZ-eggnog-restored.gbff` (~13 MB) with manifest
- `SM-NMZ_bakta/SM-NMZ-eggnog-enriched.gbff` (~14 MB) with manifest (~47 MB)
- `SM-NMZ_bakta/SM-NMZ-eggnog-restored.gbff` (~13 MB) with manifest

### Enrichment Baselines (from `integrate-eggnog-plan.md`)

These standalone counts are relative to pristine Bakta, before confidence
filtering and before cross-source reconciliation:

| Metric | C14 | SM |
|---|---:|---:|
| eggNOG annotation rows | 4,563 | 4,497 |
| Blank paired gene loci with `Preferred_name` | 763 | 737 |
| Novel GO pairs | 9,828 | 9,828 |
| Novel EC pairs | 1,079 | 1,075 |
| Novel KO pairs | 2,541 | 2,524 |
| Novel COG-ID pairs | 766 | 721 |
| Novel CAZy-family pairs | 69 | 69 |

From `merge-koala-baktfold-plan.md`:

| Metric | C14 | SM |
|---|---:|---:|
| Baktfold matched features | 9,586 | 9,464 |
| Baktfold new genes | 800 | 780 |
| Baktfold new ECs | 178 | 168 |
| Baktfold new structural xrefs | 58 | 98 |
| KofamScan leading-`*` hits | 3,397 | 3,333 |
| KofamScan hit CDSs | 3,181 | 3,120 |
| KofamScan new gene-KO pairs | 1,916 | 1,900 |
| KofamScan already-present Bakta KOs | 1,481 | 1,433 |

Cross-source overlaps (from `integrate-eggnog-plan.md`):

| Overlap | C14 | SM |
|---|---:|---:|
| Baktfold/eggNOG EC duplicates | 95 | 90 |
| KofamScan/eggNOG novel-KO duplicates | 1,609 | 1,600 |
| Baktfold/eggNOG gene symbol conflicts | 154 | 152 |

> **Note:** These numbers are sourced from the plan documents, not independently
> re-verified against the data files during this review.

---

## Verification & Quality Status

### Test Suite

All **21 pytest items** pass (20 test functions, one parametrized with 2 cases):

```
tests/test_merge_pipeline.py .....................   [100%]
21 passed in 0.99s
```

Coverage areas:
- Structural parity rejection (shifted coordinates, swapped records, duplicate
  locus tags, mismatched feature counts)
- Baktfold prefix whitelisting (rejects `pdbx:`, `cathode:`,
  `afdb_v6_extra:` while accepting `pdb:`, `cath:`, `afdb_v6:`)
- Ambiguous Baktfold gene symbol rejection
- KofamScan multi-hit retention and existing-KO deduplication
- CRLF and LF preservation; legacy Latin-1 byte passthrough
- Empty-addition byte identity and final-newline style preservation
- Output path collision detection and missing query ID rejection
- Numeric validation (KO syntax, threshold, E-value) and qualifier allowlist
- GenBank qualifier wrapping and Biopython round-trip parse
- eggNOG header-driven parsing, gene suffix cleaning, and inference provenance
- PFAM note emission and eggNOG OG notes; manifest `feature_note` evidence class
- PFAM confidence using position 12; OG notes unscored
- eggNOG field validation (malformed GO, version mismatch)
- XLSX version requirement (fails before optional `openpyxl` import)
- Confidence filtering and existing-value deduplication; COG no-confidence exception
- Cross-source reconciliation (exact duplicate collapse, gene conflict policies)
- Unified KofamScan + eggNOG merge with duplicate KO collapse
- Pseudogene translation restoration (limited to translationless `/pseudogene` CDSs)

### Static Analysis

```
Ruff check:   All checks passed (7 files)
Ruff format:  7 files already formatted
Mypy:         Success, no issues found in 6 source files
```

No explicit ruff or mypy configuration files exist in the repository; all three
tools run with their default settings.

---

## CLI Quick Reference

All CLIs support `--json` for machine-readable output. Error exit code is 2.

```bash
# Unified multi-source enrichment
python enrich_bakta.py \
  --bakta BAKTA.gbff --faa BAKTA.faa \
  --baktfold BAKTFOLD.gbff \
  --kofamscan KofamKOALA.txt --kofamscan-version 1.3.0 \
  --eggnog query.emapper.annotations --eggnog-version 3.0.0-beta6 \
  --min-eggnog-confidence low \
  --gene-conflict-policy skip \
  --output ENRICHED.gbff --manifest ENRICHED.manifest.json \
  --context-report ENRICHED.context.json

# Standalone Baktfold graft
python graft_baktfold_additions.py \
  BAKTA.gbff BAKTFOLD.gbff OUTPUT.gbff \
  --manifest OUTPUT.manifest.json

# Standalone KofamScan merge (with optional --baktfold)
python merge_kofamscan_bakta.py \
  BAKTA.gbff BAKTA.faa KofamKOALA.txt OUTPUT.gbff \
  --baktfold BAKTFOLD.gbff \
  --kofamscan-version 1.3.0 \
  --manifest OUTPUT.manifest.json

# Standalone eggNOG merge
python merge_eggnog_bakta.py \
  BAKTA.gbff BAKTA.faa query.emapper.annotations OUTPUT.gbff \
  --eggnog-version 3.0.0-beta6 \
  --min-eggnog-confidence low \
  --clean-gene-suffix \
  --manifest OUTPUT.manifest.json \
  --context-report OUTPUT.context.json

# Pseudogene translation restoration (dry-run, then write)
python restore_bakta_translations.py \
  BAKTA.gbff BAKTA.faa query.emapper.annotations RESTORED.gbff \
  --eggnog-version 3.0.0-beta6 --dry-run
python restore_bakta_translations.py \
  BAKTA.gbff BAKTA.faa query.emapper.annotations RESTORED.gbff \
  --eggnog-version 3.0.0-beta6 --manifest RESTORED.manifest.json
```

Shared flags across CLIs:
- `--no-comment-note` — omit record-level `COMMENT` provenance blocks.
- `--no-feature-provenance` — omit per-feature `/inference` and `/note`
  provenance qualifiers.
- `--merge-timestamp VALUE` — inject an explicit timestamp (otherwise omitted
  for deterministic output).

---

## Scientific Guardrails

The pipeline embeds two explicit caveats:

1. **Genomic evidence ≠ biological phenotype.** KO assignments, Foldseek
   structural hits, and orthology-based functional transfers are computational
   predictions. They do not demonstrate in vivo expression, enzyme activity,
   pathway completeness, or phenotype. This caveat is injected into every
   enriched file's `COMMENT` block.

2. **Enrichment conventions ≠ INSDC submission compliance.** The generated
   GenBank files parse cleanly via Biopython, but database names in `/db_xref`
   qualifiers (e.g., `KEGG:`, `CAZy:`) follow this project's enrichment
   conventions. Formal submission to NCBI/ENA/DDBJ requires the current
   [INSDC db_xref controlled vocabulary](https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/);
   non-INSDC xrefs would need conversion to `/note` qualifiers.

---

## Observations & Potential Improvements

### Strengths

- The byte-preservation approach with cryptographic reversibility is well
  engineered and avoids the entire class of round-trip serialization bugs.
- The pre-write validation strategy means failed runs never leave partial
  outputs.
- Deterministic output (no auto-timestamps) enables reproducible builds and
  idempotent reruns.
- The manifest records every candidate annotation including suppressed and
  filtered values, providing a complete audit trail.
- Clean separation between direct annotation (GBFF) and higher-order context
  (JSON sidecar) avoids qualifier pollution.

### Potential Improvements

- **No `pyproject.toml` or packaging.** The codebase operates as loose scripts
  with direct module imports. Tests must be run via `python -m pytest` (not
  plain `pytest`) to get the root directory on `sys.path`. A minimal
  `pyproject.toml` would formalize dependencies (`biopython`, optional
  `openpyxl`), enable `pip install -e .`, and fix the import path issue.
- **Single test file.** All 21 tests live in one 612-line file. As the suite
  grows, splitting by module (e.g., `test_baktfold.py`, `test_kofam.py`,
  `test_eggnog.py`, `test_reconciliation.py`) would improve navigability.
- **No coverage measurement.** Adding `--cov` to the pytest invocation would
  quantify which code paths lack test coverage.
- **GFF3 enrichment.** The engine operates exclusively on GenBank flat files.
  Extending to Bakta's `.gff3` output would support downstream visualization
  in genome browsers.
- **Batch processing.** A wrapper to automatically discover and process
  sample-paired inputs across a cohort directory would reduce manual CLI
  invocations.
