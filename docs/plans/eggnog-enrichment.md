# Audited Implementation Plan: eggNOG-mapper Annotation Enrichment

## Audit verdict

The original plan has the right high-level direction—graft evidence onto pristine
Bakta bytes through the existing insertion engine—but it is not safe to implement
as written. This revision corrects the following issues found against the current
code and the C14/SM inputs:

1. An eggNOG table contains query identifiers but no query protein sequences.
   Accepting only `BAKTA.gbff + *.emapper.annotations` cannot prove that the table
   was generated from the proteins in that GBFF. EggNOG enrichment must require the
   matching Bakta FAA and reuse the existing exact FAA/GBFF translation check.
2. The proposed independent planners all compare only with pristine Bakta. In a
   combined run this would duplicate 95/90 EC assignments and 1,609/1,600 KO
   assignments in C14/SM, respectively.
3. Baktfold and eggNOG propose different `/gene` values at 154 C14 and 152 SM
   blank-gene CDSs. A source order must not silently choose one or write two gene
   qualifiers.
4. The proposed `/inference` value puts the software version where the INSDC
   evidence basis belongs. The evidence basis should identify the seed ortholog;
   software and data versions belong in COMMENT and the manifest.
5. `COG_category` is not consistently a COG identifier in eggNOG-mapper 3.0.0-beta6:
   the real files contain both `COGdddd` values and functional-category letters.
   Only validated COG IDs may become `COG:COGdddd` notes.
6. CAZy cells contain values such as
   `GT2|Glycosyltransferase Family 2.`. A prefix regex such as
   `^([A-Z]+[0-9]+)` is too permissive; the family must be split from its
   description and then fully validated.
7. Gene suffix cleaning is policy, not parsing. Mutating `Preferred_name` inside
   the parser loses the original evidence and makes the manifest ambiguous.
8. The original real-data KO totals are raw eggNOG assignments, not novel
   Bakta-relative assignments. The SM COG baseline is also incorrect.
9. The proposed XLSX route does not account for loss of eggNOG header metadata,
   formulas, numeric coercion, or the currently absent `openpyxl` dependency.
10. The verification plan omits parser-negative tests, field-confidence behavior,
    cross-source conflict tests, mypy, CLI failure cases, and semantic comparison of
    the XLSX/TSV inputs.

## Current repository and data baseline

Audit snapshot: 2026-08-18, branch `main` at `3ec54ea`, matching
`origin/main`.

- Existing tracked implementation: `merge_engine.py`,
  `graft_baktfold_additions.py`, `merge_kofamscan_bakta.py`, and
  `tests/test_merge_pipeline.py`.
- Existing validation passes: 12 pytest tests, Ruff check, Ruff format check, and
  mypy.
- The eggNOG TSV/XLSX/intermediate files and this plan are currently untracked.
  They are inputs or planning artifacts, not implementation files to stage
  accidentally.
- C14 Bakta is one record with 4,636 CDSs; the eggNOG TSV contains 4,563 unique
  annotation rows, all mapping to both the FAA and GBFF. Seventy-three FAA proteins
  have no emitted annotation row.
- SM Bakta is two records with 4,576 CDSs; the eggNOG TSV contains 4,497 unique
  annotation rows, all mapping to both the FAA and GBFF. Seventy-nine FAA proteins
  have no emitted annotation row.
- Both XLSX workbooks contain one `annotations` sheet and are currently cell-for-cell
  identical to the corresponding TSV data rows and header. They do not preserve the
  TSV's `## emapper-3.0.0-beta6` metadata.

### Correct standalone eggNOG baselines

Counts below are novel CDS–value pairs relative to pristine Bakta, before
confidence filtering and before reconciliation with Baktfold/KofamScan.

| Metric | C14 | SM |
|---|---:|---:|
| Blank paired gene/CDS loci with `Preferred_name` | 763 | 737 |
| Novel GO pairs | 9,828 | 9,828 |
| Novel EC pairs | 1,079 | 1,075 |
| Novel KO pairs | 2,541 | 2,524 |
| Novel COG-ID pairs | 766 | 721 |
| Novel CAZy-family pairs | 69 | 69 |

The raw KO totals are 3,996 and 3,943; 1,455 and 1,419 of those pairs are already
represented in Bakta. Gene additions affect both the paired `gene` and `CDS`
features, so 763 loci correspond to 1,526 qualifier insertions in a standalone C14
run.

## Required invariants

1. Pristine Bakta GBFF bytes remain authoritative. Only allowlisted insertions may
   differ in output.
2. Every eggNOG query used for enrichment must occur uniquely in the supplied FAA
   and as one GBFF CDS `locus_tag`, and its normalized FAA sequence must equal the
   CDS `/translation`.
3. Missing eggNOG rows for some FAA proteins are valid. Unknown, duplicate, or
   sequence-mismatched query IDs are fatal.
4. Existing Bakta `/gene`, `/product`, `/translation`, `/protein_id`,
   `/EC_number`, `/db_xref`, `/note`, and `/inference` values are never replaced.
5. EggNOG-derived gene names are written to both the paired `gene` and `CDS` only
   when both features exist, share the same locus/location relationship, and neither
   carries a non-empty Bakta gene name.
6. All sources are parsed and validated before output or manifest promotion.
7. The final insertion list is reconciled across sources before
   `finalize_merge()`. Exact duplicates are emitted once; provenance from every
   supporting source remains in evidence rows.
8. Different proposed gene symbols for the same blank locus are skipped by default
   and reported as conflicts. An explicit CLI policy may prefer one source.
9. Output remains deterministic unless the caller explicitly supplies
   `--merge-timestamp`.
10. Re-running the same source set against its output is byte-idempotent.
11. Biopython parse success proves GenBank syntax, not INSDC submission acceptance.
    KEGG and CAZy xrefs follow this project's enrichment convention; the README
    must not describe the output as submission-validator compliant.

## Annotation and provenance policy

### Supported output fields

| eggNOG field | Validation and normalization | GBFF behavior |
|---|---|---|
| `Preferred_name` | Preserve raw value; optionally strip one terminal `_digits` in planning only | Add `/gene` to the paired gene and CDS only when both are blank |
| `GOs` | Comma-split; require `GO:` plus seven digits | Add novel `/db_xref="GO:ddddddd"` |
| `EC` | Comma-split; strip optional `ec:`; require a fully specified four-part EC | Add novel `/EC_number="X.Y.Z.W"` |
| `KEGG_ko` | Comma-split; strip optional `ko:`; require `Kddddd` | Add novel `/db_xref="KEGG:Kddddd"`, deduplicating Bakta notes/xrefs and other planners |
| `COG_category` | Accept validated `COGdddd` identifiers or category-letter metadata; never reinterpret letters as IDs | Add only novel `/note="COG:COGdddd"` |
| `CAZy` | Split each token once at `|`; require a complete GH/GT/PL/CE/AA/CBM family identifier, including an optional subfamily suffix | Add novel `/db_xref="CAZy:FAMILY"` |
| `PFAMs` | Comma-split; require non-empty whitespace-free tokens and preserve the source token | Add novel CDS `/note="PFAM:TOKEN"`; apply confidence index 12 |
| `eggNOG_OGs` | Comma-split; require non-empty whitespace-free ortholog-group tokens and preserve the source token | Add novel CDS `/note="eggNOG_OG:TOKEN"`; no confidence position |

Do not promote eggNOG pathways, modules, reactions, BRITE, TC, BiGG,
`tax_ceiling`, or donor-lineage fields to feature qualifiers. Preserve them in
`eggnog_context` manifest rows and, when `--context-report` is supplied, a
schema-versioned JSON sidecar grouped by query/CDS.

### Confidence

The parser must validate the 13-character `annotation_confidence` vector and its
documented field order. Map confidence by field:

- `Preferred_name`: index 0
- `GOs`: index 1
- `EC`: index 2
- `KEGG_ko`: index 3
- `KEGG_Pathway`: index 4
- `KEGG_Module`: index 5
- `KEGG_Reaction`: index 6
- `KEGG_rclass`: index 7
- `BRITE`: index 8
- `KEGG_TC`: index 9
- `CAZy`: index 10
- `BiGG_Reaction`: index 11
- `PFAMs`: index 12

`COG_category`, `eggNOG_OGs`, and transfer-provenance fields have no confidence
position and must not be assigned one. Add `--min-eggnog-confidence
{low,medium,high}`, defaulting to `low` so the verified baseline above remains
the default. Every candidate's raw code and emitted/existing/filtered/conflict
status must be recorded in the manifest; context-only values record
`sidecar_only` status and their confidence status.

### Per-feature and record provenance

- When at least one eggNOG value is emitted for a CDS, add one deduplicated
  `/inference`:

  `DESCRIPTION:similar to AA sequence:eggNOG:<seed_ortholog>`

- Keep eggNOG-mapper/data version, table hash, FAA hash, filters, row number,
  E-value, score, raw/normalized value, and field confidence in the manifest.
- Add one source COMMENT per record, keyed by input and FAA hash, containing the
  detected/declared eggNOG version and an evidence-not-phenotype caveat.
- A caller-supplied `--eggnog-version` must match a TSV header version. XLSX input
  has no version header, so `--eggnog-version` is required for XLSX.
- The existing `--no-comment-note` and `--no-feature-provenance` controls remain
  available and do not remove manifest evidence.

The structured inference format follows the current
[NCBI evidence-qualifier guidance](https://www.ncbi.nlm.nih.gov/genbank/evidence/).
The project should also document that `/db_xref` database names are governed by
the [INSDC controlled vocabulary](https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/);
project-specific xrefs may need conversion to notes for formal submission.

## Implementation phases

### Phase 0 — Data hygiene and regression baseline

Modify `.gitignore` to cover the actual private/intermediate inputs:

```gitignore
*-emapper_annotations.xlsx
*.emapper.annotations
*.emapper.hits
*.emapper.seed_orthologs*
*-timing.json
```

Do not delete, move, stage, or commit any existing genomic/eggNOG input. Record the
current test/lint/type-check baseline and keep all generated validation outputs
under `.test-output/`.

Gate:

```powershell
python -m pytest -q --basetemp=".test-output\phase0-tmp" -o cache_dir=".test-output\phase0-cache"
ruff check merge_engine.py graft_baktfold_additions.py merge_kofamscan_bakta.py tests
ruff format --check merge_engine.py graft_baktfold_additions.py merge_kofamscan_bakta.py tests
python -m mypy merge_engine.py graft_baktfold_additions.py merge_kofamscan_bakta.py
git diff --check
```

### Phase 1 — Strict eggNOG input model and parsers

Add `merge_eggnog_bakta.py` with:

- an immutable `EggnogHit` model using tuples rather than mutable lists;
- a strict TSV parser driven by the declared `#query` header, not hard-coded
  positions;
- explicit numeric validation for E-value and score;
- exact validation for unique query IDs, row width, identifier syntax, confidence
  length/codes/order, required fields, UTF-8, and non-empty data;
- separate raw parsing and policy normalization, so suffix cleaning never changes
  `EggnogHit.preferred_name`;
- a path-level XLSX adapter using optional `openpyxl` with
  `read_only=True`, exactly one `annotations` sheet, the same required columns,
  string-preserving cells, and formula rejection;
- a descriptive XLSX dependency error that leaves TSV operation dependency-free;
- version extraction from `## emapper-*` and mismatch detection against the CLI;
- FAA parsing and FAA/GBFF sequence validation reused from a source-neutral helper
  refactored out of `merge_kofamscan_bakta.py`.

Refactor the existing Kofam path to use the shared protein-validation helper without
changing its CLI or output.

Tests:

- minimal valid TSV; CRLF; reordered columns; footer lines;
- malformed/missing/duplicate headers and rows;
- duplicate/unknown query IDs and protein mismatches;
- invalid GO/EC/KO/COG/CAZy/confidence/numeric fields;
- zero annotation rows versus valid partial FAA coverage;
- version detection, explicit-version mismatch, and XLSX version requirement;
- TSV/XLSX semantic equivalence, guarded with the XLSX optional dependency;
- suffix cleaning retains both raw and normalized names.

Gate: existing tests plus the new parser tests, Ruff, format check, and mypy.

### Phase 2 — Standalone eggNOG planner and CLI

Implement `plan_eggnog_additions()` and a standalone CLI:

```text
python merge_eggnog_bakta.py BAKTA.gbff BAKTA.faa \
  query.emapper.annotations OUTPUT.gbff \
  --eggnog-version 3.0.0-beta6 \
  --min-eggnog-confidence low \
  --manifest OUTPUT.manifest.json \
  --context-report OUTPUT.context.json
```

Planner requirements:

- index CDS and paired gene features uniquely and record-locally;
- deduplicate GO, KO, EC, COG, and CAZy values against all relevant Bakta
  representations;
- emit paired gene qualifiers atomically or neither;
- apply optional suffix cleanup only after raw evidence is captured;
- add one seed-specific inference per CDS with at least one emitted value;
- create one evidence record per candidate annotation, including suppression reason;
- add PFAM and eggNOG ortholog-group notes to CDS features with one detailed
  manifest record per normalized value;
- write higher-order context to `enrich-bakta.eggnog-context.v1` JSON when
  `--context-report` is supplied, with one grouped entry per query/CDS;
- record source table/FAA hashes, parser format, version, row counts, mapping counts,
  confidence counts, and addition counts;
- use the existing comment, qualifier, byte-splice, collision, atomic-write, and
  reverse-reconstruction machinery;
- fail before writing if any validation fails.

Tests:

- every supported qualifier and all deduplication representations;
- existing or inconsistent gene qualifiers are protected;
- paired gene/CDS insertion behavior;
- confidence thresholds and COG's no-confidence exception;
- CAZy description stripping with anchored family validation;
- inference uses seed ortholog, not software version;
- PFAM confidence uses the final annotation-confidence position while OG notes
  remain unscored;
- higher-order context is retained in the manifest/sidecar without changing
  the GBFF when the sidecar is omitted;
- no-addition, CRLF, legacy-byte, path-collision, manifest, and idempotence cases;
- reverse-splice byte identity and Biopython parsing.

Real-data dry-run gate:

- reproduce the corrected standalone baselines in the table above;
- assert every annotation query maps to the supplied FAA and GBFF translation;
- independently recount manifest statuses rather than trusting only CLI totals;
- write outputs only beneath `.test-output/`;
- validate outputs with both the merge engine and `gbparse validate`.

### Phase 3 — Cross-source reconciliation and unified orchestrator

Add a reconciliation function, preferably in `merge_engine.py`, that receives all
planned insertions plus source evidence before `finalize_merge()`.

Rules:

1. An exact duplicate qualifier/value for the same feature is emitted once.
2. All supporting sources remain visible in evidence rows and reconciliation
   metadata.
3. Different `/gene` candidates for the same blank locus are handled as a paired
   gene/CDS conflict. Default `--gene-conflict-policy skip` emits neither candidate
   and reports a warning/count. Optional `prefer-eggnog` and `prefer-baktfold`
   policies must be explicit.
4. Different EC or KO values may coexist; only identical values are collapsed.
5. Reconciliation is deterministic and independent of dictionary/set iteration.
6. Source comments and distinct evidence inferences remain source-specific.

Add `enrich_bakta.py` as the canonical multi-source CLI while retaining existing
entry points:

```text
python enrich_bakta.py --bakta BAKTA.gbff --faa BAKTA.faa \
  --baktfold BAKTFOLD.gbff \
  --kofamscan KofamKOALA.txt --kofamscan-version 1.3.0 \
  --eggnog query.emapper.annotations --eggnog-version 3.0.0-beta6 \
  --output ENRICHED.gbff --manifest ENRICHED.manifest.json \
  --context-report ENRICHED.context.json
```

CLI validation:

- require at least one evidence source;
- require `--faa` whenever KofamScan or eggNOG is supplied;
- reject a version flag without its corresponding source;
- validate every input and all cross-source relationships before finalization;
- preserve existing `merge_kofamscan_bakta.py` behavior for backward compatibility;
- give every planner a deterministic, non-overlapping order range, then reconcile.

Synthetic tests must cover all source subsets, exact duplicate evidence,
multi-valued EC/KO coexistence, same-name gene collapse, different-name gene
conflict policies, comment/inference deduplication, manifest source retention, and
combined idempotence.

Real combined-source gate:

- C14: explicitly account for 95 Baktfold/eggNOG EC overlaps, 1,609
  Kofam/eggNOG novel-KO overlaps, and 154 raw gene conflicts;
- SM: explicitly account for 90 EC overlaps, 1,600 KO overlaps, and 152 raw gene
  conflicts;
- prove that no reconciled feature contains duplicate qualifier/value pairs;
- prove reverse-splice reconstruction to pristine Bakta bytes;
- validate C14 and two-record SM outputs with Biopython and `gbparse validate`;
- rerun each combined command against its output and require byte-identical output
  with zero new insertions.

### Phase 4 — Documentation and final validation

Update `README.md` with:

- eggNOG TSV and optional XLSX usage;
- mandatory FAA identity validation;
- supported fields and confidence threshold semantics;
- gene pairing, cleanup, and combined conflict policy;
- xref/submission caveat;
- evidence-level manifest fields and genomic-evidence caveat;
- PFAM/ortholog-group notes and higher-order context sidecar;
- unified and backward-compatible CLI examples.

Run:

```powershell
python -m pytest -q --basetemp=".test-output\final-tmp" -o cache_dir=".test-output\final-cache"
ruff check merge_engine.py graft_baktfold_additions.py merge_kofamscan_bakta.py merge_eggnog_bakta.py enrich_bakta.py tests
ruff format --check merge_engine.py graft_baktfold_additions.py merge_kofamscan_bakta.py merge_eggnog_bakta.py enrich_bakta.py tests
python -m mypy merge_engine.py graft_baktfold_additions.py merge_kofamscan_bakta.py merge_eggnog_bakta.py enrich_bakta.py
python merge_eggnog_bakta.py --help
python enrich_bakta.py --help
git diff --check
git status --short
```

Inspect `git status` and stage only approved code, tests, documentation, and
`.gitignore`. Keep every genomic source, eggNOG output, Kofam table, generated
GBFF, manifest, cache, and `.test-output` artifact untracked/ignored.

## Completion criteria

Implementation is complete only when:

- all existing behavior remains backward compatible;
- both eggNOG formats produce semantically identical plans for the current paired
  inputs;
- every emitted annotation has query/sequence identity and manifest provenance;
- corrected standalone counts are independently reproduced;
- combined overlaps/conflicts are reconciled without duplicate or ambiguous genes;
- C14 and multi-record SM pass syntax, semantic, byte-reconstruction, and
  idempotence checks;
- all unit, integration, Ruff, format, mypy, CLI, and whitespace gates pass.

## Out of scope

- Replacing Bakta products or existing gene names;
- promoting eggNOG pathway/module/reaction/BRITE/TC/BiGG/PFAM fields;
- inferring EC numbers from free-text definitions;
- interpreting annotation as expression, enzyme activity, pathway completeness, or
  phenotype;
- claiming formal GenBank submission validity without an INSDC submission
  validator.
