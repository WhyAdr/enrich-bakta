# Changelog

## 0.4.0 - 2026-10-04

- Added InterProScan annotation enrichment source supporting verified producer
  version 5.59-91.0 with profile-driven TSV layout validation (`ipr-go-pathways`,
  `ipr-go`, `ipr-pathways`, `ipr-only`).
- Implemented single-pass streaming TSV parser enforcing a 1 MiB line ceiling,
  captured input SHA-256, strict placeholder coherence, AntiFam QC tracking,
  and deterministic pathway and support digests.
- Added strict query sequence verification against Bakta matched FAA and GenBank
  CDS features with exact MD5 and length cross-validation.
- Implemented member-database note promotion (`Pfam`, `TIGRFAM`) with version-aware
  witness matching against existing annotations (e.g. bare `PF02566` matches `PF02566.20`).
- Integrated InterProScan functional proposals (`InterPro`, `GO`) and substantive
  member evidence into candidate ledger v3 and merge manifest v2.
- Added multi-source context envelope (`enrich-bakta.context.v1`) combining
  eggNOG and InterProScan sidecars, embedding finalized sidecar SHA-256 into manifest metadata.
- Added standalone compatibility launcher `merge_interproscan_bakta.py` and
  integrated unified `--interproscan` workflow options.
- Validated scientific acceptance on published C14 and SM datasets, failing closed
  on missing translations and legacy unverified provenance.

## 0.3.1 - 2026-10-04

- Hardened complete-CDS translation validation and preserved imported/restored
  translation lineage through every standalone and unified workflow.
- Added exhaustive Baktfold decisions, typed evidence roles, class-specific
  provenance support, per-node reconciliation reasons, and final-output checks.
- Restored published eggNOG XLSX compatibility and retained substantive Kofam
  hit evidence for existing KO annotations.
- Added a checked-in merge-manifest schema, dependency-floor coverage, and
  separate dataset-integrity and scientific-rerun CI gates.
- Synchronized source-entry projections with final decisions, retained planned
  reasons separately, and counted every suppressed candidate class.
- Added a reviewed C14/SM scientific baseline that fails closed on input,
  producer, policy, output, count, translation-lineage, or manifest drift.

## 0.3.0 - 2026-10-03

- Reorganized the implementation into the `enrich_bakta_lib` package while
  retaining tested root compatibility launchers and the legacy normalizer.
- Added semantic GenBank identity checks, translation-evidence provenance,
  schema-aware eggNOG parsing, typed candidate decisions, and staged artifacts.
- Added curated C14, SM, and BK71A publication manifests and Git LFS policy.
- Added the GPL-3.0-or-later code license and release validation documentation.

Functional outputs that use `import-faa` remain explicitly imported protein
evidence; they are not claims of independent genomic translation validation.
