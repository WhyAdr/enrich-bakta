# Changelog

## Unreleased

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
