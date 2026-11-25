# Deferred Features & Future Development Notes

- **Date:** 2026-11-25
- **Status:** Parked for future releases (post architectural refactoring)
- **Context:** Explicitly deferred out of the [architectural refactoring plan](architecture-refactoring-plan.md) so core module cleanup and workflow unification stay zero-risk and behavior-preserving.

---

## 1. Enrichment Summary Report (`--summary-report`)

Generate a compact, human-readable Markdown or self-contained HTML report derived strictly from the finalized merge manifest (`enrich-bakta.merge-manifest.v2`) and optional context sidecar (`enrich-bakta.context.v1` / source-specific context):

- **Per-source annotation gain:** Added vs. already-present (`supported_existing`) counts broken down by qualifier (`gene`, `EC_number`, `db_xref`, `note`, `inference`) and source (`Baktfold`, `KofamScan`, `eggNOG`, `InterProScan`).
- **Cross-source overlap summary:** Exact shared-support counts (e.g., GO terms corroborated by both eggNOG and InterProScan, KO hits corroborated by KofamScan and eggNOG), with explicit wording noting shared upstream databases are not independent experimental proof.
- **QC & conflict highlights:**
  - AntiFam QC hits (query ID, accession, score, coordinates) from InterProScan context.
  - Suppressed or resolved `gene` symbol conflicts and authoritative base-name collisions.
  - Skipped partial ECs and historical invalid EC tokens (`--baktfold-invalid-ec-policy skip`).
- **Translation lineage table:** Exact restored/imported CDSs carried forward, their origins (`imported_faa` vs `genomically_validated`), and parent manifest bindings.
- **Design rule:** Build this as a pure post-processor consuming manifest + context JSON so it works both inline (`--summary-report PATH`) and offline from previously generated manifests without mutating GBFF bytes.

---

## 2. InterProScan Cross-Format Audit Tool (`audit-interproscan`)

Optional standalone/subcommand audit utility that compares a runtime InterProScan TSV against companion JSON and/or GFF3 exports (automating the manual cross-format check performed during the `0.4.0` audit):

- **Streaming execution:** Stream JSON one protein object at a time (published JSON files reach ~800 MB–1 GB) and stream GFF3/TSV line by line so peak memory remains bounded.
- **Coverage & parity checks:**
  - Confirm every TSV `(query, normalized_analysis, signature_accession, start, end)` occurrence exists in JSON/GFF3.
  - Report analyses present in JSON/GFF3 but omitted by native TSV (e.g., PIRSR in 5.59-91.0) and queries present in JSON with zero TSV hits.
  - Cross-check producer version strings reported in JSON/GFF3 headers against the caller's asserted `--interproscan-version`.
  - Verify recalculated sequence MD5s and flag JSON protein objects shared across multiple query cross-references.
- **Boundary rule:** Remains strictly an audit/QC utility — JSON/GFF3 must never become required runtime inputs or silently supplement TSV enrichment.

---

## 3. Locus-by-Locus Diff Tool (`enrich-bakta diff`)

Machine-readable and human-readable GenBank/manifest diff tool for biological review before INSDC/ENA/NCBI submission:

- **Modes:**
  - **Manifest mode:** Read an `enrich-bakta.merge-manifest.v2` JSON and render a grouped per-locus summary of emitted qualifiers, provenance links (`supporting_candidate_ids`), and suppressed candidates with reason codes.
  - **Pairwise GBFF mode:** Given `BASE.gbff` and `ENRICHED.gbff`, verify record/feature alignment via `core/genbank.py` and emit per-CDS qualifier additions/removals (supporting JSON, TSV, and concise terminal diff views).
- **Use case:** Lets biologists quickly inspect all newly proposed `/gene`, `/EC_number`, `/db_xref`, and `/note` additions on a genome or filter to specific loci/sources without grepping raw multi-megabyte flatfiles.

---

## 4. Future Evidence Sources & Producer Profiles

Once `SourcePlan` / `SourceAdapter` and `workflows/pipeline.py` are in place:

- **Updated producer profiles:**
  - InterProScan `5.6x` / `5.7x` profiles (including tested support for source-decorated GO syntax such as `GO:0005524(InterPro)` if needed).
  - Additional eggNOG-mapper schema versions (e.g., `2.1.x`).
- **New modular bioinformatics adapters:**
  - `dbCAN` (CAZyme substrate/family annotations), `AMRFinderPlus` (antimicrobial resistance determinants), `DefenseFinder` (anti-phage systems), or `antiSMASH` (secondary metabolite clusters) — each plugging into the same byte-preserving splicer, candidate ledger v3, and FAA/GBFF parity checks.
- **Cross-file promotion hardening:**
  - Optional write-ahead journal / recovery protocol for `stage_artifacts()` if atomic multi-file rollback across GBFF + manifest + context sidecar becomes a requirement.
