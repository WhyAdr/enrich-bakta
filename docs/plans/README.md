# Plans and roadmap

## Active plans and deferred notes

- [Architectural refactoring implementation plan](architecture-refactoring-plan.md) — ready-to-execute handoff covering `core/merge_engine.py` decomposition, redundant Biopython parse removal, uniform `SourcePlan` / `SourceAdapter` contracts, and workflow pipeline unification.
- [Deferred features and future development notes](deferred-features-notes.md) — parked ideas (`--summary-report`, InterProScan TSV/JSON/GFF3 cross-audit tool, `enrich-bakta diff`, and future evidence adapters).

## Historical plans

These source-integration and release plans are retained for provenance. Their completed
parts are reflected in the current package and release documentation; they are
not an alternate execution path.

- [Implementation plan (2026-09-30)](enrich-bakta-implementation-plan-2026-09-30.md)
- [Refined handoff (2026-10-01)](enrich-bakta-refined-handoff-2026-10-01.md)
- [eggNOG enrichment plan](eggnog-enrichment.md)
- [Kofam/Baktfold merge plan](kofam-baktfold-merge.md)
- [Post-push follow-through hardening plan](enrich-bakta-followthrough-hardening-2026-10-04.md)
- [InterProScan enrichment plan](interproscan-enrichment.md)
- [InterProScan Luna handoff (2026-10-04)](enrich-bakta-interproscan-luna-handoff-2026-10-04.md)

The October 4 follow-through plan reviewed `cc98e4f` and supplied the design
and test evidence implemented by commits `a69b4e3` through `84fdc49`. Its
remaining residual items were closed by the October 4 residual review work
through `907cccb`; it is retained as a historical execution record, not as an
open task list.

The InterProScan enrichment plan (`interproscan-enrichment.md`) outlines the
design, baselines, and implementation roadmap for integrating InterProScan 5.59-91.0
evidence datasets onto Bakta GBFF records.
