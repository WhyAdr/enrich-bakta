# Plans and roadmap

## Active plans and future notes

- [Deferred features and future development notes](deferred-features-notes.md) — higher-level future backlog and onboarding scope for new external evidence sources (`dbCAN`, `AMRFinderPlus`, `DefenseFinder`, `antiSMASH`, future InterProScan/eggNOG versions).

## Completed and historical plans

These architecture, source-integration, and release plans are retained for provenance.
Their completed parts are reflected in the current package and release documentation;
they are not an alternate execution path.

- [Post-refactor features implementation plan](post-refactor-features-plan.md) — **completed** for `0.5.0` with enrichment summaries, InterProScan cross-format audit, locus diff, producer registry, and opt-in promotion journal.
- [Architectural refactoring implementation plan](architecture-refactoring-plan.md) — **completed** in commits `a5eead4` through `5ebee37` (`core/merge_engine.py` split, `6 -> 2` GenBank parse reduction, `MergeContext`/`SourcePlan`/`SourceAdapter` standardization, and `workflows/pipeline.py` unification).
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
