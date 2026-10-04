# Historical audits

The September 29 audit is preserved verbatim in this folder as
[`enrich-bakta-main-audit-2026-09-29.md`](enrich-bakta-main-audit-2026-09-29.md).
The earlier independent review is preserved as
[`independent-audit-2026-08-25.md`](independent-audit-2026-08-25.md).
The final October 4 reviewer handoff is preserved as
[`enrich-bakta-residual-review-2026-10-04.md`](enrich-bakta-residual-review-2026-10-04.md).

Applicability note: that audit reviewed commit
`5d6754d306e80cddcad7da5fc12b155ab44d2e85` before the local C14, SM, and BK71A
files were inventoried. Its findings remain the historical rationale for the
implementation changes in this checkout; its statements about absent datasets,
old paths, and pre-hardening behavior are not current-state claims. The current
dataset inventory and publication status are documented in
[`../data/README.md`](../data/README.md) and
[`../data/MANIFEST.json`](../data/MANIFEST.json).

The October 4 residual review assessed `84fdc49`. Its three manifest defects,
scientific-baseline requirements, exact imported-query checks, manual
scientific CI run, and version synchronization were resolved by commits
`3eee306`, `e4e79f7`, and `907cccb`. The review remains useful provenance, but
its action items are no longer open against current `main`.
