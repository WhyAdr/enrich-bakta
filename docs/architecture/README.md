# Architecture and design history

The [shared pipeline](../../src/enrich_bakta_lib/workflows/pipeline.py) captures
one input snapshot, runs source adapters, reconciles candidate decisions, and
passes finalized output to the artifact writer. The [adapter onboarding guide](adapter-onboarding.md)
describes the evidence, fixture, and validation contracts for new sources and
producer profiles.

The [report, audit, and diff commands](../../README.md#enrichment-reports-and-audits)
give callers separate views of enrichment output, producer-format parity, and
locus-level changes.

## Historical review

These documents record design and implementation decisions from before the
`enrich-bakta` package reorganization. They remain useful historical context,
but the canonical interfaces are the modules under `src/enrich_bakta_lib/` and
the current root `README.md`.

- [Bakta integration implementation review](implementation-review-2026-08-25.md)
