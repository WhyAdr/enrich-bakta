# Adding an evidence source

Enrichment adapters turn a versioned producer export into deterministic
GenBank insertion candidates and traceable evidence rows. Producer-specific
parsing and planning live in [`sources/`](../../src/enrich_bakta_lib/sources).
Use a small native fixture for normal tests and review real exports before
registering a producer version or offering a new runtime adapter.

## Start with the producer export

1. Record the producer and database versions, export options, native layout,
   and fixture origin. Keep the source bytes or an authorized test fixture
   with its SHA-256. If only a reformatted export is available, establish
   its own format contract and identify it as such.
2. Establish which identifiers map to a unique existing feature, which
   coordinates describe the evidence, and which fields can support a
   qualifier. Check multiple rows per query, repeated values, partial rows,
   placeholders, and conflicting annotations on real exports.
3. Validate the declared layout and reject unverified versions before
   parsing. Hash the bytes that were actually parsed. Bound line and
   aggregate sizes, distinguish invalid tokens from valid missing values,
   and report errors with row or query information.

InterProScan's `INTERPROSCAN_PROFILES` contains the verified `5.59-91.0`
native TSV profile. Each profile identifies its accepted layouts. Add a
new registry entry only after real exports verify the version, layout,
identifier behavior, and allowable field syntax. Version-specific syntax
such as source-decorated GO terms requires fixture-backed normalization and
negative tests. JSON and GFF3 belong to the standalone audit, while
InterProScan enrichment consumes TSV. A producer version change in eggNOG
also requires fixtures for its own supported export layouts.

## Implement the adapter protocol

Implement `SourceAdapter.plan(ctx, *, starting_order=0)` from
[`sources/base.py`](../../src/enrich_bakta_lib/sources/base.py). Use the base
bytes, parsed features, and available FAA and translation evidence from
`MergeContext`. Use exact query identifiers and the existing feature UID
helper. Reject ambiguous target mappings, and retain imported translation
lineage and its validation status.

Return a `SourcePlan` containing the stable source identifier, source hash,
candidate insertions, evidence rows, counts, source-specific manifest
metadata, and every input path. An adapter that provides contextual
evidence may also return a context key and payload. Order insertion proposals
from `starting_order` and sort evidence collections wherever the input order
does not carry meaning. Keep the planner free of file writes.

Use qualifier grammars in
[`sources/value_rules.py`](../../src/enrich_bakta_lib/sources/value_rules.py).
Keep accession and qualifier token validation ASCII-strict. Add grammars
only for producer fields whose semantics and canonical form have been
verified. Preserve the original field in evidence rows when the promoted
value is normalized.

## Fit the candidate ledger and schema

Review [`core/decisions.py`](../../src/enrich_bakta_lib/core/decisions.py)
when adding a source identifier or evidence row kind. Ledger candidates
need exact target UIDs, source and input hash, normalized qualifiers and
values, an explicit `candidate_role`, and a compatible `evidence_class`.
Identify producer provenance through the existing insertion role machinery.
When an annotation supports a value already in the GenBank file, preserve its
support witness. Record skipped or suppressed evidence with a reason and
regression-test that it produces no unexpected insertion.

Update the appropriate conditional `allOf` entry definition in the
[v2 merge manifest schema](../../schemas/merge-manifest.v2.schema.json)
for a new entry type. Verify emitted and non-emitted evidence rows and
keep fields and meanings expected by existing manifest consumers.

## Integrate through the shared pipeline

Wire the adapter and producer-specific options in
[`workflows/enrich.py`](../../src/enrich_bakta_lib/workflows/enrich.py).
If a source-specific command is needed, route it through the same adapter
and [`workflows/pipeline.py`](../../src/enrich_bakta_lib/workflows/pipeline.py)
for reconciliation, ledger validation, semantic checks, and artifact
staging. Keep the existing ordering of adapters stable and test deliberate
conflicts with sources that may emit the same qualifier.

Potential sources such as dbCAN, AMRFinderPlus, DefenseFinder, and antiSMASH
require their own fixture and evidence review before any adapter is
registered. A new adapter should pass tests for valid and malformed exports,
unverified versions, duplicate and unmatched identifiers, supported and
suppressed candidates, and single-source and unified runs. Run the
published scientific gates and compare existing GBFF hashes and insertion
counts. Describe the scope and limits of every new annotation claim in
the user-facing documentation.
