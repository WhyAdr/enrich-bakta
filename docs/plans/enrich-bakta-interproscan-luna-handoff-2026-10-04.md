# Sol → Luna: InterProScan feature execution brief

Prepared for Wahyu on 2026-10-04 following an independent audit of remote main
at e5e0499fc602c58be55264d75d8ed8e92a972636 and all nine published InterProScan
exports. The local repository plan was refined at
`docs/plans/interproscan-enrichment.md`; no feature implementation, remote push,
release tag or CI dispatch was performed during this audit.

The execution specification below is self-contained. Luna can begin from the
pinned main and apply its requirements directly. The critical corrections are
explicit native TSV layouts, shared planner/ledger Pfam support, bounded
snapshot parsing, exact query identity, known-legacy lineage protection, and
correct finalized context hashing. Validation completed: 146 tests, quality
checks, 42 published artifacts, all nine export hashes and 22 existing
scientific acceptance checks. New InterProScan outputs still require review
before their scientific baselines are pinned.

---

# InterProScan enrichment: audited implementation plan and Luna handoff

- **Audit date:** 2026-10-04
- **Repository:** https://github.com/WhyAdr/enrich-bakta
- **Audited main:** `e5e0499fc602c58be55264d75d8ed8e92a972636`
- **Package baseline:** 0.3.1; merge-manifest v2; candidate-ledger v3
- **Status:** implementation specification; feature implementation and new output goldens pending
- **Execution owner:** Luna

This replaces the earlier plan at this path. Implement against the audited
main, preserving its reconciled candidate ledger, strict translation-evidence
lifecycle, byte-preserving merge engine, XLSX header aliases, and scientific
acceptance gates. Earlier Opus/follow-through records describe older revisions;
do not apply their patches verbatim. The tracked historical records remain
useful regression context.

## 1. Audit verdict and corrections

The TSV-first approach and conservative qualifier policy are sound. The
following corrections are necessary before implementation:

| Finding | Required correction |
|---|---|
| Native 14-column TSV is ambiguous | Select an explicit version/layout profile; column 14 can be GO **or pathways**, depending on producer flags. Never guess from width. |
| Versioned Bakta Pfam notes suppress a bare proposal, but the current ledger validator requires exact note equality | Share a narrowly defined member-note equivalence function between planning and semantic validation, with an actual base-value witness. |
| Published TSVs are 132–157 MB; a single line reaches 93,639 bytes and pathway occurrences reach millions | Stream captured binary bytes, hash the same snapshot, bound aggregation and row size, and summarize pathway context. |
| BK71A has legacy restoration provenance, not the current translation-evidence marker | Add an explicit known-legacy provenance preflight for InterProScan. Existing marker detection alone cannot enforce this case. |
| Context is bound late in finalization | Serialize and hash the finalized sidecar before serializing the manifest; do not hash its preliminary payload. |
| JSON merges some identical proteins into multiple query cross-references | Identity is exact query-to-CDS UID, never MD5-to-one-CDS. |
| The prior acceptance commands included the 86-file local inventory | Clean-clone acceptance uses the 42-file PUBLISHED-MANIFEST; the local inventory is optional when all its assets exist. |
| The plan named a nonexistent core/value_rules.py | The current shared module is src/enrich_bakta_lib/sources/value_rules.py. |

These are implementation requirements, not claims that these capabilities
already exist. No InterProScan runtime module or CLI exists at the audited main.

## 2. Review evidence and limits

Checks independently completed on the audited revision:

- 146 tests passed; Ruff check, Ruff format (38 files), mypy (14 source
  modules), and git diff --check passed.
- Existing real-data acceptance scripts passed all 22 checks: 13 in
  run_followthrough_real_checks.py and 9 in run_followthrough_chain_checks.py.
  Pinned first-pass hashes, TSV/XLSX equivalence, no-op outputs, final candidate
  projections, imported-origin opt-in and full origin carry-forward passed.
- Git LFS materialized the published inputs; all 42 published artifacts passed
  tools/validate_dataset_manifest.py, including all nine InterProScan exports.
- All TSVs were streamed for row shape, query identity, coordinates,
  MD5/length consistency, token grammars, pair counts, and source SHA-256.
- GFF3 and JSON were independently inspected, with JSON streamed one protein
  object at a time. Every TSV match occurrence was found in both other formats
  using (query, normalized analysis name, signature accession, start, end).
  This establishes identity/location coverage, not equality of every
  annotation field across formats.
- The native 5.59-91.0 TSV writer was inspected to resolve optional-column
  behavior. GFF3 and JSON independently report version 5.59-91.0.
- Canonical genbank-parser 0.9.3 at
  09bdef65771ae292b7cf49f4b39f648719cc10df was used for GenBank review.
  Pristine C14 and SM each have one informational
  EXCEPTIONAL_TRANSLATION_ANNOTATION diagnostic, at OJMJKD_000205 and
  NNLMDO_000200 respectively; BK71A's historical target has none. These
  diagnostics are not proof of independent genomic translation validation.
- Current main's push CI succeeded (run 37204044368). The manual/tag-only
  dataset-bytes and scientific-reruns jobs were not executed by this review.

The prior follow-through fix is already present, as are durable scientific
baseline assertions and the 0.3.1 release changes. Keep these gates. No new
InterProScan output hashes or insertion counts can be accepted before the
feature exists and its outputs are reviewed.

### 2.1 Corpus baseline

Counts of annotations below use unique (query ID, canonical accession) pairs,
unless explicitly labelled occurrences. They are parser/identity baselines,
not predicted insertion totals.

| Metric | C14 | SM | BK71A |
|---|---:|---:|---:|
| Target CDS features | 4,636 | 4,576 | 4,131 |
| TSV rows | 40,234 | 39,684 | 32,791 |
| 13-column rows, without integrated InterPro entry | 12,814 | 12,640 | 10,431 |
| 15-column rows, with integrated InterPro entry | 27,420 | 27,044 | 22,360 |
| Unique TSV query CDSs | 4,457 | 4,381 | 3,834 |
| Target CDSs unrepresented in TSV | 179 | 195 | 297 |
| Query–InterPro pairs | 16,160 | 15,930 | 12,994 |
| Query–GO pairs | 9,407 | 9,332 | 6,880 |
| Query–Pfam pairs | 6,219 | 6,088 | 5,131 |
| Query–TIGRFAM pairs | 1,766 | 1,758 | 1,454 |
| Query–CDD pairs, context only | 2,591 | 2,571 | 2,183 |
| Pfam pairs already in base Bakta notes, ignoring accession version | 5 | 15 | 9 |
| GO pairs already in target Bakta db_xref | 3,858 | 3,857 | 3,680 |
| GO pairs absent from target Bakta db_xref | 5,549 | 5,475 | 3,200 |
| GO pairs also present in eggNOG | 2,063 | 2,054 | not established |
| Reactome token occurrences | 4,534,168 | 4,416,062 | 3,717,389 |
| MetaCyc token occurrences | 2,441,941 | 2,410,944 | 2,177,198 |
| Maximum raw TSV line bytes | 93,639 | 93,639 | 93,638 |
| Additional GFF3/JSON match occurrences beyond TSV | 1,944 | 1,850 | 1,968 |
| JSON query cross-references, including unmatched queries | 4,636 | 4,576 | 4,131 |
| JSON protein result objects | 4,634 | 4,573 | 4,122 |
| JSON protein objects with multiple query cross-references | 2 | 3 | 7 |

Native TSV deliberately omits some analyses, including PIRSR in this producer
version. Do not promise complete JSON/GFF3 match coverage from TSV or silently
supplement the runtime TSV with either format. Published JSON can include
queries with no exported TSV match.

Zero TSV identity/location occurrences were missing from GFF3/JSON; JSON
sequence MD5 recalculation found zero mismatches. C14 and SM have zero TSV-to-FAA
MD5 or length mismatches. Among query CDSs carrying exactly one target
translation, TSV-to-translation comparisons also found zero mismatches.
Missing translations are unresolved cases, not zero-mismatch validations.

TSV SHA-256 values (full file sizes/hashes for all formats remain in
docs/data/PUBLISHED-MANIFEST.json):

| Sample | Bytes | SHA-256 |
|---|---:|---|
| C14 | 157031166 | 5de249ac97ed468bdae75c5ab483d91f5d2f613e7dbbc2e922ee263ec441ac6b |
| SM | 153607485 | 3cb5f0e8b6723faae26f450c4cc288526452cd8e9483e7aa0920029f2d3f29fe |
| BK71A | 131729054 | 0b44aab2a511ebec53f1fd78ae004674ad060c3acb44a2a31864bdbfd974aa21 |

Pathway counts are repeated token occurrences, not unique pathways or evidence
of genome-level pathway completeness. Labels may refer to other taxa; do not
interpret a Reactome taxon label as a bacterial pathway assignment.

## 3. Runtime input and parser contract

### 3.1 Explicit producer profile

Initial verified support is native InterProScan **5.59-91.0**, with InterPro
lookup enabled. Require --interproscan-version; reject unverified versions
until a separately tested profile is added. Record the supplied version as a
producer assertion: TSV alone cannot authenticate the producing executable.

Add --interproscan-tsv-layout with the following profiles. Default to
ipr-go-pathways for the published corpus, and record the selected layout.

| Layout | Integrated row | Unintegrated row | Optional interpretation |
|---|---|---|---|
| ipr-go-pathways | 15 fields | 13 fields | 14=GO; 15=pathways |
| ipr-go | 14 fields | 13 fields | 14=GO |
| ipr-pathways | 14 fields | 13 fields | 14=pathways |
| ipr-only | 13 fields | 13 fields | no GO or pathway field |

Unintegrated native rows end with two "-" InterPro placeholders. Require
integrated rows to have the selected expanded shape, including empty optional
fields if the writer emitted them. Validate InterPro accession/description
placeholder coherence. A valid native 11-column export with lookup disabled is
outside this first-release profile: diagnose unsupported layout, rather than
calling the producer output intrinsically malformed. Do not infer a layout
from the first row or from whether a token happens to resemble a GO term.

The first eleven fields are query ID, MD5, protein length, analysis, signature
accession, signature description, match start, match end, score, status, and
date. Fields 12–13 are InterPro accession and description. TSV MD5 is column 2.

### 3.2 Streaming validation and source binding

- Read binary lines with an explicit limit; use a documented 1 MiB raw-line
  ceiling initially, with a contextual error if exceeded. Do not call
  read_input_bytes() for the TSV or build a raw-row list.
- Hash exactly the bytes being parsed. Use a captured/spooled snapshot or a
  single binary pass with bounded aggregation. Two independent reads for
  hashing and parsing without a snapshot contract are insufficient.
- Decode UTF-8 strictly; support LF, CRLF, and a final line without newline.
  Strip only line endings; preserve trailing tab-delimited empty fields.
  Reject embedded NUL/control characters incompatible with the field grammar.
- Require nonempty query/analysis/signature fields, 32 hexadecimal MD5,
  positive protein length, integer inclusive
  1 <= start <= end <= protein length, and repeated-query MD5/length agreement.
- Missing score "-" is valid. Preserve method-specific scores and status/date
  as audit values; do not impose one E-value threshold on heterogeneous methods.
  A producer "T" is not experimental validation.
- For the verified native profile, accept bare GO:ddddddd tokens separated by
  "|", and absent/empty fields. Source-decorated GO syntax is deferred until
  its producer profile and complete grammar are tested. Never strip arbitrary
  suffixes to make an invalid token look valid.
- Full-match IPRdddddd, PFddddd (with optional positive numeric version for
  comparison), and TIGRddddd grammars only. CDD and other context-only
  signatures must be retained as raw context without imposing an unapproved
  promotion grammar on otherwise valid rows.
- A malformed row fails before any artifact is committed. Report source,
  physical row, query if known, field, and stable reason code. Avoid printing
  a 90 KB pathway field in an error.
- Aggregate by exact query, output field, and canonical value. Retain source
  occurrence count, first supporting row, compact row/location samples, and a
  deterministic support digest; document sample truncation. Do not create
  millions of identical GO candidates or retain every repeated raw pathway.
- Bound query, candidate, and diagnostic aggregation with documented limits;
  fail explicitly on overflow rather than silently truncate annotations.
  Pathway context uses streaming counters and order-defined digests. Unique
  pathway counts are optional only with a declared bounded/disk-backed method.

Include TSV in other_inputs for collision protection. This does **not** hash
it: the adapter/workflow must explicitly put its captured SHA-256 into
source_hashes, InterProScan metadata, COMMENT when enabled, and context. Use the same source
snapshot throughout planning. JSON/GFF3 are audit assets, not runtime inputs;
do not load a near-GB JSON file just to obtain the producer version.

## 4. Identity and translation lineage

Require the matched FAA on every InterProScan functional merge. Use existing
parse_faa(), normalize_protein(), parse_genbank_bytes(), cds_by_locus(), and
the translation-evidence loader rather than a competing identity parser.

Each TSV query must map to exactly one CDS locus tag and one stable target
feature UID. Every candidate carries that UID; disable query/value fallback
matching for this source. Reject unknown/ambiguous locus tags, missing FAA
queries, inconsistent repeated identities, and mismatched MD5 or length.

Check normalized FAA and encoded /translation against both TSV MD5 and length.
Use the existing normalization policy, including its terminal-stop handling;
do not silently alter internal residues to force equality. Identical sequence
MD5s at different loci remain separate targets. C14's OJMJKD_001145 and
OJMJKD_022100 are a real regression example.

Equality with an encoded /translation is protein parity. It does not by itself
establish a complete genomic CDS model or prove pristine annotation history.
Reuse complete-CDS validated-only restoration semantics (cds=True), including
failures for partial, fuzzy, remote, exceptional, out-of-bounds, and unknown
translation-table models. Imported FAA equality must never be labelled genomic
validation.

### 4.1 Corpus-specific acceptance routes

| Sample | Functional route |
|---|---|
| C14 | Pristine data/C14/bakta/C14-NMZ.gbff and matched .faa are eligible for direct InterProScan acceptance. |
| SM | Pristine GBFF must fail on five queried CDSs without translations. Restore using the existing import-faa workflow and matched eggNOG input; accept subsequent InterProScan enrichment only with the bound manifest and explicit imported-origin opt-in. |
| BK71A | Historical restored GBFF is parser/audit material only. FAA and a bound restoration ledger are absent; functional acceptance remains blocked. |

SM unresolved TSV queries are exactly:
NNLMDO_002705, NNLMDO_007000, NNLMDO_007005, NNLMDO_013725,
NNLMDO_023610. The existing restoration baseline carries **six** SM imported
origins overall. Retain all six through InterProScan, including NNLMDO_011850,
which is not represented in its queried set.

BK71A has 20 queried CDSs without encoded translations:
BCCDEH_00094, BCCDEH_00095, BCCDEH_00135, BCCDEH_00142,
BCCDEH_01028, BCCDEH_01029, BCCDEH_01567, BCCDEH_01601,
BCCDEH_01613, BCCDEH_02181, BCCDEH_02307, BCCDEH_03097,
BCCDEH_03101, BCCDEH_03102, BCCDEH_03191, BCCDEH_03458,
BCCDEH_03690, BCCDEH_03692, BCCDEH_03700, BCCDEH_03701.

Its COMMENT records normalization by normalize_baktfold.py with restored Bakta
provenance and preserved Baktfold v0.1.0 history. It does **not** contain the
current ##enrich-bakta:translation-evidence:v1## marker.
has_translation_evidence_marker() therefore does not protect this historical
case. Add a narrow known-legacy restoration classifier to the InterProScan
preflight, with synthetic tests, and fail closed when it cannot establish the
existing bound lineage contract. Do not use filenames as the classifier or
claim to detect arbitrary unmarked historical edits.

Do not invent BK71A FAA/ledger from JSON sequences or GBFF translations.
Recovery requires genuine appropriate source inputs and the established
restoration process; an invalid/missing ledger cannot be waived.

### 4.2 Existing lifecycle to preserve

Use --translation-evidence-manifest and --allow-imported-translations.
No --allow-unverified-lineage flag or second restoration ledger.

Carry every verified origin entry through every downstream merge, including
untouched CDSs, and bind entries to actual output SHA-256. Hash the parent
manifest bytes actually loaded. Require its parent evidence on restoration
reruns, validate strict enums/booleans/digests/duplicate metadata and entry
agreement, and preserve original origins on no-op reruns. Marked restored input
must never finalize with an empty/dropped ledger. Imported origins remain
unavailable to functional sources without explicit opt-in.

## 5. Scientific promotion policy

| Evidence | Default qualifier / handling | Candidate role and evidence class |
|---|---|---|
| InterPro accession | db_xref=InterPro:IPRdddddd | functional_proposal; InterPro |
| GO term | db_xref=GO:ddddddd | functional_proposal; GO |
| Pfam signature | note=PFAM:PFddddd, canonical base accession | substantive_evidence; member:Pfam |
| TIGRFAM signature | note=TIGRFAM:TIGRddddd | substantive_evidence; member:TIGRFAM |
| CDD and other member signatures | context only; no initial promotion switch | no functional candidate |
| Reactome/MetaCyc and other pathway labels | context sidecar only | no functional candidate |
| AntiFam | explicit QC context only | no functional candidate |
| InterProScan inference | inference=protein motif:InterProScan:5.59-91.0 | producer_provenance; InterProScan:producer |
| Source COMMENT | deterministic source/version/hash policy marker | context |

Both Pfam and TIGRFAM are enabled by default through
--interproscan-member-dbs Pfam,TIGRFAM. An explicit empty allowlist disables
member-note promotion; reject unsupported requested names. This allowlist
controls member notes, not which otherwise valid integrated analyses can
contribute InterPro or GO evidence. CDD promotion is deferred, not an
undocumented opt-in.

Do not change product, gene, EC_number, pseudo, coordinates, sequence, or
translation from this source. Do not add GO ancestors, ontology inference,
live database validation, or pathway completeness calls. Retain scores,
descriptions, and producer dates as context rather than qualifier prose.
Shared eggNOG/InterPro support can share upstream databases; it is not
independent experimental confirmation.

AntiFam ANF00260 occurs at OJMJKD_003040 and NNLMDO_003025 (score 3.4E-6).
Expose these as gene-model QC flags. Do not automatically mark pseudo, delete
the CDS, or veto other analyses. Any future veto requires a separate reviewed
policy.

### 5.1 Pfam semantic support: planner AND validator

Current validate_candidate_ledger() uses exact note membership. A planner
marking PFAM:PF00001 as supported_existing against PFAM:PF00001.33 fails
the final-output check. Fix this as part of the feature, not as a post-release
patch.

Define one shared affirmative member-note parser/equivalence rule used by the
adapter, ledger base support, and ledger output support. Preserve original
Bakta bytes; record the actual existing value as a support witness. Accept
only explicitly structured standalone PFAM notes and their valid versions.
Reject prefix collisions, negative/free-text mentions, and malformed suffixes.
Do not loosen all note membership or insertion-value checks.

First-release semantic equivalence applies to standalone member notes.
Composite eggNOG PFAMs: notes remain under their existing behavior; do not
pretend generic reconciliation already provides equivalence to those notes.
Tests must document this scope. Any later extension to composite notes needs
an exact affirmative grammar and shared planner/validator tests.

## 6. Adapter, ledger and reconciliation

Create src/enrich_bakta_lib/sources/interproscan.py and a typed
InterProScanPlan analogous to EggnogPlan, containing insertions, evidence_rows,
stats, context_report, and the captured source hash/profile. A hit iterator
consumes the validated bounded binary snapshot; immutable hit records must
not leak a raw giant field into each retained candidate.

Source adapters return plans and rows, never construct CandidateDecision
directly. Emit interproscan_candidate source rows for proposed and
supported-existing outcomes; retain reasons and compact supporting facts for
suppressed/context outcomes in their appropriate schema.

Reuse the real insertion APIs:

```python
qualifier_insertion(
    base.data, feature, qualifier, value, source, source_value, order,
    candidate_role=..., evidence_class=...
)
comment_insertion(base, record, marker, lines, source, order)
```

Use monotonic source insertion order, reconcile_insertions(),
build_candidate_ledger(), and finalize_merge(). Do not duplicate raw-byte
formatting or splice logic. Restrict additions to eligible CDS qualifiers and
deterministic COMMENT; all other base bytes remain intact.

Central integration must:

- Normalize source names interproscan/interpro to InterProScan, recognize new
  entry types, and map fields to final qualifier/value without description
  or prefix guessing.
- Require exactly one adapter-supplied target UID on every InterProScan
  functional/member/provenance candidate; no fallback resolution.
- Assign stable candidate IDs per query/field/canonical value and preserve
  planned_status/planned_reason_code separately from final outcomes.
- Reconcile exact same-target, same-qualifier, same-value GO proposals across
  eggNOG and InterProScan; retain both candidate identities and exact support
  links to the one surviving output. Test the existing generic path before
  adding special cases.
- Synchronize each source projection with candidate_id, final_status, status,
  reason_code, insertion_ids, and supporting_candidate_ids. Current
  synchronization does not project the last two arrays; explicitly extend
  it and its tests rather than assume they exist.
- Count every final suppressed_* status in suppressed_candidate_count.
- Preserve merge-manifest v2 / ledger v3 unless a deliberate breaking schema
  decision is reviewed. Add conditional typed InterProScan entry/context
  definitions; do not globally reject existing legacy source entry shapes.

Emit at most one InterProScan producer inference per CDS **only when a new
InterProScan feature annotation survives**. Supported-existing evidence alone
does not create a new inference. Its supporting_candidate_ids must reference
eligible accepted InterProScan candidates on the exact same target.

Map InterProScan:producer narrowly to InterPro, GO, member:Pfam, and
member:TIGRFAM. Apply it in both construction and validation. Preserve Kofam
substantive evidence semantics. Same source/locus coincidence is not support.
If all linked newly emitted annotations disappear during reconciliation, prune
the inference, update its final reason/projection, and validate no orphan
provenance remains. If useful accepted support survives, exact links remain;
do not keep provenance alive merely because an unrelated source row survived.

## 7. Context, digest and transaction contract

Reuse --context-report. Change the current eggNOG-only API guard/help to allow
InterProScan and reject context output when neither relevant source is supplied.

- eggNOG-only preserves enrich-bakta.eggnog-context.v1 and its established
  semantics.
- InterProScan-only uses enrich-bakta.interproscan-context.v1.
- Both sources use enrich-bakta.context.v1, preserving each source payload
  under a named section and having top-level metadata for finalizer bindings.

InterProScan context contains captured source hash, asserted version, layout,
row shapes/counts, matched query counts, identity diagnostics, member promotion
policy, QC flags, pathway occurrence counters globally/per query, and
deterministic compact support/pathway digests. Define digest inputs, delimiters,
and order precisely and test them. No raw TSV duplication and no inferred
unique pathway total. Context-only rows do not enter the functional ledger.

Serialization order in finalize_merge() must be explicit:

1. Compute and semantically validate final GBFF bytes and candidate ledger.
2. Bind output/base hashes, insertion counts, tool version, and verified
   translation evidence into operation metadata and the sidecar.
3. Serialize finalized sidecar bytes exactly once, then hash those bytes.
4. Put context_report_sha256 in manifest metadata, then serialize manifest.
5. Stage GBFF, manifest, and sidecar before promotion with stage_artifacts().

Do not create a manifest–sidecar hash cycle or a sidecar self-hash. Collision
checks include all source inputs, FAA, parent manifest, output and sidecars.
Validation and staging failures before promotion must leave existing outputs
untouched. Current stage_artifacts() promotes files sequentially with atomic
per-file replacements; it explicitly does not provide a cross-file transaction.
Preserve and test that contract, including clear reporting if promotion fails
after one file has been replaced. Do not promise multi-file rollback. A true
cross-file recovery protocol would be a separate hardening change.

Distinguish two determinism assertions:

- Identical inputs/policies/tool version in two fresh runs produce identical
  GBFF/sidecar bytes (with automatic timestamps disabled).
- Output-as-input reruns produce unchanged GBFF and zero insertions, but
  sidecar operation base_sha256/insertions legitimately change. Require
  unchanged source-context digest and correct new operation bindings;
  do not require a first-pass sidecar hash on a no-op operation.

Default COMMENT markers must be deterministic and duplicate-safe. Use the
existing --merge-timestamp only when explicitly supplied; no wall-clock date.

## 8. CLI and package surfaces

Unified enrich() API and CLI add:

```text
--interproscan PATH
--interproscan-version 5.59-91.0       required with --interproscan
--interproscan-tsv-layout LAYOUT      default: ipr-go-pathways
--interproscan-member-dbs Pfam,TIGRFAM
```

Require --faa with this source. Reuse --manifest, --context-report,
--translation-evidence-manifest, --allow-imported-translations,
--no-comment-note, --no-feature-provenance, and --merge-timestamp.
No separate --interproscan-context-report.

Add a thin merge_interproscan_bakta.py compatibility shim and console entry
enrich-bakta-interproscan pointing to the canonical implementation. Document
standalone positional arguments consistently with existing source workflows;
no duplicate parser or merge logic. Update pyproject.toml py-modules,
public import aliases where required, package tests, wheel/sdist contents, and
installed entrypoint help/version tests outside the checkout.

## 9. Ordered execution for Luna

Work on an isolated implementation branch/worktree. Preserve every published
input and existing unrelated local file. Keep generated outputs in
.test-output. Before implementation, verify the audited SHA or review the delta
if main advanced; reproduce the 146-test baseline and 42-artifact validation.

| Phase | Files / work | Exit criterion |
|---|---|---|
| A: fixtures and shared semantics | tests; sources/value_rules.py; core/decisions.py | Native optional layouts and versioned Pfam support are specified with positive/negative fixtures; semantic validator accepts the actual witness without weakening unrelated sources. |
| B: stream/identity adapter | sources/interproscan.py | Snapshot SHA, bounded parsing/aggregation, exact UID and MD5/length checks; malformed or unresolved input has stable diagnostics and no writes. |
| C: ledger integration | core/decisions.py; core/merge_engine.py | Correct roles/classes, exact support, final projections, cross-source GO reconciliation, provenance pruning and all suppressed_* counts. |
| D: workflow/context/package | workflows/enrich.py; shim; pyproject.toml; schemas | Existing lineage gates plus legacy preflight; finalized sidecar hash; unified/standalone/API parity; typed compatible schemas and installed wheel behavior. |
| E: acceptance and release docs | tools; tests; CI; README; CHANGELOG; docs/data | Corpus baselines, first/no-op reruns, combined regressions and reviewed output goldens; existing scientific gates remain enforced. |

Focused fixtures must cover:

1. All four layout profiles; integrated/unintegrated rows; ambiguous column 14;
   native lookup-off unsupported diagnostics; truncated/extra columns.
2. LF/CRLF, trailing empty field, missing final newline, strict UTF-8,
   oversized row and aggregate limit failures; correct captured source hash.
3. Full accession/GO grammar, missing scores, mixed member methods, bounded
   source samples and pathway occurrences; no pathway feature qualifiers.
4. Exact query identities; identical proteins at two CDSs; repeated
   MD5/length conflict; duplicate/missing loci; FAA and translation mismatch.
5. Pfam version witness survives planner plus base/output ledger validation;
   false prefix, malformed version and negative prose do not count as support.
6. Supported-existing vs emitted behavior; deterministic ordering, stable
   IDs, combined exact GO support and final source projections.
7. One producer inference only for surviving new annotations; support on the
   wrong target/class/source rejected; exact orphan pruning; Kofam semantics.
8. Pristine translationless failure; bound import-faa opt-in; full origin
   carry-forward, parent/output hashes, no-op ledger preservation, strict
   complete-CDS validated-only failures.
9. Legacy restoration COMMENT with no current marker fails InterProScan
   preflight even when a synthetic FAA is supplied; no fabricated ledger.
10. eggNOG-only context compatibility; InterPro-only and combined envelopes;
    finalized sidecar SHA binding; fresh determinism vs no-op operation changes;
    collision/staging failures and the documented sequential-promotion limit.
11. All source combinations, provenance/comment switches, XLSX header alias,
    public shim/API identity, --help/--version and installed wheel execution.
12. Manifest/context schemas, all final suppressed_* totals, changed scientific
    baseline rejection, and release-version baseline migration.

Do not write production expected counts/hashes from the corpus. Keep them in
reviewable baseline fixtures/tools. Reject baseline drift; do not silently
refresh a golden to make a test pass.

## 10. Acceptance and scientific gates

Run from a clean checkout with published LFS inputs available:

```bash
python -m pip install -e '.[dev,xlsx]'
git lfs pull
python -m pytest -q -p no:cacheprovider
python -m ruff check .
python -m ruff format --check .
python -m mypy src
git diff --check
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
python tools/run_followthrough_real_checks.py
python tools/run_followthrough_chain_checks.py
```

docs/data/MANIFEST.json is the wider 86-artifact local inventory. It is not a
required clean-public-clone gate.

After implementation, C14 standalone-source smoke via unified CLI:

```bash
python enrich_bakta.py \
  --bakta data/C14/bakta/C14-NMZ.gbff \
  --faa data/C14/bakta/C14-NMZ.faa \
  --interproscan data/C14/evidence/interproscan/C14-NMZ.interproscan.tsv \
  --interproscan-version 5.59-91.0 \
  --interproscan-tsv-layout ipr-go-pathways \
  --output .test-output/C14-interpro-enriched.gbff \
  --manifest .test-output/C14-interpro-enriched.manifest.json \
  --context-report .test-output/C14-interpro-context.json
```

For SM, assert a deterministic nonzero pristine failure naming exactly the
five unresolved queried CDSs. Then use the established restore-translations
workflow with matched FAA/eggNOG and --translation-policy import-faa. Supply
that actual bound restoration manifest plus --allow-imported-translations to
the accepted InterProScan run. Confirm all six imported origins survive.

For BK71A, recompute parser/context baselines and unresolved identities, and
assert missing-input/known-legacy lineage failures. Synthetic fixtures must
separately prove the legacy preflight; failure merely because --faa was absent
does not prove lineage enforcement. No successful real functional baseline
for this historical target is claimed.

For C14 and accepted restored SM, independently inspect:

- Parsed qualifier sets, byte preservation, source IDs/roles/classes, candidate
  support, planned/final reasons, source projections and suppression totals.
- Output-as-input GBFF equality and zero new insertions; rerun with parent
  evidence where required. Compare source-context digests and actual sidecar
  hash bindings using the determinism rules above.
- Canonical gbparse validation compared with the original's known diagnostics;
  do not turn an existing informational exception into a false "genomic
  validation passed" claim.
- Combined Baktfold/KofamScan/eggNOG/InterProScan output, imported origins, TSV
  and XLSX eggNOG paths, and current three-source acceptance compatibility.
- Baseline JSON containing exact input/source hashes and producer assertions,
  profile/policies, output and first/no-op sidecar hashes, annotation sets,
  role/status/insertion counts, unresolved/imported IDs, projections and final
  reason counts. Full sidecar hashes and stable source digests are different
  fields.

Extend a durable InterProScan comparison tool and scientific CI invocation.
Review the first outputs before pinning the new baseline. A no-op alone does
not establish scientific correctness.

The current CI already runs dataset-bytes/scientific-reruns on manual dispatch
**or v* tags**. Verify both gates on the exact completed revision and retain
reports. Ordinary push success alone is not acceptance of the new real-data
workflow. No remote dispatch or publishing is performed by this document.

## 11. Release and completion handoff

Recommended feature version is 0.4.0. Before committing the release bump,
update pyproject.toml, TOOL_VERSION, README/CHANGELOG, CLI/version checks and
scientific baselines together. Existing
docs/data/followthrough-scientific-baseline.json pins 0.3.1 and exact output
hashes; versioned COMMENT changes can change those hashes even without a
scientific difference. Preserve the old accepted baseline/history, explicitly
review semantic before/after differences, and create a new compatible release
baseline. Do not replace three-source regression coverage with only the new
four-source cases.

Luna's completion report must identify the exact implemented revision, files
changed, tests/quality checks, installed package checks, accepted C14/SM
real-data counts/hashes, expected SM/BK71A failures, preserved origin ledgers,
sidecar hash bindings, reviewed baseline changes, and scientific CI reports.
If genuine BK71A recovery inputs remain absent, report that case blocked while
completing the supported feature.

**Execution prompt for Luna:**

> Implement this audited plan against e5e0499 (or first review any newer-main
> delta). Follow the phase gates and preserve existing source behavior and
> input bytes. Begin with explicit native TSV layouts, shared Pfam semantic
> support, exact protein/CDS identity, and bound translation lineage. Complete
> the adapter, central ledger integration, finalized context hashing, CLI,
> schemas, packaging and scientific acceptance. Keep BK71A functional use
> fail-closed until genuine lineage inputs exist. Produce a reviewable patch
> and evidence report; do not publish, tag, or push without the owner's
> instruction. Record real blockers rather than inventing inputs or goldens.

## 12. Primary references

- Native v5 output fields, score semantics, matched-sequence output:
  https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html
- Exact reviewed 5.59-91.0 TSV writer, including optional fields/PIRSR omission:
  https://github.com/ebi-pf-team/interproscan/blob/5.59-91.0/core/io/src/main/java/uk/ac/ebi/interpro/scan/io/match/writer/ProteinMatchesTSVResultWriter.java
- INSDC case-sensitive db_xref vocabulary:
  https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/
- ENA inference qualifier grammar:
  https://www.ebi.ac.uk/ena/WebFeat/qualifiers/inference.html
- InterPro entry/GO/pathway interpretation:
  https://interpro-documentation.readthedocs.io/en/latest/entries_info.html
- AntiFam scope and limitations:
  https://interpro-documentation.readthedocs.io/en/latest/antifam.html

The corpus-specific counts and lineage findings above come from the pinned
repository inputs, not from these external documentation pages.
