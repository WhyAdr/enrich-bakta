# enrich-bakta: validated implementation plan and patch handoff

Prepared for Wahyu and Luna, 2026-10-01 (Asia/Jakarta).

**Disposition:** approve the direction after the corrections below. The supplied
draft is a sound inventory, but it is not yet an executable release plan. Its
package naming, legacy placement, restored-protein evidence lifecycle, final
reconciliation ledger, and definition of artifact atomicity need explicit decisions.
This document supplies those decisions, tested foundation diffs, and acceptance
contracts for the remaining implementation.

**Exact baseline:** `WhyAdr/enrich-bakta`, `main`,
`5d6754d306e80cddcad7da5fc12b155ab44d2e85`, package `0.2.0`.
The remote and a fresh clone agreed at review time. No remote changes, data
publication, release, or changes to Wahyu's Windows checkout were performed.

The uploaded September 30 draft is the source of A01–A12 and the stated local
dataset inventory. The September 29 Astra audit itself is absent from this clone
and was not separately supplied. I independently checked the draft's claims
against code, tests, history, and reproductions; I do not claim to have read that
missing audit or inspected its private C14/SM inputs.

## 1. What was actually checked

The review covered all six supported Python modules, the historical normalizer,
the test suite, README, packaging, CI, ignore rules, the recent commit sequence,
and the existing August audit with its maintainer dispositions.

| Check | Result |
|---|---|
| Existing baseline tests | **31 passed** |
| Baseline Ruff, formatting, mypy, diff whitespace | Passed; mypy checked six supported modules |
| Baseline wheel and sdist | Built successfully |
| New adversarial checks against baseline | **19 failed**, as expected: they assert corrected behavior |
| Package-layout prototype | **43 passed**: original 31 plus 12 compatibility checks |
| Semantic/XLSX/staging prototype | **45 passed**: original 31 plus 14 foundation checks |
| Foundation at declared Biopython floor 1.83 | **45 passed** |
| Package + foundation patches applied in supplied order to a fresh worktree | **57 passed**; Ruff, formatting, mypy, diff checks passed |
| Installed layout wheel, outside checkout | All five commands passed help/version; all six old module aliases and both module entry points worked; `pip check` passed |
| Real C14, SM, BK71A reruns | **Not performed: the datasets are not tracked in this clone** |
| Windows behavior | Not executed here; remains an explicit release gate |

Runtime details are in `review-environment.json`. The full baseline failure
output and the independently runnable adversarial probes are included in the
bundle. These probes are review evidence, not a finished replacement for the
repository's test suite: some encode a particular skip-versus-error choice that
Luna must adapt to the finalized policy. Keep their biological counterexamples.

The official eggNOG-mapper v3 usage documentation confirms the 13-field order
already used by this repository and explicitly documents comment-free outputs.
Its public fixture also uses a git-describe producer string,
`v3.0.0-beta5-2-ga3fdb5b`. Thus schema identity and producer-version text must be
separate concepts. The correct fix is neither inventing v2 confidence nor
rejecting legitimate comment-free v3 exports indiscriminately.

Sources checked on this date:

- [Official v3 output/confidence documentation](https://github.com/eggnogdb/eggnog-mapper/blob/main/USAGE.md)
- [Official v7 fixture with producer metadata](https://github.com/eggnogdb/eggnog-mapper/blob/main/tests/fixtures_v7/test_diamond.emapper.annotations)
- [Repository baseline](https://github.com/WhyAdr/enrich-bakta/tree/5d6754d306e80cddcad7da5fc12b155ab44d2e85)

Before committing a schema registry, record the exact upstream source commit and
fixture SHA-256 in its documentation; the moving `main` links above are discovery
references, not permanent producer-contract pins.

## 2. Corrections to the original plan

1. **Resolve the package name now.** Use `src/enrich_bakta_lib/` for canonical
   implementation and retain the existing root module names as compatibility
   aliases. A root `enrich_bakta.py` shadows a `src/enrich_bakta/` package under
   ordinary source-tree import precedence. This was reproduced, not hypothetical.
   The distribution remains `enrich-bakta`; command names remain unchanged.
2. **Keep legacy outside the package.** Place the historical implementation at
   `legacy/normalize_baktfold.py`, with a source-only root launcher. Exclude both
   from the supported wheel and legacy from quality checks. The draft's
   `src/enrich_bakta/legacy/` contradicts its own exclusion requirement.
3. **Preserve the core import too.** The shim list must include `merge_engine.py`.
   Re-exporting only selected functions can also break exception/class identity
   and monkeypatches. Module aliases preserve these better than wildcard imports.
4. **Build tooling belongs in the gate's environment.** Add `build` to a
   dedicated build installation step or development extra. It is absent from the
   current `dev` extra; the draft's `python -m build` assumes an undeclared tool.
5. **Separate no-op moves from behavior changes.** Apply the two migration patches
   together, pass the original suite and installed-wheel checks, then commit.
   Apply semantic foundation separately. Never commit the intermediate move-only
   tree: its compatibility modules have not yet been recreated.
6. **Make restoration lineage survive downstream merges.** A new flag on the
   restoration command alone does not prevent a later Kofam/eggNOG merge from
   treating the imported FAA as ordinary GBFF-encoded evidence. Carry and verify
   the origin ledger in every FAA-dependent workflow.
7. **Make reconciliation operate on decisions, not display strings.** Candidate
   IDs, feature identities, support links, and final counts must be typed. Parsing
   `emitted_qualifiers` text afterward cannot reliably reconstruct decisions.
8. **Define A11's guarantee precisely.** Preflight and staging failures preserve
   prior artifact contents. A failure during the second or later promotion can
   leave earlier files replaced. Neither multiple `os.replace` calls nor
   best-effort rollback is a multi-file transaction.
9. **Keep scientific policy changes visible.** Baktfold currently promotes
   syntactically valid partial ECs, and a test requires it. eggNOG intentionally
   records partial ECs without promotion. Shared validation must not silently
   erase this source-specific distinction.
10. **Do not turn absent local data into a release claim.** The stated 0.5 GB and
    largest 47 MB artifact were not verifiable here. Inventory the actual Windows
    files before choosing storage or promising clean-clone real-data validation.
11. **Do not make the directory rename a correctness gate.** Renaming
    `D:\W\baktfold` is an optional local maintenance step after virtualenv/editable
    reinstall. It is independent of remote scientific release readiness.
12. **Use a release branch/PR by default.** Drafting does not authorize publishing.
    Review code, data permissions, and real-input results before the separate
    publication step. Do not have this handoff silently push to `main`.

## 3. Validated findings and required outcomes

| ID | Current-code evidence | Required result |
|---|---|---|
| A01 | `join(1..3,` followed by `7..9)` versus `4..6)` passes parity; equivalent one-line wrapping is rejected | Complete semantic locations, with ordered parts, strand, operator, fuzzy boundaries and remote refs retained |
| A02 | Multiline `test` + `protein` becomes `testprotein` in raw qualifiers | One authoritative semantic qualifier view; emitted values and exact target features verified |
| A03 | A 14-aa unrelated FAA entry restores into a 9-nt pseudogene; empty `/translation` receives another value | Explicit import/validation policy, raw qualifier cardinality, provenance lifecycle, plausible sequence checks |
| A04 | Version 99.0.0, version 2.1.15 with a v3-shaped row, repeated versions, empty legend all parse | Explicit schema registry, no empty-legend fallback, no silent version overwrite |
| A05 | Skipped conflicting eggNOG name still has candidate `status=emitted`; producer provenance was already planned | Final ledger and provenance derive from surviving annotated targets |
| A06 | Empty gene feature gets `abc` while its CDS retains authoritative `original` | Treat paired naming as one decision; preserve and diagnose authoritative/source inconsistencies |
| A07 | `pdb:` and EC `nonsense` are transferred; `_123` becomes empty `/gene` with suffix cleaning | Field-specific validation and explicit non-promoted invalid/missing outcomes |
| A08 | `notGO:0000001x` and `KEGG:K000010` suppress shorter valid annotations | Complete structured-token matching; prose is not automatically affirmative evidence |
| A09 | Duplicate identical CDS locus tags pass FAA validation; a `##bad` query row disappears | Raw cardinality validation before indexing and schema-aware comment/data classification |
| A10 | Captured XLSX says `original`, replaced disk path says `changed`, parser selects `changed` | Parse the captured bytes; close workbook on every branch; report malformed workbook errors consistently |
| A11 | Manifest path is a directory; primary output is overwritten before error | Validate all destinations, serialize all payloads, stage all files before the first promotion |
| A12 | Baktfold translation mismatch appears in returned parity stats but not the written manifest | Typed complete diagnostics, policy/schema flags, final support/counts and no-op outcomes in artifacts |

These are genuine release concerns even with green CI. The test suite covers
existing happy paths and selected failure paths; it does not currently cover
these counterexamples.

## 4. Implementation order and review boundaries

### P0 — Preserve baseline and recover evidence

- Record exact HEAD, Python/dependencies, original source hashes, and all current
  direct/installed command outputs. Preserve the private audit unchanged when
  Luna locates it; add it to `docs/audits/` with an applicability banner.
- Inventory ignored datasets with byte sizes/checksums before moving anything.
  Include original Bakta GBFF and FAA, every producer table, and all restoration
  manifests, not just convenient enriched outputs.
- Run all baseline gates in the exact execution environment. Preserve unrelated
  working files and inspect tracked/untracked status separately.

### P1 — Package migration, zero intended biological changes

Apply `01a-package-moves.patch`, then `01b-package-compatibility.patch` against the
exact baseline. These two are one review/commit unit. They implement:

```text
src/enrich_bakta_lib/
  core/merge_engine.py
  sources/baktfold.py
  sources/kofam.py
  sources/eggnog.py
  workflows/enrich.py
  workflows/restore_translations.py
  __main__.py
legacy/normalize_baktfold.py
```

All six supported root imports are aliases to canonical modules; five old direct
script commands invoke canonical `main()`. Both `python -m enrich_bakta` and
`python -m enrich_bakta_lib` work once installed; old source-tree direct commands
also work without relying on a globally installed editable package. No empty
`compat/` directory is needed. Add `compat/README.md` only if it explains an
actual compatibility contract.

Verify wheel contents explicitly. Merely building the wheel is insufficient:
an early prototype built successfully while omitting every root shim because
setuptools inferred the wrong module root. The corrected package-dir mapping is
included in the supplied patch and the corrected wheel was tested.

Gate: original tests, 12 compatibility tests, Ruff, format, mypy, build, installed
wheel checks from a temporary directory, all five help/version commands, all six
module aliases, legacy exclusion, and `pip check`. Linux passed locally. Add
Windows to CI for path, newline, handle-cleanup and launcher behavior.

### P2 — Semantic core, identity, XLSX snapshot and artifact staging

Apply `03-foundation-after-package-move.patch` after P1. The alternative
`02-foundation-flat-baseline-reference.patch` is the same implementation against
the old flat layout; **do not apply both**.

This tested foundation implements A01/A02, the CDS-cardinality part of A09, A10,
the staging/preflight part of A11, and central sidecar tool-version stamping:

- Raw bytes, spans, and concatenated raw location text remain available.
- Biopython records/features are aligned by record/feature ordinal; counts,
  feature types and non-null locations are checked. Qualifiers come from the
  semantic parse, including normal free-text whitespace and translation rules.
- `RawFeature.location_key` compares complete parsed locations. It does not
  collapse a join into its outermost coordinates or erase fuzzy markers.
- `feature_key()` rejects repeated identical/conflicting and blank locus values
  before CDS indexing; targeted missing IDs still fail in FAA validation.
- Finalization verifies exact byte reversibility **and** the intended qualifier
  values on the intended features, preserving original qualifiers and their
  list order. A correctly reversible payload with the wrong semantic value now
  fails before output publication.
- XLSX uses `BytesIO(snapshot)`, `data_only=False`, and closes the workbook in
  `finally`. The existing mocked XLSX test is updated to supply bytes and a close
  method, reflecting the real API lifecycle rather than skipping disk access.
- Manifest/context serialization completes before any primary replacement.
  Destination regular-file/parent/ancestry/alias checks precede staging and are
  repeated before promotion. Every temporary file is staged and flushed first.
- An injected second-staging failure preserves both old files. An injected
  second-promotion failure demonstrates and documents partial publication.

Foundation limitations to keep explicit:

- This is not the full A03–A12 patch. In particular it does not implement schema
  selection, candidate finalization, name-pair policy, or restoration lineage.
- Canonical location strings are parser-derived identity keys, not a serialized
  long-term interchange schema. Test the declared Biopython floor and selected
  current version whenever changing that dependency.
- Hand-built `RawFeature` objects without a semantic location retain a raw-key
  fallback for API compatibility. Supported file-based workflows use the
  authoritative parser. Require semantic locations for any future evidence
  importer that accepts manually constructed documents.
- Parsing now validates semantics inside `parse_genbank_bytes`; callers still
  have redundant explicit semantic checks for compatibility. Consolidate them
  later with measured performance and golden outputs, not in the migration.
- Destination checks do not serialize concurrent writers or prevent a hostile
  concurrent filesystem mutation. The contract assumes one writer per output.
- Temp-file fsync does not promise parent-directory persistence across sudden
  power loss on every platform. Do not advertise crash-durable transactions.

Gate: original 31 tests, 14 new foundation tests, compatibility tests, semantic
round-trip, zero-insertion rerun, LF/CRLF and mixed newline fixtures, complete
compound/fuzzy/remote locations, stage/promote failure injection and no temp-file
leaks. The combined patch tree passed 57 tests and the standard quality gates.

### P3 — A04 schema registry and A07/A08 source-value policy

Implement this next, before restoration or final candidate ledger integration.
Keep producer metadata separate from the selected parser schema.

**Schema selection matrix:**

| Input | Decision |
|---|---|
| Declared known producer + valid legend | Select supported registry entry and verify legend/header/vector |
| Declared known producer + absent legend | Use only that entry's documented fallback; record `version_registry` |
| Declared producer + empty/repeated/conflicting legend | Reject; empty is not absent |
| Missing producer metadata, explicit `--eggnog-schema` | Use registered schema with confidence source `explicit_schema`; record producer as unknown or user-supplied version |
| Missing metadata, no schema selection | Reject |
| Unknown declared producer/version | Reject automatic selection; explicit compatible schema is an auditable caller assertion, only allowed if documented policy enables it, and never labeled producer-verified |
| v2 input | Reject for this v3 adapter; no synthetic all-passing confidence |
| XLSX | Same registry and vector/header checks; explicit schema required when trustworthy producer metadata is unavailable |

Start automatic matching with the actual supported `3.0.0-beta6` producer. Do not
claim all versions starting with `3` have been tested. Retain raw version text;
recognize an optional leading `v` and a git-describe suffix without conflating
them with a guarantee that that revision implements the registered contract.
Test an official pinned fixture before extending the automatic matcher.

Thread `eggnog_schema` through TSV/XLSX/path parsing, standalone eggNOG,
unified enrichment, `prepare_restoration`, `restore`, and every corresponding CLI.
Add defaulted dataclass fields at the end to preserve existing `EggnogTable`
positional construction where feasible. Reject duplicate version declarations
even if identical, rather than silently taking the last one.

Example integration diff below is a **contract sketch**, not a separately
validated/apply-ready patch. Names should follow the final registry module:

```diff
--- a/src/enrich_bakta_lib/sources/eggnog.py
+++ b/src/enrich_bakta_lib/sources/eggnog.py
@@ EggnogTable
     confidence_contract_source: str
+    schema_id: str = ""
+    schema_selection_source: str = ""
+    producer_version_raw: str | None = None
@@ parse_eggnog_tsv signature
-    data: bytes, *, expected_version: str | None = None
+    data: bytes, *, expected_version: str | None = None,
+    schema_id: str | None = None
@@ version declaration handling
             if match:
+                if version is not None:
+                    raise MergeError("duplicate eggNOG producer-version declaration")
                 version = match.group(1)
@@ confidence selection
-    confidence_order = declared_confidence_order or CONFIDENCE_SCORED_FIELDS
+    contract = resolve_eggnog_contract(
+        producer_version=version,
+        requested_version=expected_version,
+        requested_schema=schema_id,
+        declared_order=declared_confidence_order,
+        columns=tuple(columns),
+    )
+    confidence_order = contract.confidence_field_order
```

Header order and confidence order are related but not identical. A supported
projection can omit annotation columns while retaining the full 13-character
vector; document that projection explicitly. Do not reject a table solely
because it lacks unused higher-order columns that the current parser permits.
Validate confidence by field name/index, not by the header position after
projection. Add all 13 one-hot checks at low/medium/high and all `-` positions.

**Value policies:** centralize grammar and normalization, but keep promotion
decisions source-specific. Return a result containing raw value, normalized
value, syntax status, and reason instead of dropping failures invisibly.

| Field | Default behavior |
|---|---|
| Full EC | Validate four numeric components; no claim of database membership |
| Partial EC | Validate only documented four-component trailing-dash forms; preserve existing Baktfold promotion unless intentionally changed; eggNOG remains non-promoted |
| `afdb_v6:`, `cath:`, `pdb:` | Exact prefix, nonempty identifier, no whitespace/control-only suffix; syntax alone does not validate the structural assignment |
| Preferred name | Check missing sentinels before suffix removal, then reject empty/missing result afterward |
| GO/KO/CAZy | Full tokens including numeric suffixes/subfamilies; no prefix matching |

For existing `db_xref` values use exact full matches. For notes accept only
documented structured token-list formats or a specifically recognized producer
evidence-note format. Arbitrary prose such as `no KEGG:K00001 hit` must not count
as affirmative annotation. A boundary-aware regex alone cannot understand
negation. This deliberately conservative policy may add an explicit xref next to
incidental prose; that is safer than incorrectly suppressing valid evidence.

Illustrative shared extraction seam:

```diff
--- a/src/enrich_bakta_lib/sources/eggnog.py
+++ b/src/enrich_bakta_lib/sources/eggnog.py
@@ _existing_values
-        values.update(match.group(0) for note in feature.values("note")
-                      for match in re.finditer(TOKEN_PATTERN, note))
+        values.update(
+            token for note in feature.values("note")
+            for token in structured_note_tokens(note)
+        )
--- a/src/enrich_bakta_lib/sources/kofam.py
+++ b/src/enrich_bakta_lib/sources/kofam.py
@@ _existing_kos
-        result.update(match.group(1) for match in KEGG_RE.finditer(value))
+        result.update(
+            token.removeprefix("KEGG:")
+            for token in structured_note_tokens(value)
+            if token.startswith("KEGG:")
+        )
```

The helper must full-match the complete recognized note grammar before returning
tokens. A useful lexical form is
`(?<![A-Za-z0-9_:])(?:GO:\d{7}|KEGG:K\d{5}|CAZy:(?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?)(?![A-Za-z0-9_])`,
but use it only within the approved affirmative grammar. Add punctuation,
numeric suffix, namespace-prefix, `GH10` versus `GH10_2`, and negative-prose
regressions shared by both adapters.

### P4 — A03 explicit translation evidence, with downstream enforcement

The core distinction is **artifact equality versus independent genomic
validation**. FAA matching a `/translation` value proves those two strings agree.
It does not itself prove either sequence derives correctly from the nucleotide
CDS, even for an original non-restored GBFF.

Implement two explicit restoration policies:

- `validated-only`: require independent genomic translation validation under a
  supported CDS model. Unsupported partial/exceptional/disrupted cases remain
  unresolved and do not receive translations under this policy.
- `import-faa`: explicitly import sequence evidence tied to an acknowledged
  source lineage; record it as imported. This never becomes independent genomic
  confirmation merely because a later GBFF contains the imported string.

Require `--translation-policy` for restoration writes. Dry-run can list
unresolved targets and the prospective policy decisions without pretending the
FAA source has been independently validated. Carry the option through planner,
prepare, restore, CLI, metadata, and updated documentation/examples.

Required checks before a restoration insertion:

1. Raw `/translation` cardinality: absence can be imported; empty, duplicate
   identical or conflicting values fail. Do not filter blank values first.
2. Query IDs/locus cardinality, record/feature identity, location and FAA lineage
   must agree. Preserve `/pseudo` and `/pseudogene`; accept either recognized
   marker by **key presence**, because `/pseudo` has an empty flag value.
3. Validate amino-acid symbols explicitly, including any deliberately supported
   ambiguity/selenocysteine/pyrrolysine symbols. Normalize a terminal stop only
   according to documented policy. Do not accept gaps, control characters or
   arbitrary punctuation. Internal stops require a supported exceptional model
   or an unresolved/rejected outcome, not automatic identity approval.
4. Validate `/codon_start`, `/transl_table`, fuzzy/partial locations and exceptions
   before genomic translation. Use `SeqFeature.extract` on the complete parsed
   location; do not slice min/max genomic bounds or lose strand/part order.
5. For complete ordinary CDSs, use the appropriate genetic code and CDS-start
   rules, including valid alternative initiators; normalize terminal stop
   consistently and compare the independently derived protein exactly.
6. For partial, recoded, frameshifted or disrupted pseudogenes, do not force an
   intact-CDS model. Report the unsupported model and policy decision. Simple
   length plausibility is a screen, not identity proof; source-repaired models
   need their own documented bounds and provenance.
7. Record no-op requested targets as well as imported, validated and unresolved
   targets. Existing translated targets remain FAA-equality checked.

Suggested explicit evidence object:

```diff
--- /dev/null
+++ b/src/enrich_bakta_lib/core/translation_evidence.py
@@ proposed typed contract
+@dataclass(frozen=True)
+class TranslationEvidence:
+    feature_uid: str
+    protein_sha256: str
+    origin: str  # gbff_encoded | imported_faa | genomically_validated
+    artifact_sequence_match: bool
+    genomic_validation_status: str  # passed | unresolved | not_attempted
+    validation_reason: str
+    original_base_sha256: str
+    faa_sha256: str
+    restoration_policy: str | None
+    producer_lineage: dict[str, object]
```

This is a design sketch; add imports, validators, schema, serialization, and
callers together. It is not a drop-in new module.

**Lineage must survive enrichment.** Require a restoration/translation-evidence
manifest when enriching a self-identified restored input, verify that it binds
to the input bytes and feature/protein identities, and require explicit
`--allow-imported-translations` for use of imported proteins in every FAA-based
functional adapter. Thread this through Kofam, eggNOG and unified CLI/API, not
only eggNOG. A fresh merge manifest must carry the verified origin ledger forward
and bind it to its own output hash, so reruns can verify the new parent artifact.

Hash binding detects disagreement and accidental substitution. It does not
authenticate a fabricated producer history. Without independent evidence,
same-tag equality or a user assertion stays an import/attestation. Do not promise
to recover biological provenance after someone deliberately strips all markers
and ledgers. Change README language from generic genomic validation to the exact
artifact/genomic distinction above.

Suggested workflow integration seams:

```diff
--- a/src/enrich_bakta_lib/workflows/restore_translations.py
+++ b/src/enrich_bakta_lib/workflows/restore_translations.py
@@ plan_translation_restoration signature
     faa_data: bytes,
+    translation_policy: str,
+    source_lineage: dict[str, object] | None = None,
@@ translation presence
-        translations = [v for v in feature.values("translation") if v.strip()]
+        translations = feature.values("translation")
+        if len(translations) > 1 or (translations and not translations[0].strip()):
+            raise MergeError("empty or duplicate translation qualifier")
@@ pseudogene eligibility
-        if not feature.values("pseudogene"):
+        if not ({"pseudo", "pseudogene"} & feature.qualifiers.keys()):
             raise MergeError(...)
@@ planner evidence classification
-                "source": "matched Bakta FAA",
-                "status": "restored_pseudogene_translation",
+                "translation_origin": translation_evidence.origin,
+                "genomic_validation_status": translation_evidence.genomic_validation_status,
+                "status": "imported_translation" if translation_evidence.origin == "imported_faa"
+                          else "validated_translation",
--- a/src/enrich_bakta_lib/core/merge_engine.py
+++ b/src/enrich_bakta_lib/core/merge_engine.py
@@ validate_faa_gbff signature
     source_name: str,
+    translation_evidence: Mapping[str, TranslationEvidence] | None = None,
+    allow_imported_translations: bool = False,
@@ after exact FAA/qualifier equality
+    enforce_translation_origins(
+        base, unique_ids, translation_evidence,
+        allow_imported=allow_imported_translations,
+    )
```

The `...` and helper names above explicitly identify integration work still to be
implemented; do not feed these illustrative blocks to `git apply`.

Gate: unrelated same-tag proteins, same-length wrong proteins, impossible
ordinary-CDS lengths, invalid symbols, internal stops, legitimate partial and
alternative-start cases, supported/unsupported exceptions, both pseudogene
markers, empty/duplicate translations, no-op targets, forged/mismatched ledger
hashes, and a full restore → eggNOG/Kofam → rerun chain retaining import status.
Validate real C14/SM/BK71A exceptions before finalizing compatibility behavior.

### P5 — A05/A06/A12 final candidate ledger and paired naming

This is the largest integration change. It requires all adapters, not a patch
only to `reconcile_insertions()`.

Define a stable `feature_uid` from record/feature ordinal plus source identity
and location key. Keep readable record ID/locus tag separately. Candidate IDs
should derive deterministically from source snapshot hash, row/feature ordinal,
field, source token and target set. Use canonical JSON bytes plus SHA-256, not
Python's randomized `hash()`. IDs describe a given input snapshot; do not claim
they stay unchanged after the input file itself changes.

Each candidate contains typed target identities, normalized value, raw source,
planned decision, final decision, supporting insertion IDs, confidence/identity
diagnostics and policy reason. Paired names have two explicit targets; one
candidate must not silently become a half-pair.

Suggested model and orchestration changes:

```diff
--- /dev/null
+++ b/src/enrich_bakta_lib/core/decisions.py
@@ proposed decision model
+@dataclass
+class CandidateDecision:
+    candidate_id: str
+    source_id: str
+    source_sha256: str
+    target_feature_uids: tuple[str, ...]
+    field: str
+    qualifier: str
+    raw_value: str
+    normalized_value: str
+    planned_status: str
+    final_status: str | None = None
+    insertion_ids: tuple[str, ...] = ()
+    reason: str = ""
--- a/src/enrich_bakta_lib/workflows/enrich.py
+++ b/src/enrich_bakta_lib/workflows/enrich.py
@@ reconciliation pipeline
-    insertions, reconciliation_rows, reconciliation = reconcile_insertions(
-        insertions, gene_conflict_policy=gene_conflict_policy
-    )
-    evidence_rows.extend(reconciliation_rows)
+    reconciled = reconcile_candidates(
+        base=base, source_plans=source_plans,
+        gene_conflict_policy=gene_conflict_policy,
+    )
+    provenance = plan_final_provenance(base, reconciled, source_metadata)
+    insertions = [*reconciled.insertions, *provenance.insertions]
+    evidence_rows = serialize_final_decisions(reconciled)
+    metadata.update(reconciled.final_counts)
@@ finalize payload
+    metadata["policies"] = effective_policy_values
+    metadata["source_diagnostics"] = source_diagnostics
+    metadata["translation_evidence"] = verified_translation_evidence
```

This sketch describes the desired dependency direction. `source_plans`, helpers
and schemas must be implemented and tested as one integration unit.

**Reconciliation rules:**

- Validate authoritative base gene/CDS pairs before proposing source names. If
  either paired feature has a nonempty authoritative name, conservatively avoid
  adding a new source name to the other half. Record base conflicts separately
  from source conflicts; do not silently repair them.
- Detect internally inconsistent Baktfold paired names before cross-source
  preference policy. Reject or suppress that naming candidate with an explicit
  reason while considering independent permitted annotations separately.
- Keep strict same-record, same-locus and complete-location matching. Do not
  broaden pairing to overlapping features to rescue a candidate.
- Distinguish unsupported/unpaired protein-coding gene/CDS pairs from legitimate
  noncoding features. Preserve currently supported feature types through an
  explicit per-feature policy; do not ban every standalone RNA-associated gene.
- Collapse identical values by exact target feature identity. Shared support
  keeps all contributing candidate IDs even when only one insertion is needed.
- Default conflicting blank-pair names to suppressed. An explicit preference
  selects a value and both targets as a unit. If the requested source provides
  no valid candidate, report unresolved/suppressed rather than manufacturing one.
- Preserve raw `planned_status`; compute `final_status` after reconciliation.
  Suggested final values are `emitted`, `supported_existing`, `shared_support`,
  `suppressed_conflict`, `filtered_confidence`, `invalid_value`, `unpaired_gene`,
  `unresolved_identity`, and `not_applicable`. Distinguish rejected whole-file
  validation from skipped candidates.
- Build feature provenance from surviving or explicitly supported-existing
  evidence according to a documented rule. A conflict-only gene source must not
  emit a gene-evidence note/inference claiming a transfer occurred.
- Retain Kofam hit score/threshold/E-value notes under the established policy:
  `--no-feature-provenance` suppresses producer inference, not substantive hit
  evidence. Model those evidence notes explicitly; do not remove them merely
  because an identical KO was already present.
- File-level comments may say a source was consulted with zero new additions.
  If so, identify that clearly. Their counts must be final, record-specific and
  separate from feature-level transferred-evidence claims.

**Ledger invariants before publication:**

1. Every `emitted`/`shared_support` decision references actual insertion IDs on
   exact output features with the intended semantic value.
2. Every functional insertion has at least one accepted supporting candidate;
   provenance insertions state their supported candidate IDs.
3. Every suppressed decision has no insertion for that rejected value/target.
4. All source contributions are retained when duplicate values collapse.
5. Unique inserted qualifiers, accepted source candidates, emitting CDSs and
   provenance insertions are separate counts. A paired name inserts two
   qualifiers but represents one naming decision.
6. No-op candidates and existing values retain relevant identity diagnostics.

Central finalization must validate these invariants against the semantic output,
not only trust a planner's counters.

**Baktfold translation-mismatch decision:** default to recording the complete
parity diagnostics and suppressing protein-derived functional/name/structural
transfers for a mismatched CDS and its paired naming candidate. If a maintainer
chooses to permit them, expose a deliberate policy and mark the weaker protein
identity in every affected candidate. Structural genomic parity alone does not
validate a different protein. Apply the policy identically in standalone,
Kofam+Baktfold and unified workflows; never drop `stats['parity']` at the wrapper.

Version the merge manifest schema to v2 when changing semantics/typed structure.
Keep an explicit v1 reader/migration note if previous sidecars are read. JSON is
the reconstruction authority; TSV is a deterministic display projection with
defined escaping/JSON encoding for list/dict fields. Do not stringify nested
objects into Python repr and call that a lossless ledger.

Record tool version, original/input/output hashes, adapter schema/format,
effective flags including no-comment/no-feature-provenance/clean suffix,
confidence contract selection, restoration policy/origin and verified parent
ledger identity in every applicable artifact. Record the real support identities,
not merely source names joined with ` | `.

Gate: skip and both preferences; equal names; multiple sources supporting one
value; conflict-only sources; authoritative half-pairs; inconsistent base/source
pairs; standalone legitimate noncoding features; mismatched Baktfold proteins;
every no-op restoration target; and manifest-to-output assertions over all rows.

### P6 — Complete A09 parser/error consistency and packaging gates

The foundation fixes CDS cardinality, but comment/data classification remains.
After the TSV header, distinguish schema-recognized producer comments/footer
lines from tabular data-shaped `#` and `##` rows. Never unconditionally skip all
double-hash lines. Recognize legitimate producer command metadata that may contain
tabs; classification must use the selected producer grammar and row shape.
Header/metadata declarations after data begins should be rejected if the schema
does not allow them. Add valid footer/command fixtures and malformed reserved-ID
rows together to prevent exchanging a silent skip bug for false rejection.

Review Kofam separately: its parser accepts only leading-asterisk hit rows and
legitimately ignores below-threshold rows. Keep that producer behavior. Do not
globally forbid negative profile scores/thresholds; the existing audit correctly
rejected that assumption. Validate selected hit rows and report malformed
header/comment/data boundaries explicitly.

Reject ambiguous target IDs before indexing, including identical duplicates,
blank locus values, cross-record duplicate CDS IDs and ambiguous paired genes.
Do not insist every unreferenced CDS has a locus identifier if the supported
input contract does not require it. Include record ordinal in internal identity
or explicitly reject duplicate record IDs where record-ID keys would alias.

Normalize expected file/serialization/parsing failures to contextual `MergeError`
at public workflow boundaries and CLI exit 2. Keep unexpected programming errors
visible to tests. Avoid `except Exception` around whole CLIs; the XLSX boundary
can wrap its library's multiple parser exceptions as in the foundation.

Add a build job and wheel/sdist artifact inspection. Install the wheel from a
temporary directory with the checkout absent from import paths. Test base-only
installation and the optional XLSX extra independently. Add a targeted floor
dependency job and Linux/Windows supported-Python matrix; do not rely on an
editable install to prove packaging compatibility.

Single-source and unified outputs need separate checks. Sequential enrichment
can have different provenance/order from one-pass merging; require documented
functional/decision equivalence, not accidental byte identity between different
execution histories. Repeating the same defined workflow on its own enriched
output should yield zero new insertions and identical GBFF bytes. Fresh duplicate
runs on identical inputs should yield identical output and artifacts when paths
and explicit metadata are held constant. A rerun manifest may truthfully differ
because its input/output roles differ; do not require a first-run insertion ledger
to equal a second-run no-op ledger.

### P7 — Documentation and datasets

Keep the original draft's historical-doc mapping, but preserve original authors,
dates, exact reviewed commits and any superseding disposition. Replace hard-coded
Windows links with relative links. Add a `docs/README.md` navigation index.

Use the original `data/C14`, `data/SM`, `data/BK71A`, `data/fixtures` layout after
an actual inventory. Place the original Bakta FAA next to its original GBFF;
restored/enriched outputs go in explicitly labeled derived directories. Never
use derived outputs as pristine identity references.

Recommended manifest is a versioned canonical JSON document plus generated TSV:

```diff
--- /dev/null
+++ b/docs/data/MANIFEST.example.json
@@ minimum conceptual structure; implement separate full JSON Schema and validator
+{
+  "schema": "enrich-bakta.dataset-manifest.v1",
+  "artifacts": [
+    {
+      "path": "C14/bakta/original.gbff",
+      "sample": "C14",
+      "role": "upstream_input",
+      "size_bytes": 0,
+      "sha256": "<64 lowercase hex characters>",
+      "lineage": [],
+      "producer": {"name": "Bakta", "version": "<verified>", "database": "<verified>"},
+      "status": "historical_unverified",
+      "redistribution": {"license": "<verified>", "permission": "<verified>"},
+      "validation": {"code_commit": null, "checks": []}
+    }
+  ]
+}
```

The example above is a data-shape sketch. Implement separate schema and data files; no
placeholder metadata or zero sizes may enter a real published manifest.

The validator must reject absolute paths, `..` traversal, duplicate normalized
paths, symlink escapes outside the dataset root, invalid hash/size values and
missing lineage references. Verify path containment after resolution, file type,
actual size/hash, permitted roles/statuses, acyclic lineage and required producer
metadata. For published datasets, permission/sensitivity fields must be resolved,
not `unknown`. The same validator should support ordinary small CI fixtures.

Existing `.gitignore` patterns ignore annotation files at any depth. Moving them
under `data/` does not make them publishable. Choose a curated tracked-subtree
policy or storage-specific manifest; verify every intended file with
`git check-ignore -v` and `git ls-files`. Prefer explicit scoped patterns over a
broad `git add -f data/` that bypasses protections accidentally.

Measure ordinary Git, LFS and release-artifact behavior before selecting storage.
For LFS/release assets, a clean clone must explicitly materialize bytes and verify
checksums; an LFS pointer is not a GBFF fixture. Keep small synthetic fixtures in
ordinary Git and full biological datasets behind a manually invoked validation
job if runtime/storage warrants it. Inspect inherited producer/database
redistribution conditions as the original plan requires. The current clone also
lacks a repository `LICENSE`; select a code license deliberately before describing
the repository as licensed open source, independently of data permission.

Dataset acceptance requires commands, versions, input hashes and exact expected
functional/provenance outcomes. Feature counts alone cannot validate annotation
transfer. Historical derived files may be preserved with unverified status, but
they must not become scientifically validated merely by passing format parsing.

### P8 — Release review and optional local rename

Gate all twelve findings with the tests above, both package-install modes, a
clean-clone fixture validation, real C14/SM/BK71A reruns where applicable and
independent review of the final decision ledger. Build from the committed tree,
verify wheel/sdist contents, source links and exact artifact hashes.

Prepare a review branch/PR containing small concerns in this order: preservation,
migration, foundation, source contracts/values, translation evidence lifecycle,
decision ledger/paired naming, remaining parser/CI hardening, docs, curated data.
Do not tag or publish a release claiming conservative transfer before the
restoration and reconciliation contracts pass.

Only after code/data release review and the user's publication instruction,
publish the authorized commits/artifacts and record the exact remote SHA. The
plan's permission gates apply to that future publication; they do not prevent
Luna implementing and testing the local review branch now.

Rename the local checkout last if still desired. Recreate/reinstall editable
environments after moving it; old `.pth` paths and command launchers may retain
the previous path. Check IDE tasks, notebooks and local configs, not just tracked
README strings. Keep the move reversible. Release readiness does not depend on
the basename of Wahyu's folder.

## 5. Luna execution instructions

1. Verify exact HEAD and working-tree state. If `main` moved, rebase/revalidate
   every supplied patch; do not assume audit conclusions or hunk context persist.
2. Apply the tested migration pair and pass its gates before committing it.
3. Apply the tested canonical-path foundation patch and pass combined gates.
4. Implement P3–P6 as separate reviewed units. Illustrative diff sketches specify
   contracts; they are not complete implementation patches or validated APIs.
5. Preserve all 19 biological/operational counterexamples as regression coverage.
   Convert probes to production tests with finalized failure/skip decisions,
   source schemas, legitimate producer metadata and real positive examples.
6. Report every changed policy, exact tests and unresolved real-data case. Do not
   replace a difficult counterexample with a weaker test merely to get green CI.
7. Complete local code, artifact and dataset checks before asking Wahyu for a
   concrete publication decision. Keep unrelated local data and audits intact.

## 6. Tested diff blocks

The following appendices contain the complete tested patches. They are included
verbatim in this handoff and as standalone `.patch` files in the bundle. Apply
01a → 01b → 03. Do not apply the optional flat-baseline reference as well.

Validation of these patches is deliberately narrower than release approval:
57 passing combined tests prove the migration and foundation, not completion of
the remaining restoration/schema/reconciliation policies or real-data reruns.

### Appendix A — package moves

```diff
diff --git a/normalize_baktfold.py b/legacy/normalize_baktfold.py
similarity index 100%
rename from normalize_baktfold.py
rename to legacy/normalize_baktfold.py
diff --git a/src/enrich_bakta_lib/__init__.py b/src/enrich_bakta_lib/__init__.py
new file mode 100644
index 0000000..32258c6
--- /dev/null
+++ b/src/enrich_bakta_lib/__init__.py
@@ -0,0 +1 @@
+"""enrich-bakta implementation package."""
diff --git a/src/enrich_bakta_lib/__main__.py b/src/enrich_bakta_lib/__main__.py
new file mode 100644
index 0000000..883eeff
--- /dev/null
+++ b/src/enrich_bakta_lib/__main__.py
@@ -0,0 +1,3 @@
+from .workflows.enrich import main
+
+raise SystemExit(main())
diff --git a/src/enrich_bakta_lib/core/__init__.py b/src/enrich_bakta_lib/core/__init__.py
new file mode 100644
index 0000000..32258c6
--- /dev/null
+++ b/src/enrich_bakta_lib/core/__init__.py
@@ -0,0 +1 @@
+"""enrich-bakta implementation package."""
diff --git a/merge_engine.py b/src/enrich_bakta_lib/core/merge_engine.py
similarity index 100%
rename from merge_engine.py
rename to src/enrich_bakta_lib/core/merge_engine.py
diff --git a/src/enrich_bakta_lib/sources/__init__.py b/src/enrich_bakta_lib/sources/__init__.py
new file mode 100644
index 0000000..32258c6
--- /dev/null
+++ b/src/enrich_bakta_lib/sources/__init__.py
@@ -0,0 +1 @@
+"""enrich-bakta implementation package."""
diff --git a/graft_baktfold_additions.py b/src/enrich_bakta_lib/sources/baktfold.py
similarity index 99%
rename from graft_baktfold_additions.py
rename to src/enrich_bakta_lib/sources/baktfold.py
index 9066c21..7ae24e8 100644
--- a/graft_baktfold_additions.py
+++ b/src/enrich_bakta_lib/sources/baktfold.py
@@ -17,7 +17,7 @@ import sys
 from pathlib import Path
 from typing import Any

-from merge_engine import (
+from enrich_bakta_lib.core.merge_engine import (
     TOOL_VERSION,
     Insertion,
     MergeError,
diff --git a/merge_eggnog_bakta.py b/src/enrich_bakta_lib/sources/eggnog.py
similarity index 99%
rename from merge_eggnog_bakta.py
rename to src/enrich_bakta_lib/sources/eggnog.py
index d231898..59907b1 100644
--- a/merge_eggnog_bakta.py
+++ b/src/enrich_bakta_lib/sources/eggnog.py
@@ -12,7 +12,7 @@ from decimal import Decimal, InvalidOperation
 from pathlib import Path
 from typing import Any

-from merge_engine import (
+from enrich_bakta_lib.core.merge_engine import (
     TOOL_VERSION,
     Insertion,
     MergeError,
diff --git a/merge_kofamscan_bakta.py b/src/enrich_bakta_lib/sources/kofam.py
similarity index 98%
rename from merge_kofamscan_bakta.py
rename to src/enrich_bakta_lib/sources/kofam.py
index 378b0fa..ba0f04a 100644
--- a/merge_kofamscan_bakta.py
+++ b/src/enrich_bakta_lib/sources/kofam.py
@@ -12,8 +12,7 @@ from decimal import Decimal, InvalidOperation
 from pathlib import Path
 from typing import Any

-from graft_baktfold_additions import plan_baktfold_additions
-from merge_engine import (
+from enrich_bakta_lib.core.merge_engine import (
     TOOL_VERSION,
     Insertion,
     MergeError,
@@ -28,9 +27,10 @@ from merge_engine import (
     sha256_bytes,
     validate_genbank_semantics,
 )
-from merge_engine import (
+from enrich_bakta_lib.core.merge_engine import (
     validate_faa_gbff as _validate_faa_gbff,
 )
+from enrich_bakta_lib.sources.baktfold import plan_baktfold_additions

 KO_RE = re.compile(r"K\d{5}\Z")
 KEGG_RE = re.compile(r"(?:^|(?<=[^A-Za-z0-9]))KEGG:(K\d{5})\b")
diff --git a/src/enrich_bakta_lib/workflows/__init__.py b/src/enrich_bakta_lib/workflows/__init__.py
new file mode 100644
index 0000000..32258c6
--- /dev/null
+++ b/src/enrich_bakta_lib/workflows/__init__.py
@@ -0,0 +1 @@
+"""enrich-bakta implementation package."""
diff --git a/enrich_bakta.py b/src/enrich_bakta_lib/workflows/enrich.py
similarity index 96%
rename from enrich_bakta.py
rename to src/enrich_bakta_lib/workflows/enrich.py
index 62c3013..b8ed269 100644
--- a/enrich_bakta.py
+++ b/src/enrich_bakta_lib/workflows/enrich.py
@@ -9,9 +9,7 @@ import sys
 from pathlib import Path
 from typing import Any

-from graft_baktfold_additions import plan_baktfold_additions
-from merge_eggnog_bakta import parse_eggnog_path, plan_eggnog_additions
-from merge_engine import (
+from enrich_bakta_lib.core.merge_engine import (
     TOOL_VERSION,
     MergeError,
     finalize_merge,
@@ -21,7 +19,9 @@ from merge_engine import (
     sha256_bytes,
     validate_genbank_semantics,
 )
-from merge_kofamscan_bakta import parse_kofam_table, plan_kofam_additions
+from enrich_bakta_lib.sources.baktfold import plan_baktfold_additions
+from enrich_bakta_lib.sources.eggnog import parse_eggnog_path, plan_eggnog_additions
+from enrich_bakta_lib.sources.kofam import parse_kofam_table, plan_kofam_additions


 def enrich(
diff --git a/restore_bakta_translations.py b/src/enrich_bakta_lib/workflows/restore_translations.py
similarity index 98%
rename from restore_bakta_translations.py
rename to src/enrich_bakta_lib/workflows/restore_translations.py
index d86fb6c..78aee08 100644
--- a/restore_bakta_translations.py
+++ b/src/enrich_bakta_lib/workflows/restore_translations.py
@@ -9,8 +9,7 @@ import sys
 from pathlib import Path
 from typing import Any

-from merge_eggnog_bakta import EggnogTable, parse_eggnog_path
-from merge_engine import (
+from enrich_bakta_lib.core.merge_engine import (
     TOOL_VERSION,
     Insertion,
     MergeError,
@@ -25,6 +24,7 @@ from merge_engine import (
     validate_faa_gbff,
     validate_genbank_semantics,
 )
+from enrich_bakta_lib.sources.eggnog import EggnogTable, parse_eggnog_path


 def plan_translation_restoration(
```

### Appendix B — compatibility and packaging

```diff
diff --git a/.gitignore b/.gitignore
index 64cd8cd..0034fb9 100644
--- a/.gitignore
+++ b/.gitignore
@@ -18,3 +18,9 @@ BK71A-restored.gbff
 *-timing.json
 /.test-output/
 .pytest_cache/
+
+.mypy_cache/
+.ruff_cache/
+*.egg-info/
+build/
+dist/
diff --git a/enrich_bakta.py b/enrich_bakta.py
new file mode 100644
index 0000000..248d029
--- /dev/null
+++ b/enrich_bakta.py
@@ -0,0 +1,18 @@
+"""Compatibility entry point for enrich_bakta_lib.workflows.enrich."""
+
+from __future__ import annotations
+
+import importlib
+import sys
+from pathlib import Path
+
+_source_root = Path(__file__).resolve().parent / "src"
+if (_source_root / "enrich_bakta_lib").is_dir():
+    sys.path.insert(0, str(_source_root))
+
+_implementation = importlib.import_module("enrich_bakta_lib.workflows.enrich")
+
+if __name__ == "__main__":
+    raise SystemExit(_implementation.main())
+else:
+    sys.modules[__name__] = _implementation
diff --git a/graft_baktfold_additions.py b/graft_baktfold_additions.py
new file mode 100644
index 0000000..33947cd
--- /dev/null
+++ b/graft_baktfold_additions.py
@@ -0,0 +1,18 @@
+"""Compatibility entry point for enrich_bakta_lib.sources.baktfold."""
+
+from __future__ import annotations
+
+import importlib
+import sys
+from pathlib import Path
+
+_source_root = Path(__file__).resolve().parent / "src"
+if (_source_root / "enrich_bakta_lib").is_dir():
+    sys.path.insert(0, str(_source_root))
+
+_implementation = importlib.import_module("enrich_bakta_lib.sources.baktfold")
+
+if __name__ == "__main__":
+    raise SystemExit(_implementation.main())
+else:
+    sys.modules[__name__] = _implementation
diff --git a/merge_eggnog_bakta.py b/merge_eggnog_bakta.py
new file mode 100644
index 0000000..ff389e9
--- /dev/null
+++ b/merge_eggnog_bakta.py
@@ -0,0 +1,18 @@
+"""Compatibility entry point for enrich_bakta_lib.sources.eggnog."""
+
+from __future__ import annotations
+
+import importlib
+import sys
+from pathlib import Path
+
+_source_root = Path(__file__).resolve().parent / "src"
+if (_source_root / "enrich_bakta_lib").is_dir():
+    sys.path.insert(0, str(_source_root))
+
+_implementation = importlib.import_module("enrich_bakta_lib.sources.eggnog")
+
+if __name__ == "__main__":
+    raise SystemExit(_implementation.main())
+else:
+    sys.modules[__name__] = _implementation
diff --git a/merge_engine.py b/merge_engine.py
new file mode 100644
index 0000000..deee600
--- /dev/null
+++ b/merge_engine.py
@@ -0,0 +1,18 @@
+"""Compatibility entry point for enrich_bakta_lib.core.merge_engine."""
+
+from __future__ import annotations
+
+import importlib
+import sys
+from pathlib import Path
+
+_source_root = Path(__file__).resolve().parent / "src"
+if (_source_root / "enrich_bakta_lib").is_dir():
+    sys.path.insert(0, str(_source_root))
+
+_implementation = importlib.import_module("enrich_bakta_lib.core.merge_engine")
+
+if __name__ == "__main__":
+    pass
+else:
+    sys.modules[__name__] = _implementation
diff --git a/merge_kofamscan_bakta.py b/merge_kofamscan_bakta.py
new file mode 100644
index 0000000..bfdc1de
--- /dev/null
+++ b/merge_kofamscan_bakta.py
@@ -0,0 +1,18 @@
+"""Compatibility entry point for enrich_bakta_lib.sources.kofam."""
+
+from __future__ import annotations
+
+import importlib
+import sys
+from pathlib import Path
+
+_source_root = Path(__file__).resolve().parent / "src"
+if (_source_root / "enrich_bakta_lib").is_dir():
+    sys.path.insert(0, str(_source_root))
+
+_implementation = importlib.import_module("enrich_bakta_lib.sources.kofam")
+
+if __name__ == "__main__":
+    raise SystemExit(_implementation.main())
+else:
+    sys.modules[__name__] = _implementation
diff --git a/normalize_baktfold.py b/normalize_baktfold.py
new file mode 100644
index 0000000..05e995a
--- /dev/null
+++ b/normalize_baktfold.py
@@ -0,0 +1,6 @@
+"""Source-only compatibility launcher for the historical normalizer."""
+from pathlib import Path
+import runpy
+
+if __name__ == "__main__":
+    runpy.run_path(str(Path(__file__).resolve().parent / "legacy/normalize_baktfold.py"), run_name="__main__")
diff --git a/pyproject.toml b/pyproject.toml
index e9a29ad..0859e71 100644
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -19,11 +19,11 @@ dev = [
 ]

 [project.scripts]
-enrich-bakta = "enrich_bakta:main"
-enrich-bakta-baktfold = "graft_baktfold_additions:main"
-enrich-bakta-kofam = "merge_kofamscan_bakta:main"
-enrich-bakta-eggnog = "merge_eggnog_bakta:main"
-enrich-bakta-restore-translations = "restore_bakta_translations:main"
+enrich-bakta = "enrich_bakta_lib.workflows.enrich:main"
+enrich-bakta-baktfold = "enrich_bakta_lib.sources.baktfold:main"
+enrich-bakta-kofam = "enrich_bakta_lib.sources.kofam:main"
+enrich-bakta-eggnog = "enrich_bakta_lib.sources.eggnog:main"
+enrich-bakta-restore-translations = "enrich_bakta_lib.workflows.restore_translations:main"

 [tool.setuptools]
 py-modules = [
@@ -35,22 +35,23 @@ py-modules = [
   "restore_bakta_translations",
 ]

+[tool.setuptools.package-dir]
+"" = "."
+enrich_bakta_lib = "src/enrich_bakta_lib"
+
+[tool.setuptools.packages.find]
+where = ["src"]
+include = ["enrich_bakta_lib*"]
+
 [tool.pytest.ini_options]
 testpaths = ["tests"]

 [tool.ruff]
-extend-exclude = ["*.md", "normalize_baktfold.py"]
+extend-exclude = ["*.md", "normalize_baktfold.py", "legacy/"]

 [tool.ruff.lint]
 select = ["E4", "E7", "E9", "F", "I"]

 [tool.mypy]
 ignore_missing_imports = true
-files = [
-  "enrich_bakta.py",
-  "graft_baktfold_additions.py",
-  "merge_eggnog_bakta.py",
-  "merge_engine.py",
-  "merge_kofamscan_bakta.py",
-  "restore_bakta_translations.py",
-]
+files = ["src/enrich_bakta_lib"]
diff --git a/restore_bakta_translations.py b/restore_bakta_translations.py
new file mode 100644
index 0000000..6308dbc
--- /dev/null
+++ b/restore_bakta_translations.py
@@ -0,0 +1,20 @@
+"""Compatibility entry point for enrich_bakta_lib.workflows.restore_translations."""
+
+from __future__ import annotations
+
+import importlib
+import sys
+from pathlib import Path
+
+_source_root = Path(__file__).resolve().parent / "src"
+if (_source_root / "enrich_bakta_lib").is_dir():
+    sys.path.insert(0, str(_source_root))
+
+_implementation = importlib.import_module(
+    "enrich_bakta_lib.workflows.restore_translations"
+)
+
+if __name__ == "__main__":
+    raise SystemExit(_implementation.main())
+else:
+    sys.modules[__name__] = _implementation
diff --git a/tests/test_compatibility.py b/tests/test_compatibility.py
new file mode 100644
index 0000000..884ae41
--- /dev/null
+++ b/tests/test_compatibility.py
@@ -0,0 +1,48 @@
+from __future__ import annotations
+
+import importlib
+import subprocess
+import sys
+from pathlib import Path
+
+import pytest
+
+PAIRS = {
+    "merge_engine": "core.merge_engine",
+    "graft_baktfold_additions": "sources.baktfold",
+    "merge_kofamscan_bakta": "sources.kofam",
+    "merge_eggnog_bakta": "sources.eggnog",
+    "enrich_bakta": "workflows.enrich",
+    "restore_bakta_translations": "workflows.restore_translations",
+}
+
+
+@pytest.mark.parametrize("legacy, canonical", PAIRS.items())
+def test_alias_preserves_module_and_class_identity(legacy, canonical):
+    old = importlib.import_module(legacy)
+    new = importlib.import_module("enrich_bakta_lib." + canonical)
+    assert old is new
+
+
+@pytest.mark.parametrize("legacy", [name for name in PAIRS if name != "merge_engine"])
+def test_direct_launcher_outside_checkout(tmp_path, legacy):
+    root = Path(__file__).resolve().parents[1]
+    result = subprocess.run(
+        [sys.executable, str(root / (legacy + ".py")), "--help"],
+        cwd=tmp_path,
+        capture_output=True,
+        text=True,
+    )
+    assert result.returncode == 0, result.stderr
+
+
+def test_module_entry_point(tmp_path):
+    root = Path(__file__).resolve().parents[1]
+    result = subprocess.run(
+        [sys.executable, "-m", "enrich_bakta", "--version"],
+        cwd=root,
+        capture_output=True,
+        text=True,
+    )
+    assert result.returncode == 0, result.stderr
+    assert "0.2.0" in result.stdout
```

### Appendix C — semantic, XLSX and staging foundation

```diff
diff --git a/src/enrich_bakta_lib/sources/eggnog.py b/src/enrich_bakta_lib/sources/eggnog.py
index d231898..35a6789 100644
--- a/src/enrich_bakta_lib/sources/eggnog.py
+++ b/src/enrich_bakta_lib/sources/eggnog.py
@@ -4,6 +4,7 @@
 from __future__ import annotations

 import argparse
+import io
 import json
 import re
 import sys
@@ -366,7 +367,7 @@ def parse_eggnog_tsv(


 def parse_eggnog_xlsx(
-    path: Path, *, expected_version: str | None = None
+    path: Path, *, expected_version: str | None = None, data: bytes | None = None
 ) -> EggnogTable:
     if expected_version is None:
         raise MergeError(
@@ -378,72 +379,84 @@ def parse_eggnog_xlsx(
         raise MergeError(
             "XLSX input requires optional dependency openpyxl; TSV needs no extra dependency"
         ) from exc
-    workbook = openpyxl.load_workbook(path, read_only=True, data_only=False)
-    if workbook.sheetnames != ["annotations"]:
-        raise MergeError(
-            "eggNOG XLSX must contain exactly one sheet named 'annotations'"
-        )
-    worksheet = workbook["annotations"]
-    rows = worksheet.iter_rows()
+    snapshot = data if data is not None else path.read_bytes()
+    workbook = None
     try:
-        header_cells = next(rows)
-    except StopIteration as exc:
-        raise MergeError("eggNOG XLSX is empty") from exc
-    columns = [
-        str(cell.value) if cell.value is not None else "" for cell in header_cells
-    ]
-    if len(columns) != len(set(columns)):
-        raise MergeError("eggNOG XLSX header contains duplicate columns")
-    if any(not column.strip() for column in columns):
-        raise MergeError("eggNOG XLSX header contains empty column names")
-    hits: list[EggnogHit] = []
-    query_ids: set[str] = set()
-    for row_number, cells in enumerate(rows, start=2):
-        if all(cell.value is None for cell in cells):
-            continue
-        if len(cells) != len(columns):
-            raise MergeError(f"eggNOG XLSX row {row_number}: wrong cell count")
-        values: list[str] = []
-        for cell in cells:
-            if cell.data_type == "f":
+        workbook = openpyxl.load_workbook(
+            io.BytesIO(snapshot), read_only=True, data_only=False
+        )
+        if workbook.sheetnames != ["annotations"]:
+            raise MergeError(
+                "eggNOG XLSX must contain exactly one sheet named 'annotations'"
+            )
+        worksheet = workbook["annotations"]
+        rows = worksheet.iter_rows()
+        try:
+            header_cells = next(rows)
+        except StopIteration as exc:
+            raise MergeError("eggNOG XLSX is empty") from exc
+        columns = [
+            str(cell.value) if cell.value is not None else "" for cell in header_cells
+        ]
+        if len(columns) != len(set(columns)):
+            raise MergeError("eggNOG XLSX header contains duplicate columns")
+        if any(not column.strip() for column in columns):
+            raise MergeError("eggNOG XLSX header contains empty column names")
+        hits: list[EggnogHit] = []
+        query_ids: set[str] = set()
+        for row_number, cells in enumerate(rows, start=2):
+            if all(cell.value is None for cell in cells):
+                continue
+            if len(cells) != len(columns):
+                raise MergeError(f"eggNOG XLSX row {row_number}: wrong cell count")
+            values: list[str] = []
+            for cell in cells:
+                if cell.data_type == "f":
+                    raise MergeError(
+                        f"eggNOG XLSX row {row_number}: formulas are not accepted"
+                    )
+                values.append("" if cell.value is None else str(cell.value))
+            hit = _build_hit(dict(zip(columns, values, strict=True)), row_number)
+            if hit.query_id in query_ids:
                 raise MergeError(
-                    f"eggNOG XLSX row {row_number}: formulas are not accepted"
+                    f"eggNOG XLSX row {row_number}: duplicate query ID {hit.query_id!r}"
                 )
-            values.append("" if cell.value is None else str(cell.value))
-        hit = _build_hit(dict(zip(columns, values, strict=True)), row_number)
-        if hit.query_id in query_ids:
+            query_ids.add(hit.query_id)
+            hits.append(hit)
+        if not hits:
+            raise MergeError("eggNOG XLSX contains no annotation rows")
+        header_order = tuple(
+            column for column in columns if column in CONFIDENCE_FIELD_INDEX
+        )
+        expected_header_order = tuple(
+            column for column in CONFIDENCE_SCORED_FIELDS if column in header_order
+        )
+        if header_order != expected_header_order:
             raise MergeError(
-                f"eggNOG XLSX row {row_number}: duplicate query ID {hit.query_id!r}"
+                "eggNOG scored columns conflict with the confidence field-order contract"
             )
-        query_ids.add(hit.query_id)
-        hits.append(hit)
-    if not hits:
-        raise MergeError("eggNOG XLSX contains no annotation rows")
-    header_order = tuple(
-        column for column in columns if column in CONFIDENCE_FIELD_INDEX
-    )
-    expected_header_order = tuple(
-        column for column in CONFIDENCE_SCORED_FIELDS if column in header_order
-    )
-    if header_order != expected_header_order:
-        raise MergeError(
-            "eggNOG scored columns conflict with the confidence field-order contract"
+        return EggnogTable(
+            tuple(hits),
+            tuple(columns),
+            expected_version,
+            "xlsx",
+            CONFIDENCE_SCORED_FIELDS,
+            "documented_default",
         )
-    return EggnogTable(
-        tuple(hits),
-        tuple(columns),
-        expected_version,
-        "xlsx",
-        CONFIDENCE_SCORED_FIELDS,
-        "documented_default",
-    )
+    except MergeError:
+        raise
+    except Exception as exc:
+        raise MergeError(f"eggNOG XLSX parsing failed: {exc}") from exc
+    finally:
+        if workbook is not None:
+            workbook.close()


 def parse_eggnog_path(
     path: Path, *, expected_version: str | None = None, data: bytes | None = None
 ) -> EggnogTable:
     if path.suffix.lower() == ".xlsx":
-        return parse_eggnog_xlsx(path, expected_version=expected_version)
+        return parse_eggnog_xlsx(path, expected_version=expected_version, data=data)
     return parse_eggnog_tsv(
         data if data is not None else path.read_bytes(),
         expected_version=expected_version,
@@ -475,7 +488,7 @@ def _index_paired_genes(
     for feature in base.features:
         if feature.feature_type != "gene" or not feature.locus_tag:
             continue
-        key = (feature.record_index, feature.locus_tag, feature.location)
+        key = (feature.record_index, feature.locus_tag, feature.location_key)
         index[key] = None if key in index else feature
     return index

@@ -483,7 +496,7 @@ def _index_paired_genes(
 def _paired_gene(
     index: dict[tuple[int, str, str], RawFeature | None], cds: RawFeature
 ) -> RawFeature | None:
-    return index.get((cds.record_index, cds.locus_tag or "", cds.location))
+    return index.get((cds.record_index, cds.locus_tag or "", cds.location_key))


 def _normalized_gene(raw: str, clean_suffix: bool) -> str:
diff --git a/src/enrich_bakta_lib/core/merge_engine.py b/src/enrich_bakta_lib/core/merge_engine.py
index a096d05..a1c060f 100644
--- a/src/enrich_bakta_lib/core/merge_engine.py
+++ b/src/enrich_bakta_lib/core/merge_engine.py
@@ -18,6 +18,8 @@ from pathlib import Path
 from typing import Any

 from Bio import BiopythonParserWarning, SeqIO
+from Bio.SeqFeature import CompoundLocation, SimpleLocation
+from Bio.SeqRecord import SeqRecord

 QUALIFIER_INDENT = b" " * 21
 TOOL_VERSION = "0.2.0"
@@ -40,6 +42,16 @@ class RawFeature:
     end: int
     qualifiers: dict[str, list[str]] = field(default_factory=dict)

+    semantic_location: SimpleLocation | CompoundLocation | None = None
+
+    @property
+    def location_key(self) -> str:
+        return (
+            str(self.semantic_location)
+            if self.semantic_location is not None
+            else self.location
+        )
+
     def values(self, key: str) -> list[str]:
         return list(self.qualifiers.get(key, []))

@@ -49,8 +61,8 @@ class RawFeature:

     @property
     def locus_tag(self) -> str | None:
-        values = list(dict.fromkeys(self.locus_tags))
-        return values[0] if len(values) == 1 else None
+        values = self.values("locus_tag")
+        return values[0] if len(values) == 1 and values[0] else None


 @dataclass
@@ -133,38 +145,6 @@ def _decode(value: bytes, label: str) -> str:
         raise MergeError(f"{label} is not valid UTF-8 at byte {exc.start}") from exc


-def _parse_qualifiers(block: bytes) -> dict[str, list[str]]:
-    qualifiers: dict[str, list[str]] = collections.defaultdict(list)
-    current_key: str | None = None
-    pieces: list[str] = []
-
-    def flush() -> None:
-        nonlocal current_key, pieces
-        if current_key is None:
-            return
-        raw = "".join(piece.strip() for piece in pieces)
-        if len(raw) >= 2 and raw.startswith('"') and raw.endswith('"'):
-            raw = raw[1:-1].replace('""', '"')
-        qualifiers[current_key].append(raw)
-        current_key = None
-        pieces = []
-
-    for raw_line in block.splitlines()[1:]:
-        match = _QUALIFIER_RE.match(raw_line)
-        if match:
-            flush()
-            current_key = _decode(match.group(1), "qualifier name")
-            pieces = [_decode(match.group(2) or b"", f"/{current_key} value")]
-        elif current_key is not None and raw_line.startswith(QUALIFIER_INDENT):
-            pieces.append(
-                _decode(
-                    raw_line[len(QUALIFIER_INDENT) :], f"/{current_key} continuation"
-                )
-            )
-    flush()
-    return dict(qualifiers)
-
-
 def parse_genbank_bytes(data: bytes, label: str = "GenBank input") -> RawDocument:
     """Locate records/features without altering the original bytes."""
     lines = _split_raw_lines(data)
@@ -177,7 +157,6 @@ def parse_genbank_bytes(data: bytes, label: str = "GenBank input") -> RawDocumen
         nonlocal current_feature
         if current_feature is None or current_record is None:
             return
-        block = data[current_feature["start"] : end]
         feature = RawFeature(
             record_index=current_record["index"],
             record_id=current_record["record_id"],
@@ -185,7 +164,6 @@ def parse_genbank_bytes(data: bytes, label: str = "GenBank input") -> RawDocumen
             location=current_feature["location"],
             start=current_feature["start"],
             end=end,
-            qualifiers=_parse_qualifiers(block),
         )
         current_record["features"].append(feature)
         current_feature = None
@@ -270,7 +248,17 @@ def parse_genbank_bytes(data: bytes, label: str = "GenBank input") -> RawDocumen
                 "feature_type": _decode(feature_match.group(1), "feature type"),
                 "location": _decode(feature_match.group(2).strip(), "feature location"),
                 "start": start,
+                "has_qualifier": False,
             }
+            continue
+        if current_feature is not None and body.startswith(QUALIFIER_INDENT):
+            continuation = body[len(QUALIFIER_INDENT) :].strip()
+            if continuation.startswith(b"/"):
+                current_feature["has_qualifier"] = True
+            elif not current_feature["has_qualifier"]:
+                current_feature["location"] += _decode(
+                    continuation, "feature location continuation"
+                )

     flush_record(len(data))
     if not records:
@@ -288,11 +276,28 @@ def parse_genbank_bytes(data: bytes, label: str = "GenBank input") -> RawDocumen
                 f"{label}: record {record.record_id!r} declares {record.declared_length} bp "
                 f"but ORIGIN contains {len(record.sequence)} bp"
             )
+    semantic_records = _semantic_records(data, label)
+    if len(records) != len(semantic_records):
+        raise MergeError(f"{label}: raw/semantic record counts differ")
+    for raw_record, semantic_record in zip(records, semantic_records, strict=True):
+        if len(raw_record.features) != len(semantic_record.features):
+            raise MergeError(f"{label}: raw/semantic feature counts differ")
+        for raw, semantic in zip(
+            raw_record.features, semantic_record.features, strict=True
+        ):
+            if raw.feature_type != semantic.type or semantic.location is None:
+                raise MergeError(
+                    f"{label}: raw/semantic feature alignment failed at byte {raw.start}"
+                )
+            raw.semantic_location = semantic.location
+            raw.qualifiers = {
+                key: list(values) for key, values in semantic.qualifiers.items()
+            }
     return RawDocument(data=data, records=records)


-def validate_genbank_semantics(data: bytes, label: str = "GenBank input") -> int:
-    """Require canonical Biopython parsing and return the record count."""
+def _semantic_records(data: bytes, label: str) -> list[SeqRecord]:
+    """Parse one authoritative semantic view, keeping raw bytes separately."""
     try:
         text = data.decode("utf-8")
     except UnicodeDecodeError:
@@ -309,19 +314,23 @@ def validate_genbank_semantics(data: bytes, label: str = "GenBank input") -> int
         raise MergeError(f"{label}: Biopython GenBank parsing failed: {exc}") from exc
     if not records:
         raise MergeError(f"{label}: Biopython found no GenBank records")
-    return len(records)
+    return records
+
+
+def validate_genbank_semantics(data: bytes, label: str = "GenBank input") -> int:
+    return len(_semantic_records(data, label))


 def feature_key(feature: RawFeature) -> tuple[str, str, str] | None:
-    if not feature.locus_tags:
+    values = feature.values("locus_tag")
+    if not values:
         return None
-    unique = list(dict.fromkeys(feature.locus_tags))
-    if len(unique) != 1 or len(feature.locus_tags) != 1:
+    if len(values) != 1 or not values[0].strip():
         raise MergeError(
             f"feature {feature.record_id}/{feature.feature_type} at byte {feature.start} "
             f"has ambiguous locus_tag values: {feature.locus_tags!r}"
         )
-    return feature.record_id, feature.feature_type, unique[0]
+    return feature.record_id, feature.feature_type, values[0]


 def index_unique_features(
@@ -389,12 +398,12 @@ def strict_parity_check(base: RawDocument, source: RawDocument) -> dict[str, Any
             base_tuple = (
                 base_feature.feature_type,
                 tuple(base_feature.locus_tags),
-                base_feature.location,
+                base_feature.location_key,
             )
             source_tuple = (
                 source_feature.feature_type,
                 tuple(source_feature.locus_tags),
-                source_feature.location,
+                source_feature.location_key,
             )
             if base_tuple != source_tuple:
                 problems.append(
@@ -469,15 +478,19 @@ def cds_by_locus(base: RawDocument) -> dict[str, RawFeature]:
     """Build a global CDS locus index, rejecting ambiguous record membership."""
     result: dict[str, RawFeature] = {}
     for feature in base.features:
-        if feature.feature_type != "CDS" or not feature.locus_tag:
+        if feature.feature_type != "CDS":
+            continue
+        key = feature_key(feature)
+        if key is None:
             continue
-        if feature.locus_tag in result:
-            other = result[feature.locus_tag]
+        locus_tag = key[2]
+        if locus_tag in result:
+            other = result[locus_tag]
             raise MergeError(
                 f"GBFF duplicate CDS locus_tag {feature.locus_tag!r} in "
                 f"{other.record_id!r} and {feature.record_id!r}"
             )
-        result[feature.locus_tag] = feature
+        result[locus_tag] = feature
     return result


@@ -556,7 +569,7 @@ def newline_for_offset(data: bytes, offset: int) -> bytes:
 def format_qualifier(key: str, value: str, newline: bytes) -> bytes:
     if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
         raise MergeError(f"invalid qualifier name: {key!r}")
-    if any(character in value for character in "\r\n\t"):
+    if any(ord(character) < 32 or ord(character) == 127 for character in value):
         raise MergeError(f"qualifier /{key} value contains a control character")
     escaped = value.replace('"', '""')
     payload = f'/{key}="{escaped}"'
@@ -711,6 +724,48 @@ def apply_insertions(
     return result, applied


+def verify_insertion_semantics(
+    base: RawDocument, output: RawDocument, applied: list[AppliedInsertion]
+) -> None:
+    """Prove semantic values and targets as well as exact byte reversibility."""
+    if len(base.records) != len(output.records):
+        raise MergeError("semantic round-trip changed record count")
+    targets = {feature.end: feature for feature in base.features}
+    additions: dict[int, list[Insertion]] = collections.defaultdict(list)
+    for item in applied:
+        insertion = item.insertion
+        if insertion.qualifier == "COMMENT":
+            continue
+        target = targets.get(insertion.offset)
+        if target is None or (
+            target.record_id,
+            target.feature_type,
+            target.locus_tag or "",
+        ) != (insertion.record, insertion.feature_type, insertion.locus_tag):
+            raise MergeError("insertion does not identify its exact target feature")
+        additions[target.start].append(insertion)
+    for before_record, after_record in zip(base.records, output.records, strict=True):
+        if before_record.sequence != after_record.sequence or len(
+            before_record.features
+        ) != len(after_record.features):
+            raise MergeError("semantic round-trip changed sequence or feature count")
+        for before, after in zip(
+            before_record.features, after_record.features, strict=True
+        ):
+            if (before.feature_type, before.location_key) != (
+                after.feature_type,
+                after.location_key,
+            ):
+                raise MergeError("semantic round-trip changed feature location/type")
+            expected = {key: list(values) for key, values in before.qualifiers.items()}
+            for insertion in additions.get(before.start, []):
+                expected.setdefault(insertion.qualifier, []).append(insertion.value)
+            if expected != after.qualifiers:
+                raise MergeError(
+                    f"semantic qualifier round-trip failed for {before.record_id}/{before.locus_tag}"
+                )
+
+
 def paths_collide(output: Path, inputs: Iterable[Path]) -> bool:
     output_resolved = output.resolve()
     for input_path in inputs:
@@ -862,9 +917,9 @@ def reconcile_insertions(
     )


-def write_manifest(
+def manifest_bytes(
     path: Path, metadata: dict[str, Any], rows: list[dict[str, Any]]
-) -> None:
+) -> bytes:
     if path.suffix.lower() == ".json":
         payload = (
             json.dumps(
@@ -879,8 +934,7 @@ def write_manifest(
             ).encode("utf-8")
             + b"\n"
         )
-        atomic_write(path, payload)
-        return
+        return payload
     preferred_columns = [
         "entry_type",
         "query_id",
@@ -915,17 +969,77 @@ def write_manifest(
                 for column in columns
             )
         )
-    atomic_write(path, ("\n".join(output) + "\n").encode("utf-8"))
+    return ("\n".join(output) + "\n").encode("utf-8")


-def write_json_sidecar(path: Path, payload: dict[str, Any]) -> None:
-    data = (
+def write_manifest(
+    path: Path, metadata: dict[str, Any], rows: list[dict[str, Any]]
+) -> None:
+    atomic_write(path, manifest_bytes(path, metadata, rows))
+
+
+def json_sidecar_bytes(payload: dict[str, Any]) -> bytes:
+    return (
         json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False).encode(
             "utf-8"
         )
         + b"\n"
     )
-    atomic_write(path, data)
+
+
+def write_json_sidecar(path: Path, payload: dict[str, Any]) -> None:
+    atomic_write(path, json_sidecar_bytes(payload))
+
+
+def preflight_artifact_paths(outputs: Iterable[Path], inputs: Iterable[Path]) -> None:
+    paths, input_paths = list(outputs), list(inputs)
+    for i, path in enumerate(paths):
+        if paths_collide(path, [*input_paths, *paths[:i]]):
+            raise MergeError(f"artifact destination aliases an input/output: {path}")
+        if path.exists() and not path.is_file():
+            raise MergeError(f"artifact destination is not a regular file: {path}")
+        for parent in path.parents:
+            if parent.exists() and not parent.is_dir():
+                raise MergeError(f"artifact parent is not a directory: {parent}")
+        resolved = path.resolve()
+        for other in [*input_paths, *paths[:i]]:
+            other_resolved = other.resolve()
+            if resolved in other_resolved.parents or other_resolved in resolved.parents:
+                raise MergeError(
+                    f"ancestor/descendant artifact path relationship: {path} / {other}"
+                )
+
+
+def stage_artifacts(
+    artifacts: list[tuple[Path, bytes]], inputs: Iterable[Path]
+) -> None:
+    """Stage every artifact before the first replacement; promotions are sequential."""
+    input_paths = list(inputs)
+    paths = [path for path, _ in artifacts]
+    preflight_artifact_paths(paths, input_paths)
+    staged: list[tuple[Path, Path]] = []
+    try:
+        for path, data in artifacts:
+            path.parent.mkdir(parents=True, exist_ok=True)
+            fd, name = tempfile.mkstemp(
+                prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
+            )
+            temporary = Path(name)
+            staged.append((temporary, path))
+            with os.fdopen(fd, "wb") as handle:
+                handle.write(data)
+                handle.flush()
+                os.fsync(handle.fileno())
+        preflight_artifact_paths(paths, input_paths)
+        for temporary, path in staged:
+            os.replace(temporary, path)
+    except OSError as exc:
+        raise MergeError(
+            f"artifact staging/promotion failed; promotions are not a cross-file transaction: {exc}"
+        ) from exc
+    finally:
+        for temporary, _ in staged:
+            temporary.unlink(missing_ok=True)


 def finalize_merge(
@@ -967,6 +1081,9 @@ def finalize_merge(
     )
     validate_genbank_semantics(merged, "merged output")
     parsed_output = parse_genbank_bytes(merged, "merged output")
+    verify_insertion_semantics(
+        parse_genbank_bytes(base_data, "Bakta input"), parsed_output, applied
+    )
     metadata = {
         **metadata,
         "tool_version": TOOL_VERSION,
@@ -984,16 +1101,18 @@ def finalize_merge(
                     "output_sha256": metadata["output_sha256"],
                     "insertions": metadata["insertions"],
                     "records": metadata["records"],
+                    "tool_version": TOOL_VERSION,
                 }
             )
-    atomic_write(output_path, merged)
+    artifacts = [(output_path, merged)]
     if manifest_path is not None:
         rows = insertion_rows(applied)
         if evidence_rows is not None:
             rows = [*evidence_rows, *rows]
-        write_manifest(manifest_path, metadata, rows)
+        artifacts.append((manifest_path, manifest_bytes(manifest_path, metadata, rows)))
     if sidecar_path is not None and sidecar_payload is not None:
-        write_json_sidecar(sidecar_path, sidecar_payload)
+        artifacts.append((sidecar_path, json_sidecar_bytes(sidecar_payload)))
+    stage_artifacts(artifacts, inputs)
     return {
         "output_sha256": metadata["output_sha256"],
         "base_sha256": metadata["base_sha256"],
diff --git a/tests/test_hardening_foundation.py b/tests/test_hardening_foundation.py
new file mode 100644
index 0000000..5db9172
--- /dev/null
+++ b/tests/test_hardening_foundation.py
@@ -0,0 +1,95 @@
+from __future__ import annotations
+
+from pathlib import Path
+
+import pytest
+from test_merge_pipeline import record_bytes, write
+
+import merge_engine as core
+
+
+@pytest.mark.parametrize(
+    "location",
+    [
+        b"join(1..3,7..9)",
+        b"complement(join(1..3,7..9))",
+        b"join(7..9,1..3)",
+        b"<1..>9",
+        b"OTHER.1:1..9",
+    ],
+)
+def test_complete_location_roundtrip(location: bytes) -> None:
+    base = record_bytes("TEST", "T_0001").replace(
+        b"     CDS             1..9", b"     CDS             " + location
+    )
+    wrapped = base.replace(b"1..3,7..9", b"1..3,\n                     7..9")
+    core.strict_parity_check(
+        core.parse_genbank_bytes(base), core.parse_genbank_bytes(wrapped)
+    )
+
+
+def test_semantic_roundtrip_rejects_incorrect_payload(tmp_path: Path) -> None:
+    data = record_bytes("TEST", "T_0001")
+    path = write(tmp_path / "base.gbff", data)
+    feature = core.parse_genbank_bytes(data).features[2]
+    insertion = core.qualifier_insertion(
+        data, feature, "note", "expected", "test", "test", 0
+    )
+    insertion.payload = insertion.payload.replace(b"expected", b"different")
+    with pytest.raises(core.MergeError, match="semantic qualifier round-trip"):
+        core.finalize_merge(
+            base_path=path,
+            base_data=data,
+            output_path=tmp_path / "out",
+            other_inputs=[],
+            insertions=[insertion],
+            manifest_path=None,
+            metadata={},
+        )
+    assert not (tmp_path / "out").exists()
+
+
+def test_all_artifacts_staged_before_replacement(
+    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
+) -> None:
+    first = write(tmp_path / "first", b"old first")
+    second = write(tmp_path / "second", b"old second")
+    original = core.tempfile.mkstemp
+    calls = 0
+
+    def fail_second(*args, **kwargs):
+        nonlocal calls
+        calls += 1
+        if calls == 2:
+            raise OSError("injected staging failure")
+        return original(*args, **kwargs)
+
+    monkeypatch.setattr(core.tempfile, "mkstemp", fail_second)
+    with pytest.raises(core.MergeError, match="injected staging failure"):
+        core.stage_artifacts([(first, b"new first"), (second, b"new second")], [])
+    assert first.read_bytes() == b"old first"
+    assert second.read_bytes() == b"old second"
+    assert not list(tmp_path.glob("*.tmp"))
+
+
+def test_second_promotion_failure_has_documented_boundary(
+    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
+) -> None:
+    first = write(tmp_path / "first", b"old first")
+    second = write(tmp_path / "second", b"old second")
+    original = core.os.replace
+    calls = 0
+
+    def fail_second(*args, **kwargs):
+        nonlocal calls
+        calls += 1
+        if calls == 2:
+            raise OSError("injected promotion failure")
+        return original(*args, **kwargs)
+
+    monkeypatch.setattr(core.os, "replace", fail_second)
+    with pytest.raises(core.MergeError, match="not a cross-file transaction"):
+        core.stage_artifacts([(first, b"new first"), (second, b"new second")], [])
+    assert first.read_bytes() == b"new first"
+    assert second.read_bytes() == b"old second"
+    assert not list(tmp_path.glob("*.tmp"))
diff --git a/tests/test_hardening_regressions.py b/tests/test_hardening_regressions.py
new file mode 100644
index 0000000..151f5b1
--- /dev/null
+++ b/tests/test_hardening_regressions.py
@@ -0,0 +1,99 @@
+from __future__ import annotations
+
+import io
+
+import openpyxl
+import pytest
+from Bio import SeqIO
+from test_merge_pipeline import eggnog_bytes, record_bytes, write
+
+import merge_eggnog_bakta as egg
+import merge_engine as core
+
+ROW = "T_0001\tseed.1\t1e-20\t100\t-\tabc\t-\t-\t-\t-\thhhhhhhhhhhhh"
+
+
+def compound(last):
+    return record_bytes("TEST", "T_0001").replace(
+        b"     CDS             1..9",
+        b"     CDS             join(1..3,\n                     " + last + b")",
+    )
+
+
+def test_A01_changed_location_continuation_rejected():
+    left, right = compound(b"7..9"), compound(b"4..6")
+    core.validate_genbank_semantics(left)
+    core.validate_genbank_semantics(right)
+    with pytest.raises(core.MergeError):
+        core.strict_parity_check(
+            core.parse_genbank_bytes(left), core.parse_genbank_bytes(right)
+        )
+
+
+def test_A01_equivalent_location_wrapping_accepted():
+    left = compound(b"7..9")
+    right = left.replace(b"join(1..3,\n                     7..9)", b"join(1..3,7..9)")
+    core.strict_parity_check(
+        core.parse_genbank_bytes(left), core.parse_genbank_bytes(right)
+    )
+
+
+def test_A02_free_text_uses_biopython_semantics():
+    data = record_bytes("TEST", "T_0001").replace(
+        b'/product="test protein"', b'/product="test\n                     protein"'
+    )
+    bio = list(SeqIO.parse(io.StringIO(data.decode()), "genbank"))[0].features[2]
+    raw = core.parse_genbank_bytes(data).features[2]
+    assert raw.values("product") == bio.qualifiers["product"]
+
+
+def test_A09_duplicate_identical_locus_rejected():
+    data = record_bytes("TEST", "T_0001").replace(
+        b'/protein_id="gnl|Bakta|test"',
+        b'/locus_tag="T_0001"\n                     /protein_id="gnl|Bakta|test"',
+    )
+    with pytest.raises(core.MergeError):
+        core.validate_faa_gbff(
+            core.parse_genbank_bytes(data),
+            {"T_0001": "MK"},
+            ["T_0001"],
+            source_name="test",
+        )
+
+
+def workbook_bytes(name):
+    workbook = openpyxl.Workbook()
+    worksheet = workbook.active
+    worksheet.title = "annotations"
+    lines = eggnog_bytes(ROW.replace("\tabc\t", f"\t{name}\t")).decode().splitlines()
+    worksheet.append(lines[1].split("\t"))
+    worksheet.append(lines[2].split("\t"))
+    stream = io.BytesIO()
+    workbook.save(stream)
+    workbook.close()
+    return stream.getvalue()
+
+
+def test_A10_xlsx_uses_captured_snapshot(tmp_path):
+    path = write(tmp_path / "table.xlsx", workbook_bytes("changed"))
+    captured = workbook_bytes("original")
+    table = egg.parse_eggnog_path(path, expected_version="3.0.0-beta6", data=captured)
+    assert table.hits[0].preferred_name == "original"
+
+
+def test_A11_manifest_directory_failure_preserves_output(tmp_path):
+    base = write(tmp_path / "base.gbff", record_bytes("TEST", "T_0001"))
+    output = write(tmp_path / "out.gbff", b"previous output")
+    manifest = tmp_path / "manifest.json"
+    manifest.mkdir()
+    with pytest.raises((core.MergeError, OSError)):
+        core.finalize_merge(
+            base_path=base,
+            base_data=base.read_bytes(),
+            output_path=output,
+            other_inputs=[],
+            insertions=[],
+            manifest_path=manifest,
+            metadata={},
+        )
+    assert output.read_bytes() == b"previous output"
diff --git a/tests/test_merge_pipeline.py b/tests/test_merge_pipeline.py
index a7359d5..4d95e39 100644
--- a/tests/test_merge_pipeline.py
+++ b/tests/test_merge_pipeline.py
@@ -635,6 +635,9 @@ def test_eggnog_xlsx_rejects_an_empty_header(
     class Workbook:
         sheetnames = ["annotations"]

+        def close(self):
+            pass
+
         def __getitem__(self, key: str):
             assert key == "annotations"
             return worksheet
@@ -642,7 +645,9 @@ def test_eggnog_xlsx_rejects_an_empty_header(
     fake_openpyxl = SimpleNamespace(load_workbook=lambda *_args, **_kwargs: Workbook())
     monkeypatch.setitem(sys.modules, "openpyxl", fake_openpyxl)
     with pytest.raises(MergeError, match="empty column names"):
-        parse_eggnog_path(tmp_path / "annotations.xlsx", expected_version="3.0.0-beta6")
+        parse_eggnog_path(
+            tmp_path / "annotations.xlsx", expected_version="3.0.0-beta6", data=b""
+        )


 def test_eggnog_confidence_and_existing_values_are_not_reinserted(
```
