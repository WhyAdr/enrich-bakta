# enrich-bakta: independent technical and scientific audit

Date: 2026-09-29 (Asia/Jakarta). Audited repository: https://github.com/WhyAdr/enrich-bakta

Pinned main: `5d6754d306e80cddcad7da5fc12b155ab44d2e85`; package version 0.2.0.

## Verdict

Keep the architecture. Do not release the current tree unchanged under its strong identity, conservative-transfer, and complete-provenance claims. The insertion engine preserves bytes correctly in the exercised workflows, but the semantic model used to decide what to insert is not yet reliable for all valid GenBank layouts. In particular, multiline locations are truncated, and restoration can manufacture the very translation identity subsequently used to authorize functional enrichment.

This is an audit, not a patch. No tracked source file or remote branch was changed. Findings below distinguish demonstrated defects, incomplete safety policies, and deliberate limitations. Severity is contextual: **High** means potentially wrong feature/function assignment or a central safety guarantee bypass; **Medium** means incorrect output, provenance, or predictable failed-run behavior with a narrower trigger. None of the tests establishes that the user's existing C14/SM results are wrong: those untracked datasets are absent from this checkout.

## Scope, architecture, and verification

Read all six supported implementation modules, the historical normalization script's role, README, packaging/CI configuration, the complete current test file, both prior reviews and the relevant plans. Examined history through the four hardening commits: `8f2ca84` (byte handling), `a164f7c` (eggNOG contract), `6f2a4df` (provenance/scaling), and `5d6754d` (API/packaging).

Architecture is appropriately small: source parsers -> source planners -> optional reconciliation -> byte splice/reversibility -> Biopython validation -> artifact writes. The main weakness is that raw byte spans and independently hand-decoded biology share one model; Biopython parsing is currently a pass/fail check, rather than the authoritative semantic view attached to those spans.

Executed in a clean Python 3.12.14 environment with Biopython 1.88, pytest 9.1.1, Ruff 0.16.9, mypy 1.20.2 and openpyxl 3.1.5:

- `python -m pytest -q`: **31 passed**.
- `python -m ruff check .`: **passed**.
- `python -m ruff format --check .`: **7 files already formatted**.
- `python -m mypy`: **no issues in 6 source files**.
- `git diff --check`: **passed**.
- Standard isolated wheel build and `pip check` passed; the six supported modules imported from an isolated wheel target and the legacy normalizer was absent. Editable installation succeeded; all five installed commands accepted `--help` and `--version` from outside the checkout.
- Additional observation-based reproduction scripts exercised multiline locations, wrapping, confidence fallback, restoration, conflicts, duplicate identifiers, path failures, malformed transfers, XLSX snapshots, and ordinary successful merges. They are supplied with their captured results, not represented as existing pytest coverage.
- Three-source merge with duplicate Kofam rows and overlapping KO/EC/name evidence: 17 initial insertions; rerun against enriched base: **0 insertions, byte-identical**. Two runs from the same pristine base were byte-identical.
- Same-path, symlink and hard-link collision probes all rejected aliasing.
- An upstream Bakta fixture from `gbouras13/baktfold` commit `827c7201df76350e312bd76dbe9be0aa79c9cb49`, `tests/test_data/assembly_bakta_output/assembly.gbff`, parsed successfully: 205,886 bytes, 2 records, 190 features, 93 FAA entries; all 92 CDSs with translations matched the FAA. Its raw-parser/Biopython qualifier comparison exposed 23 differing values, including products and inferences. Existing bytes are not rewritten by these differences.
- A deliberately lightweight 5,000-locus synthetic planning exercise (15,000 features, 30,000 insertions, 25,000 manifest rows) took approximately 0.70 seconds and process peak RSS was approximately 117 MiB. This is a planning sanity check using short proteins and many records, not a representative large-genome benchmark or a memory guarantee.

## Findings

### A01 — High: multiline locations bypass identity and gene-pair validation

**Affected code:** `merge_engine.py:168–291` (`parse_genbank_bytes`, particularly feature creation at 260–268), `strict_parity_check:343–434`; `merge_eggnog_bakta.py:471–486`.

**Reproduced.** The raw parser stores only the location text on the feature's first line; continuation lines are never appended. These valid CDS locations both become `join(1..3,`:

```text
CDS             join(1..3,        CDS             join(1..3,
                7..9)                           4..6)
```

Both inputs pass Biopython validation. Baktfold parity returns `status=passed`, and an EC is transferred despite different complete CDS locations. Separately, a gene at the first location and CDS at the second are considered a location-matched pair by the eggNOG planner, which writes the preferred name onto both.

**Why it matters:** this violates a central documented identity check, including valid compound/circular-boundary feature layouts. It can also reject equivalent locations if one producer wraps a location and another does not. Simple single-line coordinate tests cannot detect it.

**Least-disruptive fix:** keep byte offsets/splicing, but attach the complete location to each raw feature. Prefer matching each raw record/feature ordinal to the already parsed Biopython record/feature and assert counts/types agree; compare a full location representation preserving strand, part order, fuzzy boundaries and remote references. Do not reduce compound locations to min/max spans. A smaller interim repair can accumulate location continuations before the first qualifier and remove only formatting whitespace.

**Regression gate:** differing continuation parts must fail before any output; equivalent line wrapping must succeed; wrong-location gene/CDS pairs must never receive paired names. Include complement(join(...)), origin crossing and fuzzy boundaries.

### A02 — Medium: raw qualifier decoding removes meaningful spaces

**Affected code:** `merge_engine.py:136–165`, particularly `"".join(piece.strip() for piece in pieces)` at line 145; planners consuming `RawFeature.values()`.

**Reproduced, including a real upstream Bakta fixture.** A wrapped value understood by Biopython as `similar to AA sequence:UniRef:...` becomes `similar to AAsequence:UniRef:...` in the raw model. A synthetic wrapped Baktfold gene value ending one line with `iota` and starting the next with `kappa` was transferred as `iotakappa`, even though the source value was `iota kappa`.

**Why it matters:** byte preservation leaves existing source qualifiers intact, but source-derived transferred values, duplicate detection and provenance comparisons can use the wrong semantic value. The current formatter round-trip tests ask Biopython whether the new text is correct; they do not test whether the merger's own parser reads it identically. This is not proof that ordinary short gene symbols are currently corrupted.

**Least-disruptive fix:** obtain semantic qualifier values from Biopython while preserving raw spans for output. Do not indiscriminately join every continuation with a space: translations and some sequence-like tokens require different handling. Assert semantic equivalence for newly emitted qualifiers, separately from byte reversibility.

**Regression gate:** multiline free text, escaped quotes, sequences, wrapped inferences and existing structured tokens; compare both semantic views and rerun behavior.

### A03 — High, workflow safety: restoration has no independent sequence-identity proof

**Affected code:** `restore_bakta_translations.py:30–117`, especially 61–94; `merge_engine.py:484–531`; README translationless-pseudogene workflow.

**Reproduced.** A translationless `/pseudogene="unitary"` CDS spanning nine nucleotides accepted an FAA entry containing 100 tryptophans under the same locus tag. The restored GBFF then passed ordinary eggNOG FAA/GBFF matching and acquired an EC and KO. No nucleotide compatibility check, input-lineage proof or FAA alphabet validation prevented this.

The regular FAA/GBFF validator is doing its stated comparison correctly. The problem is circular validation: the restoration step imported the reference sequence used by the subsequent comparison. A locus label plus `/pseudogene` does not establish that a sequence came from the same Bakta run or is supported by that genomic locus. Preserving `/pseudogene` is good, but does not remedy an unrelated `/translation` value.

**Policy versus bug:** copying an explicitly selected FAA into a separate research GBFF is already documented as a deliberate repair policy. The defect is the unqualified “matched FAA” safety interpretation and absence of a guard against obviously unrelated sequences. Do not claim that conceptual translation of every pseudogene should pass an intact-ORF `cds=True` test; genuine disruptions would fail it.

**Least-disruptive fix:** make the exception an explicit, recorded evidence-import policy. Prefer accepting FAA evidence for a pseudogene through a dedicated planner/sidecar without inventing an ordinary translation. If the restored-GBFF workflow is retained, require source-run provenance where available (matched Bakta JSON plus nucleotide/feature hashes), validate supported alphabet and independently assess the CDS-extracted sequence with its genetic code/offset/exceptions. Unresolved frameshift cases need an explicit unverified classification, not an ordinary identity-passed label. Add an in-file import note and mandatory restoration provenance. Preserve pseudo status and avoid interpreting EC/KO transfers as intact enzymatic capacity.

**Additional reproduced edge:** an existing `/translation=""` is treated as absent and receives a second translation qualifier. Restoration succeeds, but the next enrichment rejects `['', 'MK']` as ambiguous. In an insertion-only design, fail closed on empty/duplicate translation qualifiers rather than appending another.

**Regression gate:** unrelated same-tag FAA, impossible protein length, invalid symbols/internal stops, correct partial or exceptional conceptual translations, empty/duplicate qualifiers, `/pseudo` versus `/pseudogene`, and explicit manifest proof that import and independent validation are different states.

### A04 — High hardening gap: unknown eggNOG versions silently inherit a positional contract

**Affected code:** `merge_eggnog_bakta.py:277–365`, 368–439 and 494–500.

**Reproduced acceptance, prospective misassignment risk.** TSVs declaring `99.0.0`, `2.1.12`, or no version are accepted if their columns/vector look superficially compatible. Each uses `documented_default` when the legend is absent. XLSX accepts any supplied version and likewise assumes the same order. `--eggnog-version` compares strings for TSV; it does not select or validate a supported schema. An empty declared legend also falls through the truthiness-based fallback.

**What is sound:** the currently documented v3 order is correct; the 13-character alphabet and length checks are correct; a present nonempty reordered legend is rejected; known scored-column relative order is checked. This is not the old review's incorrect assertion that the v3 confidence field does not exist. Upstream documents that confidence describes donor-cascade tier, not an experimentally calibrated correctness probability.

**Why it matters:** a future producer could retain 13 positions while changing their meaning, especially in a comments-stripped table or XLSX. The parser currently has no basis to claim the old default applies. These reproductions demonstrate acceptance without a supported contract, not an observed error in today's beta6 mapping.

**Least-disruptive fix:** register tested version/schema contracts. Allow fallback only for a documented known version or an explicit recorded schema selection when comments are absent. Reject an explicitly empty legend and conflicting/repeated version declarations. For unfamiliar versions, require a compatible explicit legend or refuse; do not silently assign confidence. Keep v2 unsupported unless a separately designed policy exists. Check confidence/field inconsistencies with useful diagnostics without inventing numeric confidence probabilities.

**Regression gate:** unknown major versions; missing metadata with and without explicit schema; empty, duplicated and reordered legends; repeated inconsistent version lines; TSV/XLSX parity; all 13 positions using one-hot confidence vectors.

### A05 — Medium: reconciliation leaves misleading final provenance and candidate statuses

**Affected code:** `merge_engine.py:761–862`; `enrich_bakta.py:143–168`; `merge_eggnog_bakta.py:708–817`; Baktfold provenance planning at `graft_baktfold_additions.py:166–204`.

**Reproduced.** A blank locus with Baktfold name `abc` and eggNOG name `xyz`, default `skip`, emits no gene names. Nevertheless, its eggNOG candidate says `status="emitted"` and lists both planned gene qualifiers. Baktfold's gene-symbol note and eggNOG inference remain even when the contested name was each source's only functional addition.

**Why it matters:** insertion rows correctly describe the actual byte changes, and a reconciliation row does record the conflict, so a careful reader can infer the outcome. But the candidate row describes pre-reconciliation planning as final emission, and source-dependent provenance is not reconciled with the annotations it was meant to support. This undermines the requested machine-readable “why emitted/skipped/suppressed” trail.

**Least-disruptive fix:** give candidates stable identities and separate `planned_status` from `final_status`, linking surviving insertion IDs. Reconcile functional candidates before deriving producer provenance and count comments. If provenance is intentionally retained for a suppressed hit, label it as suppressed/conflicting evidence. Record which actual feature(s) each source supported; do not remove legitimate support for a surviving duplicate.

**Regression gate:** skip and both preference policies end-to-end; conflict-only sources; shared annotation support; manifest-to-parsed-output consistency; final counts and comments.

### A06 — Medium: Baktfold gene transfer can contradict an authoritative paired feature

**Affected code:** `graft_baktfold_additions.py:100–145`; `merge_engine.py:777–806`.

**Reproduced.** With Bakta `/gene="authoritative"` on a gene feature and a blank CDS at the same locus, Baktfold name `other` is inserted into the CDS. The original qualifier is not overwritten, but the pair now has conflicting primary names. Baktfold planning checks blankness separately per feature; it lacks eggNOG's paired-gene guard.

**Why it matters:** “never overwrite Bakta” is literally upheld but does not prevent creating a contradictory gene/CDS interpretation. The current unified conflict logic only sees proposed additions, not authoritative existing values. It also needs to distinguish an internally inconsistent Baktfold pair from disagreement between sources: preferring Baktfold cannot resolve two different Baktfold names within one locus.

**Least-disruptive fix:** collect authoritative gene symbols across a location-validated pair before planning Baktfold names. If a source name disagrees, skip and report it. Do not repair existing biology silently. Reject or explicitly suppress inconsistent source-pair names before applying cross-source preferences.

**Regression gate:** name present only on gene, only on CDS, conflicting base pair, inconsistent source pair, missing pair, and each conflict policy.

### A07 — Medium: value validation is inconsistent and can emit bogus or empty annotations

**Affected code:** `graft_baktfold_additions.py:44–62,120–144`; `merge_eggnog_bakta.py:489–491,672–683,713–745`; `merge_engine.py:442–465` for FAA validation scope.

**Reproduced.** Baktfold `/EC_number="not-an-EC"` and `pdb:` with an empty identifier are accepted, transferred and pass output parsing. eggNOG `Preferred_name=NA` becomes `/gene="NA"`, despite `NA` being treated as missing in other annotation fields. `_123` with suffix cleaning becomes an empty `/gene`; a second run adds two more empty gene qualifiers.

**Why it matters:** a namespace allowlist does not validate identifier syntax, and Biopython does not perform biological qualifier validation. The empty-name case directly breaks idempotence. An upstream malformed annotation should not acquire apparent legitimacy merely by passing through this merger.

**Least-disruptive fix:** share field-specific validators across sources. At minimum validate full/partial EC syntax and nonempty structural identifiers; explicitly choose and report partial-EC policy for Baktfold rather than accidentally inheriting arbitrary strings. Apply missing-value handling to names before cleaning and reject/skip empty results afterward. Preserve raw values and reasons. Do not impose a simplistic three-letter bacterial-name regex or require online ontology lookup for every merge.

**Regression gate:** the reproduced values; partial EC versus malformed EC; missing sentinels; suffix-only names; structural prefixes with empty or whitespace-only identifiers; supported protein symbols. Tighten syntax without pretending that syntactic validity proves current ontology membership or enzymatic specificity.

### A08 — Medium: eggNOG note deduplication uses partial identifier matches

**Affected code:** `merge_eggnog_bakta.py:453–468`.

**Reproduced.** Existing `/note="notGO:00000010"` suppresses candidate `GO:0000001` as `existing`. The regex matches substrings without left/right identifier boundaries. Related prefix problems apply to KO and CAZy tokens (including a family matching the beginning of a longer family/subfamily token).

**Why it matters:** a valid new annotation is silently lost and its manifest incorrectly claims it already exists. Kofam's KO-note regex has better boundary handling, so behavior is source-dependent.

**Least-disruptive fix:** extract complete structured tokens with shared explicit boundaries or parse recognized note conventions. Distinguish an identifier occurrence from a negated/free-text claim if broad prose scanning is retained; do not treat arbitrary contextual mention as experimentally supported annotation.

**Regression gate:** GO/KO numeric suffixes, preceding letters, punctuation delimiters, CAZy family/subfamily cases, and matching behavior across Kofam/eggNOG.

### A09 — Medium: input identity validation is not consistently fail-closed

**Affected code:** `merge_engine.py:315–339,468–481`; `merge_eggnog_bakta.py:289–315`.

**Reproduced.** Duplicate identical `/locus_tag` qualifiers on a CDS are accepted by eggNOG/FAA mapping, whereas Baktfold's `feature_key()` rejects the same multiplicity. `cds_by_locus()` uses the property that deduplicates values rather than validating cardinality. Ambiguous unrelated CDSs can also be skipped instead of diagnosed. These are malformed-input gaps, not evidence that ordinary globally unique Bakta locus tags are unsafe.

A separate parser edge survives the reserved-ID patch: a data-shaped TSV row whose query starts `##` is discarded as a comment before query validation, while a `#` row is rejected. Reproduced a table containing a valid row and `##bad\t...`; only the valid row remained, without a skipped-row diagnostic.

**Least-disruptive fix:** centralize feature-ID cardinality validation for all paths before building indexes. Detect data-shaped reserved-marker rows after the table header while continuing to allow recognized producer footer comments. Preserve the existing rejection of duplicate CDS IDs across records.

**Regression gate:** identical and differing duplicate qualifiers, blank tags, duplicate CDS IDs across records, missing tags on targeted versus untargeted features, `#`, `##` and `*` query prefixes, valid producer comment/footer lines.

### A10 — Medium: XLSX parsing can use different bytes from its recorded hash

**Affected code:** `merge_eggnog_bakta.py:368–447`, plus XLSX callers that first read `eggnog_data`.

**Reproduced at the parsing API.** Read XLSX bytes containing `K00001`, replace the on-disk workbook with `K00002`, then call `parse_eggnog_path(path, data=captured_bytes, expected_version=...)`: it parses `K00002`. The XLSX branch ignores supplied bytes and reopens the path, while upstream callers retain the original byte buffer for hashing.

**Why it matters:** immutable inputs behave correctly, but a concurrent writer can make the manifest fingerprint describe a workbook other than the one used for annotations. The earlier snapshot-binding fix covers TSV, not XLSX.

**Least-disruptive fix:** load openpyxl from `io.BytesIO()` of the same captured bytes used for hashing; close the workbook in a `finally` block. Normalize malformed workbook exceptions to a clear `MergeError`. Add header formula validation if the promise remains “no formulas anywhere.”

**Regression gate:** deliberately alter the path after snapshot capture; prove snapshot annotations and hash stay bound; malformed archives, formulas, required headers, unknown versions and workbook closure.

### A11 — Medium, known transaction limitation: a failed sidecar write can replace the output

**Affected code:** `merge_engine.py:931–1004`, especially writes at 989 onward; `paths_collide:714–725`.

**Reproduced.** With an existing output and a manifest path that is a directory, enrichment first replaces the output, then raises `IsADirectoryError` while writing the manifest. The old output is gone despite overall failure. Output/manifest ancestor-descendant paths likewise are not rejected by equality/alias collision checks and can trigger this class of late failure.

**What is sound:** each individual write uses a sibling temporary, file flush/fsync and atomic replace. The earlier review explicitly deferred multi-artifact transactions. README says atomic output, not a crash-proof cross-file transaction. This is therefore not a broken `os.replace` implementation or justification for claiming source corruption.

**Least-disruptive fix:** serialize all artifacts and preflight destination types plus ancestor/descendant conflicts before promoting anything. Stage every artifact before promotion and give deterministic operational errors. State the remaining cross-file crash boundary. A generation directory plus a final pointer/completion record is a later design if true transaction semantics are needed; multiple renames alone are not an atomic transaction.

**Regression gate:** existing manifest directory, output parent that is a file, ancestor paths, sidecar failure after staging, injected replace errors, retention of previous output on preflight failure. Directory fsync for crash durability is separate non-blocking work.

### A12 — Medium: manifests are useful, but not a complete decision reconstruction

**Affected code:** `graft_baktfold_additions.py:65–272`; `merge_engine.py:395–434,865–918`; `merge_eggnog_bakta.py:719–729,847–876,933–941`; unified metadata selection.

**Reproduced/inspection.** Baktfold emits insertion rows but no complete candidate decision ledger for existing or suppressed source values. Translation disagreements are deliberately diagnostic in the plan (`merge-koala-baktfold-plan.md:68–69`), yet those diagnostics are absent from written manifests: a source translation `WW` versus base `MK` still contributes an EC, `parity.status` is `passed`, and the manifest contains no mismatch diagnostic. Standalone graft returns the diagnostic in stats; unified enrichment does not retain it in its result/manifest.

Additional precise gaps:

- Missing/ambiguous paired genes get `existing_gene`, because the sentinel `["missing paired gene"]` makes the next `unpaired_gene` branch unreachable. Reproduced.
- Standalone eggNOG manifest metadata omits `clean_gene_suffix`; provenance flags are not uniformly recorded. Columns/format are in the optional context sidecar but not uniformly in manifest metadata. The context sidecar itself lacks the tool version stamped onto manifests.
- TSV nested values are Python string representations rather than typed JSON; this is convenient for display but unsuitable as the canonical interchange format for full reconstruction.
- Restoration records imported loci, but not the requested already-translated loci that passed validation.

**Why it matters:** hashes establish artifact identity, not why a rejected proposal disappeared. Candidate states must also describe the final result (A05). Persisting a source hash is valuable but requires retaining the corresponding source; it does not replace a decision ledger.

**Least-disruptive fix:** write Baktfold candidate outcomes, parity diagnostics, complete policy flags and schema/format information. Check missing pairing before existing genes. Include tool version in every standalone sidecar. Make JSON the documented canonical provenance format; encode nested TSV cells as JSON if machine reconstruction is promised. Record validated/no-op restoration targets.

**Policy recommendation, not a rediscovered implementation bug:** Baktfold translation mismatches were deliberately allowed by the original plan. Consider suppressing protein-derived transfers for mismatching CDSs or requiring an explicit recorded permissive policy. Same genomic coordinates alone do not establish that the source hit was generated from the same protein. This recommendation changes policy; the present loss of diagnostics is the definite provenance defect.

## Scientific interpretation and intentional limitations

The current v3 confidence mapping is supported by upstream documentation, including PFAM at position 12. Filtering unscored OG identifiers independently from scored annotations is intentional. `--min-eggnog-confidence high` does not mean every emitted field has high scored confidence; README already explains the exceptions. Keep this behavior explicit rather than fabricating confidence for COG/OG fields.

EggNOG transfers annotations through orthology. The seed ortholog identifies the seed search relationship; it need not be the exact donor for every annotation field in a cascade. The current inference text describes sequence similarity and is not necessarily false, but it should not be advertised as full per-field donor provenance. Retain producer version/database/options and donor information where available. An optional producer query MD5, if present, is an opportunity to check that the annotation table was generated from this FAA; matching FAA to GBFF alone cannot prove that table provenance.

GO/KO/EC/CAZy values are computational assignments. Different ECs or KOs can legitimately coexist (multifunctionality, overlap or alternative hypotheses); do not resolve them by “highest score wins” across tools, whose scores are not directly comparable. A non-blocking disagreement/context report would help users distinguish multisource support from alternatives. Correct syntax is not verification of a current ontology record or catalytic specificity.

Kofam retaining starred hits, including multiple KOs, is sensible. Definition text is not promoted to product or EC, and score/e-value notes remain even when producer inference is disabled. Do not revive the old proposed prohibition on negative HMM scores/cutoffs: those values are not inherently invalid. A lower numeric score than a displayed threshold is not by itself proof of corrupt output because producer threshold scaling may be configured. The parser is fixed-position after recognizing a header, not a general arbitrary-column-order parser; document that scope or validate the exact expected order.

Baktfold structural xrefs are references to structural evidence, not proof of biochemical activity. The merger trusts the upstream Baktfold annotation policy and does not independently examine Foldseek coverage, alignment, database version or enzyme active-site conservation. That is an intentional scope boundary; a GBFF alone may not contain the information required to re-audit every upstream transfer. Accepting arbitrary EC strings, however, is a correctable validation gap (A07).

`protein structure similarity:Baktfold Foldseek:...` is a useful research description but is not among the inference types listed in the INSDC feature-table reference examined. Kofam's `profile` and eggNOG's `similar to AA sequence` use recognized types. The README already disclaims formal submission compliance for xrefs; extend this warning to inference vocabulary and restored pseudogene translations as appropriate. Do not convert structural evidence into a misleading sequence-homology statement solely to satisfy syntax. Preserve research provenance in notes/sidecars, or design a separate submission-export mode.

Translationless `/pseudo` is not accepted by restoration; only `/pseudogene` is. This is an explicit documented limitation, not a hidden bug. Likewise, v2 eggNOG files lacking the required confidence vector are intentionally unsupported, exact record/feature ordering is intentionally strict, and these tools do not promise to retranslate every ordinary CDS from DNA. Global CDS locus-tag uniqueness is appropriate for FAA IDs without record namespaces.

The broad README sentence that “pseudogene evidence is rejected” should instead say that evidence without a matching translation is rejected by the ordinary path: a translated pseudogene is currently accepted, as the documented restored workflow itself requires.

## What is solid and should not be churned

1. Preserve original bytes and splice only planned additions; do not replace the writer with whole-record Biopython serialization.
2. Keep the same captured base bytes from planning through splicing. The old base reread corruption is fixed; its regression test passes.
3. Keep byte reconstruction and hash equality as a separate preservation invariant. Strengthen semantic checking alongside it, not instead of it.
4. Keep fatal Biopython parser warnings, samefile-aware collision checks, record sequence hashes, ordered parity and no force bypass.
5. Keep source planners separated from reconciliation and finalization. Improve the data contract between them instead of adding a new framework.
6. Keep current v3 confidence positions, partial-EC nonpromotion for eggNOG, one raw row per input query, context-only higher-order annotations and exact duplicate collapse.
7. Keep no implicit timestamp, stable input-order behavior, multisource support and tested ordinary idempotence. Determinism means identical inputs/options give identical output; arbitrary source-row permutations need not be byte-identical because row provenance and source hashes legitimately change.
8. Keep the paired-gene index and ordinary in-memory implementation. No measured need for streaming reconstruction, multiprocessing or a wholesale scale rewrite emerged.
9. Packaging/CLI separation is sensible. Historical `normalize_baktfold.py` is explicitly excluded and should not be put back into the supported enrichment path.

## Release decision and prioritized implementation handoff

**Release-blocking correctness work:** A01; A02 for the semantic model feeding transfers; A03 unless the restoration command is explicitly restricted/disabled pending safe policy; A04 for unsupported confidence contracts; A05 for final decision correctness. The very small A07 empty-gene/idempotence fix should accompany these. If malformed-source robustness and manifest completeness remain advertised central guarantees, finish A06–A10 and A12 in the same hardening release rather than claiming they are already handled.

A11 can remain a documented limitation for true cross-file crash atomicity, but its predictable destination preflight failures are inexpensive to fix now. Policy changes such as rejecting all Baktfold protein mismatches should be separately explicit and tested against legitimate pseudogene differences rather than slipped into refactoring.

| Priority | Patch scope | Concrete acceptance criteria |
|---|---|---|
| 1 | Semantic input model (A01/A02/A09 cardinality) | Complete locations and authoritative qualifiers align to raw spans; adversarial pairs cannot pass; original bytes remain identical outside insertions. |
| 2 | Restoration evidence policy (A03) | Imported sequence is never mislabeled as independently validated; impossible/malformed inputs fail; empty/duplicate translations fail before output; pseudo status and provenance survive. |
| 3 | eggNOG schema + value gates (A04/A07/A08/A09 reserved rows/A10) | Unknown schema cannot inherit confidence silently; no empty names or partial-token suppression; XLSX consumes the hashed snapshot. |
| 4 | Pairwise names + reconciliation ledger (A05/A06/A12) | Existing gene/CDS authority is respected; candidate final states agree with parsed output; suppressed-only evidence is labeled or removed; every policy/diagnostic is persisted. |
| 5 | Artifact staging and operational errors (A11) | Invalid output destinations are detected before old output is replaced; per-file atomicity and remaining transaction limits are accurately documented. |
| 6 | Integration/release gate | Run full existing gates plus new regressions; exercise known C14/SM datasets locally; inspect scientific annotation deltas, not just counts; verify deterministic reruns and packaged CLI behavior. |

Use small commits corresponding to these phases. Introduce failing regression tests before fixes. Avoid changing source transfer policy and byte formatting in the same opaque patch. The source tree's existing 31 tests should remain green except where a test fixture needs to state its schema or restoration policy explicitly.

## Missing regression coverage to prioritize

- Full compound locations, line wrapping, complement, origin crossing, fuzzy boundaries and semantic/raw parser agreement.
- Every confidence position, unknown/absent version contracts, malformed legends and real producer comments/footers.
- End-to-end gene conflicts under all policies, authoritative half-populated pairs, inconsistent source pairs and final candidate statuses.
- Invalid Baktfold EC/structural IDs, missing gene sentinels, suffix-only names and rerun idempotence.
- Wrong-run restoration inputs, independent nucleotide/provenance checks, empty/duplicate translations and legitimate disrupted pseudogenes.
- Duplicate qualifier cardinality and multi-record FAA ambiguity across every CLI.
- XLSX successful round trips, snapshots, formulas, malformed files and resource closure (current tests predominantly cover preliminary XLSX rejection paths).
- Full-token deduplication for GO/KO/CAZy and wrapped qualifiers.
- Preflight/staging failures with a pre-existing output, symlink/hard-link/ancestor paths and injected write/replace errors.
- Real installed-wheel imports and execution; declared minimum Biopython version as well as current version. The current CI installs editable mode and does not prove wheel contents in an isolated consumer environment.

## Source grounding

Repository code references above are pinned to the audited commit, available at:
https://github.com/WhyAdr/enrich-bakta/tree/5d6754d306e80cddcad7da5fc12b155ab44d2e85

Primary references consulted on 2026-09-29:

- eggNOG-mapper v3 usage and confidence legend: https://github.com/eggnogdb/eggnog-mapper/blob/main/USAGE.md
- Real v3 producer fixture: https://github.com/eggnogdb/eggnog-mapper/blob/main/tests/fixtures_v7/test_diamond.emapper.annotations
- INSDC feature-table definition, version 11.4 (April 2026), inference/locus-tag/translation sections: https://www.insdc.org/submitting-standards/feature-table/
- Baktfold implementation and input fixture: https://github.com/gbouras13/baktfold/tree/827c7201df76350e312bd76dbe9be0aa79c9cb49
- KofamScan usage and output conventions: https://github.com/takaram/kofam_scan/blob/master/README.md

The supplied scripts preserve observations at this commit; they are not production patches, biological datasets, or a replacement for focused pytest assertions. No current-user dataset was available to revalidate historical isolate-level counts or quantify affected real loci.
