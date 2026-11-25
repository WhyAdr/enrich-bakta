# InterProScan integration audit and Gemini remediation handoff

- Date: 2026-10-05 (Asia/Jakarta)
- Repository: [WhyAdr/enrich-bakta](https://github.com/WhyAdr/enrich-bakta)
- Audited main: a34a2076da453961bbe29eec6664df7a2f7eb7cb
- Package: 0.4.0; merge-manifest v2; candidate-ledger v3
- Review scope: six commits after e5e0499, including the adopted plan and Phases A–E.
- Execution recipient: Gemini. This is an audit and repair specification, not an applied implementation patch.

## Verdict

The integration works on the supplied C14 and restored-SM corpus. Its main
scientific policy is conservative: InterPro/GO feature references and
allowlisted member notes are added, pathway context stays out of feature
qualifiers, and the tested imported origins survive unified merges.

It is **not ready to call complete and robust**. There are reproducible bugs in
the standalone restored-input path, legacy restoration protection, and FAA
snapshot consistency; scalability defects; and unfinished validation,
documentation and release gates. Successful corpus runs do not cover these
failure modes. No incorrect published biological results were established by
this audit; distinguish the passing supplied-data cases from the defects below.

P1 means close before declaring this integration complete. P2 means required
hardening, contract completion or usability work. These are engineering
priorities, not security vulnerability scores.

| ID | Priority | Finding | Evidence |
|---|---|---|---|
| F1 | P1 | Standalone restored-input manifest serialization crashes | Reproduced TypeError with a valid bound imported-origin ledger |
| F2 | P1 | Legacy protection misses late comments and accepts mixed unverified legacy/current-marker inputs | Two independent successful-bypass reproductions |
| F3 | P1 | Unified workflow can bind one FAA hash while the adapter validates different FAA bytes | Reproduced manifest/context SHA disagreement |
| F4 | P2 | Parser retains all hit/pathway objects; candidate limit is a hit limit; producer support lookup repeatedly scans the full ledger | Measured linear duplicate-data memory growth and inspected the quadratic lookup |
| F5 | P2 | Required producer version is silently assumed | Omitted flag accepted and 5.59-91.0 recorded |
| F6 | P2 | Malformed enabled member accessions silently disappear | Invalid PF000010 accepted, member proposal and invalid-value diagnostic absent |
| F7 | P2 | Exact one-CDS UID contract is not enforced centrally | Missing UID accepted; conflicting query/UID yields a two-target accepted decision |
| F8 | P1, unfinished gate | New scientific checks are not wired into CI or comprehensive goldens | Workflow and actual job results inspected |
| F9 | P2 | README standalone command fails; version switch and several documented semantics are wrong | Positional invocation and standalone --version both exit 2 |
| F10 | P2, unfinished contract | Source/context schemas and compact match support are incomplete | Corrupted source projection passes JSON Schema; no context schemas; support loss visible in code |

## What was independently verified

The checkout was isolated from the prior uncommitted plan edit. Source imports
and installed package metadata were confirmed to resolve to the audited 0.4.0
implementation. Previously materialized, unchanged corpus bytes were reused
and revalidated against the audited PUBLISHED-MANIFEST. Unmaterialized LFS
pointers in a new worktree are not implementation defects.

- 196 tests passed in 2.65 seconds.
- Ruff check, Ruff format (46 files), mypy (15 source files), and git diff
  --check passed.
- All 42 published artifacts passed size/SHA validation, including nine
  InterProScan exports. No published input was rewritten by this review.
- Existing real-data scripts passed all 22 checks: 13 first-pass/compatibility
  checks and nine lineage/no-op chain checks.
- All four new InterProScan real-data script checks passed: C14 acceptance,
  strict SM rejection, restored SM acceptance, and BK71A legacy rejection.
- Additional InterProScan-only C14/SM output-as-input runs had zero insertions,
  identical GBFF bytes and unchanged source support/pathway digests. SM retained
  all six imported origins.
- Canonical genbank-parser 0.9.3 output QC matched the prior inputs' diagnostics:
  one informational exceptional-translation annotation in each genome,
  OJMJKD_000205 and NNLMDO_000200 respectively. No new QC diagnostic appeared.
  This is not independent proof of every encoded CDS translation.
- All-four-source first-pass and no-op runs passed on both C14 and SM. The
  reruns were byte-identical with zero insertions and retained all four C14
  and six SM imported origins. Across the old/new scripts and these extra
  reruns, 32 real-data checks passed in total.
- GitHub push CI [run 37217008696](https://github.com/WhyAdr/enrich-bakta/actions/runs/37217008696)
  succeeded. Its dataset-bytes and scientific-reruns jobs were both skipped.
  The observed pinned-head run was a push run; no manual/tag scientific run
  was observed for this revision.

Current independently observed InterProScan-only first-pass outputs:

| Metric | C14 | Restored SM |
|---|---:|---:|
| Insertions, including COMMENT | 34,080 | 33,542 |
| Candidate decisions | 37,942 | 37,412 |
| Emitted candidates | 34,079 | 33,540 |
| Supported-existing candidates | 3,863 | 3,872 |
| Producer-provenance candidates | 4,390 | 4,304 |
| Imported origins carried | 0 on pristine input | 6 |

C14 GBFF SHA-256: 472b0bd60d3990dda6307c89de7153d2ab3fbfd9f9324180e8889cb1b47bc9e1  
C14 context SHA-256: 2113c50f83e35df548f6cbc0087750d925b30153edfb7aa04b64b093e4a9adc9  
SM GBFF SHA-256: 0df5dce4a2aa615e313eaf858d948e8ff8b8698cc11ffc4bca444a2860d4e471  
SM context SHA-256: 36de2cf21aefbffdb3b721eb8fb703f50b317584349d026884579047bc3598ae

All-four-source observations, using the freshly reproduced import-faa restored
inputs and the existing scientific policy flags:

| Metric | C14 | SM |
|---|---:|---:|
| First-pass insertions | 74,910 | 73,889 |
| Candidate decisions | 98,688 | 97,420 |
| First-pass seconds, this audit environment | 207.79 | 210.86 |
| No-op seconds | 31.67 | 32.42 |
| No-op insertions | 0 | 0 |
| Imported origins retained | 4 | 6 |

C14 four-source GBFF SHA-256: dbe840ad739a66f6259544ec65071fc70b40f931f2a7e54b57618e9e32a4d908

SM four-source GBFF SHA-256: 0cbfa457be320422f02a88e2ce49707cd142e0a8ccc403f7682d93d228714bea

These are audit observations against 0.4.0. They are not permission to
automatically adopt new goldens after a repair. Full-file hashes, annotation
sets, candidate semantics and producer/policy metadata must be reviewed together.

## F1 — Standalone restored-input serialization

Code: sources/interproscan.py, merge(), especially lines 917–925; finalizer
metadata binding and manifest_bytes() in core/merge_engine.py.

merge() assigns the TranslationEvidenceLedger directly to
metadata["translation_evidence"]. Its values are TranslationEvidence dataclass
objects. The finalizer's binding logic expects dictionaries, and JSON
serialization cannot encode those objects.

A genuine import-faa restoration, its valid bound parent manifest, explicit
allow-imported-translations and a matched TSV produce:
“Object of type TranslationEvidence is not JSON serializable”.
No GBFF or manifest was committed in this reproduction, which is the correct
failure containment; the legitimate advertised workflow still cannot complete.
The unified path works because it serializes each evidence value with as_dict().

Repair:

1. Serialize every verified origin through as_dict(), retaining the full ledger,
   not only InterProScan's queried subset.
2. Preserve the captured parent manifest SHA, imported-origin policy, FAA hash
   and source/profile metadata consistently with unified enrichment.
3. Prefer a common workflow helper or delegation over maintaining two divergent
   lineage implementations.
4. Add standalone/API/console tests for restored imported and validated origins,
   positive and negative opt-in, untouched-origin carry-forward, correctly bound
   output/parent hashes and output-as-input reruns.
5. Errors should use the normal contextual MergeError path; do not merely teach
   JSON to serialize arbitrary dataclasses without fixing evidence binding.

## F2 — Legacy lineage protection has two bypasses

Code: core/merge_engine.py lines 741–748 and sources/interproscan.py lines
468–476.

The detector searches only data[:100_000]. A valid synthetic GenBank COMMENT
containing the known legacy signature at byte 130,537 was not detected.
Enrichment succeeded without a legacy lineage ledger.

Separately, the guard is conditioned on *absence of any current marker*.
A two-record input with a genuine current restoration marker/ledger for
T_0001 and a legacy record for OJMJKD_006675 was accepted. The bound ledger
contained only T_0001; OJMJKD_006675 was unverified. Both the legacy signature
and current marker were detected. This reproduces a scope failure, not an
invalid/tampered parent hash.

Repair:

- Inspect COMMENT/header provenance for every record using the existing parsed
  record boundaries. Avoid an arbitrary file-prefix cutoff and regex matching
  across unrelated records or feature notes.
- Establish eligibility for the queried legacy record/protein. A file-wide
  current marker plus a partial valid ledger does not certify other historical
  translations.
- Fail closed where known legacy provenance cannot be resolved to genuine
  appropriate inputs and the established evidence contract. Do not fabricate
  retrospective origins or use filenames as proof.
- Test late comments, multi-record order, mixed marked/unmarked legacy records,
  valid partial ledgers, unrelated mentions and negative prose.
- Preserve the real BK71A failure. Its absent genuine FAA/verified history is a
  data limitation, not a feature implementation task Gemini should “repair”
  by synthesizing an FAA from JSON or encoded GBFF translations.

## F3 — FAA identity and provenance use different snapshots

Code: workflows/enrich.py first FAA read/hash near lines 109–118, followed by
plan_interproscan(..., faa_path, ...) at lines 272–275; the adapter reopens FAA
at lines 478 onward.

The unified workflow reads FAA once for its metadata/other adapters and the
InterProScan planner reads the path again. I changed the actual FAA on disk
immediately before the second read. The merge passed protein parity but wrote:

- Manifest FAA SHA: b101dc75fceab35a60f1632830a2cee3439521d03de1d98d5e7d7b3556bddcdd
- Context and actually validated FAA SHA:
  ca0b124a60add5abf78f743dfbbe917d349fb5215b10368a18c90b8ac7311a95

This is a controlled input-change reproduction. It does not show that the
supplied stable corpus was affected.

Repair: capture FAA bytes/proteins and digest once per workflow and pass that
same snapshot to every source planner and every metadata writer. Standalone
must use the same capture contract. Collision checking does not bind bytes.
Either reject subsequent input mutation or consistently bind the captured
snapshot; never claim a digest for bytes a different adapter actually used.

Test real file mutation between phases and assert both failure containment and
manifest/context identity agreement. Reuse the captured parent-manifest hash,
too; do not introduce another load/hash race while repairing F1.

## F4 — Streaming and scalability remain unfinished

Code: stream_interproscan_tsv() retains hits at lines 400–423 and returns the
entire list; plan_interproscan() groups that list, then builds all per-query
proposals before deduplication. Each hit retains its pathway tuple.
MAX_AGGREGATED_CANDIDATES limits hit rows, not unique semantic candidates.
core/decisions.py lines 599–605 uses any(...) over the complete candidate list
for every producer-support candidate ID.

Parser-only duplicate-row probe, fixed to one query and the same semantic
annotations:

| Rows | Retained hit objects | Retained pathway tokens | Peak Python allocations |
|---:|---:|---:|---:|
| 500 | 500 | 50,000 | 3,684,366 bytes |
| 5,000 | 5,000 | 500,000 | 36,810,450 bytes |

The full four-check InterProScan-only process peaked at 902,148 KiB
(approximately 881 MiB). This includes merge/manifest work, not solely parsing.
C14 has nearly seven million pathway token occurrences. C14 four-source first
pass took 207.79 seconds in this environment versus 57.59 seconds for the
InterProScan-only acceptance run; these are observations, not portable timing
thresholds.

Repair:

- Aggregate during the one hashed binary pass. Retain bounded per-query,
  field/canonical-value facts, occurrence counts, compact location/row samples
  and order-defined support digests. Do not retain every pathway token.
- Count/digest pathway context while streaming, with bounded global and
  per-query counters. Bound diagnostic/QC lists and counter key namespaces.
- Apply an actual unique-candidate ceiling, separate from any hit/query ceiling.
- Replace repeated full-ledger lookup with a candidate-ID/status index.
- Add repeated-row/pathway tests demonstrating state grows with unique semantic
  facts and documented sample bounds, not every duplicate occurrence.
- Record wall time and peak memory on published data without fragile CI timing
  assertions. Preserve annotation sets and byte semantics through the refactor.

The current parser does hash the same TSV bytes it parses; preserve that useful
property. The defect is retained state and repeated lookup, not a separate TSV
hash/parse race.

## F5 — Explicit producer-version assertion is missing

Code: workflows/enrich.py line 276; plan_interproscan()/merge() default
arguments; standalone and unified CLI defaults.

Omitting interproscan_version is accepted and records 5.59-91.0 in provenance.
The audited plan required an explicit assertion because the TSV has no version
header that authenticates its producer. Compatible-shaped input from another
release can therefore acquire a false asserted version.

Require the version on unified CLI/API and standalone CLI/API when this source
is enabled. Keep unsupported versions rejected. Update current tests and the
real-data script, which currently endorse omission. An explicitly supplied
version remains a caller assertion, not executable authentication; document
that distinction.

## F6 — Invalid promoted-member evidence is silently lost

Code: sources/interproscan.py lines 609–640.

For an enabled Pfam row with malformed PF000010, validation returns None.
Planning silently omits its member note, does not record invalid_value, and
does not preserve the malformed member in context. Its integrated InterPro
proposal and inference still proceed. TIGRFAM follows the same pattern.

Use full-match validation at the appropriate row/promotion boundary, with a
stable contextual error before writes, or an explicitly documented and typed
invalid-value policy. Do not silently drop malformed enabled evidence.
Preserve context-only signatures without promoting unapproved grammars.

Also complete negative tests for control characters, ASCII accession digits,
exact line-ceiling boundaries and declared native score/status/date grammar.
Do not turn these checks into a universal E-value cutoff across heterogeneous
methods. Scores/descriptions remain context.

## F7 — Central exact-CDS identity invariant is missing

Code: core/decisions.py lines 337–365 and validate_candidate_ledger().

The adapter normally supplies the correct one-CDS UID. The central builder
still adds query-derived UIDs, accepts additional supplied UIDs, and falls
back to record/value resolution.

A source row claiming query Q1 but explicit UID Q2 produced an accepted
two-target supported-existing decision when both CDSs already carried the same
accession. A row with no explicit UID also passed. This is a reproduced
defensive-contract failure; no misrouting of published TSV queries was observed.

For InterProScan specifically, require exactly one registered CDS UID supplied
by the adapter, agreeing with its exact query identity. Reject missing, unknown,
wrong, multiple or non-CDS UIDs. Do not resolve this source through fallback.
Enforce it in the builder, runtime ledger validator and conditional schema.
Keep existing paired gene/CDS and legacy source contracts intact.

Retain the narrow producer-class support mapping, exact accepted candidate
links and final source projections. Avoid relaxing all sources to make an
InterProScan exception pass.

## F8 — Scientific release gates are not closed

Code: .github/workflows/ci.yml scientific-reruns job; tests/test_interproscan_acceptance.py;
tools/run_interproscan_real_checks.py.

Even a manual/tag scientific job runs only the two pre-existing follow-through
scripts. It never invokes the new InterProScan script or uploads its reports.
The new “acceptance” unit tests inspect constant lengths/counts and return
success if no generated report exists. They do not execute corpus acceptance.

The script pins a partial C14 golden. Restored SM has no corresponding
output/context/candidate golden, and it checks the expected six imported IDs
as a subset rather than asserting the complete bound output ledger. The
checked-in script does not exercise InterProScan no-op, fresh determinism,
standalone restored input or all-four-source acceptance.

Close this as an explicit release work item:

1. Wire InterProScan real checks into the manual/tag scientific job, retaining
   dataset-bytes dependency and every existing follow-through check.
2. Upload new reports and failure tracebacks. Required scientific evidence
   must not silently pass because a report is absent.
3. Put comprehensive versioned baselines in data/docs or fixtures, rather than
   relying only on production-script constants. Bind exact inputs/producers,
   version/layout/member policy, annotation sets, first/no-op hashes,
   candidate roles/statuses/reasons, exact origins and projections.
4. Add C14, restored SM, strict pristine SM, real BK71A rejection,
   mixed-legacy negatives, standalone restored runs, two fresh identical runs,
   output-as-input runs, and four-source C14/SM cases.
5. Reload final manifests with the existing evidence loader to verify output
   bindings; inspect exact origin sets, not only expected-membership checks.
6. Execute the final required CI jobs on the exact corrected revision and keep
   the reports. Local passes here do not substitute for wiring these gates.

No scientific finding in this report authorizes changing annotation policy or
refreshing a golden until the resulting semantic differences are reviewed.

## F9 — CLI and README mismatch

Reproduced: the README's four positional inputs exit 2; the parser requires
--bakta/--faa/--interproscan/--output. Standalone --version also exits 2.

Other discrepancies:

- ipr-only is documented as 11–13 columns; runtime rejects the native
  lookup-disabled 11-column layout and accepts 13 in this profile.
- The actual note prefix is PFAM:, not the README's Pfam:.
- Actual inference is protein motif:InterProScan:5.59-91.0, not the README's
  DESCRIPTION:similar to AA sequence:InterPro:MEMBER_DB:ACCESSION.
- Unified InterProScan-only uses its own context shape; the combined envelope
  is used only when eggNOG context is also present.
- A lineage manifest does not waive missing FAA or missing target translation;
  the restoration must actually establish eligible encoded proteins.

Choose one documented standalone convention, preserve existing working named
forms, and make every example executable. Add --version consistently and test
the shim and installed console entry outside the checkout. Align README and
CHANGELOG with the actual bounded-context and provenance semantics.

ENA's inference grammar permits an algorithm with version as evidence basis:
[official qualifier definition](https://www.ebi.ac.uk/ena/WebFeat/qualifiers/inference.html).
The discrepancy above is README/runtime agreement, not a claim that the
implemented protein-motif producer provenance is intrinsically invalid.

## F10 — Typed audit contracts and match support are incomplete

Only context_report_sha256 was added to the manifest schema. New source
metadata and interproscan_candidate entries remain arbitrary; no InterProScan
or combined-context schema exists.

I corrupted an InterProScan source entry to have a numeric normalized_value,
empty target UID list and invented final_status, and changed metadata's version
to a number. The checked-in JSON Schema still accepted the manifest. The
candidate decisions themselves were unchanged and remain under the existing
typed schema; do not confuse those guarantees with source projection typing.

The adapter also retains only the first supporting row after deduplication.
There is no occurrence count, bounded supporting match/location sample or
per-candidate evidence digest. support_digest hashes unique query/field/value
keys, not their contributing matches. Context signatures are silently cut to
50 labels, with no truncation count/digest, and per-query pathways have only a
total rather than database occurrence counters.

Repair in conjunction with F4:

- Conditional InterProScan metadata/source-entry schemas and versioned context
  schemas; preserve eggNOG-only and legacy source compatibility.
- Exact source projection-to-decision agreement, including support/insertion
  references and the one-CDS rule.
- Bounded, documented occurrence/location support summaries and hashes.
  Distinguish an annotation-set digest from a match-evidence digest.
- Report truncation explicitly and define digest ordering/delimiters.
- Keep raw TSV as the authoritative detailed evidence. Do not copy the whole
  export or millions of pathway labels into sidecars.
- Add schema-negative tests and evidence-aggregation tests that alter one
  supporting match without merely changing a hard-coded expected constant.

## Ordered Gemini execution plan

| Phase | Work | Exit evidence |
|---|---|---|
| 0 | Reproduce against the pinned main; preserve all inputs and old baselines | Finding disposition log: reproduced, corrected, qualified or deferred |
| 1 | F1, F2, F3 and F5: lineage, serialization, FAA capture, version assertion | Standalone restored success; both legacy bypasses rejected; hashes agree; omission rejected |
| 2 | F4 and F6: streaming aggregation, candidate index/limits, invalid-member diagnostics | Duplicate-state bound; valid corpus semantics preserved; malformed rows fail before writes |
| 3 | F7 and F10: exact UID guards, typed source/context schemas, bounded match support | Negative ledger/schema tests and final projection checks pass |
| 4 | F8 and F9: full acceptance, CLI/docs, installed-package tests and CI wiring | Complete local reports and required scientific CI evidence on corrected head |

Use an isolated branch/worktree and .test-output for generated artifacts.
Complete source changes and focused regression tests before broad corpus gates.
The expected patch target is a 0.4.x hardening release; choose/document the exact
version with the owner’s release convention. Preserve 0.4.0 observations and
older accepted baselines as history. Tool-version/context changes can alter
sidecar and manifest hashes; review those differences explicitly.

Do not push, tag, publish, fabricate recovery inputs or send messages to other
agents from this handoff. A reviewable implementation patch and evidence report
are the requested completion artifacts.

### Required closeout commands

Use a consistent selected Python environment and a clean checkout with
published LFS data materialized. The larger local MANIFEST is not a required
public-clone gate.


```bash
python -m pip install -e '.[dev,xlsx]'
git lfs pull
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
python -m pytest -q -p no:cacheprovider
python -m ruff check .
python -m ruff format --check .
python -m mypy src
git diff --check
python tools/run_followthrough_real_checks.py
python tools/run_followthrough_chain_checks.py
python tools/run_interproscan_real_checks.py
python -m build
```

The repaired InterProScan script must execute the expanded cases described
above; its current four passing cases are insufficient. For wheel verification,
install into a clean environment outside the source checkout, check the package
version/import origin, and execute the new console --help, --version and a
small standalone restored-input merge.

Final evidence report must identify the corrected commit, changed files,
finding dispositions, tests, exact input/output/context hashes, final candidate
projections/reasons, exact origins and parent bindings, resource observations,
reviewed baseline changes, installed CLI checks and scientific CI reports.

## Reproduction appendix

Run the following scripts from the pinned repository root after installing
.[dev,xlsx]. They use temporary synthetic fixtures and do not modify published
inputs. They are **observational audit probes**, not the final regression suite:
successful bug observations can produce exit code zero. Gemini should turn the
observations into tests asserting the repaired contract.

Save the first block as reproduce_interproscan_gaps.py and run it with the
selected Python environment. It covers F1/F2/F4/F5/F6/F7.

```python
from __future__ import annotations
import copy, hashlib, json, sys, tempfile, tracemalloc
from pathlib import Path
ROOT=Path.cwd().resolve()
sys.path.insert(0,str(ROOT/"tests"))
from test_interproscan_adapter import record_bytes, faa_bytes, DUMMY_QUERY, DUMMY_SEQ, DUMMY_MD5
from test_followthrough import restored_inputs
from test_merge_pipeline import record_bytes as short_record
from enrich_bakta_lib.sources.interproscan import merge, stream_interproscan_tsv, plan_interproscan
from enrich_bakta_lib.workflows.enrich import enrich
from enrich_bakta_lib.core.merge_engine import parse_genbank_bytes,cds_by_locus,feature_uid,has_legacy_restoration_comment
from enrich_bakta_lib.core.decisions import build_candidate_ledger,validate_candidate_ledger

OUT=ROOT/".test-output"/"sol-interproscan-audit"
OUT.mkdir(parents=True,exist_ok=True)
rows=[]
def check(name, action):
    try: result=action()
    except Exception as e: result={"exception_type":type(e).__name__,"error":str(e)}
    row={"probe":name,**result};rows.append(row);print(json.dumps(row),flush=True)
def inputs(folder):
    folder.mkdir()
    b=folder/"base.gbff"; f=folder/"base.faa";t=folder/"hits.tsv"
    b.write_bytes(record_bytes("REC1",DUMMY_QUERY)); f.write_bytes(faa_bytes())
    t.write_text((ROOT/"tests/fixtures/interproscan/valid_ipr_go_pathways.tsv").read_text())
    return b,f,t
with tempfile.TemporaryDirectory(prefix="interpro-review-probes-",dir=OUT) as temp:
    tmp=Path(temp)
    def standalone_import():
        folder=tmp/"import";folder.mkdir()
        base,faa,egg,parent=restored_inputs(folder)
        t=folder/"hits.tsv"
        md5=hashlib.md5(b"MM").hexdigest()
        t.write_text(f"T_0001\t{md5}\t2\tPfam\tPF00001\tDomain\t1\t2\t1e-5\tT\t25-08-2026\tIPR000001\tDomain\n")
        output=folder/"out.gbff";manifest=folder/"out.json"
        try:
            result=merge(base,faa,t,output,manifest_path=manifest,interproscan_version="5.59-91.0",interproscan_tsv_layout="ipr-only",
                         translation_evidence_manifest=parent,allow_imported_translations=True)
            return {"unexpected_success":True,**result}
        except Exception as e:
            return {"exception_type":type(e).__name__,"error":str(e),"gbff_written":output.exists(),"manifest_written":manifest.exists()}
    check("standalone_imported_manifest",standalone_import)
    def assumed_version():
        b,f,t=inputs(tmp/"version");m=b.parent/"manifest.json"
        enrich(bakta_path=b,faa_path=f,interproscan_path=t,output_path=b.parent/"out.gbff",manifest_path=m)
        return {"version_omitted_but_accepted":True,"recorded_version":json.loads(m.read_text())["metadata"]["interproscan_version"]}
    check("missing_producer_version",assumed_version)
    def malformed_member():
        b,f,t=inputs(tmp/"member")
        line=t.read_text().splitlines()[1].split("\t");line[3]="Pfam";line[4]="PF000010"
        t.write_text("\t".join(line)+"\n")
        p=plan_interproscan(parse_genbank_bytes(b.read_bytes()),f,t,version="5.59-91.0")
        return {"invalid_pfam_accepted":True,"fields":[r["field"] for r in p.evidence_rows],
                "context_entries":p.context_report["entries"]}
    check("malformed_promoted_member",malformed_member)
    def late_legacy():
        b,f,t=inputs(tmp/"legacy")
        data=b.read_bytes()
        marker=b"            Post-processed by normalize_baktfold.py to restore Bakta provenance\n"
        padding=(b"            "+b"padding "*10+b"\n")*1400
        data=data.replace(b"FEATURES             Location/Qualifiers\n",padding+marker+b"FEATURES             Location/Qualifiers\n")
        b.write_bytes(data)
        detected=has_legacy_restoration_comment(data)
        result=enrich(bakta_path=b,faa_path=f,interproscan_path=t,
                      interproscan_version="5.59-91.0",output_path=b.parent/"out.gbff")
        return {"legacy_marker_byte_offset":data.index(b"normalize_baktfold.py"),
                "legacy_detected":detected,"merge_accepted":result["self_check"]}
    check("legacy_after_100k",late_legacy)
    def target_contract():
        data=short_record("R1","Q1",db_xrefs=("InterPro:IPR000001",))+short_record("R2","Q2",db_xrefs=("InterPro:IPR000001",))
        base=parse_genbank_bytes(data); cds=cds_by_locus(base)
        evidence={"entry_type":"interproscan_candidate","query_id":"Q1",
                  "target_feature_uids":[feature_uid(cds["Q2"])],
                  "field":"InterPro","qualifier":"db_xref","raw_value":"IPR000001",
                  "normalized_value":"InterPro:IPR000001","status":"supported_existing",
                  "emitted_qualifiers":"","candidate_role":"functional_proposal","evidence_class":"InterPro",
                  "reason_code":"value_already_present","source":"InterProScan"}
        dec,_=build_candidate_ledger(base,[],[],[evidence],source_hashes={"InterProScan":"a"*64})
        validate_candidate_ledger(dec,[],base=base,output=base)
        evidence2=copy.deepcopy(evidence);evidence2.pop("target_feature_uids")
        dec2,_=build_candidate_ledger(base,[],[],[evidence2],source_hashes={"InterProScan":"a"*64})
        validate_candidate_ledger(dec2,[],base=base,output=base)
        return {"wrong_explicit_uid_accepted_target_count":len(dec[0]["target_feature_uids"]),
                "missing_explicit_uid_accepted":len(dec2)==1}
    check("exact_uid_contract",target_contract)
    def retained_pathways():
        folder=tmp/"memory";folder.mkdir()
        results=[]
        for count in (500,5000):
            t=folder/f"repeated-{count}.tsv"
            fields=(ROOT/"tests/fixtures/interproscan/valid_ipr_go_pathways.tsv").read_text().splitlines()[1].split("\t")
            fields[-1]="|".join(f"Reactome:R-HSA-{n}" for n in range(100))
            line="\t".join(fields)+"\n";t.write_text(line*count)
            tracemalloc.start()
            sha,hits,diag=stream_interproscan_tsv(t,"ipr-go-pathways")
            current,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
            results.append({"rows":count,"unique_queries":len(diag["query_metadata"]),
                            "retained_hit_objects":len(hits),"retained_pathway_tokens":sum(len(h.pathways) for h in hits),
                            "peak_python_bytes":peak})
            del hits,diag
        return {"measurements":results}
    check("repeated_rows_memory",retained_pathways)
(OUT/"reproduction-results.json").write_text(json.dumps(rows,indent=2)+"\n")
```

The following extra probe covers FAA snapshot disagreement, source schema
typing and CLI drift (F3/F9/F10).

```python
from __future__ import annotations
import copy,json,sys,tempfile,hashlib,subprocess
from pathlib import Path
from unittest.mock import patch
ROOT=Path.cwd().resolve()
sys.path.insert(0,str(ROOT/"tests"))
from test_interproscan_adapter import record_bytes,faa_bytes,DUMMY_QUERY
import enrich_bakta_lib.sources.interproscan as ips
from enrich_bakta_lib.workflows.enrich import enrich
from jsonschema import Draft202012Validator
OUT=ROOT/".test-output"/"sol-interproscan-audit"
OUT.mkdir(parents=True,exist_ok=True)
results=[]
with tempfile.TemporaryDirectory(dir=OUT,prefix="extra-probes-") as temp:
    folder=Path(temp);b=folder/"base.gbff";f=folder/"base.faa";t=folder/"hits.tsv"
    b.write_bytes(record_bytes("REC1",DUMMY_QUERY));f.write_bytes(faa_bytes(DUMMY_QUERY,"MM"))
    t.write_text((ROOT/"tests/fixtures/interproscan/valid_ipr_go_pathways.tsv").read_text())
    original=ips.read_input_bytes
    def change_before_second_faa_read(path,label):
        if label=="FAA":f.write_bytes(faa_bytes())
        return original(path,label)
    m=folder/"manifest.json";c=folder/"context.json"
    with patch.object(ips,"read_input_bytes",side_effect=change_before_second_faa_read):
        res=enrich(bakta_path=b,faa_path=f,interproscan_path=t,interproscan_version="5.59-91.0",
                   output_path=folder/"out.gbff",manifest_path=m,context_report_path=c)
    man=json.loads(m.read_text());ctx=json.loads(c.read_text())
    results.append({"probe":"faa_snapshot_consistency","accepted":res["self_check"],
                    "manifest_faa_sha256":man["metadata"]["faa_sha256"],
                    "context_faa_sha256":ctx["metadata"]["faa_sha256"],
                    "actual_faa_sha256":hashlib.sha256(f.read_bytes()).hexdigest()})
    bad=copy.deepcopy(man)
    row=next(r for r in bad["entries"] if r["entry_type"]=="interproscan_candidate")
    row["normalized_value"]=123;row["target_feature_uids"]=[];row["final_status"]="invented_status"
    bad["metadata"]["interproscan_version"]=123
    schema=json.loads((ROOT/"schemas/merge-manifest.v2.schema.json").read_text())
    errors=list(Draft202012Validator(schema).iter_errors(bad))
    results.append({"probe":"interproscan_schema_contract","corrupted_source_entry_and_metadata_accepted":not errors})
    positional=subprocess.run([sys.executable,str(ROOT/"merge_interproscan_bakta.py"),
                    str(b),str(f),str(t),str(folder/"cli-out.gbff")],
                    capture_output=True,text=True)
    version=subprocess.run([sys.executable,str(ROOT/"merge_interproscan_bakta.py"),"--version"],
                    capture_output=True,text=True)
    results.append({"probe":"documented_cli","positional_exit_code":positional.returncode,
                    "positional_error":positional.stderr.splitlines()[-1],
                    "version_exit_code":version.returncode})
for row in results:print(json.dumps(row),flush=True)
(OUT/"extra-reproduction-results.json").write_text(json.dumps(results,indent=2)+"\n")
```

Mixed current-marker/legacy-record proof for F2:

```python
from pathlib import Path
import sys, tempfile
ROOT=Path.cwd().resolve()
sys.path.insert(0,str(ROOT/"tests"))
from test_followthrough import inputs
from test_interproscan_adapter import record_bytes,faa_bytes
from enrich_bakta_lib.workflows.restore_translations import restore
from enrich_bakta_lib.workflows.enrich import enrich
from enrich_bakta_lib.core.merge_engine import load_translation_evidence,parse_genbank_bytes
with tempfile.TemporaryDirectory() as temp:
    folder=Path(temp)
    base,faa,egg=inputs(folder,protein=b"MM")
    base.write_bytes(base.read_bytes()+record_bytes(
        "LEGACY","OJMJKD_006675",
        comment_header="COMMENT     normalize_baktfold.py: restore Bakta provenance"))
    faa.write_bytes(faa.read_bytes()+faa_bytes())
    restored=folder/"restored.gbff";parent=folder/"parent.json"
    restore(base,faa,egg,restored,manifest_path=parent,translation_policy="import-faa")
    ledger=load_translation_evidence(parent,parse_genbank_bytes(restored.read_bytes()))
    assert sorted(ledger)==["T_0001"]
    # At the audited revision this succeeds despite the legacy query lacking evidence.
    result=enrich(
        bakta_path=restored,faa_path=faa,
        interproscan_path=ROOT/"tests/fixtures/interproscan/valid_ipr_go_pathways.tsv",
        interproscan_version="5.59-91.0",output_path=folder/"out.gbff",
        manifest_path=folder/"out.json",translation_evidence_manifest=parent,
        allow_imported_translations=True)
    print({"legacy_query_in_ledger":"OJMJKD_006675" in ledger,
           "functional_merge_accepted":result["self_check"]})
```

## Gemini execution prompt

Independently reproduce this audit against
a34a2076da453961bbe29eec6664df7a2f7eb7cb, or first review any newer-main delta.
Resolve F1–F10 in the ordered phases above, preserving the original inputs and
accepted scientific policy. Keep exact protein/CDS identity, conservative
qualifier promotion, complete translation-origin carry-forward and byte
preservation. Do not weaken lineage checks to make historical BK71A succeed.

Complete focused regressions, source/context schemas, the documented CLI,
installed-package verification, comprehensive C14/SM scientific baselines and
manual/tag CI wiring. Review all changed hashes against actual semantic
differences; never refresh goldens simply to pass. Deliver a reviewable patch,
finding-resolution log and complete evidence report. Do not push, tag or
publish without the owner's instruction.
