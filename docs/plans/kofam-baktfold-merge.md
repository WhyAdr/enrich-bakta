# Plan: merge KofamScan/KOALA and Baktfold annotations onto Bakta GBFF

Status: implemented (2026-08-18)

## Goal

Create a reproducible, provenance-preserving annotation merge pipeline whose
authoritative base remains the original Bakta `.gbff`. The pipeline should:

1. retain the existing Baktfold graft behavior for genuine Baktfold additions;
2. add KofamScan KO assignments generated from the matching Bakta `.faa`;
3. preserve Bakta's original feature blocks, sequence, metadata, and
   qualifiers byte-for-byte after removing explicitly recorded insertions; and
4. make it impossible for a mismatched annotation, FASTA, or KofamScan file to
   silently enrich the wrong feature.

The output is annotation-supported evidence. A KO assignment, an EC string in
a KO definition, or a Baktfold structural hit is not proof of enzyme activity,
pathway completeness, expression, or phenotype.

## Local audit baseline

The current `graft_baktfold_additions.py` compiles and runs successfully on the
two real pairs in this workspace:

| Pair | Bakta feature entries | Baktfold feature entries | Matched keyed features | New genes | New ECs | New structural xrefs |
|---|---:|---:|---:|---:|---:|---:|
| C14 | 9,644 | 9,644 | 9,586 | 800 | 178 | 58 |
| SM | 9,528 | 9,528 | 9,464 | 780 | 168 | 98 |

For both pairs, the ordered `(record, feature type, locus_tag)` sequence and
feature location headers are identical, and the script's self-verification
passes. Baktfold does rewrite many existing qualifier values and orders, which
supports grafting onto pristine Bakta rather than restoring Bakta from the
Baktfold file.

The KofamScan files are hit-only KOALA tables with the standard header and
leading `*` hit marker:

| Input | Hit rows | Unique query IDs | CDSs with one or more hits | Existing Bakta KEGG pairs | New gene–KO pairs |
|---|---:|---:|---:|---:|---:|
| `C14-KofamKOALA.txt` | 3,397 | 3,181 | 3,181 | 1,481 | 1,916 |
| `SM-KofamKOALA.txt` | 3,333 | 3,120 | 3,120 | 1,433 | 1,900 |

Every observed Kofam query ID matches a CDS locus tag in both the Bakta FAA and
GBFF. Multiple KO hits occur for individual CDSs, so a one-value-per-gene map
would be incorrect. Existing Bakta KEGG annotations are normally
`/note="KEGG:Kxxxxx"`; the Baktfold output represents the same family of
identifiers as `/db_xref="KEGG:Kxxxxx"`.

## Findings and patches for `graft_baktfold_additions.py`

### P0: strengthen structural parity before any graft

The current `structural_parity_check()` compares global feature-type counts and
sets of locus tags. That is useful, but it can miss a swapped feature, a
duplicate `(locus_tag, feature_type)`, a record-local mismatch, or a mismatch
among features without locus tags. The current index also silently overwrites
duplicate Baktfold keys.

Replace the permissive checks with a record-aware validation report that
compares, in order:

- record IDs/accessions, record count, lengths, topology, and source sequence
  checksums;
- the ordered feature sequence and per-record feature-type counts;
- exact feature-key multiplicity, with duplicate keys treated as an error;
- `(record, feature type, locus_tag)` for keyed features;
- raw location/header text, or a normalized Biopython location equivalent; and
- CDS-to-CDS translation identity only as a diagnostic, not as a reason to copy
  Baktfold sequence data.

Do not let `--force` turn a failed identity check into a normal enrichment run.
If retained for exceptional forensic use, it should require an explicit
override, print a machine-readable mismatch report, and mark the output as
unsafe. The default path should fail before writing an output file.

### P0: make the insertion and preservation proof stronger

The current self-check proves that Bakta feature blocks occur as prefixes of
the output blocks, but arbitrary appended lines would still pass. It also uses
text decoding with `errors="replace"`, so its “byte-for-byte” claim is only a
line-level claim for UTF-8/LF inputs.

Implement a raw-byte preservation audit:

- read and splice bytes without normalizing newline style or invalid bytes;
- refuse an output path equal to either input path;
- write to a temporary sibling file and atomically replace the requested new
  output only after validation succeeds;
- record every inserted line/block with record, feature key, qualifier, and
  source value; and
- reconstruct the original bytes by removing only the recorded insertions and
  the exact comment addendum, then compare SHA-256 and byte content.

Add a second allowlist check that every inserted qualifier is one of the
computed Baktfold additions or the requested provenance text. This separates
“the original survived” from “only authorized content was added.” Make the
merge idempotent: a second run must not add duplicate provenance lines or
duplicate qualifiers.

### P1: tighten qualifier semantics

The existing qualifier scanner is intentionally minimal, but its boundaries
and value handling should be tested against Biopython and tricky fixtures.
Keep a raw scanner only for locating verbatim feature blocks; use the canonical
GenBank/Biopython parser for semantic validation where available.

Specific changes:

- require structural xrefs to have an exact allowed prefix followed by `:`;
  the current `startswith()` accepts values such as `pdbx:` and `cathode:`;
- treat blank `/gene` values as absent, but never overwrite a non-empty Bakta
  gene;
- preserve repeated qualifiers; if Baktfold supplies multiple distinct gene
  values where Bakta has none, report the ambiguity rather than silently taking
  only the first value;
- recognize ECs from both Baktfold `/EC_number` and
  `/db_xref="EC:..."`, normalize only for comparison, and emit the qualifier
  form selected by the explicit compatibility policy;
- de-duplicate exact values while retaining wildcard EC strings as annotations
  rather than upgrading their biological meaning; and
- quote/escape and wrap newly emitted qualifier values using a single GenBank
  formatter.

### P1: correct provenance and reproducibility

`/inference="protein structure similarity:Baktfold:..."` is currently added
for every type of addition. Separate provenance by evidence class:

- structural xrefs: Baktfold Foldseek/structure evidence;
- gene symbols: Baktfold-derived gene-name evidence; and
- EC additions: Baktfold-derived functional annotation evidence.

Avoid claiming that an EC or gene symbol is a protein-structure result unless
that is confirmed by the Baktfold output contract. Add the Baktfold input hash,
detected version, and counts to a deterministic manifest or comment addendum.
Make the timestamp optional or explicitly supplied so identical inputs can
produce identical output.

### P1: tests for the existing graft

Add small synthetic GenBank fixtures and real-file smoke tests covering:

- one record and multiple records;
- duplicate keys, reordered features, shifted locations, missing records, and
  missing locus tags;
- repeated and blank qualifiers, wrapped strings, and escaped quotes;
- exact xref-prefix filtering and EC normalization;
- an output with no additions;
- idempotent reruns;
- preservation of LF/CRLF and final-newline style; and
- failure before output creation when parity or allowlist validation fails.

Keep the C14 and SM runs as integration checks, with their current counts as
regression baselines unless an intentional policy change updates them.

## New KofamScan merger

Add a focused script, tentatively named `merge_kofamscan_bakta.py`, using the
same byte-preserving graft engine. A first version should support one explicit
input set:

```text
python merge_kofamscan_bakta.py \
  BAKTA.gbff BAKTA.faa KofamKOALA.txt OUTPUT.gbff \
  [--manifest OUTPUT.tsv] [--no-comment-note] [--kofamscan-version VERSION]
```

Batch discovery can follow after the single-input behavior is stable; the
current filename pairs are `C14-NMZ_bakta/C14-NMZ.*` plus
`C14-KofamKOALA.txt`, and the equivalent SM files.

### Kofam parser and validation

Parse the header-driven table rather than relying on fixed column widths:

- skip comment/header lines;
- recognize the leading `*` as a hit marker;
- parse query ID, KO, threshold, score, E-value, and the remaining definition;
- validate KO syntax (`K` plus five digits), numeric threshold/score, and
  scientific-notation or zero E-values; and
- retain every hit row, including multiple KOs for one query.

Before writing, require:

- unique FAA record IDs;
- every Kofam query ID to occur exactly in the FAA and in one CDS `locus_tag` in
  the GBFF;
- no duplicate GBFF CDS keys; and
- a normalized protein-sequence checksum match between the Bakta FAA and the
  corresponding GBFF CDS translation, with mismatches reported by locus tag.

FAA records without a Kofam hit are expected. Kofam hit IDs missing from either
the FAA or GBFF are not silently dropped.

### Emission policy

Use the Bakta GBFF as the only base. Modify CDS feature blocks only; do not
mirror KO calls onto paired `gene` features and do not replace Bakta
`/product`, `/gene`, `/translation`, `/protein_id`, `/EC_number`, or existing
inference lines.

For each Kofam hit, append deterministic, de-duplicated evidence:

1. `/db_xref="KEGG:Kxxxxx"` when the KO is not already represented in either
   Bakta `/db_xref` or `/note`.
2. A compact provenance `/note` containing the KO, threshold, score, and
   E-value, so an existing Bakta KO can still be distinguished as independently
   supported by KofamScan.
3. One `/inference="profile:KofamScan[:VERSION]"` per CDS with at least one
   Kofam hit, unless that exact inference already exists.

The exact `KEGG:` xref spelling should be tested against downstream consumers;
it matches the representation already present in the Baktfold output and the
identifier syntax used by Bakta's KEGG notes. Do not parse `[EC:...]` from the
KO definition into `/EC_number`: that is KO-level descriptive text, not an
independent EC assignment for this protein. Do not overwrite `/product` with
the Kofam definition; retain the full definition in the optional manifest.

The optional TSV/JSON manifest should include input hashes, query ID, locus tag,
record, KO, threshold, score, E-value, definition, whether the KO was already
present, and the exact qualifiers emitted. It is the authoritative place for
the detailed Kofam evidence while the GBFF remains standards-compatible and
readable.

### Combined Baktfold + Kofam operation

The preferred final architecture is one shared merge engine that validates all
sources first and applies both addition sets to pristine Bakta in one pass.
This avoids repeated comment rewrites and makes one preservation proof cover
the final artifact.

An acceptable intermediate implementation is:

1. pristine Bakta + Baktfold -> Baktfold-grafted Bakta;
2. that output + Kofam/FAA -> final KOALA-enriched GBFF; and
3. a final preservation audit against pristine Bakta with both insertion
   manifests.

The two-source implementation must retain independent provenance and must not
let a Baktfold-added qualifier suppress a Kofam evidence record, or vice versa.

## Acceptance gates

### Existing Baktfold graft

- C14 and SM pass strict parity, semantic validation, and byte-preservation
  checks.
- Existing additions remain present with no Bakta qualifier overwritten.
- False-prefix xrefs, duplicate keys, shifted locations, and malformed
  qualifiers fail safely.
- A second identical run is byte-identical to the first output.

### Kofam merger

- C14: 3,397 parsed hits, 3,181 hit CDSs, 1,916 new gene–KO pairs, 1,481
  already represented, and zero missing FAA/GBFF query IDs.
- SM: 3,333 parsed hits, 3,120 hit CDSs, 1,900 new gene–KO pairs, 1,433
  already represented, and zero missing FAA/GBFF query IDs.
- Multiple KOs on one CDS are all retained and no Kofam hit is attached to a
  different locus tag.
- Bakta translations, products, gene names, ECs, protein IDs, comments, and
  existing inference lines are unchanged except for explicitly recorded
  additions.
- Output parses as valid GenBank/Biopython input and passes the repository's
  available GenBank validation/CLI checks.
- The manifest and output are deterministic, idempotent, and contain input
  hashes and source provenance.

## Suggested implementation order

1. Extract or define the shared raw-splice, semantic-validation, qualifier-
   formatting, insertion-manifest, and byte-audit components.
2. Patch and regression-test `graft_baktfold_additions.py` against C14, SM, and
   synthetic failure fixtures.
3. Implement and test `merge_kofamscan_bakta.py` against both KOALA files,
   starting with structured KO xrefs plus compact per-hit provenance.
4. Add the one-pass combined command only after both independent mergers are
   idempotent and fully validated.
5. Keep large/private GBFF, FAA, and Kofam inputs out of any future commit;
   commit only scripts, small fixtures, tests, and this plan when explicitly
   requested.

## Implementation result

The plan is implemented in `merge_engine.py`,
`graft_baktfold_additions.py`, and `merge_kofamscan_bakta.py` with synthetic
coverage in `tests/test_merge_pipeline.py`. The Kofam command supports both the
focused four-input interface and the preferred one-pass combined operation via
`--baktfold`.

One audited input nuance is now explicit in the validation policy: Bakta FAA
files can contain translated pseudogene candidates for which Bakta omits a
GBFF `/translation`. Non-hit FAA entries are allowed, while every Kofam hit is
still required to have an exact normalized FAA-to-GBFF translation match.
