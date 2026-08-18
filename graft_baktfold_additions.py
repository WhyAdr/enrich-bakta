#!/usr/bin/env python3
"""
graft_baktfold_additions.py
============================

Start from the ORIGINAL, PRISTINE Bakta .gbff file and graft on ONLY the
genuine new information contributed by Baktfold, leaving every single byte
of Bakta's own output physically untouched. Nothing is reconstructed or
reverse-engineered -- every line of the output is either (a) copied
verbatim from the Bakta input, or (b) a newly appended qualifier line
whose provenance is traceable to a specific Baktfold field.

Rationale (why graft-onto-pristine instead of restore-from-corrupted):
  Baktfold v0.1.0's GenBank writer has a known, reproducible bug: it turns
  every /transl_table=11 into a non-standard /translation_table=False,
  overwrites the Bakta/DB version stamps in COMMENT, and rewrites the
  LOCUS date to a hardcoded placeholder (01-JAN-1980). A "restore" script
  that starts from that file has to hand-reconstruct exactly what Bakta
  would have written -- including internal accounting quirks like Bakta's
  COMMENT "CDSs" count, which is not a pure function of feature counts and
  has no way to be independently verified once the original is discarded.
  Grafting starts from ground truth instead and never has to guess.

Genuine Baktfold additions grafted onto matching Bakta features:
  1. /gene=              on features where Bakta assigned none
                          (rRNA, ncRNA, tRNA, CDS, and their paired `gene`
                          features -- see NOTE on Bakta's gene/CDS mirroring
                          below)
  2. /EC_number=          new EC numbers Baktfold assigned that Bakta did
                          not (written as /EC_number=, matching Bakta's own
                          qualifier type -- NOT /db_xref="EC:...", which
                          would be silently dropped by any /EC_number-
                          priority consumer, e.g. genbank_parser.py's
                          extract_xrefs())
  3. /db_xref=            new structural cross-references with prefixes
                          {afdb_v6, cath, pdb} -- AlphaFold DB, CATH, and
                          PDB hits from Baktfold's Foldseek search, which
                          have no Bakta equivalent
  4. /inference=          (optional, on by default) one additional
                          provenance line on any feature that received an
                          addition above, documenting that it came from
                          Baktfold's structural-similarity pipeline rather
                          than Bakta's own evidence

NOTE on gene/CDS mirroring: Bakta emits a `gene` feature and a specific
feature (CDS/tRNA/rRNA/ncRNA/tmRNA) as a pair sharing one locus_tag; both
carry the same /gene= symbol when one is assigned. Baktfold preserves this
convention (confirmed: its added-gene counts split exactly as expected
across `gene` feature-type entries and their paired specific-type
entries). Because this script matches on (locus_tag, feature_type)
independently for every feature, both halves of the pair are grafted
correctly without any special-casing.

WHAT THIS SCRIPT NEVER TOUCHES:
  - Any line inside a Bakta feature block that isn't a newly appended
    qualifier (verified automatically -- see "Self-verification" below)
  - /transl_table=11, LOCUS date, COMMENT provenance (Software/Database/
    DOI/URL/Annotation Date), /protein_id= namespace, /inference=
    strings Bakta already wrote, or the alphabetical qualifier order
  - Any feature/locus_tag Baktfold doesn't have a match for (passed
    through unchanged)

WHAT IS APPENDED OUTSIDE THE FEATURE TABLE:
  A short addendum inserted at the end of each record's COMMENT block
  (before FEATURES), documenting that Baktfold-derived annotations were
  merged in, with the Baktfold version and merge timestamp. Every existing
  byte of Bakta's COMMENT block (Software/Database/DOI/URL/Annotation
  Date/summary counts) is preserved verbatim above it. Disable with
  --no-comment-note.

Self-verification: after writing the output, the script re-reads it,
removes exactly the lines it inserted, and asserts the result is
byte-for-byte identical to the original Bakta input. This is not a
sanity heuristic -- it is a direct proof that nothing else changed.

USAGE:
  python graft_baktfold_additions.py BAKTA.gbff BAKTFOLD.gbff OUTPUT.gbff
      [--no-comment-note] [--no-inference-provenance] [--force]

EXAMPLE:
  python graft_baktfold_additions.py \\
      BK71A-bakta.gbff BK71A-baktfold.gbff BK71A-merged.gbff
"""
from __future__ import annotations

import argparse
import bisect
import collections
import re
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Minimal embedded parser
#
# Adapted from the genbank-feature-parser skill's canonical genbank_parser.py
# (regexes validated against independent grep/awk/md5sum ground truth during
# the BK71A Bakta-vs-Baktfold audit). Trimmed to what this task needs: we
# don't care about genomic coordinates or join()/order() locations here, only
# (a) precise feature-block line boundaries for verbatim copying, and
# (b) qualifier values for matching and diffing. This is deliberately
# self-contained (stdlib only) so it runs anywhere, independent of the skill.
# ---------------------------------------------------------------------------

_locus_re = re.compile(r'^LOCUS\s+(\S+)')
_feature_re = re.compile(r'^     (\w+)\s+\S')          # any feature header line
_qualifier_re = re.compile(r'^\s{19,22}/(\w+)(?:="?(.*?)"?)?\s*$')
_continuation_re = re.compile(r'^\s{19,22}(?!/)\S')


def parse_feature_blocks(lines: list[str]) -> list[dict]:
    """Parse feature blocks with exact [start_line, end_line) boundaries.

    end_line is exclusive and always equals the index of whatever comes
    next (next feature header, ORIGIN, or //), so lines[start:end] is the
    complete, exact block for that feature -- safe to copy verbatim.
    """
    features: list[dict] = []
    current: dict | None = None
    current_contig = None
    in_origin = False

    def flush(end_line: int) -> None:
        nonlocal current
        if current is not None:
            current['end_line'] = end_line
            features.append(current)
            current = None

    for i, line in enumerate(lines):
        if line.startswith('ORIGIN'):
            in_origin = True
            flush(i)
            continue
        if line.startswith('//'):
            in_origin = False
            flush(i)
            continue
        if in_origin:
            continue

        lm = _locus_re.match(line)
        if lm:
            flush(i)
            current_contig = lm.group(1)
            continue

        fm = _feature_re.match(line)
        if fm:
            flush(i)
            current = {
                'type': fm.group(1),
                'contig': current_contig,
                'start_line': i,
                'end_line': None,
                'qualifiers': collections.defaultdict(list),
            }
            continue

        qm = _qualifier_re.match(line)
        if qm and current is not None:
            key, val = qm.group(1), qm.group(2) or ''
            current['qualifiers'][key].append(val)
            continue

        if current is not None and _continuation_re.match(line) and current['qualifiers']:
            last_key = list(current['qualifiers'].keys())[-1]
            current['qualifiers'][last_key][-1] += line.strip()

    flush(len(lines))
    return features


def get1(feature: dict, key: str) -> str | None:
    vals = feature['qualifiers'].get(key, [])
    return vals[0] if vals else None


def get_all(feature: dict, key: str) -> list[str]:
    return list(feature['qualifiers'].get(key, []))


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

QUALIFIER_INDENT = ' ' * 21  # GenBank standard qualifier indent
NEW_XREF_PREFIXES = ('afdb_v6', 'cath', 'pdb')  # Baktfold structural DB hits


# ---------------------------------------------------------------------------
# Structural parity pre-check
# ---------------------------------------------------------------------------

def structural_parity_check(bakta_feats: list[dict], baktfold_feats: list[dict]) -> list[str]:
    """Return a list of problems (empty = safe to graft)."""
    problems = []

    def type_counts(feats):
        c = collections.Counter(f['type'] for f in feats)
        return c

    bc, fc = type_counts(bakta_feats), type_counts(baktfold_feats)
    if bc != fc:
        problems.append(f"feature-type counts differ: bakta={dict(bc)} baktfold={dict(fc)}")

    def locus_tags(feats):
        return {get1(f, 'locus_tag') for f in feats if get1(f, 'locus_tag')}

    bl, fl = locus_tags(bakta_feats), locus_tags(baktfold_feats)
    if bl != fl:
        missing = bl - fl
        extra = fl - bl
        problems.append(
            f"locus_tag sets differ: {len(missing)} in bakta only, "
            f"{len(extra)} in baktfold only"
        )
    return problems


# ---------------------------------------------------------------------------
# Compute additions for one Bakta feature, given its Baktfold counterpart
# ---------------------------------------------------------------------------

def compute_additions(bakta_f: dict, baktfold_f: dict) -> tuple[list[str], dict]:
    """Return (new_qualifier_strings, tally) to append to this Bakta feature.

    Order: new /db_xref= (sorted) -> new /EC_number= (sorted) -> new /gene=
    This mirrors Bakta's own observed convention of placing /gene= as the
    last qualifier in a block (confirmed against Bakta's own CDS features
    that already carry both /EC_number= and /gene=).
    """
    tally = {'gene': 0, 'ec': 0, 'xref': 0}
    additions: list[str] = []

    # 1. New structural db_xrefs (afdb_v6 / cath / pdb) -- Bakta has none of these.
    bakta_xrefs = set(get_all(bakta_f, 'db_xref'))
    new_xrefs = sorted({
        x for x in get_all(baktfold_f, 'db_xref')
        if x.startswith(NEW_XREF_PREFIXES) and x not in bakta_xrefs
    })
    for x in new_xrefs:
        additions.append(f'/db_xref="{x}"')
        tally['xref'] += 1

    # 2. New EC numbers. Baktfold stores ALL ECs (old + new) as db_xref="EC:...";
    #    Bakta stores them as /EC_number=. Emit as /EC_number= to match Bakta's
    #    qualifier type and stay compatible with /EC_number-priority consumers.
    bakta_ecs = set(get_all(bakta_f, 'EC_number'))
    baktfold_ecs = {x[3:] for x in get_all(baktfold_f, 'db_xref') if x.startswith('EC:')}
    new_ecs = sorted(baktfold_ecs - bakta_ecs)
    for ec in new_ecs:
        additions.append(f'/EC_number="{ec}"')
        tally['ec'] += 1

    # 3. New /gene= -- only if Bakta assigned none. Bakta's own value is always
    #    kept as source of truth if present.
    if not get_all(bakta_f, 'gene'):
        baktfold_genes = get_all(baktfold_f, 'gene')
        if baktfold_genes:
            additions.append(f'/gene="{baktfold_genes[0]}"')
            tally['gene'] += 1

    return additions, tally


# ---------------------------------------------------------------------------
# COMMENT-block addendum (appended, never overwrites Bakta's own text)
# ---------------------------------------------------------------------------

def build_addendum(baktfold_version: str) -> list[str]:
    ts = datetime.now().strftime('%Y-%m-%d %H:%M')
    return [
        '            ',
        '            ##Post-processing:##',
        f'            Merged with Baktfold v{baktfold_version} annotations',
        '            (https://github.com/gbouras13/baktfold) using',
        f'            graft_baktfold_additions.py on {ts}.',
        '            Bakta provenance above is unmodified; only new /gene=,',
        '            /EC_number=, and /db_xref=(afdb_v6|cath|pdb) qualifiers',
        '            were appended to existing features.',
    ]


def detect_baktfold_version(baktfold_lines: list[str]) -> str:
    for line in baktfold_lines[:60]:
        m = re.search(r'Software:\s*v?([\d.]+)', line)
        if m:
            return m.group(1)
    return 'unknown'


def insert_comment_addenda(out_lines: list[str], baktfold_version: str) -> tuple[list[str], set[int]]:
    """Insert the addendum before each record's FEATURES header.

    Returns (new_lines, inserted_indices) where inserted_indices are the
    indices (in the returned list) of every line this function added, for
    use in self-verification.
    """
    addendum = build_addendum(baktfold_version)
    result: list[str] = []
    inserted: set[int] = set()
    in_comment = False

    for line in out_lines:
        if line.startswith('COMMENT'):
            in_comment = True
            result.append(line)
            continue
        if in_comment and line.startswith('FEATURES'):
            for a in addendum:
                inserted.add(len(result))
                result.append(a)
            in_comment = False
            result.append(line)
            continue
        result.append(line)

    return result, inserted


# ---------------------------------------------------------------------------
# Main graft
# ---------------------------------------------------------------------------

def graft(
    bakta_path: Path,
    baktfold_path: Path,
    output_path: Path,
    add_comment_note: bool = True,
    add_inference_provenance: bool = True,
    force: bool = False,
) -> dict:
    bakta_lines = bakta_path.read_text(encoding='utf-8', errors='replace').split('\n')
    baktfold_lines = baktfold_path.read_text(encoding='utf-8', errors='replace').split('\n')

    bakta_feats = parse_feature_blocks(bakta_lines)
    baktfold_feats = parse_feature_blocks(baktfold_lines)

    problems = structural_parity_check(bakta_feats, baktfold_feats)
    if problems:
        msg = "Structural parity check FAILED:\n  " + "\n  ".join(problems)
        if not force:
            sys.exit(
                msg + "\n\nThese two files don't look like they came from the same "
                "annotation run -- grafting would silently mismatch features.\n"
                "Re-run with --force only if you've confirmed this is expected."
            )
        else:
            print("WARNING: " + msg + "\n(continuing because --force was given)", file=sys.stderr)

    # Index Baktfold features by (locus_tag, feature_type)
    baktfold_index: dict[tuple[str, str], dict] = {}
    for f in baktfold_feats:
        lt = get1(f, 'locus_tag')
        if lt:
            baktfold_index[(lt, f['type'])] = f

    stats = {
        'features_processed': 0,
        'features_matched': 0,
        'gene_added': 0,
        'ec_added': 0,
        'xref_added': 0,
        'features_with_any_addition': 0,
        'gene_by_type': collections.Counter(),
        'ec_by_type': collections.Counter(),
        'xref_by_type': collections.Counter(),
    }

    baktfold_version = detect_baktfold_version(baktfold_lines)

    out: list[str] = []
    inserted_indices: set[int] = set()
    cursor = 0

    for f in bakta_feats:
        # Copy any inter-feature lines (headers, ORIGIN/sequence, //, next LOCUS...)
        out.extend(bakta_lines[cursor:f['start_line']])

        # Copy this feature's block verbatim
        out.extend(bakta_lines[f['start_line']:f['end_line']])
        stats['features_processed'] += 1

        lt = get1(f, 'locus_tag')
        additions: list[str] = []
        tally = {'gene': 0, 'ec': 0, 'xref': 0}
        if lt:
            key = (lt, f['type'])
            bf = baktfold_index.get(key)
            if bf is not None:
                stats['features_matched'] += 1
                additions, tally = compute_additions(f, bf)

        if additions:
            stats['features_with_any_addition'] += 1
            for a in additions:
                out_index = len(out)
                out.append(QUALIFIER_INDENT + a)
                inserted_indices.add(out_index)
            if add_inference_provenance:
                prov = f'/inference="protein structure similarity:Baktfold:{baktfold_version}"'
                out_index = len(out)
                out.append(QUALIFIER_INDENT + prov)
                inserted_indices.add(out_index)

            stats['gene_added'] += tally['gene']
            stats['ec_added'] += tally['ec']
            stats['xref_added'] += tally['xref']
            if tally['gene']:
                stats['gene_by_type'][f['type']] += 1
            if tally['ec']:
                stats['ec_by_type'][f['type']] += 1
            if tally['xref']:
                stats['xref_by_type'][f['type']] += 1

        cursor = f['end_line']

    # Trailing lines after the last feature
    out.extend(bakta_lines[cursor:])

    # COMMENT addendum (separate pass, own tracked insertions)
    if add_comment_note:
        out, comment_inserted = insert_comment_addenda(out, baktfold_version)
        # Re-map: comment_inserted indices are relative to the NEW list; since
        # insert_comment_addenda rebuilds the list, shift stored qualifier
        # insertion indices are no longer valid against `out` -- so we verify
        # feature-qualifier insertions separately below via content, not index.
    else:
        comment_inserted = set()

    output_path.write_text('\n'.join(out), encoding='utf-8')

    stats['self_check'] = self_verify(
        bakta_lines, output_path, add_comment_note
    )
    return stats


# ---------------------------------------------------------------------------
# Self-verification: prove nothing else changed
# ---------------------------------------------------------------------------

def self_verify(bakta_lines: list[str], output_path: Path, comment_note_added: bool) -> bool:
    """Strip everything this script could have inserted out of the output and
    confirm the remainder is byte-identical to the original Bakta input.

    Strategy: rather than trusting recorded indices (fragile across the
    two-pass COMMENT insertion), re-derive it structurally: re-parse the
    output, and for every feature, its ORIGINAL Bakta qualifier lines must
    appear as a contiguous prefix of its output qualifier lines (everything
    after that prefix is, by construction, either absent or an appended
    line). We reconstruct the Bakta-equivalent text by dropping (a) any
    qualifier line not present at the same relative position in the
    matching Bakta feature, and (b) any COMMENT addendum block, then diff
    against the true original.
    """
    output_lines = output_path.read_text(encoding='utf-8', errors='replace').split('\n')

    bakta_feats = parse_feature_blocks(bakta_lines)
    output_feats = parse_feature_blocks(output_lines)

    if len(bakta_feats) != len(output_feats):
        print(f"SELF-CHECK FAIL: feature count changed "
              f"({len(bakta_feats)} -> {len(output_feats)})", file=sys.stderr)
        return False

    reconstructed: list[str] = []
    cursor = 0
    ok = True

    for bf, of in zip(bakta_feats, output_feats):
        # Copy inter-feature lines from OUTPUT, but strip a COMMENT addendum
        # block if present (recognizable by its literal marker line).
        between = output_lines[cursor:of['start_line']]
        if comment_note_added:
            between = _strip_addendum(between)
        reconstructed.extend(between)

        expected_block = bakta_lines[bf['start_line']:bf['end_line']]
        actual_block = output_lines[of['start_line']:of['end_line']]
        n = len(expected_block)
        if actual_block[:n] != expected_block:
            print(f"SELF-CHECK FAIL: feature at bakta line {bf['start_line']+1} "
                  f"was modified in place (not just appended to)", file=sys.stderr)
            ok = False
        reconstructed.extend(actual_block[:n])
        cursor = of['end_line']

    tail = output_lines[cursor:]
    if comment_note_added:
        tail = _strip_addendum(tail)
    reconstructed.extend(tail)

    if reconstructed != bakta_lines:
        print("SELF-CHECK FAIL: reconstructed text does not match original "
              "Bakta input byte-for-byte.", file=sys.stderr)
        # Show first differing line for debugging
        for i, (a, b) in enumerate(zip(reconstructed, bakta_lines)):
            if a != b:
                print(f"  first diff at line {i+1}:\n"
                      f"    original:      {b!r}\n"
                      f"    reconstructed: {a!r}", file=sys.stderr)
                break
        return False

    return True


def _strip_addendum(lines: list[str]) -> list[str]:
    marker = '            ##Post-processing:##'
    if marker not in lines:
        return lines
    idx = lines.index(marker)
    # Addendum starts one line above the marker (a blank continuation line)
    # and runs through the fixed-length block built in build_addendum().
    # We strip contiguous "            " (12-space) continuation lines that
    # follow the blank separator immediately preceding the marker.
    start = idx - 1 if idx > 0 and lines[idx - 1].strip() == '' else idx
    end = idx + 1
    while end < len(lines) and lines[end].startswith('            ') and lines[end].strip() != '':
        end += 1
    return lines[:start] + lines[end:]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Graft Baktfold's genuine additions onto the original, pristine Bakta .gbff.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('bakta', type=Path, help='Original Bakta .gbff file')
    parser.add_argument('baktfold', type=Path, help='Baktfold-processed .gbff file')
    parser.add_argument('output', type=Path, help='Output merged .gbff file')
    parser.add_argument('--no-comment-note', action='store_true',
                         help="Don't append the merge provenance note to COMMENT")
    parser.add_argument('--no-inference-provenance', action='store_true',
                         help="Don't append a per-feature /inference= line for grafted additions")
    parser.add_argument('--force', action='store_true',
                         help='Proceed even if the structural parity check fails')
    args = parser.parse_args()

    if not args.bakta.exists():
        sys.exit(f"Error: Bakta file not found: {args.bakta}")
    if not args.baktfold.exists():
        sys.exit(f"Error: Baktfold file not found: {args.baktfold}")

    print(f"Bakta input:    {args.bakta}")
    print(f"Baktfold input: {args.baktfold}")
    print(f"Output:         {args.output}")
    print()

    stats = graft(
        args.bakta, args.baktfold, args.output,
        add_comment_note=not args.no_comment_note,
        add_inference_provenance=not args.no_inference_provenance,
        force=args.force,
    )

    print("Grafting summary:")
    print(f"  Bakta features processed:          {stats['features_processed']:,}")
    print(f"  Matched in Baktfold:                {stats['features_matched']:,}")
    print(f"  Features receiving any addition:    {stats['features_with_any_addition']:,}")
    print(f"  /gene= added:                        {stats['gene_added']:,}  {dict(stats['gene_by_type'])}")
    print(f"  /EC_number= added:                   {stats['ec_added']:,}  {dict(stats['ec_by_type'])}")
    print(f"  /db_xref= added (afdb/cath/pdb):     {stats['xref_added']:,}  {dict(stats['xref_by_type'])}")
    print()
    check = stats['self_check']
    print(f"Self-verification (output minus insertions == original Bakta): "
          f"{'PASS' if check else 'FAIL'}")
    if not check:
        print("WARNING: output written but did not pass self-verification. "
              "Inspect before using.", file=sys.stderr)
        return 1
    print()
    print(f"Done. Output written to: {args.output}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
