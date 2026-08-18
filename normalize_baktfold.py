#!/usr/bin/env python3
"""
normalize_baktfold.py
=====================

Restore Bakta v1.12.0 GenBank qualifier format and provenance that Baktfold
v0.1.0 overwrote or dropped, while preserving Baktfold's genuine new
annotations (gene symbols for ncRNAs/rRNAs/tRNAs, AFDB/CATH/PDB cross-
references, and additional EC numbers).

WHAT IS RESTORED (Baktfold -> Bakta):
  1. /translation_table=False      ->  /transl_table=11            (every CDS)
  2. /db_xref="EC:X.Y.Z.W"         ->  /EC_number="X.Y.Z.W"
  3. /db_xref="<prefix>:<value>"   ->  /note="<prefix>:<value>"
     for prefix in {UniRef, UniParc, RefSeq, KEGG, BlastRules, COG,
                    PFAM, IS, NCBIFam, NCBIProtein, SO}
  4. /db_xref="GO:..." and /db_xref="RFAM:..."   ->  kept as /db_xref=
     (these were /db_xref= in Bakta already)
  5. /db_xref="afdb_v6:..." / "cath:..." / "pdb:..."  ->  kept
     (new cross-references added by Baktfold)
  6. /protein_id="gnl|Baktfold|XXX"  ->  /protein_id="gnl|Bakta|XXX"
  7. /inference="ab initio prediction:Bakta:0.1"
        ->  /inference="ab initio prediction:Bakta:1.12"
  8. /inference="profile:INFERNAL:1.1.5"   (on tmRNA)
        ->  /inference="profile:aragorn:1.2"
  9. /inference="profile:tRNAscan-SE:2.0.12"  (on tRNAs)
        ->  /inference="profile:tRNAscan:2.0"
 10. LOCUS date "01-JAN-1980"  ->  user-supplied date (default 08-JUL-2026)
 11. COMMENT block  ->  restored to Bakta v1.12.0 provenance, with a short
     post-processing note appended documenting that Baktfold annotations
     were merged in.

WHAT IS PRESERVED (Baktfold additions kept):
  - /gene= additions on rRNA (30), ncRNA (36), tRNA (3), CDS (292), gene (352)
  - /db_xref="afdb_v6:..." (51), /db_xref="cath:..." (16), /db_xref="pdb:..." (14)
  - Additional /EC_number= entries Baktfold introduced (101 new EC assignments)

WHAT IS NOT RESTORED:
  - Per-feature qualifier ordering. Baktfold alphabetically sorts qualifiers
    within each feature; Bakta uses a feature-type-specific order. The output
    remains valid GenBank either way - downstream tools should not depend on
    qualifier order.

USAGE:
  python normalize_baktfold.py INPUT.gbff OUTPUT.gbff
      [--date 08-JUL-2026]
      [--annotation-date "07/08/2026, 11:50:32"]
      [--bakta-version 1.12.0]
      [--bakta-db "v6.0, full"]

EXAMPLE:
  python /home/z/my-project/scripts/normalize_baktfold.py \\
      /home/z/my-project/upload/BK71A-baktfold.txt \\
      /home/z/my-project/download/BK71A-restored.gbff
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Transform tables
# ---------------------------------------------------------------------------

# Prefixes that Bakta originally stored as /note= (Baktfold moved them to
# /db_xref=). Converting back restores Bakta's qualifier format.
NOTE_PREFIXES = {
    "UniRef",
    "UniParc",
    "RefSeq",
    "KEGG",
    "BlastRules",
    "COG",
    "PFAM",
    "IS",
    "NCBIFam",
    "NCBIProtein",
    "SO",
}

# Inference string replacements: Baktfold -> Bakta.
INFERENCE_REPLACEMENTS = {
    "ab initio prediction:Bakta:0.1": "ab initio prediction:Bakta:1.12",
    "profile:INFERNAL:1.1.5":         "profile:aragorn:1.2",
    "profile:tRNAscan-SE:2.0.12":     "profile:tRNAscan:2.0",
}


# ---------------------------------------------------------------------------
# Bakta COMMENT block template (matching the exact whitespace alignment
# Bakta v1.12.0 emits).
# ---------------------------------------------------------------------------

def build_bakta_comment(
    bakta_version: str,
    bakta_db: str,
    annotation_date: str,
    cds_count: int,
    trna_count: int,
    tmrna_count: int,
    rrna_count: int,
    ncrna_count: int,
    reg_count: int,
    crispr_count: int,
    ori_count: int,
    orit_count: int,
    gaps_count: int,
    pseudogene_count: int,
) -> list[str]:
    """Construct a Bakta-style COMMENT block with exact alignment.

    Bakta left-justifies each label to width 31, immediately followed by
    ":: " and the value (no extra space before the colons).
    """
    def row(label: str, value) -> str:
        return f"            {label:<31}:: {value}"

    return [
        "COMMENT     Annotated with Bakta",
        f"            Software: v{bakta_version}",
        f"            Database: {bakta_db}",
        "            DOI: 10.1099/mgen.0.000685",
        "            URL: github.com/oschwengers/bakta",
        "            ",
        "            ##Genome Annotation Summary:##",
        row("Annotation Date", annotation_date),
        row("CDSs", f"{cds_count:>5,}"),
        row("tRNAs", f"{trna_count:>5,}"),
        row("tmRNAs", f"{tmrna_count:>5,}"),
        row("rRNAs", f"{rrna_count:>5,}"),
        row("ncRNAs", f"{ncrna_count:>5,}"),
        row("regulatory ncRNAs", f"{reg_count:>5,}"),
        row("CRISPR Arrays", f"{crispr_count:>5,}"),
        row("oriCs/oriVs", f"{ori_count:>5,}"),
        row("oriTs", f"{orit_count:>5,}"),
        row("gaps", f"{gaps_count:>5,}"),
        row("pseudogenes", f"{pseudogene_count:>5,}"),
        "            ",
        "            Post-processed by normalize_baktfold.py to restore Bakta",
        "            provenance and qualifier format. Additional /gene=,",
        "            /EC_number=, and /db_xref=(afdb_v6|cath|pdb) annotations",
        "            contributed by Baktfold v0.1.0",
        "            (https://github.com/gbouras13/baktfold) are preserved.",
    ]


# ---------------------------------------------------------------------------
# First pass: count features and pseudogenes from the input file
# ---------------------------------------------------------------------------

def count_features(input_path: Path) -> dict:
    """First pass: tally feature types and pseudogene count.

    The COMMENT block in Bakta reports:
      - CDSs           = CDS features minus pseudogene CDS
      - tRNAs/tmRNAs/... = direct feature counts
      - pseudogenes    = (features with /pseudogene=) / 2  (gene + CDS pair)
    """
    counts: dict[str, int] = {}
    pseudogene_qualifier_count = 0
    feat_re = re.compile(r"^     (\w+)\s+\S")

    with input_path.open("r") as f:
        for line in f:
            m = feat_re.match(line)
            if m:
                ftype = m.group(1)
                counts[ftype] = counts.get(ftype, 0) + 1
            if "/pseudogene=" in line:
                pseudogene_qualifier_count += 1

    pseudogene_count = pseudogene_qualifier_count // 2  # gene + CDS each have it
    cds_total = counts.get("CDS", 0)
    cds_comment = cds_total - pseudogene_count  # Bakta's COMMENT convention

    return {
        "CDSs": cds_comment,
        "tRNAs": counts.get("tRNA", 0),
        "tmRNAs": counts.get("tmRNA", 0),
        "rRNAs": counts.get("rRNA", 0),
        "ncRNAs": counts.get("ncRNA", 0),
        "regulatory": counts.get("regulatory", 0),
        "CRISPR": counts.get("CRISPR", 0),
        "oriCs/oriVs": counts.get("rep_origin", 0),
        "oriTs": 0,  # Bakta reports oriTs separately; rep_origin covers oriC/oriV
        "gaps": 0,
        "pseudogenes": pseudogene_count,
    }


# ---------------------------------------------------------------------------
# Per-line qualifier restoration
# ---------------------------------------------------------------------------

DBXREF_TO_NOTE_RE = re.compile(r'/db_xref="([^:"]+):([^"]+)"')
DBXREF_EC_RE = re.compile(r'/db_xref="EC:([^"]+)"')


def restore_qualifier_line(line: str) -> tuple[str, str]:
    """Restore a single /qualifier line.

    Returns (new_line, change_type) where change_type is one of:
      'none', 'transl_table', 'ec_number', 'note_from_dbxref',
      'protein_id', 'inference'
    """
    # 1. /translation_table=False  ->  /transl_table=11
    if "/translation_table=False" in line:
        return line.replace("/translation_table=False", "/transl_table=11"), "transl_table"

    # 2. /db_xref="EC:X"  ->  /EC_number="X"
    m = DBXREF_EC_RE.search(line)
    if m:
        ec_value = m.group(1)
        old = f'/db_xref="EC:{ec_value}"'
        new = f'/EC_number="{ec_value}"'
        return line.replace(old, new), "ec_number"

    # 3. /db_xref="<NOTE_PREFIX>:..."  ->  /note="<NOTE_PREFIX>:..."
    m = DBXREF_TO_NOTE_RE.search(line)
    if m and m.group(1) in NOTE_PREFIXES:
        prefix, value = m.group(1), m.group(2)
        old = f'/db_xref="{prefix}:{value}"'
        new = f'/note="{prefix}:{value}"'
        return line.replace(old, new), "note_from_dbxref"

    # 4. protein_id="gnl|Baktfold|..."  ->  protein_id="gnl|Bakta|..."
    if 'gnl|Baktfold|' in line:
        return line.replace('gnl|Baktfold|', 'gnl|Bakta|'), "protein_id"

    # 5. Inference string replacements
    for old, new in INFERENCE_REPLACEMENTS.items():
        if old in line:
            return line.replace(old, new), "inference"

    return line, "none"


# ---------------------------------------------------------------------------
# Main normalization
# ---------------------------------------------------------------------------

def normalize(input_path: Path, output_path: Path, args, feat_counts: dict) -> dict:
    """Stream the input file, restoring Bakta format line by line."""
    stats = {
        "locus_date_restored": 0,
        "comment_restored": 0,
        "transl_table_restored": 0,
        "ec_number_restored": 0,
        "note_restored_from_dbxref": 0,
        "protein_id_restored": 0,
        "inference_restored": 0,
        "total_lines": 0,
    }

    bakta_comment_lines = build_bakta_comment(
        bakta_version=args.bakta_version,
        bakta_db=args.bakta_db,
        annotation_date=args.annotation_date,
        cds_count=feat_counts["CDSs"],
        trna_count=feat_counts["tRNAs"],
        tmrna_count=feat_counts["tmRNAs"],
        rrna_count=feat_counts["rRNAs"],
        ncrna_count=feat_counts["ncRNAs"],
        reg_count=feat_counts["regulatory"],
        crispr_count=feat_counts["CRISPR"],
        ori_count=feat_counts["oriCs/oriVs"],
        orit_count=feat_counts["oriTs"],
        gaps_count=feat_counts["gaps"],
        pseudogene_count=feat_counts["pseudogenes"],
    )

    in_comment = False
    output_lines: list[str] = []

    with input_path.open("r") as f:
        for line in f:
            stats["total_lines"] += 1
            line_no_nl = line.rstrip("\n")

            # 1. LOCUS date restoration
            if line.startswith("LOCUS") and "01-JAN-1980" in line:
                output_lines.append(line_no_nl.replace("01-JAN-1980", args.date))
                stats["locus_date_restored"] += 1
                continue

            # 2. COMMENT block: skip Baktfold's COMMENT, emit restored Bakta one
            if line.startswith("COMMENT"):
                in_comment = True
                continue
            if in_comment:
                if line.startswith("FEATURES"):
                    output_lines.extend(bakta_comment_lines)
                    output_lines.append(line_no_nl)
                    stats["comment_restored"] += 1
                    in_comment = False
                # else: still inside COMMENT, skip the Baktfold version
                continue

            # 3. Qualifier lines: apply restorations
            stripped = line_no_nl.lstrip()
            if stripped.startswith("/"):
                new_line, change_type = restore_qualifier_line(line_no_nl)
                if change_type != "none":
                    stats_key = {
                        "transl_table": "transl_table_restored",
                        "ec_number": "ec_number_restored",
                        "note_from_dbxref": "note_restored_from_dbxref",
                        "protein_id": "protein_id_restored",
                        "inference": "inference_restored",
                    }[change_type]
                    stats[stats_key] += 1
                output_lines.append(new_line)
            else:
                output_lines.append(line_no_nl)

    # Write the output
    with output_path.open("w") as f:
        for line in output_lines:
            f.write(line + "\n")

    return stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Restore Bakta v1.12.0 format from a Baktfold-processed .gbff file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input", type=Path, help="Input Baktfold .gbff file")
    parser.add_argument("output", type=Path, help="Output restored .gbff file")
    parser.add_argument(
        "--date",
        default="08-JUL-2026",
        help="LOCUS date stamp (default: 08-JUL-2026)",
    )
    parser.add_argument(
        "--annotation-date",
        default="07/08/2026, 11:50:32",
        help='Annotation Date value in COMMENT (default: "07/08/2026, 11:50:32")',
    )
    parser.add_argument(
        "--bakta-version",
        default="1.12.0",
        help="Bakta software version string (default: 1.12.0)",
    )
    parser.add_argument(
        "--bakta-db",
        default="v6.0, full",
        help='Bakta database version string (default: "v6.0, full")',
    )
    # Override the dynamically computed COMMENT counts. Bakta's internal
    # "CDSs" count (4,074) differs from a naive (total CDS - pseudogenes)
    # because Bakta classifies some CDS differently. Pass these to match
    # the original Bakta COMMENT exactly.
    parser.add_argument(
        "--cds-count", type=int, default=None,
        help="Override CDSs count in COMMENT (default: CDS features - pseudogenes)",
    )
    parser.add_argument(
        "--pseudogene-count", type=int, default=None,
        help="Override pseudogenes count in COMMENT (default: /pseudogene= qualifiers / 2)",
    )
    parser.add_argument("--trna-count", type=int, default=None, help="Override tRNAs count")
    parser.add_argument("--tmrna-count", type=int, default=None, help="Override tmRNAs count")
    parser.add_argument("--rrna-count", type=int, default=None, help="Override rRNAs count")
    parser.add_argument("--ncrna-count", type=int, default=None, help="Override ncRNAs count")
    parser.add_argument("--regulatory-count", type=int, default=None, help="Override regulatory ncRNAs count")
    parser.add_argument("--crispr-count", type=int, default=None, help="Override CRISPR Arrays count")
    parser.add_argument("--ori-count", type=int, default=None, help="Override oriCs/oriVs count")
    parser.add_argument("--orit-count", type=int, default=None, help="Override oriTs count")
    parser.add_argument("--gaps-count", type=int, default=None, help="Override gaps count")
    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(f"Error: input file not found: {args.input}")

    print(f"Normalizing: {args.input}")
    print(f"Output:      {args.output}")
    print(f"LOCUS date:  {args.date}")
    print(f"Bakta ver:   v{args.bakta_version}, DB {args.bakta_db}")
    print()

    stats, feat_counts = normalize_with_overrides(args)

    print("Feature counts used in COMMENT:")
    for k, v in feat_counts.items():
        print(f"  {k:<20} {v}")
    print()
    print("Restoration summary:")
    print(f"  LOCUS dates restored:              {stats['locus_date_restored']}")
    print(f"  COMMENT blocks restored:           {stats['comment_restored']}")
    print(f"  /transl_table=11 restored:         {stats['transl_table_restored']}")
    print(f"  /EC_number= restored:              {stats['ec_number_restored']}")
    print(f"  /note= restored from /db_xref=:    {stats['note_restored_from_dbxref']}")
    print(f"  protein_id namespace restored:     {stats['protein_id_restored']}")
    print(f"  inference strings restored:        {stats['inference_restored']}")
    print(f"  Total lines processed:             {stats['total_lines']:,}")
    print()
    print(f"Done. Output written to: {args.output}")
    return 0


def normalize_with_overrides(args) -> tuple[dict, dict]:
    """Wrap normalize() to apply CLI count overrides."""
    feat_counts = count_features(args.input)
    if args.cds_count is not None:
        feat_counts["CDSs"] = args.cds_count
    if args.pseudogene_count is not None:
        feat_counts["pseudogenes"] = args.pseudogene_count
    if args.trna_count is not None:
        feat_counts["tRNAs"] = args.trna_count
    if args.tmrna_count is not None:
        feat_counts["tmRNAs"] = args.tmrna_count
    if args.rrna_count is not None:
        feat_counts["rRNAs"] = args.rrna_count
    if args.ncrna_count is not None:
        feat_counts["ncRNAs"] = args.ncrna_count
    if args.regulatory_count is not None:
        feat_counts["regulatory"] = args.regulatory_count
    if args.crispr_count is not None:
        feat_counts["CRISPR"] = args.crispr_count
    if args.ori_count is not None:
        feat_counts["oriCs/oriVs"] = args.ori_count
    if args.orit_count is not None:
        feat_counts["oriTs"] = args.orit_count
    if args.gaps_count is not None:
        feat_counts["gaps"] = args.gaps_count
    stats = normalize(args.input, args.output, args, feat_counts)
    return stats, feat_counts


if __name__ == "__main__":
    sys.exit(main())
