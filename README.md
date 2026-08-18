# enrich-bakta

Tools and workflows for enriching Bakta bacterial genome annotations (`.gbff`) with novel annotations from **Baktfold** (AlphaFold/Foldseek structural homology) and **KofamScan/KofamKOALA** (KEGG Orthology profile HMMs) while preserving Bakta's original feature blocks, sequence, metadata, and qualifiers.

---

## Architectural Principles

1. **Graft-onto-Pristine (Ground Truth Integrity)**  
   Instead of trying to "restore" files modified/corrupted by downstream annotators (e.g., Baktfold v0.1.0's altered `/transl_table`, hardcoded dates, and modified `COMMENT` blocks), `enrich-bakta` treats the original Bakta `.gbff` as the authoritative source of truth and grafts *only* genuine novel evidence.
2. **Byte-Level Preservation & Verification**  
   Every pipeline step includes verification proving that removing newly inserted qualifier lines leaves the original Bakta file byte-for-byte unmodified.
3. **Structured & Granular Provenance**  
   Additions are accompanied by distinct evidence-class inference qualifiers (e.g. Foldseek structural similarity vs. KofamScan HMM profile hits) and record-level comment addenda.

---

## Workspace Structure & Components

| File / Component | Purpose |
|---|---|
| [`graft_baktfold_additions.py`](graft_baktfold_additions.py) | Grafts Baktfold structural cross-references (`/db_xref="afdb_v6:..."`, `cath:...`, `pdb:...`), novel `/EC_number`s, and unannotated `/gene` symbols onto pristine Bakta `.gbff` files. Includes strict self-verification. |
| [`normalize_baktfold.py`](normalize_baktfold.py) | Normalization utility to fix qualifier formatting, inference strings, and provenance in Baktfold outputs back to Bakta standards. |
| [`merge-koala-baktfold-plan.md`](merge-koala-baktfold-plan.md) | Engineering specification for single-pass merging of Baktfold structural additions and KofamScan KOALA assignments with strict parity and semantic checks. |

---

## Usage

### Grafting Baktfold Additions
```bash
python graft_baktfold_additions.py BAKTA.gbff BAKTFOLD.gbff OUTPUT.gbff
```
*Options:*
- `--no-comment-note`: Skip the comment addendum recording the merge.
- `--no-inference-provenance`: Do not append `/inference="protein structure similarity:Baktfold..."` qualifiers.
- `--force`: Force writing output even if non-fatal parity warnings occur.

### Normalizing Baktfold Files
```bash
python normalize_baktfold.py BAKTFOLD_INPUT.gbff RESTORED_OUTPUT.gbff
```
