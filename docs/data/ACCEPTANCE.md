# Local dataset acceptance record

Date: 2026-10-03 (Asia/Jakarta)

This is a local, reproducible acceptance record for the materialized files in
the working tree. The complete local inventory is still separate from the
curated public bundle. The curated files are owner-authorized for publication
under CC BY 4.0, while their scientific status remains
`historical_unverified`.

Environment:

- Python 3.12.10 (Windows x86-64)
- Biopython 1.87
- pytest 9.1.1, Ruff 0.16.1, mypy 2.3.0
- Dataset inventory: [`MANIFEST.json`](MANIFEST.json), validated with
  `python tools/validate_dataset_manifest.py docs/data/MANIFEST.json`
- Curated publication inventory: [`PUBLISHED-MANIFEST.json`](PUBLISHED-MANIFEST.json),
  validated with `python tools/validate_dataset_manifest.py
  docs/data/PUBLISHED-MANIFEST.json`

All commands used the input paths and SHA-256 values recorded in that manifest.
The Baktfold commands used `--baktfold-invalid-ec-policy skip` explicitly for
the historical provisional tokens `3.5.1.n3`, `3.6.5.n1`, and `4.2.2.n1`.
Those values were recorded as `invalid_value` and were never promoted.

| Check | Result |
|---|---|
| C14 Baktfold parity/planning | Passed; parity passed, 0 translation mismatches, 8 invalid-token occurrences recorded, 800 genes, 178 ECs, 58 structural xrefs, 2,026 planned insertions |
| SM Baktfold parity/planning | Passed; parity passed, 0 translation mismatches, 8 invalid-token occurrences recorded, 780 genes, 168 ECs, 98 structural xrefs, 2,023 planned insertions |
| C14 Baktfold + Kofam one-pass | Passed with tool version `0.3.0`; 3,397 Kofam hits; output SHA-256 `64103efd00493fe725e57e1a828fcde40ef0884ecfbaf43f2291c6fbe301ad8b` |
| SM Baktfold + Kofam one-pass | Passed with tool version `0.3.0`; 3,333 Kofam hits; output SHA-256 `4e2a6e21762aa6469105c09c89c84794d645a4e643c983a592aa3097250a87ee` |
| C14 `import-faa` restoration → eggNOG | Passed; 4 restored targets, downstream required the evidence manifest and explicit imported-translation opt-in, 29,392 emitted values; output SHA-256 `3a6c4cebcfb497a9f0c0b754ae677effa292e8af909a4249a1abea675f3affd7` |
| SM `import-faa` restoration → eggNOG | Passed; 6 restored targets, downstream required the evidence manifest and explicit imported-translation opt-in, 29,086 emitted values; output SHA-256 `6513327169594fb9c09b18a5ced830aca66f06ee3baa6beae60eb1fdd955878b` |
| C14 restored unified Baktfold + Kofam + eggNOG | Passed with tool version `0.3.0`; 52,658 candidate decisions, 37,665 emitted candidates, 8 invalid-token occurrences; output SHA-256 `6c48d23dd3ef6dfd28cd4bd6416842462cc5d5af6bc03e1d3cf48b1cbdcda77a` |
| BK71A | Structural parse passed: 2 records, 8,640 features, 4,131 CDSs. No pristine Bakta FAA/evidence inputs are present, so no functional rerun is claimed. |

The default Baktfold policy remains `reject`; `skip` is an explicit historical
compatibility choice and must remain visible in any future rerun manifest.

The C14 and SM `import-faa` results are deliberately labeled as imported
translation evidence: downstream enrichment required the restoration manifest
and `--allow-imported-translations`. BK71A is structural-only because no
pristine Bakta FAA/evidence bundle is present.

Repeating the C14 and SM Baktfold+Kofam and eggNOG-only commands with identical
inputs produced byte-identical GBFF hashes. The repeat outputs and manifests
were written under the ignored `.test-output/real-gates/` directory.
