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

## Follow-through acceptance — 2026-10-04

This section is a new acceptance run; the 2026-10-03 observations above remain
historical and are not rewritten. The validated code revision was
`e5def78aca8bac14f9ecc2d213d6042d0ccd56ea` on Windows 11
(`10.0.26200`, x86-64), Python 3.12.10. The primary run used Biopython 1.87,
openpyxl 3.1.5, jsonschema 4.26.0, and setuptools 84.0.0. A separate isolated
floor environment used Biopython 1.83 and setuptools 77.0.3.

Commands and policy flags:

```text
git lfs pull
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
python enrich-bakta-followthrough-bundle-2026-10-04/run_real_checks.py
python enrich-bakta-followthrough-bundle-2026-10-04/run_chain_checks.py
python tools/run_followthrough_real_checks.py --output .test-output/committed-followthrough-e5def78
python tools/run_followthrough_chain_checks.py --output .test-output/committed-followthrough-e5def78
```

The Baktfold-containing checks used `--baktfold-invalid-ec-policy skip`.
Unified checks supplied `--kofamscan-version 1.3.0`, the latest parent
translation-evidence manifest, and `--allow-imported-translations`. The
standalone Baktfold+Kofam checks intentionally omitted a Kofam version to match
the supplied comparison. No merge timestamp was supplied.

| Check | C14 | SM |
|---|---:|---:|
| Published artifacts | 33 total artifacts passed size, SHA-256, path, and lineage validation | Same shared gate |
| TSV/XLSX normalized hits | 4,563; exact supported-field equality | 4,497; exact supported-field equality |
| Imported restoration targets retained | 4 | 6 |
| Baktfold + Kofam insertions | 10,522 | 10,380 |
| Baktfold + Kofam output SHA-256 | `7a2a8482346b60503cf6f03a0a654294b0ee23cc4af0ec22f8e0016295907862` | `63bb4e91beb0035d0b18f0a4d5538429d317940a6fe6d9c4fa78f9348976360b` |
| eggNOG insertions | 33,944 | 33,567 |
| eggNOG output SHA-256 | `3a6c4cebcfb497a9f0c0b754ae677effa292e8af909a4249a1abea675f3affd7` | `6513327169594fb9c09b18a5ced830aca66f06ee3baa6beae60eb1fdd955878b` |
| Unified insertions | 41,500 | 41,013 |
| Unified output SHA-256 | `90cdfc6897410cdc20fd4e3f83aa2b316463306831e2b30ca15566568458c117` | `ad3408c57bd5982c4e0bc9c638f290da08c7867d84f13c155ffbfe608c9fb1b6` |
| Restoration, eggNOG, unified output-as-input reruns | 0 insertions; exact GBFF byte identity; 4 imported origins retained | 0 insertions; exact GBFF byte identity; 6 imported origins retained |
| Imported-protein opt-in removal | Blocked contextually | Blocked contextually |

The isolated dependency-floor environment passed all 133 tests and built the
wheel and sdist with `python -m build --no-isolation`. A clean environment
outside the checkout passed `pip check`, preserved all six legacy module
aliases, and ran all five console commands plus both module entry points.

The supplied bundle scripts completed in 293 seconds (real checks) and 196
seconds (chain checks). The checked-in durable gates completed in 302 and 436
seconds respectively on this machine. Correctness and reproducibility passed;
these variable Windows timings are recorded as observations, not a performance
acceptance claim. GitHub's Windows/Linux jobs remain the remote execution gate.

The accepted C14/SM contract is now machine-readable in
[`followthrough-scientific-baseline.json`](followthrough-scientific-baseline.json).
The durable gates fail on drift in input hashes or recorded producers, tool or
policy versions, first-pass output hashes and counts, imported query IDs or
origins, source-entry/final-decision projections, and suppression totals.
Intentional scientific changes require a separately reviewed baseline update.
The dataset and scientific jobs run on manual dispatch and on version-tag pushes;
their compact JSON reports are retained as workflow artifacts.
