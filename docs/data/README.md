# Dataset inventory

The repository-local layout is:

```text
data/
  C14/{bakta,baktfold,evidence,derived/{restored,enriched}}
  SM/{bakta,baktfold,evidence,derived/{restored,enriched}}
  BK71A/derived/restored/
  fixtures/
```

`bakta/`, `baktfold/`, and `evidence/` contain the historical material found
in the working tree. `derived/` contains restored or enriched outputs and must
never be used as the pristine identity reference for a new merge. `fixtures/`
is reserved for small synthetic files suitable for ordinary Git and CI.

The complete local inventory is [`MANIFEST.json`](MANIFEST.json); the generated
[`MANIFEST.tsv`](MANIFEST.tsv) is a review-friendly projection and is not the
authoritative machine-readable document. The owner-authorized curated
publication bundle is described separately by
[`PUBLISHED-MANIFEST.json`](PUBLISHED-MANIFEST.json) and its TSV projection.
See [`PUBLISHING.md`](PUBLISHING.md) for the exact selection policy.
The curated data grant is recorded in
[`DATA-LICENSE-CC-BY-4.0.txt`](../../DATA-LICENSE-CC-BY-4.0.txt).

The inventory is intentionally separate from the code and records the current
local SHA-256 and byte size for every materialized artifact. Historical
artifacts remain `historical_unverified` when producer versions, database
versions, or independent lineage are not fully confirmed, even when their
redistribution is owner-authorized.

The dated local functional/provenance checks are recorded in
[`ACCEPTANCE.md`](ACCEPTANCE.md); they do not change the unverified status of
the materialized inputs.

Validate the materialized inventory from the repository root with:

```text
python tools/refresh_dataset_manifest.py docs/data/MANIFEST.json docs/data/MANIFEST.json
python tools/validate_dataset_manifest.py docs/data/MANIFEST.json
python tools/render_dataset_manifest_tsv.py docs/data/MANIFEST.json docs/data/MANIFEST.tsv
python tools/validate_dataset_manifest.py docs/data/PUBLISHED-MANIFEST.json
```

The current C14/SM Baktfold files contain the legacy tokens
`3.5.1.n3`, `3.6.5.n1`, and `4.2.2.n1`. They are outside the supported EC
grammar and are never promoted by default; historical reruns may pass
`--baktfold-invalid-ec-policy skip` to record them as explicit
`invalid_value` decisions. This does not broaden the grammar or validate the
underlying database assignment.

Large biological inputs remain ignored by default. The curated owner-authorized
bundle uses Git LFS; a clean clone must materialize the LFS bytes and pass the
published-manifest checksum check.
