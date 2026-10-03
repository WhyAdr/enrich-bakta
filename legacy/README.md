# Legacy utilities

`normalize_baktfold.py` is retained for historical reproducibility only. It is
not imported by the supported package, is excluded from Ruff and packaging
quality gates, and must not be used as the byte-preserving enrichment path.

New workflows should use the canonical modules under `src/enrich_bakta_lib/`
or the compatibility launchers at repository root.
