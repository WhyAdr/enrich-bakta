# Pinned eggNOG fixtures

## Official upstream fixture

- Source: https://github.com/eggnogdb/eggnog-mapper/blob/9aef364017956369a6216a2fbe276473e8a2f0d3/tests/fixtures_v7/test_diamond.emapper.annotations
- Upstream commit: `9aef364017956369a6216a2f0d3`
- File SHA-256: `bbb509d84d5a88d7d67056333bcb960e8367b549f3d7cc46eada0132ff4a8a5a`
- Size: 25,591 bytes, 43 lines
- Producer text: `emapper-v3.0.0-beta5-2-ga3fdb5b`

This file is retained as an upstream-format and version-parsing fixture. It is
not relicensed by enrich-bakta; preserve the upstream notice and attribution.

## Owner-controlled beta6 positive fixture

`beta6-owner-sample.annotations` is a minimal extracted parser fixture from
the owner-authorized C14 annotation artifact:

- Source artifact: `data/C14/evidence/C14-NMZ-query.emapper.annotations`
- Source artifact SHA-256: `ef1b713b307b2f2a14947edc6808d70cf0021b0e8b66e7ab235ce3ddbc957987`
- Producer text: `emapper-3.0.0-beta6`

It exercises automatic selection of the registered beta6 schema without
vendoring the complete C14 evidence table into the test package.
