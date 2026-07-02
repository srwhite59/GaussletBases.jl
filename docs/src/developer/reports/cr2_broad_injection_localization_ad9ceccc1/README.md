# Cr2 Broad Original Localization Audit (ad9ceccc1)

## Plain-Language Result

Z-localization changes the interpretation of the broad-original set. The nonlocalized SVD dropped set mostly resolves into atom-local localized modes, not one smooth bridge between Cr centers. Summed across the 14 dropped SVD directions, the localization overlap is about `12.245` atom-local, `0.807` bond-axis/midbond, `0.0` transverse-bond, and `0.944` outer-tail; these sums are out of 14 directions, not normalized to one.

This supports a possible future rule that is stricter for bond-spanning broad directions and slightly more permissive for clearly atom-local broad directions. It does not justify a simple global relaxation of the `0.99` representability cut.

The currently dropped 14 nonlocalized directions are not safely described as dangerous bond-spanning modes either. Most of their weight becomes atom-local after z localization, with one notable outer-tail-heavy dropped direction and modest bond-axis admixture in several modes.

## Decision Answers

- Dropped 14 directions: mostly atom-local after localization, with one outer-tail-heavy direction and smaller bond-axis admixtures.
- High-fake localized 0.95..0.99 directions: 20 atom-local modes appear in this band.
- `0.99` strictness: likely too strict for some atom-local modes, but still defensible for bond-spanning modes.
- Fate rule: localization helps, but the policy is now multi-criterion: fake-RDM, representability, and shape class. A scalar representability-only rule would hide this distinction.

## Dimensions And Baseline

- Artifact recipe source: `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/cr2_fixed95fec2b8_ida.jld2`
- Current source commit: `ad9ceccc1`
- Cr2 setup: `ns=7`, `lmax=2`, `cc-pV5Z`, `uncontracted=false`
- Current base dimension: `6915`
- Compact residual count: `30`
- Broad original directions: `108`
- Nonlocalized strict reduced set: `94` kept, `14` dropped
- Full broad+protected B: min `7.732092411339e-01`, counts `<0.999/<0.99/<0.9 = 40/14/4`
- Reduced B: min `9.905723760441e-01`, count `<0.99 = 0`

## Localized Counts

- Atom-local modes: `88`
- Bond-axis/midbond modes: `8`
- Transverse-bond modes: `8`
- Outer-tail modes: `4`
- Ambiguous modes: `0`
- Strict `fake>=1e-3 && repr>=0.99`: `80`
- Relaxed atom-local `fake>=1e-3 && repr>=0.95`: `20`
- High-fake bond-spanning below strict representability: `2`

## Method

The probe rebuilt the current Cr2 aspect-shell basis from the old artifact recipe, then used the source-backed ordered compact-first residual selector to build `M = [G, R_compact]`. It protected the 30 compact originals, projected them out of the remaining broad originals, Gram-cleaned the broad remainder in the exact GTO overlap, and diagonalized the exact Gaussian `z` moment matrix in that S-orthonormal broad subspace.

Representability is `||M' S z_i||` for each localized mode. Fake-RDM is the same diagnostic proxy used in the earlier Cr2 representability audit. Shape classes are diagnostics, not production selection defaults.

This report is consistent with the repo guardrails in `docs/src/developer/architecture/gausslet_algorithm_refresher.md` and `docs/src/developer/architecture/gausslet_methods_fundamentals.md`: injected or protected Gaussian directions are not residual MWG channels, and diagnostic weights are not IDA quadrature authority.

## Files

- `localized_broad_modes.tsv`: one row per z-localized broad mode
- `localized_rule_summary.tsv`: counts and fake-RDM sums by rule bucket and shape
- `dropped_direction_localization.tsv`: where the current 14 dropped SVD directions localize
- `nonlocalized_vs_localized.tsv`: overlap map between SVD-broad and z-localized modes

## Validation

- Probe elapsed: `5.489240250100e+01 s`
- Package load and `git diff --check` were run by repo-doer in the same pass.
