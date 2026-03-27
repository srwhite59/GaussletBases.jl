# Radial XGaussian Paper-Parity Contract Note

This note records how `GaussletBases` should interpret the recovered
paper-parity radial `x`-Gaussian contract summarized in the durable archive
memo at:

- `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/reports/software_reviews/gaussletbases_radial_xgaussian_paper_parity_contract_2026-03-26.md`

## Settled paper-parity reference

For the radial-gausslet manuscript line, the canonical paper-parity two-function
supplement is:

- `alpha = 0.0936`
- `alpha = 0.0236`

The stronger recovered legacy support is numerically very close:

- `0.09358986806`
- `0.02357750369`

The recovered paper-parity quadrature/evaluation contract for that manuscript
line is:

- `h = 0.001`
- `sigma = 3`
- `s0 = 6.5`
- `rmax_int = 80`

`s0 = 7.5` is not currently established as better for that paper-parity line.

## Repo interpretation

The repo default two-`xgaussian` supplement is already paper parity by design.
That default is meant to align the standard repo-facing radial supplement with
the published manuscript contract rather than with later drifting legacy cache
settings.

This does **not** automatically imply that every current repo engineering
default should adopt the same quadrature/evaluation contract. In particular, the
recovered paper-parity quadrature contract above should be treated as the
manuscript/reference line unless and until a separate repo decision adopts more
of it as a runtime default.

## Non-parity legacy drift

Later live-cache settings that drift away from the paper-parity pair belong to a
different contract and should not be treated as the manuscript reference by
default.

So for future repo comparisons:

- paper-parity comparisons should use the paper pair and the recovered
  manuscript quadrature contract explicitly
- later drifting legacy cache settings should be labeled as a different
  provenance family, not folded into the main paper-parity statement
