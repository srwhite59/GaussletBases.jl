# Cartesian Nested IDA Weight-transfer Milestone

This note records the nesting branch-point reached in commit
`5b0460b9d561b2e15d69e00f79d8897a313fb916`.

## Milestone Summary

Exact bug, in one sentence:

- the nested fixed-block path was linearly contracting
  `pgdg.pair_factor_terms` even though those PGDG terms were already
  weight-divided IDA objects, so the contraction was being applied to the
  wrong object

Exact fix, in one sentence:

- the nested fixed-block path now reconstructs the raw IDA numerator, carries
  the contracted fixed-block integral weights through the same map, and only
  then divides by the new weight outer product

Why this was a real milestone:

- before this fix, the corrected complete-shell nonrecursive line still looked
  physically wrong because `⟨Vee⟩` was grossly inflated
- after this fix, that large discrepancy closes almost entirely, so the line is
  no longer blocked by a spurious fixed-block IDA transfer bug

## Legacy Rule That Had To Be Restored

In the current PGDG bundle, the one-dimensional pair-factor terms are already
stored as weight-divided IDA objects, not as raw numerators.

In
[ordinary_mapped_backends.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/ordinary_mapped_backends.jl),
the active path does:

```julia
pair_factors = [
    pair_factor_terms_raw[term, :, :] ./ weight_outer
    for term in eachindex(exponents_value)
]
pair_factor_terms = _term_tensor(pair_factors)
```

So `pgdg_intermediate.pair_factor_terms` are already divided by the finalized
one-dimensional integral weights.

The legacy nested/PGDG line did not contract those divided objects directly.
It restored the raw numerator first, transformed that, transformed the weights,
and only then reweighted the result. The same rule appears in:

- [PureGaussianGausslet.jl](/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl#L991)
  through
  [PureGaussianGausslet.jl](/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl#L1081)
- [facedecimate3steve.jl](/Users/srw/Dropbox/GaussletModules/facedecimate3steve.jl#L346)
  through
  [facedecimate3steve.jl](/Users/srw/Dropbox/GaussletModules/facedecimate3steve.jl#L351)
- [fd4.jl](/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/fd4.jl#L463)
  through
  [fd4.jl](/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/fd4.jl#L468)

The restored rule is:

```math
M = D(w)\,V\,D(w), \qquad
M' = C^T M C, \qquad
w' = C^T w, \qquad
V' = D(1/w')\,M'\,D(1/w').
```

## Numerical Before/After

Test case used throughout:

- He
- stabilized count-17 fixed-`a`
- `a = 1/4`
- `xmax = 10`
- `s = 0.626026121152214`

### Projected-parent Fixed-only `⟨Vee⟩`

| Case | Parent comparison | Before fix | After fix | After shift vs parent |
| --- | ---: | ---: | ---: | ---: |
| Parent fixed space | `1.2486322227019246` | `1.2486322227019246` | `1.2486322227019246` | `0.0` |
| Shell-plus-core | parent anchor | `1.2486612176560734` | `1.248633734959628` | `+1.5122577034e-6` |
| Corrected complete four-shell | parent comparison line | `1.818451883926288` | `1.2487650946681053` | `+1.3287196618e-4` |

### Full Nearest/GGT `1s^2` `⟨Vee⟩`

| Case | Shell-plus-core comparison | Before fix | After fix | After shift vs shell-plus-core |
| --- | ---: | ---: | ---: | ---: |
| Shell-plus-core | physical anchor | `1.2489197930296585` | `1.2489021100178215` | `0.0` |
| Corrected complete four-shell | nonrecursive comparison line | `1.816303754935915` | `1.2489346080730335` | `+3.2498055212e-5` |

The same full nearest/GGT check gives:

- shell-plus-core: `E1 = -1.9985629071304167`
- corrected complete four-shell: `E1 = -1.9981842264017424`

So the large previous `⟨Vee⟩` discrepancy is gone. The remaining difference is
small and no longer points to a broken interaction transfer.

## Fixed-block Weight Check

The transformed fixed-block weights are finite and strictly positive in the
intended regime.

- shell-plus-core:
  - min `0.061468156103430815`
  - max `2.977482219200775`
  - nonpositive entries: `0`
- corrected complete four-shell:
  - min `0.061468156103430815`
  - max `3.2591670283995673`
  - nonpositive entries: `0`

## Validation

Validation run for this milestone:

- `julia --project=. test/runtests.jl`
- `julia --project=docs docs/make.jl`

The docs build only emitted the usual Documenter no-deploy-environment warning.

## Current Trust Status

What is now trusted:

- shell-plus-core remains the trusted physical nonrecursive anchor
- the nested fixed-block IDA transfer in the current code now follows the
  correct legacy weight-aware rule
- the claim that the corrected complete-shell `⟨Vee⟩` blow-up was caused by the
  wrong IDA transfer, not by a deeper unresolved interaction failure

What is still only diagnostic:

- the corrected complete four-shell sequence as a fully promoted nonrecursive
  anchor for the roadmap as a whole
- any broader claim that the nonrecursive nesting line is completely settled
  beyond this He branch-point

## Practical Roadmap Consequence

This is no longer "another `⟨Vee⟩` bug hunt" because the dominant interaction
error has been explained, fixed, and reduced from an `O(10^{-1})` to `O(1)`
failure to a small `10^{-5}` to `10^{-4}` difference.

The next short trust-establishment pass should decide whether the corrected
complete-shell sequence is now good enough to promote from a diagnostic line to
the first trustworthy nonrecursive comparison state beyond shell-plus-core, or
whether one more small cleanup is needed before that promotion. In practice,
that pass should judge the corrected complete-shell line against the parent and
shell-plus-core anchors on:

- projected-parent fixed-only interaction transfer
- final nearest/GGT `⟨Vee⟩`
- remaining one-body / `E1` gap
- continued finite-positive transformed weights

If that pass is favorable, the nesting roadmap can move on from interaction
repair to trust establishment, and then back toward the radial public/tryout
cleanup path and later radial-vs-Cartesian comparison.
