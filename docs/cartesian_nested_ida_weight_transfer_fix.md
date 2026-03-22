# Cartesian Nested IDA Weight-transfer Fix

This note records the first fixed-fixed interaction correction after the
coverage-fixed shell-sequence diagnostics.

The starting diagnosis was:

- shell-plus-core remained the trusted physical nonrecursive anchor
- the corrected complete shell sequence already recovered the parent low-energy
  one-body content well
- `E1` was back in the right regime
- but `⟨Vee⟩` was still far too large, and the excess was already in the
  fixed-fixed block

That made the next question very specific:

- is the nested fixed-block path transforming the wrong IDA object?

## Current Code Meaning

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

Before this fix, the nested shell packet then linearly contracted those
already-divided `pair_terms` directly.

That is the wrong transfer for an IDA object once a nontrivial contraction is
applied.

## Legacy Evidence

The legacy nested/PGDG line does the weight-aware transform instead.

In
[PureGaussianGausslet.jl](/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl#L991),
the raw one-dimensional Coulomb numerators are first assembled and transformed:

- raw primitive-to-basis numerators at
  [PureGaussianGausslet.jl](/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl#L991)
  through
  [PureGaussianGausslet.jl](/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl#L1009)
- only after the contraction are they divided by the new weights at
  [PureGaussianGausslet.jl](/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl#L1079)
  through
  [PureGaussianGausslet.jl](/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl#L1081)

In the older nested decimation line, the same rule appears even more directly:

- multiply by old weights first at
  [facedecimate3steve.jl](/Users/srw/Dropbox/GaussletModules/facedecimate3steve.jl#L346)
  and
  [fd4.jl](/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/fd4.jl#L463)
- transform the raw numerator
- transform the weights
- divide by the new weights afterward at
  [facedecimate3steve.jl](/Users/srw/Dropbox/GaussletModules/facedecimate3steve.jl#L348)
  through
  [facedecimate3steve.jl](/Users/srw/Dropbox/GaussletModules/facedecimate3steve.jl#L351)
  and
  [fd4.jl](/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/fd4.jl#L465)
  through
  [fd4.jl](/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/fd4.jl#L468)

So the legacy rule is:

```math
M = D(w)\,V\,D(w), \qquad
M' = C^T M C, \qquad
w' = C^T w, \qquad
V' = D(1/w')\,M'\,D(1/w').
```

## Implemented Fix

The nested shell packet in
[cartesian_nested_faces.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/cartesian_nested_faces.jl)
now follows that rule.

For the fixed-block interaction packet:

1. reconstruct the raw one-dimensional pair numerator from the stored PGDG
   `pair_factor_terms` and parent integral weights
2. build the three-dimensional raw support numerator on the shell support
3. contract that raw numerator through the shell coefficient matrix
4. contract the parent three-dimensional integral weights through the same map
5. divide the contracted numerator by the new fixed-block weight outer product

The nested packet and fixed-block adapter now also carry the contracted
fixed-block integral weights explicitly.

## Weight Diagnostics

On the stabilized He count-17 fixed-`a` case:

- shell-plus-core contracted weights:
  - min `0.061468156103430815`
  - max `2.977482219200775`
  - nonpositive entries: `0`
- corrected complete four-shell contracted weights:
  - min `0.061468156103430815`
  - max `3.2591670283995673`
  - nonpositive entries: `0`

So the intended nested fixed-block weights are finite and strictly positive in
this regime.

## Interaction-transfer Comparison

Same stabilized He fixed-`a` count-17 case, fixed-only projected parent-ground
diagnostic:

- parent fixed-space `⟨Vee⟩`: `1.2486322227019246`
- shell-plus-core:
  - old direct-IDA contraction: `1.2486612176560734`
  - corrected weight-aware transfer: `1.248633734959628`
- corrected complete four-shell sequence:
  - old direct-IDA contraction: `1.818451883926288`
  - corrected weight-aware transfer: `1.2487650946681053`

So the old nested path really was transforming the wrong IDA object.

## End-to-end He nearest/GGT Result After the Fix

Same final `1s^2` nearest/GGT check:

- shell-plus-core:
  - `E1 = -1.9985629071304167`
  - `⟨Vee⟩ = 1.2489021100178215`
- corrected complete four-shell sequence:
  - `E1 = -1.9981842264017424`
  - `⟨Vee⟩ = 1.2489346080730335`

The large previous discrepancy is gone.

## Conclusion

Yes: the current nested fixed-block path had been transforming the wrong IDA
object.

The legacy code transformed the raw numerator and only then reweighted by the
new contracted basis integrals.

With that correction in place, the bad corrected-complete-sequence
`⟨Vee⟩` discrepancy closes almost entirely on the stabilized He test. The
remaining difference from shell-plus-core is small, on the order of
`10^{-5}` to `10^{-4}`, not the earlier `O(10^{-1})` to `O(1)` failure.
