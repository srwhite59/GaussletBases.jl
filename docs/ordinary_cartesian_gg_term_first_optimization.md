# Ordinary Cartesian GG Term-First Optimization Note

This note records the first real optimization pass after the Qiu-White timing
diagnosis.

The main bottleneck is no longer the old mixed raw-layer mistake in the
Qiu-White constructor. It is now the shared mapped ordinary gausslet-side
one-dimensional bundle, especially the Gaussian-factor and pair-factor path.

That work is generic `gg` infrastructure, not Qiu-White-specific logic.

So the next optimization target in this pass is:

- the mapped-gausslet `gg` Gaussian-factor path
- especially the mapped-gausslet `gg` pair-factor path

The structural model is the legacy Coulomb-expansion assembly pattern used in
`PureGaussianGausslet.jl`:

- keep the Coulomb-expansion index as the short inner reduction
- prepare one-dimensional factor data once in term-first form
- assemble three-dimensional entries directly from those one-dimensional term
  arrays

The target access pattern is therefore closer to:

```julia
sum(coulco .* Fx[:,ix,jx] .* Fy[:,iy,jy] .* Fz[:,iz,jz])
```

than to:

- build one full dense 3D matrix per Coulomb term
- then sum those dense term matrices afterward

The intended result of this pass is:

- shared mapped ordinary `gg` factor data available in term-first form
- ordinary Cartesian IDA and the Qiu-White reference path both consuming that
  shared representation
- no new broad timing framework
- no new Qiu-White-specific surgery before the generic `gg` path is improved

## First timing result

On the same light paper-like He case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ`
- `interaction_treatment = :mwg`

the direct basis-space `gg` scalar/factor path did improve materially:

- `basis_representation(...)`: about `15.43 s`
- new direct basis-space `x^2` plus Gaussian-factor build: about `4.90 s`

That replaces the older separate costs of roughly:

- `basis_representation(...)`: about `15.12 s`
- `_x2_matrix(basis)`: about `18.97 s`
- `gaussian_factor_matrices(...)`: about `5.08 s`

So the term-first direct-basis treatment is already a real win for the shared
`gg` scalar/factor side.

But the dominant remaining cost is still the shared `gg` pair-factor path.
The new direct basis-space pair-factor builder still did not finish within
about `70 s` on that same case, and the full Qiu-White constructor still did
not clear phase 1 in a reasonable time.

So the next bottleneck after this pass is now even clearer:

- not the later dense 3D assembly
- not the Qiu-White-specific residual-space logic
- the shared mapped-gausslet `gg` pair-factor builder still dominates
