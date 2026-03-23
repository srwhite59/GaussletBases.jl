# Radial Trust Milestone Note

This note records a short trust-establishment sweep for the current radial
line after three coupled changes:

- runtime family tables trimmed to machine-significant tails
- public `rmax` separated from internal `build_umax` and `quadrature_umax`
- interval-sampled setup-grid construction replacing the old dense setup path

The point of this pass was not another optimization. It was to check that the
faster construction path still reproduces the expected radial scientific
behavior on representative cases.

The sweep script is:

- `tmp/work/radial_trust_sweep.jl`

## Representative checks

### Standard hydrogen front door

Recipe:

- `:G10`
- `Z = 1`
- `s = 0.2`
- `c = s / (2Z)`
- `rmax = 30`
- default two-`xgaussian` front-door supplement

Result:

- basis length `35`
- `diag.overlap_error = 8.797053967946056e-6`
- `diag.D = 1.3554300216193975e-7`
- hydrogen ground-state energy `E0 = -0.4999999972288096`
- build time about `6.03 s`

So the front-door hydrogen check remains in the expected scientific regime.

### `(l,m)` hydrogen one-body path

Recipe:

- `:G10`
- `Z = 1`
- `s = 0.2`
- `rmax = 30`
- `lmax = 2`
- no extra `xgaussians` to match the existing `(l,m)` reference fixture

Result:

- basis length `33`
- overlap error `4.0596246031296215e-7`
- lowest energy `-0.4999999945176998`
- lowest `l = 1` energy `-0.12499999988503978`
- lowest `l = 2` energy `-0.05555046143346199`
- build time about `0.39 s`

This keeps the expected hydrogenic `-1/(2n^2)` pattern through the existing
angular decomposition route.

### Smaller nondefault radial recipe

Recipe:

- `:G10`
- `Z = 1`
- `s = 0.3`
- `c = s / (2Z)`
- `rmax = 12`
- no `xgaussians`

Result:

- basis length `21`
- `diag.overlap_error = 4.0596245811904994e-7`
- `diag.D = 2.4278240367355398e-4`
- hydrogen ground-state energy `E0 = -0.49999997523710266`
- build time about `0.27 s`

So the new setup path also behaves sensibly away from the standard public
recipe.

## Runtime-table trim check

For `:G10`:

- runtime radius `75`
- preserved high-precision internal radius `132`
- max absolute value difference on a dense sample over `[-8, 8]`:
  `4.440892098500626e-16`

So the machine-significant runtime tables are materially shorter without
changing representative basis values at a scientifically meaningful level.

## Current standard timing

From the cache-study timing on the recommended two-`xgaussian`, `K = 6`
recipe:

- `build_basis(spec)`: about `6.5 s`
- `basis_diagnostics(rb)`: about `0.6 s`
- `radial_quadrature(rb)`: about `0.44 s`
- `atomic_operators(rb, grid; Z = 2, lmax = 2)`: about `1.6 s`
- full front-door sequence: about `9.4 s`

That is a large improvement over the earlier `~173 s` build path.

## Practical conclusion

The fast interval-sampled build path is now trusted for the radial line.

The present radial branch is no longer primarily a live stabilization front.
Its scientific checks are in the right regime, the public extent story is now
coherent, and the standard public path is fast enough that the next work can
move on to comparison/cleanup priorities rather than emergency radial
stabilization.
