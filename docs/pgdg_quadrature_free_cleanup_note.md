The next cleanup target after the earlier bottleneck isolation was to remove any remaining numerical quadrature from the active PGDG-mediated path.

That cleanup is now the right conceptual boundary:

- if a path is using the shared PGDG intermediate layer, its core 1D data should come from analytic Gaussian primitive integrals plus contraction
- the PGDG 1D auxiliary line should then receive its final COMX cleanup inside the bundle, so the carried overlap block is orthonormal again
- any remaining midpoint or numerical quadrature should belong only to demoted diagnostic helpers or to separate pure distorted-gausslet reference work, not to the active QW-PGDG path

The legacy model for this cleanup is the same one used in the older integral modules:

- analytic 1D Gaussian / Gaussian-factor / pair-factor primitives from `GTO1ds.jl` and `Gaucoulomb.jl`
- contraction / assembly patterns like the ones used in `PureGaussianGausslet.jl`

So the active PGDG-mediated route should now be read as:

1. analytic primitive Gaussian data
2. contraction into the PGDG auxiliary line
3. COMX cleanup/localization of that auxiliary line inside the shared bundle
4. downstream Cartesian assembly from that cleaned shared layer

That makes the active path conceptually distinct from a future pure Qiu-White distorted-gausslet reference route:

- active path now: QW-PGDG / PGDG-mediated
- later possible path: pure QW, if a true distorted-gausslet reference implementation is still wanted

The immediate validation target remains the He-like `1s^2` scalar check:

- `⟨Vee⟩ = 5/4 = 1.25`

and the timing question becomes narrower:

- after the PGDG cleanup, any remaining cost should be attributed to the clean analytic/contraction structure itself, not to leftover midpoint quadrature inside the PGDG-dependent path.
