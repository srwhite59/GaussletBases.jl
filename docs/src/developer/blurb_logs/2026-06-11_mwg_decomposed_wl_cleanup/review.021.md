Accepted.

The pass did the narrow thing requested: it audited the residual-moment phase
and removed the six separate decomposed WL moment route calls from the active
atom+GTO hot path. The new `DecomposedWLMomentMatrixSet` is a real compute
concept for this seam: six retained gausslet moment matrices needed by residual
MWG construction. It is not a giant route report replacement.

The important timing result is phase-local and clear:

```text
residual_moment_matrices cold   6.912517084s -> 1.074904167s
residual_moment_matrices warm   0.00290625s  -> 0.002948792s
```

Top-level cold route time improved more modestly:

```text
cold route   30.843591209s -> 29.012344166s
warm route    0.492626583s ->  0.575187292s
```

The warm route increase is not in the new residual-moment phase, which remains
about 0.003s. Treat it as probe noise or unrelated warm-phase variation unless
it repeats in later measurements.

Physics stayed pinned to the old nested/QW oracle:

```text
Be S+P RHF total        -14.574514244574639
old nested/QW oracle    -14.574514244574694
delta                    5.51e-14 Ha
final dimension          636
retained gausslets       615
units / pairs            131 / 8646
```

Fallback and final-basis contract flags remain clean: no full-parent CPB,
direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final
density-density, or generalized final-basis solve.

Manager validation:

- `julia --project=. test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`
- `julia --project=. test/nested/cartesian_wl_factorized_backend_equivalence_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

The same shape is still visible earlier in the atom+GTO route. The active route
builds overlap, kinetic, and electron-nuclear by-center as three separate
decomposed WL one-body route reports. The current cold phase table shows those
three phases now cost about 7.67s combined:

```text
electron_nuclear_by_center   4.311094166s
overlap                      1.940648875s
kinetic                      1.421223375s
```

The next pass should audit whether a narrow decomposed WL one-electron matrix
set can do for these phases what the moment-matrix set did for residual
moments: carry matrices through the hot atom+GTO path, leaving report-shaped
wrappers as compatibility/reference surfaces.

-- repo-manager@macmini
