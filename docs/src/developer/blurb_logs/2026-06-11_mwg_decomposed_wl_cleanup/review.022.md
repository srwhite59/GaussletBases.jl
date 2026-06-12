Accepted.

The pass followed the residual-moment pattern for the gausslet one-electron
part of the active atom+GTO route. The route now uses a narrow
`DecomposedWLOneElectronMatrixSet` instead of staging separate overlap,
kinetic, and electron-nuclear route reports before combined GTO assembly.

The important timing result:

```text
cold route elapsed    29.012344166s -> 25.465580292s
warm route elapsed     0.575187292s ->  0.495644334s

old cold overlap + kinetic + nuclear sum    7.672966416s
new cold one-electron matrix-set phase      4.432942625s
```

The by-center nuclear convention remains intact. The new by-center payload
still records center index/key/location/charge, keeps matrices separated by
center, and marks `nuclear_charge_applied = false` and `centers_summed = false`.
Combined GTO assembly still applies charges and sums centers at Hamiltonian
assembly.

Physics stayed pinned to the old nested/QW oracle:

```text
Be S+P RHF total        -14.574514244574639
old nested/QW oracle    -14.574514244574694
delta                    5.51e-14 Ha
final dimension          636
retained gausslets       615
units / pairs            131 / 8646
```

Fallback/final-basis flags remain clean: no full-parent CPB, direct Cartesian
fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or
generalized final-basis solve.

Manager validation:

- `julia --project=. test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`
- `julia --project=. test/nested/cartesian_wl_factorized_backend_equivalence_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

One caution: the private atom+GTO route result still stores the one-electron
matrix set in the legacy overlap/kinetic/nuclear payload slots. That is
acceptable as a temporary compatibility shape, but it is now the most visible
report/compute mixing to clean up if it starts showing up in specialization
audits.

Next target:

Do not immediately add another broad struct layer. First audit the new
`decomposed_wl_one_electron_matrix_set` phase itself. It is now the largest
single cold route phase at about 4.43s cold and 0.054s warm. The implementation
still calls internal report-shaped factorized helpers for overlap, kinetic, and
electron-nuclear. If the audit confirms those helper result shapes are the
remaining cost, replace only the inner hot-path construction with direct shared
factorized matrix fills while keeping existing wrapper helpers as
compatibility/reference surfaces.

-- repo-manager@macmini
