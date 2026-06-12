Accepted as a full-route attribution pass.

The controlled variants answer the immediate question. Timing sinks are not the
problem, and minimal metadata does not help:

```text
current timing, full density       24.095634292s cold route
timing disabled, full density      24.585734292s cold route
minimal metadata, full density     28.330674s cold route
```

The large split is density-density enabled versus one-electron only:

```text
current full density     24.095634292s cold route
one-electron only        10.57373825s cold route
delta                    about 13.5s
```

So parent-axis metadata, timing wrappers, and the one-electron helper internals
should not be the next source edit. The next target is the final
density-density / residual MWG route shape.

The route result type remains very large:

```text
full-density type-name length       about 44k
one-electron-only type-name length  about 25k
route property count                46
```

That supports a future lean compute/result split, but the current measurement
points first to the density-enabled portion of the route. The next pass should
attribute this region, not guess at a broad route-result refactor.

Relevant source seams:

- `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
  - `:gausslet_density_density`
  - `:final_basis_density_density`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_density_density.jl`
  - `route_global_combined_gto_final_basis_density_density_matrix(...)`
  - `_route_global_combined_gto_residual_mwg_density_density_components(...)`
  - `_route_global_combined_gto_residual_mwg_component_result(...)`
  - `_route_global_combined_gto_density_density_result(...)`

No source or test changes were made. No tests were added. This pass correctly
used ignored `tmp/work` probes.

-- repo-manager@macmini
