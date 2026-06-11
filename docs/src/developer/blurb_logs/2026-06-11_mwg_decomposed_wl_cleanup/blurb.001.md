Purpose:

Implement route-global combined GTO residual moment matrix assembly for the
final-basis MWG representation.

Why now:

The current committed blocker for real side13 He + GTO is:

```text
:missing_route_global_combined_gto_residual_moment_matrices
```

Commit `81721733` added
`route_global_residual_gto_mwg_representation(...)`, but it is blocked unless
raw combined final-basis moment matrices are supplied. Synthetic tests show
centers/widths are moment-derived from final-basis residual vectors, and raw
GTO density-density is still not accepted as final operator data.

Current state:

- `c31dd405` added the final-basis density-density readiness blocker.
- `81721733` added `route_global_residual_gto_mwg_representation(...)`.
- The MWG representation is blocked unless raw combined final-basis moment
  matrices are supplied.
- Raw GTO density-density is still not accepted as final operator data.

Exact task:

Build the real route-global combined GTO moment matrices needed by
`route_global_residual_gto_mwg_representation(...)`.

Required raw combined moment matrices:

- `position_x`
- `position_y`
- `position_z`
- `x2_x`
- `x2_y`
- `x2_z`

These should be assembled in the same raw combined gausslet+GTO basis used by
the existing combined one-electron path:

```text
[ gausslet/gausslet   gausslet/GTO
  GTO/gausslet        GTO/GTO      ]
```

Existing reusable pieces:

- mixed gausslet/GTO position/x2 CPB-local blocks already exist;
- GTO/GTO position/x2 blocks already exist;
- decomposed WL gausslet/gausslet route-global position/x2 may already exist or
  may need a narrow wrapper;
- combined layout machinery already exists;
- final-basis projection already exists.

Code surfaces to inspect first:

- `src/cartesian_pair_block_materialization/route_global_combined_gto_matrix_assembly.jl`
- `src/cartesian_pair_block_materialization/route_global_mixed_gto_blocks.jl`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_final_basis.jl`
- `src/cartesian_pair_block_materialization/one_body_global_position.jl`
- `src/cartesian_pair_block_materialization/one_body_global_x2.jl`
- `src/CartesianCPBBlockProviders.jl`

Named surfaces are starting points, not authority. If a name is stale or
inconsistent with the live repo, stop and report.

Implementation target:

Add a narrow helper, name flexible, conceptually like:

```julia
route_global_combined_gto_residual_moment_matrices(...)
```

It should return the six raw combined moment matrices in a compact object
suitable for passing to:

```julia
route_global_residual_gto_mwg_representation(
    final_basis_projection,
    combined_matrices,
    supplement;
    raw_moment_matrices = ...
)
```

Use side13 He + H cc-pVTZ `lmax = 0` as the first real fixture if a fixture is
needed.

Trust boundary:

Do not build electron-electron `V_gR` or `V_RR` yet unless the moment assembly
makes that trivial and clearly correct.

Do not run He + GTO RHF yet.

Do not use raw GTO density-density as final e-e.

Do not use generalized-overlap final solve.

Do not add a long-running side17 gate.

Do not add:

- He + GTO RHF acceptance
- raw GTO/GTO density-density as accepted final operator
- generalized-overlap final solve
- full-parent CPB fallback
- old fixed-block matrix authority
- flat `_qwrg_diatomic_*` fallback
- ordinary Cartesian IDA fallback
- PQS changes
- broad metadata tests

Validation:

- focused test/probe only;
- verify all six moment matrices are materialized;
- verify raw combined shapes match `total_combined_dimension`;
- verify `route_global_residual_gto_mwg_representation(...)` materializes for
  a real side13 He + GTO final-basis projection;
- verify residual centers/widths are finite and widths are positive;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- exact moment assembly surfaces used;
- block shapes;
- real side13 residual centers/widths summary;
- whether the MWG representation is now ready for `V_gR` / `V_RR`;
- validation run;
- deletion/shrinkage report.
