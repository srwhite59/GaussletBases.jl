Purpose:

Attribute the cold compile pressure inside the density-enabled part of the Be
atom+GTO route before changing source.

Why now:

The full-route attribution probe showed the decisive split:

```text
current full density cold route     24.095634292s
one-electron-only cold route        10.57373825s
density-enabled delta               about 13.5s
```

Timing sinks and minimal metadata did not explain the remaining cost. The next
target is the density-density / residual MWG portion of the route, not parent
metadata or timing wrappers.

Relevant source seams:

```text
src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl
    :gausslet_density_density
    :final_basis_density_density

src/cartesian_pair_block_materialization/route_global_combined_gto_density_density.jl
    route_global_combined_gto_final_basis_density_density_matrix(...)
    _route_global_combined_gto_residual_mwg_density_density_components(...)
    _route_global_combined_gto_residual_mwg_component_result(...)
    _route_global_combined_gto_density_density_result(...)
```

Exact task:

Create or update an ignored `tmp/work` developer probe that isolates the
density-enabled route phases. Do not change production source unless the probe
identifies one tiny, well-supported cleanup.

Measure at least:

1. `gausslet_density_density` alone:
   - route call time cold/warm;
   - result type-name length and property count;
   - whether the matrix and key metadata are enough for downstream use.

2. `residual_mwg_representation` alone if practical:
   - cold/warm time;
   - result type-name length and property count;
   - residual centers/widths shape.

3. `final_basis_density_density` alone with existing inputs:
   - full call cold/warm;
   - time inside `_route_global_combined_gto_residual_mwg_density_density_components(...)`
     if practical;
   - time inside `ParentGaussletBases._qwrg_final_residual_mwg_component_blocks(...)`
     if practical;
   - time in `_route_global_combined_gto_residual_mwg_component_result(...)`
     and `_route_global_combined_gto_density_density_result(...)` if practical.

4. Combined density path variants:
   - gausslet density only;
   - gausslet density + residual MWG representation;
   - full final density-density.

5. Report the concrete result shape:
   - type-name length;
   - property count;
   - whether full upstream route objects are carried when only matrices,
     centers, widths, and the final-basis projection are needed.

Interpretation required:

- If `_qwrg_final_residual_mwg_component_blocks(...)` dominates cold time, rank
  it as an old-QW/MWG analytic-kernel precompile or extraction target.
- If result construction/report shape dominates, propose the smallest compute
  object split, for example a final density interaction compute payload plus
  compact audit summary.
- If `gausslet_density_density` dominates, inspect whether it is carrying a
  broad route report where downstream only needs the matrix and density-weight
  convention.
- If no single subphase dominates, rank the remaining compile-pressure sources
  with evidence.

Trust boundary:

Measurement first. No public API/export/default changes. No PQS, ECP, Be2, H2,
high-l Be, driver defaults, acceptance fixtures, raw GTO final density-density,
generalized final solve, full-parent CPB fallback, direct Cartesian fallback,
or ordinary Cartesian IDA fallback.

Test policy:

Do not add tests. Use ignored `tmp/work` probes. If no production source
changes are made, no Julia test is required beyond a load check.

Validation:

- run the new/updated density-route attribution probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if production source changes are made, run the smallest directly affected
  focused test and the Be warm/cold probe.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed, if any;
- probe artifact paths;
- density subphase timing table;
- result shape/type-size table;
- interpretation and ranked next target;
- whether any source change was made;
- validation run;
- deletion/shrinkage report.

-- repo-manager@macmini
