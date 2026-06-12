Purpose:

Replace the measured Be S+P mixed-GTO bottleneck with a factorized/projected
mixed-GTO route, or stop with the smallest precise design blocker.

Why now:

The Be S+P decomposed/final-basis RHF route is scientifically matched to the
old nested/QW oracle:

- new RHF total: `-14.574514244574662`
- old nested/QW oracle total: `-14.574514244574694`
- delta: `3.2e-14 Ha`

The previous pass hoisted reusable GTO/GTO self blocks out of the retained-unit
loop. That was correct, but the timing shows it was not the real bottleneck.
The current mixed-GTO subphase table says:

```text
phase                                      elapsed_s          count
gto_gto_self_block_construction            0.823398625        1
per_unit_total                             169.936329462      131
unit_coefficient_construction              0.010550254        131
per_unit_provider_local_block_construction 168.066350504      131
support_coefficient_construction           0.000976289        131
retained_contraction                       0.00206471         131
row_placement_coverage                     0.000661169        131
```

So the next target is not GTO/GTO self blocks, retained contraction, or
placement. It is the repeated per-unit mixed CPB/GTO local block construction
inside:

```julia
route_global_mixed_gto_blocks_from_decomposed_units(...)
```

Known code surfaces:

- active mixed-GTO route:
  `src/cartesian_pair_block_materialization/route_global_mixed_gto_blocks.jl`
- private one-center atom plus supplement seam:
  `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
- provider-level mixed CPB/GTO blocks:
  `src/CartesianCPBBlockProviders.jl`
- factorized decomposed WL operator backend and sidecar:
  search for the factorized retained-basis backend used by the He/Be
  decomposed WL operator route
- old nested/QW fixed-block route:
  use only as computational-shape oracle and timing comparator, not as route
  authority

Exact task:

1. Inspect the existing factorized decomposed-WL retained-basis sidecar and the
   old nested/QW mixed-supplement/cross-table code shape. Do not assume the
   current per-unit CPB-local mixed block construction is the final production
   algorithm.

2. Design the narrowest active-route replacement for mixed GTO blocks in the
   shellification-backed atom plus supplement path. The intended shape is:
   - reusable parent-axis / supplement cross tables;
   - projection through the retained factorized WL basis;
   - direct route-global mixed rows for overlap, kinetic, position/x2, and
     electron-nuclear by-center;
   - no per-retained-unit provider local block materialization in the hot path.

3. If a clean implementation fits in this pass, implement it privately in CPBM
   for the active `UnitPairIndexTable` / shellification-backed path. Preserve
   the current provider/local block route as reference or compatibility unless
   it becomes clearly dead.

4. If the full replacement is too broad, implement a measured prototype for one
   representative term, preferably mixed overlap or mixed electron-nuclear, and
   stop with:
   - exact numerical comparison against the existing mixed route;
   - measured timing;
   - the smallest remaining design blocker for completing all mixed terms.

5. Preserve all physics and final-basis contracts:
   - final dimension `636`;
   - RHF total remains within tight tolerance of the old nested/QW oracle;
   - final overlap remains effectively identity;
   - no fallback flag changes from false to true;
   - raw GTO density-density is still not accepted as final operator data;
   - no generalized final-basis solve.

Trust boundary:

Keep this private/internal. Do not add public APIs, exports, route defaults,
PQS, ECP, high-l Be, Be2, Cr, H2, full-parent CPB fallback, direct Cartesian
fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or
generalized final-basis solves.

Do not add tests by default. Use `tmp/work` probes for exploratory comparison
and timing. Add a tracked test only if it protects a new live contract that is
not already covered by the Be probe or a compact module-contract test, and say
which older coverage it replaces or why it earns carrying cost.

Decision rules:

- If the factorized/projected path reproduces the old mixed-route matrices and
  materially reduces `per_unit_provider_local_block_construction`, use it for
  the active shellification-backed mixed-GTO route.
- If only one term is replaced, keep it behind the active private route only
  when the Be probe still matches the old oracle; otherwise leave it as a
  `tmp/work` prototype and report.
- If the replacement would become a broad framework, stop and report the
  smaller design with exact files/functions and timing evidence.
- After the replacement is clean, name the first old route/local provider
  surface that becomes less necessary and what would have to be true before it
  can be deleted or quarantined.

Artifacts:

- update the Be S+P probe artifact:
  `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
- update phase timing:
  `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`
- update mixed-GTO subphase timing:
  `tmp/work/be_atom_sp_mixed_gto_subphase_timings.tsv`
- add a temporary prototype artifact only if needed under `tmp/work/`

Validation:

- run the Be S+P probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- run a cheap focused existing test only if it directly covers the edited
  mixed-GTO surface; otherwise explain why the Be probe is the direct
  validation.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- whether this is a full active-route replacement or a measured prototype;
- numerical comparison to the previous mixed route and old Be oracle;
- before/after timings for total Be probe, `mixed_gto_blocks`, and the dominant
  mixed-GTO subphase;
- whether per-unit provider local block materialization remains in the hot
  path;
- validation run;
- deletion/shrinkage report.
