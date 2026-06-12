Review result:

Accepted as a useful narrow optimization and timing-attribution pass, with the
important conclusion that the easy GTO/GTO self-block hoist was not the main
runtime problem.

Main result:

- Be S+P final-basis RHF still converges through the decomposed/final-basis
  route.
- final dimension: `636`
- final overlap identity error: `1.0125523569644675e-10`
- final density-density symmetry error: `0.0`
- RHF iterations: `27`
- new RHF one-electron energy: `-19.066200470580668`
- new RHF electron-electron energy: `4.4916862260060055`
- new RHF total: `-14.574514244574662`
- old nested/QW oracle total: `-14.574514244574694`
- RHF total delta from old oracle: `3.197442310920451e-14`

Timing result:

The full probe moved from about `357.4` seconds to about `342.5` seconds. The
route elapsed time moved from about `354.8` seconds to about `339.9` seconds.
The `mixed_gto_blocks` phase moved from about `188.6` seconds to about `177.2`
seconds.

The new mixed-GTO subphase timing is the key result:

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

Interpretation:

The pass correctly removed repeated GTO/GTO self-block construction from the
active mixed-GTO route, but that was a secondary cost. The measured bottleneck
is now precise: per-unit mixed CPB/GTO local block construction accounts for
about `168` seconds over `131` retained units. The retained contraction and row
placement machinery are not the problem at this fixture size.

Deletion/shrinkage review:

The active route no longer uses the old per-unit full
`cpb_gto_supplement_local_operator_bundle(...)` path for route-global mixed-GTO
assembly. Provider-level bundle construction remains available for provider
tests and direct provider consumers. The route-local helpers were simplified
from "build mixed plus GTO/GTO self every unit" to "build GTO/GTO self once,
then build only mixed blocks per unit." No tests were added.

Remaining stale/duplicate surface:

Per-unit mixed CPB/GTO local block materialization is now measured as the
production bottleneck for larger atom plus supplement fixtures. The next pass
should make that surface less necessary by introducing a factorized/projected
mixed-GTO route for the active shellification-backed path, or stop with the
smallest precise design blocker if that replacement is too broad.

Validation reviewed:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
- artifact: `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`
- artifact: `tmp/work/be_atom_sp_mixed_gto_subphase_timings.tsv`

Next target:

Replace the per-unit mixed CPB/GTO local block construction in the active
shellification-backed mixed-GTO route with a factorized or projected route that
reuses parent-axis/supplement cross tables. Preserve the current Be S+P oracle
agreement and keep old provider/local block paths as reference or compatibility
surfaces unless they become clearly dead.
