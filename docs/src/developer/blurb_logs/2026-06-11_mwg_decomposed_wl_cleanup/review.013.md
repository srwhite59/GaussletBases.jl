Review result:

Accepted as the first complete Be S+P decomposed/final-basis RHF diagnostic.
The route now materializes final-basis density-density, runs closed-shell RHF
in the final orthonormal basis, and matches the old nested/QW oracle to
roundoff.

Main result:

- status: `:materialized_final_basis_be_sp_rhf_probe`
- route status:
  `:materialized_decomposed_atom_gto_final_basis_density_density_route`
- retained gausslet / supplement / final dimensions: `615 / 21 / 636`
- final overlap identity error: `1.0125523569644675e-10`
- final density-density symmetry error: `0.0`
- RHF converged: true
- RHF iterations: `27`
- one-electron energy: `-19.066200470580668`
- electron-electron energy: `4.4916862260060055`
- RHF total: `-14.574514244574662`
- old nested/QW oracle total: `-14.574514244574694`
- delta from old oracle: `3.197442310920451e-14`

Phase timing:

The pass added private opt-in phase timing to the Be seam. The total elapsed
time was about `357.4` seconds. The dominant phase is clearly
`mixed_gto_blocks` at about `188.6` seconds. Other notable phases are residual
moments at about `36.0` seconds, electron-nuclear at about `28.2` seconds,
overlap at about `22.3` seconds, and gausslet density-density at about `21.3`
seconds. RHF itself took about `2.6` seconds.

Validation reviewed:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
- artifact: `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`

Deletion/shrinkage review:

No H/H2+ fixture-local GTO wiring became removable yet. The phase-attributed
run identifies the next actual replacement/optimization target:
`route_global_mixed_gto_blocks_from_decomposed_units(...)`, especially repeated
mixed-GTO block materialization across 131 retained units.

Commit/push:

Pending manager commit/push of the private timing hook, numerical-contract
update, response/review log, and the next published blurb.

Next target:

Audit and optimize mixed GTO route-global block materialization. Start with
subphase timing inside `mixed_gto_blocks`; then remove repeated work where the
source is clear, especially repeated GTO/GTO self-block construction across
retained units.
