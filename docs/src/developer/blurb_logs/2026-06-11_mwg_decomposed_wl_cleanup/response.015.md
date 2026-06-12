Implemented the narrow active-route factorized mixed gausslet/GTO block path for shellification-backed decomposed WL inventories.

Changed file:

- `src/cartesian_pair_block_materialization/route_global_mixed_gto_blocks.jl`

What changed:

- Added a private factorized retained-basis projection path inside `route_global_mixed_gto_blocks_from_decomposed_units`.
- The path consumes the existing decomposed WL factorized retained-basis sidecar derived from retained-unit coefficient maps.
- It builds reusable parent-axis/GTO cross tables once through the old QW cross-table kernel shape, then projects them into the retained gausslet rows.
- It materializes mixed overlap, kinetic, position, x2, and uncharged electron-nuclear-by-center mixed blocks without looping over retained units and without constructing per-unit CPB/GTO provider local blocks.
- It reuses the CPB provider-owned `_cpb_gto_internal_supplement` conversion so the factorized path uses the same legacy `_AtomicCartesianShellOrbital3D` source shape as the established provider path.
- It is currently enabled only for the active one-center origin atomic `UnitPairIndexTable` path. Unsupported cases fall back to the existing per-unit provider path.

Numerical validation:

- Be S+P final-basis route completed.
- `rhf_total_energy = -14.574514244574639`
- old nested/QW oracle total = `-14.574514244574694`
- delta from old oracle = `5.5067062021407764e-14 Ha`
- final dimension = `636`
- retained gausslet dimension = `615`
- retained supplement count = `21`
- final overlap rank = `636`
- final overlap identity error = `1.4233342578506758e-10`
- no full-parent CPB, direct Cartesian product assembly, ordinary Cartesian IDA, raw GTO density-density final operator, or generalized final solve was used.

Timing result:

- Previous route elapsed: `339.913994709s`
- New route elapsed: `171.710878583s`
- Previous total elapsed: `342.53537725s`
- New total elapsed: `174.253094s`
- Previous `mixed_gto_blocks`: `177.196382542s`
- New `mixed_gto_blocks`: `9.177082084s`
- Previous mixed subphase `per_unit_provider_local_block_construction`: `168.066350504s` across 131 units
- New mixed subphase `factorized_projection_total`: `0.338807417s` once
- New mixed subphase `gto_gto_self_block_construction`: `0.777837292s` once

Updated artifacts:

- `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
- `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`
- `tmp/work/be_atom_sp_mixed_gto_subphase_timings.tsv`

The Be probe confirmed the BasisSets source path as:

- `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets`

Validation run:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

No separate test was added. The Be S+P probe is the direct live-route validation for this blurb, and adding another metadata-style test would duplicate the acceptance probe without protecting a new contract.

Deletion/shrinkage report:

- The active Be S+P mixed-GTO hot path no longer uses the 131-unit per-unit provider-local block materialization loop.
- That loop remains as compatibility/reference fallback for unsupported shapes, including non-`UnitPairIndexTable` inventories, multi-center cases, and off-origin centers.
- No old source file or test became fully dead in this pass because the fallback path still has live compatibility value outside the one-center atomic factorized route.
- No new test replaces older coverage; no new test was added.
- Remaining stale/duplicate surface to retire next: once factorized mixed-GTO support is generalized beyond one-center origin atomic cases, the per-unit CPB/GTO local provider route can be narrowed further to explicit diagnostic/reference status.

-- repo-doer@macmini
