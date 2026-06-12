What changed:

- Added opt-in phase timing to the private seam
  `_white_lindsey_decomposed_atom_gto_final_basis_route(...)`.
- The timing hook is private/internal and records only existing phase calls into
  a caller-supplied vector. It does not add public API, exports, driver
  defaults, global timing state, fallback routes, or route behavior.
- Updated `tmp/work/be_atom_sp_decomposed_final_basis_probe.jl` to run Be S+P
  with `build_density_density = true`, pass a phase timing sink, write a timing
  TSV, and attempt closed-shell RHF only after final-basis density-density
  materializes.
- Updated artifacts:
  - `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
  - `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`

Timing attribution:

- Timing was added to the private production seam as an opt-in sink, then used
  by the probe.
- Carrying cost: one small private timing helper plus phase wrappers around
  already-existing seam calls. No broad timing framework was added.

Be S+P density-density/RHF result:

- Route status:
  `:materialized_decomposed_atom_gto_final_basis_density_density_route`.
- Final probe status:
  `:materialized_final_basis_be_sp_rhf_probe`.
- Blocker: `nothing`.
- Fixture:
  Be, all-electron, `Z = 4`, `q/ns = 5/5`, `lmax = 1`, `d = 0.15`,
  parent side count `15`, endpoints
  `(-10.056710949734484, 10.056710949734484)`.
- Supplement:
  `legacy_atomic_gaussian_supplement("Be", "cc-pV5Z"; lmax = 1,
  basisfile = "/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets")`.
- Shellification-backed decomposed WL inventory: `true`.
- Low-order seed inventory used: `false`.
- Fallback flags:
  full-parent CPB `false`, direct Cartesian product assembly `false`,
  ordinary Cartesian IDA operators `false`, raw GTO final density-density
  `false`, generalized final solve `false`.
- Dimensions:
  retained gausslet `615`, units `131`, pairs `8646`, raw supplement `21`,
  retained supplement `21`, dropped supplement `0`, final dimension `636`.
- Old oracle dimensions:
  fixed/residual/final `615 / 21 / 636`; final dimension delta `0`.
- Final-basis diagnostics:
  overlap identity error `1.0125523569644675e-10`, overlap rank `636`,
  overlap condition estimate `1.0000000002807463`, Hamiltonian symmetry error
  `1.2108358760087867e-10`, final density-density symmetry error `0.0`.

RHF comparison to old oracle:

- New route RHF converged: `true`.
- Iterations: `27`.
- New one-electron energy: `-19.066200470580668`.
- Old oracle one-electron energy: `-19.06620047058102`.
- New electron-electron energy: `4.4916862260060055`.
- Old oracle electron-electron energy: `4.491686226006327`.
- New RHF total: `-14.574514244574662`.
- Old oracle RHF total: `-14.574514244574694`.
- RHF total delta from old oracle: `3.197442310920451e-14`.
- Density trace: `2.000000000000001`.
- Electron count: `4.000000000000002`.

Phase timing table:

```text
phase                                      elapsed_s        status
parent_axis_setup                         0.059548625      completed
shellification_decomposed_inventory       2.552077583      available_white_lindsey_decomposed_unit_pair_inventory
combined_gto_basis_layout                 1.668265625      available_route_global_combined_gto_basis_layout
overlap                                   22.281612959     materialized_route_global_overlap_matrix
kinetic                                   13.04325925      materialized_route_global_kinetic_matrix
electron_nuclear_by_center                28.230689042     materialized_route_global_electron_nuclear_by_center_matrix_set
mixed_gto_blocks                          188.576837208    materialized_route_global_mixed_gto_blocks
combined_one_electron_assembly            2.126697333      materialized_route_global_combined_gto_one_electron_matrices
final_basis_projection                    0.032196084      materialized_route_global_combined_gto_final_basis_projection
residual_moment_matrices                  36.012651167     materialized_route_global_combined_gto_residual_moment_matrices
residual_mwg_representation               0.194609834      materialized_route_global_residual_gto_mwg_representation
gausslet_density_density                  21.287211417     materialized_route_global_density_density_interaction_matrix
final_basis_density_density               8.358507167      materialized_route_global_combined_gto_final_basis_density_density_matrix
rhf_solve                                 2.621941875      converged_final_basis_be_sp_restricted_hartree_fock
```

- Route elapsed: `354.822137417` seconds.
- RHF elapsed: `2.621941875` seconds.
- Total elapsed: `357.444079292` seconds.
- Dominant phase: `mixed_gto_blocks` at about `188.58` seconds.
- Next optimization/replacement target: mixed GTO route-global block
  materialization for the Be S+P shellification-backed route. Residual moments,
  gausslet density-density, final density-density, and RHF all completed and are
  smaller than the mixed-GTO phase.

Files changed:

- `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
- `tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
- `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`
- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.013.md`

Validation:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
  passed and produced the RHF result above.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- No additional existing acceptance test was run. The production change is a
  private opt-in timing hook on the Be seam, and the Be probe is the direct
  validation for that seam. Adding a long-term test would not satisfy the
  current carrying-cost policy for this exploratory Be pass.

Deletion/shrinkage:

- No old H/H2+ fixture-local GTO wiring was deleted. Those paths still own
  their fixture-specific scientific acceptance assertions and are not made dead
  by this timing pass.
- No tests were added. The probe remains under `tmp/work` and does not add
  permanent coverage.
- No compatibility path, fallback path, or metadata layer was added.
- Nothing was deleted because this pass revealed the next real cost center
  rather than making an existing surface obsolete.
- Remaining stale/duplicate surfaces to retire next:
  the H/H2+ fixture-local construction sequence can be considered for migration
  only after the shared atom+supplement seam is reused by an endpoint without
  weakening those scientific assertions. The immediate Be-specific target is
  reducing repeated mixed-GTO block materialization work.

-- repo-doer@macmini
