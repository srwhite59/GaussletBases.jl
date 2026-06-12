What changed:

- Added a private CPBM construction seam:
  `_white_lindsey_decomposed_atom_gto_final_basis_route(...)` in
  `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`.
- The seam wires mapped parent axes, shellification-backed decomposed WL
  inventory, route-global overlap/kinetic/electron-nuclear, mixed GTO blocks,
  combined one-electron assembly, final-basis projection, and optional
  residual-MWG/final density-density materialization.
- The seam is internal/private. It does not add exports, public driver defaults,
  PQS, ECP, high-l Be, Be2, Cr, H2 work, full-parent CPB fallback, direct
  Cartesian fallback, ordinary IDA fallback, raw GTO final density-density, or a
  generalized final solve.
- Added `include_moment_blocks` to
  `route_global_mixed_gto_blocks_from_decomposed_units(...)` so one-electron-only
  routes can materialize mixed overlap/kinetic/nuclear blocks without also
  forcing position/x2 residual-moment blocks. The default remains `true`, so the
  existing H/H2+ + GTO final-basis path keeps its residual-moment behavior.
- Added a developer probe:
  `tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`.
- Wrote the probe summary artifact:
  `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`.

Be S+P current-route result:

- Status:
  `:materialized_decomposed_atom_gto_final_basis_one_electron_route`.
- Blocker: `nothing` for final-basis one-electron materialization.
- Authorized supplement source used:
  `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets`.
- Fixture:
  Be, `Z = 4`, `q/ns = 5/5`, `lmax = 1`, `d = 0.15`, parent side count `15`,
  endpoints `(-10.056710949734484, 10.056710949734484)`.
- Shellification-backed decomposed WL inventory used: `true`.
- Low-order seed inventory used: `false`.
- Fallback flags:
  full-parent CPB `false`, direct Cartesian product assembly `false`,
  ordinary Cartesian IDA operators `false`, raw GTO final density-density
  `false`, generalized final solve `false`.
- Current-route dimensions:
  retained gausslet dimension `615`, unit count `131`, pair count `8646`,
  raw supplement count `21`, retained supplement count `21`, final dimension
  `636`.
- Old oracle dimension comparison:
  fixed `615`, residual `21`, final `636`; final dimension delta `0`.
- Final-basis diagnostics:
  final overlap identity error `1.0125523569644675e-10`, final overlap rank
  `636`, final overlap condition estimate `1.0000000002807463`, final
  Hamiltonian symmetry error `1.2108358760087867e-10`.
- Density-density/RHF was not attempted in the durable rerun. The first
  uninstrumented full-density run reached the gausslet density-density stage and
  was stopped because it lacked useful phase attribution. The next blocker is
  therefore not the old driver-owned wiring blocker; it is a measured,
  phase-attributed final density-density/RHF run for Be S+P.
- Old oracle RHF values remain comparator-only:
  one-electron `-19.06620047058102`, electron-electron `4.491686226006327`,
  total `-14.574514244574694`.

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
- `src/cartesian_pair_block_materialization/route_global_mixed_gto_blocks.jl`
- `tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.012.md`

Validation:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
  passed and wrote the summary above. The Be S+P route elapsed time was about
  `148.24` seconds for the one-electron final-basis materialization.
- `julia --project=. test/nested/cartesian_wl_gausslet_h_atom_acceptance_runtests.jl`
  passed:
  H gausslet-only `78/78`, H2+ gausslet-only `75/75`, H + GTO final-basis
  `62/62`, H2+ + GTO final-basis `46/46`.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage:

- No old H/H2+ fixture-local GTO wiring was deleted in this pass. Those tests
  still own fixture-specific center records, supplement choices, scientific
  assertions, and active acceptance energies.
- The new seam makes the route sequence reusable for Be S+P and is the first
  step toward shrinking duplicated H/H2+ fixture-local wiring, but it does not
  yet replace those acceptance reports. The next removal condition is migrating
  H/H2+ GTO acceptance setup to the shared atom/supplement route seam or a
  closely related driver-facing helper while preserving their endpoint
  assertions.
- The mixed-GTO route was simplified for one-electron-only callers by bypassing
  the full CPB/GTO local operator bundle when residual moment blocks are not
  requested. That removed unnecessary position/x2 mixed block construction from
  the Be one-electron probe path.
- No new long-term tests were added. The Be script is a `tmp/work` developer
  probe, not acceptance coverage.
- Remaining stale/duplicate surfaces to retire next:
  H/H2+ fixture-local GTO construction sequence duplication, and any
  one-electron-only probe path that still forces residual moment or
  density-density setup before it is needed.

-- repo-doer@macmini
