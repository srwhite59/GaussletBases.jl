Review result:

Accepted with one important follow-up: the driver-owned Be S+P seam now
materializes the current-route final-basis one-electron path, but the one-
electron route is slow enough that the next density/RHF pass must be
phase-attributed rather than a blind long run.

What landed:

- private CPBM seam:
  `_white_lindsey_decomposed_atom_gto_final_basis_route(...)`
- seam file:
  `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
- include wiring in `CartesianPairBlockMaterialization.jl`
- `include_moment_blocks` option in
  `route_global_mixed_gto_blocks_from_decomposed_units(...)`

Be S+P current-route result:

- status:
  `:materialized_decomposed_atom_gto_final_basis_one_electron_route`
- shellification-backed decomposed inventory: true
- retained gausslet dimension: `615`
- unit/pair count: `131 / 8646`
- raw/retained supplement count: `21 / 21`
- final dimension: `636`
- final overlap identity error: `1.0125523569644675e-10`
- final overlap rank: `636`
- final overlap condition estimate: `1.0000000002807463`
- final Hamiltonian symmetry error: `1.2108358760087867e-10`
- old oracle final dimension delta: `0`

Boundary review:

The new seam does not add exports, public defaults, PQS, ECP, full-parent CPB
fallback, direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO
final density-density, or generalized final solves. The Be probe stayed under
`tmp/work`.

The default parent-axis backend is `:numerical_reference`, which is consistent
with the current He decomposed acceptance helper. This is acceptable for this
checkpoint, but future production-facing route language should keep backend
scope explicit.

Validation reviewed:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `julia --project=. test/nested/cartesian_wl_gausslet_h_atom_acceptance_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- manager reran `julia --project=. -e 'using GaussletBases; println("load ok")'`
- manager reran `git diff --check`
- artifact: `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`

Deletion/shrinkage review:

The old `:missing_driver_owned_decomposed_be_sp_fixture_wiring` blocker is
removed. No H/H2+ fixture-local GTO wiring was deleted yet, because those
acceptance reports still own endpoint-specific setup and assertions. The
mixed-GTO one-electron-only path avoids constructing residual moment blocks
when they are not requested, which is a real simplification for one-electron
probes.

Commit/push:

Pending manager commit/push of source, docs, response/review log, and the next
published blurb.

Next target:

Run a phase-attributed Be S+P density-density/RHF pass through the new seam.
Do not start with a blind full-density run; first add or use timing
attribution that can identify whether the cost is moment blocks, residual MWG,
gausslet density-density, final density assembly, or RHF.
