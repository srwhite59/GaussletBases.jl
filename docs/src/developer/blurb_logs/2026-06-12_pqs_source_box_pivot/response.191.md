Pass 191 response

No H2 implementation was done. This pass recovered the old WL/QW H2 reference, recommended the first driver-owned H2 PQS target shape, and deleted exactly one stale scaffold test.

Exact old WL/QW H2 reference found:

- Primary source: `docs/src/developer/high_order_endcap_panel_h2_chemistry_reproduction_2026-05-16.md`
- Supporting source: `docs/src/developer/high_order_mainline_import_readiness_2026-05-15.md`
- Historical handoff source: `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/blurb.009.md`

Reference fixture:

- molecule/geometry: H2 at `R = 4.0` bohr, bond-aligned on `:z`
- nuclear charges: `[1.0, 1.0]`
- centered coordinates implied for the driver target: `(0.0, 0.0, -2.0)` and `(0.0, 0.0, 2.0)`
- parent/basis constructor: `bond_aligned_homonuclear_qw_basis`
- family/backend: `family = :G10`; old-standard chemistry reproduction used ordinary/QW operators
- parent controls: `core_spacing = 0.5`, `xmax_parallel = 6.0`, `xmax_transverse = 4.0`, `bond_axis = :z`
- supplement: H/cc-pVTZ molecular Gaussian supplement with `lmax = 1` (S/P)
- nested size control: `nside = 5`
- density interaction route: finite IDA/QW density-density Hamiltonian
- residual interaction treatment: `interaction_treatment = :ggt_nearest`
- solver level: restricted closed-shell HF plus ED in the historical reproduction

Exact documented values:

- default complete-rectangular WL/QW route:
  - fixed block size `(1215, 463)`
  - final dimension `481`
  - residual count `18`
  - HF total `-0.910938264352`
  - ED total `-1.015613837691`
  - BO error `0.776415` mHa
  - ED-HF lowering `104.675573` mHa
  - shared layer types `_CartesianNestedCompleteShell3D`, `_CartesianNestedCompleteShell3D`
  - shared layer columns `(98, 114)`
  - fixed overlap error `1.38e-14`
  - operator overlap error `2.44e-12`
  - H symmetry error `0.0`
  - V symmetry error `0.0`
  - residual owner set `(1, 2)`
  - HF converged in `9` iterations
- endcap/panel `q=4,L=4` WL/QW route:
  - fixed block size `(1215, 443)`
  - final dimension `461`
  - residual count `18`
  - HF total `-0.910977315003`
  - ED total `-1.015663743783`
  - BO error `0.726509` mHa
  - ED-HF lowering `104.686429` mHa
  - shared layer types `_CartesianNestedEndcapPanelShellLayer3D`, `_CartesianNestedEndcapPanelShellLayer3D`
  - shared layer columns `(96, 96)`
  - fixed overlap error `1.38e-14`
  - operator overlap error `1.03e-12`
  - H symmetry error `0.0`
  - V symmetry error `0.0`
  - residual owner set `(1, 2)`
  - HF converged in `9` iterations

H1/J availability:

- I found no tracked H1-only or H1/J self-Coulomb scalar for the cited H2 `R = 4.0` rows.
- The tracked references record HF totals and ED totals, plus overlap/symmetry/residual diagnostics. Do not invent H1/J constants for the future driver test.

Recommended first H2 PQS driver target:

- First target should be the PQS source-box analog of the default complete-rectangular WL/QW H2 route, not the endcap/panel `q=4,L=4` route.
- Reason: the default complete-rectangular row is the baseline old-standard route with final dimension `481` and HF total `-0.910938264352`; matching it first avoids repeating the earlier mistake of comparing a PQS route to a different shellification policy.
- The endcap/panel `q=4,L=4` row should remain a secondary comparison after the default analog is driver-owned and comparable.
- Do not overinterpret `q = n_s`: the equivalence requirement is parent geometry, parent mapping/backend, shellification/retained-basis policy, residual supplement policy, and density-interaction convention. Matching a scalar `q` is not enough.

Future H2 driver input shape:

Candidate future file, not created in this pass:

```julia
route_family = :pqs_source_box
route_kind = :bond_aligned_diatomic_fixed_q_complete_core_shell

atom_symbols = ("H", "H")
nuclear_charges = (1, 1)
atom_locations = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))

parent_axis_family = :G10
bond_axis = :z
bond_length = 4.0
core_spacing = 0.5
xmax_parallel = 6.0
xmax_transverse = 4.0
reference_spacing = 1.0
tail_spacing = 10.0

q = 5
n_s = 5
fixed_source_mode_shape = true
diatomic_shellification_policy = :complete_rectangular

run_h1 = true
run_h1_j = true
run_private_rhf = true
private_rhf_electron_count = 2
private_rhf_fixture_role = :physics_endpoint
route_configured_diatomic_ham_interaction_treatment = :ggt_nearest

comparison_reference_label = "WL/QW H2 R=4.0 default complete rectangular"
wl_rhf_total = -0.910938264352
```

Important caveat: the documented WL/QW HF total includes the H/cc-pVTZ S/P residual supplement. If the first PQS driver target is gausslet-only, it should not compare directly to `-0.910938264352`. Either the future driver target must carry the same supplement/residual policy, or the comparison should first be restricted to comparable basis inventory/H1/H1-J diagnostics until the supplement path exists.

Future artifact/check shape:

- `system/*`: atom symbols, charges, centered coordinates, bond length, bond axis, electron count
- `config/*`: route family/kind, parent constructor/mapping/backend, core spacing, extents, q/n_s, shellification policy, interaction treatment, supplement policy
- `basis/final_dimension`
- `basis/fixed_block_size` or compact equivalent
- `basis/shared_layer_types`
- `basis/shared_layer_columns`
- `basis/final_overlap_identity_error`
- `density_interaction/status`
- `density_interaction/density_gauge`
- `density_interaction/raw_pair_factor_convention`
- `physics/h1_lowest` only after available for the H2 route
- `physics/h1_j_self_coulomb` only after available for the H2 route
- `private_rhf/status`
- `private_rhf/total_energy`
- `private_rhf/iteration_count`
- `private_rhf/converged`
- `private_rhf/residual`
- `comparison/reference_label`
- `comparison/wl_rhf_total`
- `comparison/delta_rhf`

Do not store final-basis self-overlap matrices as downstream working data. Keep overlap as scalar diagnostics such as identity error.

Deleted stale test/probe:

- Deleted `test/nested/cartesian_pair_block_one_body_block_set_consumption_skeleton_runtests.jl`.
- Safety verification:
  - no include in `test/nested/runtests.jl`
  - no include in `test/nested/integration_runtests.jl`
  - no top-level/group runner include found
  - no `test/src/bin` caller remained after deletion
  - not a He endpoint, not WL H/H2 acceptance/reference, not a scientific endpoint
  - only references found were inside the file itself and stale archived handoff logs
- Rationale: this was a metadata-only block-set consumption skeleton scaffold from the mixed one-body transition. The default runner already keeps the smaller `cartesian_pair_block_one_body_consumer_smoke_runtests.jl` for routine status/count validation.

Line budget:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
0	205	test/nested/cartesian_pair_block_one_body_block_set_consumption_skeleton_runtests.jl
```

Arithmetic: added `0`, deleted `205`, net `-205`.

Validation:

- `rg -n "cartesian_pair_block_one_body_block_set_consumption_skeleton|one_body_block_set_consumption_skeleton|CartesianPairBlockMaterialization one-body block-set consumption skeleton" test src bin`: no hits.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`: passed, `load ok`.
- `git diff --check`: passed.
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`: `0 added / 205 deleted`.
- `git status --short --branch`:

```text
## main...origin/main
 D test/nested/cartesian_pair_block_one_body_block_set_consumption_skeleton_runtests.jl
```

Deletion/shrinkage report:

- deleted:
  - `test/nested/cartesian_pair_block_one_body_block_set_consumption_skeleton_runtests.jl`
- simplified:
  - Removed unreferenced metadata-only mixed one-body skeleton test pressure; default small smoke remains.
- quarantined:
  - H2 remains audit-only; no input file, source code, test, RHF change, or driver behavior was added.
- not deleted because:
  - live included materializer/report integration tests remain included in `test/nested/integration_runtests.jl`.
  - H/H2 WL acceptance/reference tests are scientific endpoints and were not touched.
  - new explicit He endpoint tests were not touched.
- exact remaining caller/blocker:
  - H2 PQS implementation is blocked on a manager-approved driver pass that chooses whether to carry the H/cc-pVTZ S/P residual supplement in the driver target or to start with a gausslet-only H2 diagnostic that does not compare against the documented HF totals.

-- repo-doer@macmini
