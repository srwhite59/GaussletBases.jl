Pass 200 response - physical H2 gausslet-only target audit.

Old WL gausslet-only inventory:

- Documentation row: default complete-rectangular H2 R=4 fixed block `(1215, 463)`, supplemented final dimension `481`, residual count `18`, shared layer columns `(98, 114)`.
- Direct old-source inspection with `bond_aligned_homonuclear_qw_basis(... bond_length=4.0, core_spacing=0.5, xmax_parallel=6.0, xmax_transverse=4.0 ...)` and `bond_aligned_diatomic_nested_fixed_source(... nside=5 ...)` confirmed:
  - parent axis lengths `(9, 9, 15)`, parent size `1215`;
  - fixed block size `(1215, 463)`;
  - retained order is core/child first, then shared layers:
    - core/child columns `1:251`;
    - shared layer 1 columns `252:349`, count `98`;
    - shared layer 2 columns `350:463`, count `114`;
  - shared shell support counts `(578, 362)`;
  - one child/core sequence with support count `275` and core retained count `251`;
  - child/core support is the 5 x 5 x 11 atom-contact working box, i.e. two full 5^3 atom cores plus the 25-row contact plane as support;
  - `midpoint_slab_column_range = nothing`; there is no separate midpoint/product retained unit in the old gausslet-only 463-column fixed block.
- Interpretation: the `463` gausslet-only columns are `251` atom-contact child/core retained columns plus two complete-rectangular shared molecular shell layers `(98, 114)`. The later `481` is `463 + 18` H/cc-pVTZ S/P residual supplement columns.

Proposed physical H2 PQS target:

- Route purpose: physical gausslet-only H2 R=4 endpoint target, not the current 221-dimensional `source_box_diagnostic`.
- Source/final units:
  - `atom_contact_core`: full atom-contact core support with full 5^3 left and right atom-core interiors retained, plus contact-plane support;
  - `shared_shell_1`: PQS source-box-first retained representation of the first old shared molecular shell;
  - `shared_shell_2`: PQS source-box-first retained representation of the second old shared molecular shell.
- Support inventory to match before operator work:
  - atom-contact core support count `275`;
  - shared shell support counts `(578, 362)`;
  - full parent support count `1215`.
- Retained order should match the old fixed-block convention first: core/child retained columns, then shell layers in source order.
- Expected retained counts if the PQS target exactly matches the old gausslet-only fixed-block inventory: `(251, 98, 114)`, final dimension `463`. If a later PQS shell cleanup changes rank, that should be an explicit target-result field, not an implicit route drift.
- The midpoint/product slab should not remain as an independent retained unit for this target. The 25-row contact plane belongs inside the atom-contact core support.
- Full atom-core interiors are retained through the atom-contact core unit, not by the current boundary-only `pqs_left`/`pqs_right` source boxes.
- Shell layers should be represented as source-box-first PQS shell units over the old shared molecular shell supports, with shell/support-row contraction staying oracle/debug.

Recommended pass-201 seam:

- Introduce a new route kind rather than mutating the existing diagnostic route. Suggested name:
  - `:bond_aligned_diatomic_physical_gausslet_core_shell_pqs`
- First implementation pass should add only a compact private target/inventory payload, for example:
  - `_PQSDiatomicPhysicalGaussletCoreShellTargetPayload`
- Suggested payload fields:
  - `status`, `blocker`, `route_family`, `route_kind`, `parent_axis_counts`, `support_units`, `retained_units`, `support_counts`, `retained_counts`, `retained_order`, `expected_final_dimension`, `retained_atom_core_interiors`, `source_plan_role`, `supplement_policy`, `summary`, `metadata`.
- Reuse:
  - current bond-aligned parent-axis construction;
  - old nested H2 source/fixed-source inventory as an oracle for counts/order;
  - atom-growth construction-plan region facts where they expose the atom/contact/shared partition;
  - existing driver/report staged style.
- Do not mutate:
  - `:bond_aligned_diatomic_fixed_q_complete_core_shell`;
  - current `pqs_left/product/pqs_right` diagnostic source-plan objects;
  - current 221-dimensional H1 path;
  - H1-J/RHF/density interaction.

Future physical H2 artifact status before H1-J/RHF:

- `artifact_role = :physical_gausslet_endpoint_target`
- `physics_endpoint_ready = false`
- `retained_atom_core_interiors = true`
- `source_plan_role = :atom_contact_core_plus_pqs_shared_shells`
- `comparison_ready = false` until a gausslet-only WL/QW comparison row is explicitly produced or selected
- `supplement_policy = :none`

Stale scaffold deleted:

- Deleted `test/nested/cartesian_pair_block_one_body_plan_batch_runtests.jl`.
- Safe because it is not included by default/integration/top-level runners, has no source caller, is not a He/H2 driver endpoint, is not a WL H/H2 scientific reference, and is not needed by this audit.
- The live mixed one-body consumer smoke remains in `test/nested/cartesian_pair_block_one_body_consumer_smoke_runtests.jl` and covers the compact direct/PQS/LW status/count contract.

Validation:

- `rg -n "cartesian_pair_block_one_body_plan_batch_runtests|CartesianPairBlockMaterialization mixed one-body direct-only plan batch|CartesianPairBlockMaterialization mixed one-body direct/PQS plan batch|CartesianPairBlockMaterialization mixed one-body plan batch missing inputs|CPBMBatchDispatch" test src bin`
  - no matches; exit code 1 as expected.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed, printed `load ok`.
- `git diff --check`
  - passed.
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - `0	480	test/nested/cartesian_pair_block_one_body_plan_batch_runtests.jl`
  - scoped source/test/bin budget is net -480.
- `git status --short --branch`
  - `## main...origin/main`
  - ` D test/nested/cartesian_pair_block_one_body_plan_batch_runtests.jl`

Deletion/shrinkage report:

- deleted: `test/nested/cartesian_pair_block_one_body_plan_batch_runtests.jl` (480 lines).
- simplified: mixed one-body batch test pressure; the remaining compact smoke is the routine live-contract gate.
- quarantined: none.
- not deleted because: no route implementation was performed; current H2 221 diagnostic route objects remain until the physical target route exists.
- exact remaining caller/blocker: physical H2 PQS endpoint is blocked by the absence of a new route-owned target/inventory payload for full atom-contact core interiors plus shared shell layers; the existing `:bond_aligned_diatomic_fixed_q_complete_core_shell` route remains boundary-source-box diagnostic only.

-- repo-doer@macmini
