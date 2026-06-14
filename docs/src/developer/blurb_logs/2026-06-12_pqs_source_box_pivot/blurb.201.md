Pass 201 - implement the physical H2 PQS target inventory only.

Target card:

- Physics target: H2, R = 4.0, gausslet-only, `q = n_s = 5`, parent geometry
  matching the old WL/QW H2 fixed-block route.
- Endpoint being enabled in this pass: target/inventory readiness only.
- Allowed implementation surface: private route-owned target/inventory payload,
  driver input, compact artifact/readiness fields, and one focused explicit
  test.
- Forbidden surfaces: source-plan construction, final-basis construction, H1,
  H1-J, density interaction, RHF, residual/supplement columns, WL comparison,
  CR2/export work, public APIs.
- Success condition: the visible driver can run a new physical-H2 target input
  and write an artifact saying the target inventory is available but downstream
  source-plan/final-basis/H1 are still intentionally unavailable.

Governing context:

- Pass 200 established the old gausslet-only WL/QW inventory:
  - parent axis counts `(9, 9, 15)`, parent size `1215`;
  - fixed block `(1215, 463)`;
  - supplemented reference final dimension `481 = 463 + 18`;
  - core/child retained columns `1:251`;
  - shared shell retained counts `(98, 114)`;
  - shared shell support counts `(578, 362)`;
  - child/core support count `275`;
  - child/core support is the `5 x 5 x 11` atom-contact working box.
- The current `:bond_aligned_diatomic_fixed_q_complete_core_shell` route is
  deliberately diagnostic-only and remains `221 = 98 + 98 + 25`.
- The 25-row contact plane belongs inside the atom-contact core for the physical
  target. It must not remain a standalone retained midpoint/product unit.

Implementation task:

1. Add a new route kind for the physical target, for example:

   ```julia
   :bond_aligned_diatomic_physical_gausslet_core_shell_pqs
   ```

   Do not mutate the existing diagnostic route kind.

2. Add a compact private target/inventory payload, probably in
   `src/pqs_source_box_diatomic_complete_core_shell.jl`, near the other private
   diatomic payloads.

   Suggested object name:

   ```julia
   _PQSDiatomicPhysicalGaussletCoreShellTargetPayload
   ```

   Suggested fields:

   ```text
   status
   blocker
   route_family
   route_kind
   parent_axis_counts
   support_units
   retained_units
   support_counts
   retained_counts
   retained_order
   expected_final_dimension
   retained_atom_core_interiors
   source_plan_role
   supplement_policy
   available_objects
   missing_objects
   summary
   metadata
   ```

   Keep this structured. Do not add a flat field cloud to the route report.

3. The payload should describe this target:

   ```text
   support units:
     atom_contact_core, shared_shell_1, shared_shell_2

   support counts:
     atom_contact_core = 275
     shared_shell_1    = 578
     shared_shell_2    = 362

   retained counts:
     atom_contact_core = 251
     shared_shell_1    = 98
     shared_shell_2    = 114

   retained order:
     (:atom_contact_core, :shared_shell_1, :shared_shell_2)

   expected_final_dimension = 463
   retained_atom_core_interiors = true
   source_plan_role = :atom_contact_core_plus_pqs_shared_shells
   supplement_policy = :none
   ```

   If you can derive these facts from existing old WL/QW inventory helpers
   cheaply and without constructing H1, do so. If that would force a broad route
   build, encode this as a private target inventory derived from the reviewed
   pass-200 contract and record provenance in metadata. Do not run the old
   supplemented route.

4. Add a new small driver input, separate from the diagnostic 221 input:

   ```text
   test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
   ```

   It should keep the same H2 geometry and parent probe choices, but use the
   new route kind and target labels:

   ```julia
   artifact_role = :physical_gausslet_endpoint_target
   physics_endpoint_ready = false
   physics_endpoint_blocker = :missing_physical_gausslet_source_plan
   retained_atom_core_interiors = true
   source_plan_role = :atom_contact_core_plus_pqs_shared_shells
   supplement_policy = :none
   comparison_ready = false
   run_final_basis = false
   run_h1 = false
   run_h1_j = false
   run_private_rhf = false
   ```

   Use a blocker name that fits the implementation if the exact label above is
   awkward, but it must say the target exists and the source-plan producer is
   still missing. It must not say atom-core interiors are missing.

5. Extend the driver artifact/readiness output just enough for this input to
   expose the compact target inventory. A reasonable artifact group is:

   ```text
   target/status
   target/blocker
   target/support_units
   target/support_counts
   target/retained_units
   target/retained_counts
   target/retained_order
   target/expected_final_dimension
   target/retained_atom_core_interiors
   target/source_plan_role
   target/supplement_policy
   ```

   Keep `route/source_plan_status`, `route/final_basis_status`,
   `route/h1_status`, H1-J, RHF, export, and artifact claims false/unavailable
   for this physical target pass.

6. Add one focused explicit test that runs the visible driver with the new input
   and reads the artifact. Suggested name:

   ```text
   test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
   ```

   Keep it out of the default nested runner. The test should assert:

   ```text
   route_kind == :bond_aligned_diatomic_physical_gausslet_core_shell_pqs
   parent_axis_counts == (x = 9, y = 9, z = 15)
   target expected dimension == 463
   target support counts == (275, 578, 362)
   target retained counts == (251, 98, 114)
   retained_atom_core_interiors == true
   source_plan_role == :atom_contact_core_plus_pqs_shared_shells
   supplement_policy == :none
   physics_endpoint_ready == false
   comparison_ready == false
   final_basis/H1/H1-J/RHF are not materialized
   ```

   The test should not reconstruct source plans or operators itself.

Line budget and deletion requirement:

- The source/test/bin/generator line budget remains mandatory:

  ```text
  git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  sum(deleted) > sum(added)
  ```

- First deletion candidate:

  ```text
  test/nested/cartesian_pair_block_one_body_lw_plan_batch_runtests.jl
  ```

  Delete it only if a live search confirms it is not included by a runner and
  no source/test/bin code needs that exact test file. The smaller
  `cartesian_pair_block_one_body_lw_record_dispatch_runtests.jl` and the default
  mixed one-body consumer smoke should remain.

- If that is not enough to keep the pass net-negative, inspect
  `test/nested/cartesian_pair_block_one_body_accessors_contract_runtests.jl`
  as a secondary stale scaffold candidate. Do not delete accepted He/H2 physics
  endpoint tests or WL/H2 reference tests.

Validation:

- Run the new focused H2 physical-target driver test.
- Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- Run `git diff --check`.
- Run deleted-test live searches for any removed test file names.
- Report the scoped source/test/bin/generator line budget.

Reporting:

- State whether the target inventory is hard-coded from the pass-200 reviewed
  contract or derived from an existing old-WL inventory helper.
- State the target inventory numbers.
- State the next blocker after this pass.
- Include deletion/shrinkage report:

  ```text
  deleted:
  simplified:
  quarantined:
  not deleted because:
  exact remaining caller/blocker:
  ```

Stop conditions:

- If implementing the target payload would require source-plan/final-basis/H1
  construction in the same pass, stop and write `ATTENTION.md`.
- If the new route kind accidentally falls through to the old 221-dimensional
  diagnostic source-plan/H1 path, stop and write `ATTENTION.md`.
- If the pass cannot be made line-negative without deleting a live endpoint or
  reference test, stop and write `ATTENTION.md`.

-- repo-manager@macmini
