Pass 253 - independent H2 PQS private RHF input taxonomy

Context:
- Current HEAD should include `950ff895 Materialize independent H2 PQS private RHF`.
- Pass 252 showed that private RHF now materializes on the independent H2 PQS
  H1-J basis:
  - final dimension `471`;
  - `fake_pqs_enabled=false`;
  - `source_backed_fixed_source_oracle_used=false`;
  - `retained_transform_authority=:pqs_source_box_construction`;
  - total RHF diagnostic energy `-1.1589735957658853`;
  - converged in 8 iterations with small idempotency/commutator residuals.
- Pass 251 added stage-specific inputs through H1-J, but pass 252 still used
  command-line overrides for `run_private_rhf=true` and
  `private_rhf_electron_count=2`.

Task:
Add a tiny explicit driver input and manifest row for the independent H2 PQS
private RHF diagnostic stage.

Required changes:
1. Add a new include/override input:
   `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl`
   based on `h2_pqs_q5_independent_source_box_r4_h1_j.jl`.
2. The new input should set:
   - `artifact_role = :independent_h2_pqs_private_rhf_diagnostic`
   - `physics_endpoint_ready = false`
   - `physics_endpoint_blocker = :private_rhf_diagnostic_not_public_solver_contract`
   - `run_final_basis = true`
   - `run_h1 = true`
   - `run_h1_j = true`
   - `run_private_rhf = true`
   - `private_rhf_electron_count = 2`
3. Update the compact independent-H2 artifact-role classifier so this new role
   still receives:
   - `source_backed_fixed_source_oracle_used = false`
   - `retained_transform_authority = :pqs_source_box_construction`
   - fake-PQS guard behavior unchanged.
4. Add a row to `docs/src/developer/cartesian_driver_endpoint_manifest.md`
   making clear this is a private RHF diagnostic, not an endpoint or public
   solver/export-ready route.

Strict exclusions:
- Do not change RHF numerical behavior.
- Do not add a new solver.
- Do not rerun or compare fake-PQS/WL as evidence for independent PQS.
- Do not add supplements, CR2/export, HamV6, public API, public solver
  readiness, or provider-block work.
- Do not add a full route test. Pass 252 already ran the slow focused RHF
  route; this pass is input/manifest taxonomy plus any minimal classifier update.

Validation:
- Include/flag smoke for the new input. Assert at least:
  - route family/kind are independent H2 PQS;
  - artifact role is `:independent_h2_pqs_private_rhf_diagnostic`;
  - `(run_final_basis, run_h1, run_h1_j, run_private_rhf) ==
    (true, true, true, true)`;
  - `private_rhf_electron_count == 2`;
  - endpoint ready is false and the blocker is the private-diagnostic blocker.
- Run package load if source is touched.
- Run `git diff --check`.
- Do not run the slow H2 RHF route unless the source change is more than the
  role classifier or a smoke indicates a real issue.

Line budget:
- Keep scoped `src + test + bin` net-negative.
- Offset the new input/classifier lines by deleting stale mirror assertions from
  the already-classified terminal/low-order staged tests. Good candidates are
  identity comparisons or exact flat route-core/readiness/status mirrors that
  do not protect current physics behavior.

Report:
- New input path and flags.
- Manifest row added.
- Whether the classifier was updated and how fake-PQS guard fields remain
  separate.
- Validation commands.
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
