Pass 171 - add final-basis PQS Vee to CR2 inspection artifact, net-negative

Role: repo-doer@macmini

Task type: implementation plus cleanup/deletion.

Purpose:

CR2 confirmed the Be2 WL/PQS artifact is real inspection data, but asked the
right consumer-contract question: for a first HF-facing smoke, the minimal
object is `H1 :: N x N` plus `Vee :: N x N`. WL already has that shape. PQS has
real two-body data, but currently exposes it only as a pre-final density
interaction plus final-to-pre-final coefficients.

Add a dense final-basis PQS density-density interaction matrix to the private
CR2 artifact while preserving the existing pre-final provenance/debug fields.
Keep solver/export readiness false.

Hard line-budget rule:

- Final tracked implementation diff must be net-negative across `src`, `test`,
  and the tracked CR2 generator by:

  ```text
  git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  ```

- Require `sum(deleted) > sum(added)`.
- Do not satisfy this by deleting scientific endpoint tests or moving code to
  untracked/tmp files.
- If you cannot keep the pass net-negative, write `.agent_handoffs/ATTENTION.md`
  and stop.

Cleanup/deletion target to pay for this pass:

- Delete the stale integration test:

  ```text
  test/nested/white_lindsey_materialized_seed_runtests.jl
  ```

- Remove its include from:

  ```text
  test/nested/integration_runtests.jl
  ```

Reason: this test is a private one-center/materialized-seed development
scaffold. The current CR2 artifact uses the route-configured diatomic
atom-growth White-Lindsey path instead, and the accepted diatomic route smoke
already validates the live WL ham bundle surface. Do not delete
`src/white_lindsey_materialized_seed.jl` or the seed oracle summary in this
pass; they still have source references and can be considered in a later
cleanup.

Decision rule for deletion:

- If you find `white_lindsey_materialized_seed_runtests.jl` is still the only
  test protecting a live scientific endpoint or active route consumer, stop and
  report the exact gap instead of deleting it.
- If it is only preserving old seed route-shadow vocabulary, delete it.

Implementation surface:

- Primary source helper:

  ```text
  src/pqs_source_box_diatomic_complete_core_shell.jl
  _pqs_source_box_route_driver_be2_cr2_inspection_bundle_payload
  ```

- Focused test:

  ```text
  test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
  ```

- Generator:

  ```text
  tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  ```

  Touch the generator only if it needs readback/summary printing updates. The
  schema itself should come from the source helper above, not a generator-only
  patch.

PQS final-basis Vee contract:

- Existing objects in the source helper:

  ```julia
  pair_matrix = Matrix{Float64}(handoff.pre_final_pair_matrix)
  coefficients = Matrix{Float64}(handoff.final_to_pre_final_coefficients)
  ```

- Treat `coefficients` as the final-to-pre-final density map. Form the dense
  final-basis density-density matrix as:

  ```julia
  final_interaction_matrix = transpose(coefficients) * pair_matrix * coefficients
  ```

- If inspection of the current code shows the orientation is opposite, stop and
  report the exact evidence before implementing. Do not guess.

Artifact schema:

- Under `routes/pqs_source_box/two_body`, add:

  ```text
  interaction_matrix
  interaction_matrix_representation_kind = :final_basis_density_density_matrix
  interaction_matrix_derivation = :final_to_pre_final_density_congruence
  interaction_matrix_formula = :transpose_final_to_pre_final_times_pre_final_pair_times_final_to_pre_final
  interaction_matrix_shape
  interaction_matrix_symmetry_defect
  interaction_matrix_finite
  ```

- Keep existing fields:

  ```text
  representation_kind = :pre_final_density_interaction
  pre_final_pair_matrix
  final_to_pre_final_coefficients
  pre_final_weights
  support_weights
  support_raw_pair_numerator
  density_gauge
  raw_pair_factor_convention
  ```

- Do not remove pre-final provenance in this pass.
- Do not mark PQS solver/export ready.

Validation in the focused test:

- Assert the new PQS `interaction_matrix` exists and has shape `(221, 221)`.
- Assert it is finite and symmetric to the same practical tolerance as the
  pre-final matrix.
- Add a small convention check using a few deterministic final density vectors,
  for example basis vector(s) and a short normalized vector:

  ```text
  d_pre = final_to_pre_final_coefficients * d_final
  d_final' * interaction_matrix * d_final
    ~= d_pre' * pre_final_pair_matrix * d_pre
  ```

- Also assert the new representation/derivation labels.
- Keep the test compact. Do not add broad metadata assertions.

Generated artifact validation:

- Run:

  ```text
  julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  ```

  This may take several minutes because it builds both PQS and WL artifact
  routes. That is acceptable because the generated artifact is the deliverable.

- Read back enough JLD2 to confirm:
  - clean producer commit/dirty status if run after commit is not required yet;
    before commit, `dirty=true` is expected;
  - PQS `routes/pqs_source_box/two_body/interaction_matrix` shape `(221, 221)`;
  - PQS `interaction_matrix_representation_kind`;
  - PQS solver/export flags remain false;
  - WL route remains available and unchanged in kind.

Hard boundaries:

- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense four-index `Vee`,
  `V6`, `Vblocks`, or solver code.
- Do not convert this to a solver-ready contract.
- Do not change PQS route construction, WL route construction, lattice size, or
  physics fixture.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Validation commands:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Report:

- files changed;
- `src`/`test`/tracked-generator lines added, deleted, and net result;
- deleted test/include and why it was stale;
- final interaction matrix formula and orientation evidence;
- focused test results;
- generator run/readback summary;
- whether generated JLD2/TSV remain ignored;
- `git status --short --branch`;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
