Pass 172 - tighten CR2 HF-facing artifact convention labels, net-negative

Role: repo-doer@macmini

Task type: implementation plus cleanup/shrinkage.

Purpose:

CR2 has now consumed the Be2 WL/PQS artifact as provisional matrix data. Both
routes have real final-basis `H1` and final-basis density-density `Vee`; CR2
can form provisional `V*n` and `H1 + Diagonal(V*n)` objects. The remaining
blocker is a reviewed HF density-density Fock/energy convention, not missing
matrices.

Tighten the private artifact schema so each route declares its overlap
convention and the common final-basis `Vee` label consistently, while keeping
solver/export readiness false.

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

  Touch the generator if needed for WL convention labels or readback output.
  Count it in the line budget.

Required schema changes:

PQS:

- Add route/final-basis convention fields such as:

  ```text
  overlap_convention = :orthonormal_identity_by_contract
  overlap_matrix_stored = false
  overlap_identity_defect = 0.0
  ```

- Do not add or preserve a PQS self-overlap matrix as downstream data.
- Keep `routes/pqs_source_box/two_body/interaction_matrix`.
- Keep pre-final provenance/debug fields.

WL:

- Add convention fields such as:

  ```text
  overlap_convention = :orthonormal_identity_with_diagnostic_matrix
  overlap_matrix_stored = true
  overlap_identity_defect = norm(overlap - I, Inf)
  ```

- Keep stored `routes/white_lindsey/one_body/overlap` as diagnostic/available
  route data.
- Do not require CR2 to treat WL as a generalized-overlap HF route when the
  identity defect is small and the convention says diagnostic identity-noise.

Both routes:

- Standardize on the common final-basis Vee key:

  ```text
  two_body/interaction_matrix_representation_kind =
      :final_basis_density_density_matrix
  ```

- WL currently has `two_body/representation_kind`; keep it only if needed for
  compatibility, but add the common key.
- Add a compact blocked HF convention section or route metadata fields such as:

  ```text
  density_density_hf_convention_status =
      :missing_reviewed_density_density_hf_fock_energy_convention
  density_density_hf_convention_blocker =
      :missing_reviewed_density_density_hf_fock_energy_convention
  ```

- Keep:

  ```text
  cr2_solver_ready = false
  cr2_export_ready = false
  hfdmrg_ready = false
  rhf_ready = false
  dmrg_ready = false
  ```

Cleanup/shrinkage target:

- First look for stale or repeated assertions in
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`.
  Prefer replacing repeated status/metadata assertions with the new concise
  convention checks.
- If that is not enough for net-negative, inspect these concrete stale
  non-default candidates before choosing any broader deletion:

  ```text
  test/nested/cartesian_cpb_overlap_placement_facts_runtests.jl
  test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl
  test/nested/cartesian_cpb_reviewed_overlap_placement_plan_runtests.jl
  ```

  These are metadata/blocker-heavy overlap-placement scaffolds, not Be2/Ham
  scientific endpoint checks. Prefer shrinking repeated nonclaim/blocker
  assertions if a small cut is enough; delete only if inspection shows the file
  is preserving superseded transitional vocabulary and no active include or
  endpoint depends on it.
- Do not delete CR2 inspection generator or the focused Be2 artifact test.

Validation:

Run:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Read back the generated JLD2 enough to confirm:

- PQS overlap convention is identity-by-contract and no PQS overlap matrix is
  stored;
- WL overlap convention is diagnostic identity-noise with stored overlap and
  identity defect;
- both routes expose
  `two_body/interaction_matrix_representation_kind =
  :final_basis_density_density_matrix`;
- both routes carry the density-density HF convention blocker;
- solver/export flags remain false;
- generated outputs remain ignored.

Hard boundaries:

- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense four-index `Vee`,
  `V6`, `Vblocks`, or solver code.
- Do not mark the artifact solver-ready.
- Do not add PQS overlap as downstream data.
- Do not change route construction, WL/PQS fixture size, or physics values.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Report:

- files changed;
- tracked implementation lines added, deleted, and net result;
- schema fields added/standardized;
- PQS and WL overlap convention values;
- common Vee key used by both routes;
- HF convention blocker value;
- validation commands/results;
- generated artifact provenance/readback summary;
- ignored artifact status;
- `git status --short --branch`;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
