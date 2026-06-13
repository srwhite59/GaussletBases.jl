Pass 173 - shrink CR2 artifact by removing downstream WL overlap matrix

Role: repo-doer@macmini

Task type: generator optimization/shrinkage.

Purpose:

The Be2 WL/PQS artifact is now HF-input-readable for CR2 in the minimal matrix
sense, but the generator is slow and the ignored JLD2 bundle is large
(`~242M`). Before the next CR2 handoff, shrink the artifact contract where the
data are diagnostic only.

The first target is the White-Lindsey overlap matrix. CR2 observed that WL `S`
is identity to numerical noise (`~4e-14`). Like PQS, this should be a
GaussletBases-side diagnostic/check, not downstream working data. CR2 should
consume both routes as orthonormal final-basis `H1 + Vee` inputs, with overlap
convention labels and identity-defect diagnostics.

Hard line-budget rule:

- Final tracked implementation diff must be net-negative across `src`, `test`,
  and the tracked CR2 generator by:

  ```text
  git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  ```

- Require `sum(deleted) > sum(added)`.
- Count the tracked generator.
- Do not satisfy this by deleting scientific endpoint tests or moving code to
  untracked/tmp files.
- If the generator-only change is not net-negative, inspect these concrete
  stale non-default candidates before broader deletion:

  ```text
  test/nested/cartesian_cpb_overlap_placement_facts_runtests.jl
  test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl
  ```

Implementation surface:

- Prefer editing only:

  ```text
  tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  ```

- Do not edit `src/` unless a source helper already writes WL overlap into the
  artifact. Currently WL population is generator-local.
- Do not add new tests for this generator-only artifact shrink.

Required change:

- Keep computing/checking WL overlap enough to validate:

  ```text
  overlap_identity_defect = norm(overlap - I, Inf)
  ```

- Stop storing the dense WL overlap matrix in the artifact:

  ```text
  routes/white_lindsey/one_body/overlap
  ```

- Update WL final-basis convention fields to indicate the matrix is diagnostic
  and not stored downstream, for example:

  ```text
  overlap_convention = :orthonormal_identity_diagnostic_checked_not_stored
  overlap_matrix_stored = false
  overlap_identity_defect = <computed defect>
  ```

- Keep:

  ```text
  routes/white_lindsey/one_body/hamiltonian
  routes/white_lindsey/two_body/interaction_matrix
  routes/white_lindsey/two_body/interaction_matrix_representation_kind
  routes/white_lindsey/hf_convention/*
  ```

- PQS should remain unchanged:

  ```text
  overlap_convention = :orthonormal_identity_by_contract
  overlap_matrix_stored = false
  no PQS overlap matrix
  ```

Validation/readback:

Run the generator:

```text
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Then read back the JLD2 and report:

- producer commit and dirty flag;
- bundle size before/after if available from `ls -lh`;
- PQS overlap convention and no PQS overlap matrix;
- WL overlap convention;
- WL `overlap_matrix_stored == false`;
- absence of `routes/white_lindsey/one_body/overlap`;
- WL `overlap_identity_defect`;
- PQS/WL H1 and Vee shapes still present;
- solver/export flags still false.

Also run:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Hard boundaries:

- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense four-index `Vee`,
  `V6`, `Vblocks`, or solver code.
- Do not mark the artifact solver-ready.
- Do not add PQS overlap or WL overlap as downstream working data.
- Do not change route construction, WL/PQS fixture size, or physics values.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Report:

- files changed;
- tracked implementation lines added, deleted, and net result;
- artifact size before/after;
- generator runtime if measured cheaply with Julia-level timing or shell
  elapsed already available;
- WL overlap convention/readback;
- H1/Vee shape readback for both routes;
- validation commands/results;
- generated outputs ignored status;
- `git status --short --branch`;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
