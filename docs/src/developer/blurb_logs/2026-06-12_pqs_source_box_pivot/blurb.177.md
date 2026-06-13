Pass 177 - add WL/PQS comparison-readiness gate and WL Z audit labels

Role: repo-doer@macmini

Task type: generator schema correction / comparability guard.

Purpose:

The current Be2 CR2 inspection artifact is useful, but it must not imply that
the current WL and PQS matrices are an apples-to-apples comparison:

```text
PQS final dimension = 221
WL final dimension  = 2287
```

CR2 can inspect both routes and can run route-local smoke experiments, but the
current pair is not the intended same-retained-basis WL-vs-PQS comparison. The
artifact should say that explicitly.

There is also a suspicious WL input:

```text
nuclear_charges = (4, 4)
white_lindsey_Z = 2.0
```

Do not decide the physics meaning of `white_lindsey_Z` in this pass. Add a
compact audit label so downstream consumers do not miss the mismatch.

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
- If the generator-only change is not net-negative, do not hunt unrelated old
  tests. Instead simplify stale generator metadata in the same artifact schema.

Implementation surface:

- Edit only:

  ```text
  tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  ```

Required change:

- Add compact comparison-readiness fields to the generated artifact. Prefer a
  shared comparison group such as:

  ```text
  comparison/wl_pqs/final_dimension_match
  comparison/wl_pqs/pqs_final_dimension
  comparison/wl_pqs/wl_final_dimension
  comparison/wl_pqs/comparison_ready
  comparison/wl_pqs/comparison_blocker
  comparison/wl_pqs/comparison_role
  ```

- For the current artifact, readback should be:

  ```text
  pqs_final_dimension = 221
  wl_final_dimension = 2287
  final_dimension_match = false
  comparison_ready = false
  comparison_blocker = :wl_pqs_final_dimension_mismatch
  comparison_role = :separate_route_inspection_not_same_basis_comparison
  ```

- Add compact WL Z-audit fields. Prefer a route-local metadata group:

  ```text
  routes/white_lindsey/metadata/white_lindsey_Z
  routes/white_lindsey/metadata/nuclear_charge_tuple
  routes/white_lindsey/metadata/white_lindsey_Z_audit_status
  routes/white_lindsey/metadata/white_lindsey_Z_audit_blocker
  ```

- For the current artifact, do not claim this is correct. Use an explicit
  review-required label, for example:

  ```text
  white_lindsey_Z_audit_status = :requires_review
  white_lindsey_Z_audit_blocker = :white_lindsey_Z_differs_from_nuclear_charge
  ```

- Keep route-level read-only readiness intact:

  ```text
  routes/pqs_source_box/readiness/cr2_read_only_inspector_ready = true
  routes/white_lindsey/readiness/cr2_read_only_inspector_ready = true
  ```

- Keep solver/export false.
- Do not mark CR2/HF/HFDMRG ready.
- Do not try to find or implement a dimension-matched WL route in this pass.

Line-budget guidance:

- Avoid adding a field cloud. Keep the comparison group compact.
- If needed for net-negative budget, remove stale per-route readiness aliases
  from the generator that duplicate `cr2_solver_ready=false` /
  `cr2_export_ready=false` and are not used by the current CR2 inspector. Do
  not remove the primary `cr2_*` readiness keys.

Validation/readback:

Run:

```text
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
julia --project=. -e 'using JLD2; p="tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2"; jldopen(p,"r") do f; println("pqs_dim=", f["comparison/wl_pqs/pqs_final_dimension"]); println("wl_dim=", f["comparison/wl_pqs/wl_final_dimension"]); println("dim_match=", f["comparison/wl_pqs/final_dimension_match"]); println("comparison_ready=", f["comparison/wl_pqs/comparison_ready"]); println("comparison_blocker=", f["comparison/wl_pqs/comparison_blocker"]); println("comparison_role=", f["comparison/wl_pqs/comparison_role"]); println("wl_Z=", f["routes/white_lindsey/metadata/white_lindsey_Z"]); println("wl_Z_audit_status=", f["routes/white_lindsey/metadata/white_lindsey_Z_audit_status"]); println("wl_Z_audit_blocker=", f["routes/white_lindsey/metadata/white_lindsey_Z_audit_blocker"]); println("pqs_readonly=", f["routes/pqs_source_box/readiness/cr2_read_only_inspector_ready"]); println("wl_readonly=", f["routes/white_lindsey/readiness/cr2_read_only_inspector_ready"]); println("pqs_solver_ready=", f["routes/pqs_source_box/readiness/cr2_solver_ready"]); println("wl_solver_ready=", f["routes/white_lindsey/readiness/cr2_solver_ready"]); end'
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch --ignored tmp/work/be2_wl_pqs_cr2_inspection_artifact
```

Hard boundaries:

- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense four-index `Vee`,
  `V6`, `Vblocks`, or solver code.
- Do not mark the artifact solver-ready.
- Do not claim WL/PQS comparison-ready while dimensions differ.
- Do not change fixture size, route construction, H1/Vee values, or physics.
- Do not add final-basis overlap matrices downstream.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Report:

- files changed;
- line-budget added/deleted/net result;
- comparison-readiness readback;
- WL Z audit readback;
- H1/Vee dimensions and route read-only/solver flags;
- artifact size and generator timing if available from the run;
- generated outputs ignored status;
- validation commands/results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
