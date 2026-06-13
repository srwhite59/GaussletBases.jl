Pass 175 - shrink WL handoff to HF-minimal one-body payload

Role: repo-doer@macmini

Task type: generator artifact-size/read-speed optimization.

Purpose:

Pass 174 showed the generator runtime is dominated by route construction:

```text
generator.phase.pqs_route_build=44.31484975
generator.phase.wl_route_build=186.2394435
generator.phase.bundle_write=0.786514
generator.phase.total=232.239269167
```

The next useful optimization before another CR2 handoff is to shrink the WL
payload CR2 has to read. CR2’s minimal HF-facing contract needs final-basis
`H1` and final-basis density-density `Vee`. The WL artifact still carries dense
one-body split sidecars (`kinetic_one_body` and per-center nuclear matrices)
that are not needed for the first HF smoke.

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
- If the generator-only change somehow is not net-negative, the next concrete
  stale candidate is:

  ```text
  test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl
  ```

Implementation surface:

- Prefer editing only:

  ```text
  tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  ```

Required change:

- Stop writing these dense WL sidecar matrices to the generated JLD2 payload:

  ```text
  routes/white_lindsey/one_body/kinetic_one_body
  routes/white_lindsey/one_body/nuclear_one_body_by_center
  ```

- Keep the final-basis one-body Hamiltonian:

  ```text
  routes/white_lindsey/one_body/hamiltonian
  ```

- Keep compact split metadata/status, for example:

  ```text
  nuclear_term_storage
  kinetic_one_body_stored = false
  nuclear_one_body_by_center_count = <existing count>
  nuclear_one_body_by_center_stored = false
  one_body_split_storage = :not_stored_hf_minimal_handoff
  ```

- Do not change route construction, H1, Vee, overlap convention, solver flags,
  fixture size, or physics values.
- Do not remove PQS provenance/debug fields in this pass.

Validation/readback:

Run:

```text
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
julia --project=. -e 'using JLD2; p="tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2"; jldopen(p,"r") do f; ks=Set(String.(keys(f["routes/white_lindsey/one_body"]))); println("producer_commit=", f["producer/repo_commit"]); println("producer_dirty=", f["producer/dirty"]); println("wl_h1=", size(f["routes/white_lindsey/one_body/hamiltonian"])); println("wl_vee=", size(f["routes/white_lindsey/two_body/interaction_matrix"])); println("has_wl_kinetic=", "kinetic_one_body" in ks); println("has_wl_nuclear_by_center=", "nuclear_one_body_by_center" in ks); println("wl_kinetic_stored=", f["routes/white_lindsey/one_body/kinetic_one_body_stored"]); println("wl_nuclear_by_center_stored=", f["routes/white_lindsey/one_body/nuclear_one_body_by_center_stored"]); println("wl_solver_ready=", f["routes/white_lindsey/readiness/cr2_solver_ready"]); end'
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch --ignored tmp/work/be2_wl_pqs_cr2_inspection_artifact
```

Also report bundle size before/after if available from `ls -lh`.

Hard boundaries:

- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense four-index `Vee`,
  `V6`, `Vblocks`, or solver code.
- Do not mark the artifact solver-ready.
- Do not store PQS or WL final-basis overlap matrices downstream.
- Do not alter the density-density HF convention blocker.
- Do not add another tracked timing report.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Report:

- files changed/deleted;
- line-budget added/deleted/net result;
- artifact size before/after;
- generator timing output by top-level phase;
- readback for absent WL split matrices and retained H1/Vee;
- generated outputs ignored status;
- validation commands/results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
