Pass 174 - instrument CR2 generator top-level phases and delete stale placement-facts scaffold

Role: repo-doer@macmini

Task type: generator optimization groundwork plus test carrying-cost deletion.

Purpose:

Pass 173 shrank the Be2 WL/PQS CR2 artifact from `242M` to `202M`, but the
generator still took about `240 s` in a fresh Julia process. The route-local
printed timings are only a few seconds, so the next generator pass should expose
the uninstrumented top-level cost before changing numerical construction.

At the same time, keep the hard line-budget pressure real. The old CPB overlap
placement-facts test is metadata-only development scaffolding and is not
referenced by the default runner or source. Retire it in this pass rather than
adding instrumentation for free.

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

Implementation surface:

- Edit the tracked generator:

  ```text
  tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  ```

- Delete this stale non-default scaffold unless your inspection finds a live
  source/default-runner caller:

  ```text
  test/nested/cartesian_cpb_overlap_placement_facts_runtests.jl
  ```

- Do not delete `test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl`
  in this pass unless the line budget somehow still fails.

Required generator change:

- Add coarse top-level timing prints around these phases:

  ```text
  generator.phase.pqs_route_build
  generator.phase.pqs_payload
  generator.phase.provenance
  generator.phase.wl_route_build
  generator.phase.wl_population
  generator.phase.bundle_write
  generator.phase.total
  ```

- Use Julia-level timing such as `@elapsed`.
- Keep the output compact and human-readable.
- Do not add a timing dependency, logging framework, docs, or another tracked
  report file.
- Do not change fixture size, route semantics, H1/Vee values, overlap policy,
  solver flags, or artifact schema except adding timing lines to stdout.

Validation:

Run:

```text
rg -n "cartesian_cpb_overlap_placement_facts_runtests" test src docs/src/developer
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
julia --project=. -e 'using JLD2; p="tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2"; jldopen(p,"r") do f; println("producer_commit=", f["producer/repo_commit"]); println("producer_dirty=", f["producer/dirty"]); println("pqs_h1=", size(f["routes/pqs_source_box/one_body/hamiltonian"])); println("pqs_vee=", size(f["routes/pqs_source_box/two_body/interaction_matrix"])); println("wl_h1=", size(f["routes/white_lindsey/one_body/hamiltonian"])); println("wl_vee=", size(f["routes/white_lindsey/two_body/interaction_matrix"])); println("pqs_solver_ready=", f["routes/pqs_source_box/readiness/cr2_solver_ready"]); println("wl_solver_ready=", f["routes/white_lindsey/readiness/cr2_solver_ready"]); end'
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch --ignored tmp/work/be2_wl_pqs_cr2_inspection_artifact
```

Hard boundaries:

- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense four-index `Vee`,
  `V6`, `Vblocks`, or solver code.
- Do not mark the artifact solver-ready.
- Do not store PQS or WL final-basis overlap matrices downstream.
- Do not add another tracked timing report.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Report:

- files changed/deleted;
- line-budget added/deleted/net result;
- generator timing output by top-level phase;
- artifact readback for commit/dirty, H1/Vee shapes, solver flags;
- generated outputs ignored status;
- validation commands/results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
