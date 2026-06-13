Pass 176 - instrument WL generator subphases and delete stale placement pilot

Role: repo-doer@macmini

Task type: generator runtime optimization groundwork plus stale test deletion.

Purpose:

Pass 175 shrank the generated Be2 WL/PQS CR2 bundle from `202M` to `82M`, but
the generator still takes about `233 s`. Top-level timing shows:

```text
generator.phase.pqs_route_build=45.027609625
generator.phase.wl_route_build=186.269930834
generator.phase.wl_population=0.169682791
generator.phase.bundle_write=0.687189083
generator.phase.total=232.861505792
```

The bottleneck is now specifically `wl_route_build`. Before changing route
construction, instrument that helper's internal phases so the next pass can
remove or bypass the right work.

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

- Delete the remaining stale non-default CPB overlap placement scaffold unless
  your inspection finds a live source/default-runner caller:

  ```text
  test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl
  ```

Required generator change:

- Add compact timing prints inside `_be2_white_lindsey_atom_growth_route()` for
  the sequential build phases:

  ```text
  generator.phase.wl.system
  generator.phase.wl.recipe
  generator.phase.wl.parent
  generator.phase.wl.shells
  generator.phase.wl.units
  generator.phase.wl.transforms
  generator.phase.wl.pairs
  generator.phase.wl.assembly
  generator.phase.wl.report
  generator.phase.wl.materialization
  ```

- Use Julia-level timing such as `@elapsed`.
- Keep output compact and human-readable.
- Do not add a timing dependency, logging framework, docs, or another tracked
  report file.
- Do not change fixture size, route semantics, H1/Vee values, overlap policy,
  solver flags, or artifact schema.

Validation:

Run:

```text
rg -n "cartesian_cpb_overlap_placement_pilot_runtests" test src docs/src/developer
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
julia --project=. -e 'using JLD2; p="tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2"; jldopen(p,"r") do f; println("producer_commit=", f["producer/repo_commit"]); println("producer_dirty=", f["producer/dirty"]); println("pqs_h1=", size(f["routes/pqs_source_box/one_body/hamiltonian"])); println("pqs_vee=", size(f["routes/pqs_source_box/two_body/interaction_matrix"])); println("wl_h1=", size(f["routes/white_lindsey/one_body/hamiltonian"])); println("wl_vee=", size(f["routes/white_lindsey/two_body/interaction_matrix"])); println("wl_solver_ready=", f["routes/white_lindsey/readiness/cr2_solver_ready"]); end'
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
- Do not change the density-density HF convention blocker.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Report:

- files changed/deleted;
- line-budget added/deleted/net result;
- WL subphase timing output;
- artifact size;
- readback for H1/Vee shapes and solver flags;
- generated outputs ignored status;
- validation commands/results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
