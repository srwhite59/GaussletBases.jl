Pass 217 - implement the physical H2 private RHF execution adapter

Role:
You are `repo-doer@macmini` implementing one bounded private-RHF execution
adapter for GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`,
`BlurbStyle.md`, and the unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `ec18a623 Block H2 physical RHF execution adapter`
- The H2 physical q5 gausslet-only route has:
  - source plan available;
  - final basis available;
  - H1 available;
  - density/H1-J seam available;
  - private RHF input contract available;
  - explicit `private_rhf_electron_count = 2`;
  - private RHF execution payload blocked on
    `:missing_physical_gausslet_rhf_execution_adapter`.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. H2 221 remains
source-box diagnostic only. He 419 is the already validated atom endpoint.

Purpose:
Remove the blocker named in pass 216 by implementing the physical-H2 private
RHF execution adapter. This should attempt private RHF execution on the 463
physical H2 route using reviewed existing RHF math, not a new solver product.

Task:
Make `_pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload`
actually attempt execution when the physical-H2 RHF input contract is available.

Allowed approaches:
1. Preferred: refactor the existing private RHF helper accessors/predicates in
   `src/pqs_multilayer_complete_core_shell_rhf.jl` so they accept both the old
   complete-core/shell object family and the physical-H2 object family by
   explicit reviewed object-kind/status predicates.
2. Acceptable: add a small physical-H2 wrapper that reuses the same lower math
   for:
   - H1-Aufbau initial density;
   - one-step Fock/energy from final density;
   - SCF iteration and final recomputed diagnostics.

Do not fake old object kinds and do not duplicate the SCF algorithm wholesale.
If you need more than a small wrapper, stop with a precise blocker explaining
which old helper must be generalized.

Required execution contract:
- closed-shell RHF only;
- electron count comes only from explicit `private_rhf_electron_count`;
- nocc = 1 for this H2 fixture;
- initial density may be H1-Aufbau;
- use existing private SCF controls only;
- default controls may be fixed-point or existing Fock-DIIS, but do not tune DIIS
  in this pass;
- final returned density and final one-step/Fock diagnostics must correspond to
  the same density;
- convergence must be gated on final recomputed diagnostics, not just
  iteration-input checks;
- if SCF does not converge, return a blocked/nonconverged payload with residual
  facts instead of changing controls ad hoc.

Expected success transition if converged:

```text
private_rhf/executed = true
private_rhf/materialized = true
private_rhf/converged = true
private_rhf/electron_count = 2
private_rhf/occupation_nocc = 1
private_rhf/total_energy = finite
private_rhf/iteration_count > 0
private_rhf/final_density_one_step_consistency_status = reviewed/recomputed
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_h2_gausslet_only_reference_comparison
comparison/ready = false
```

Expected blocked/nonconverged transition if execution runs but fails to
converge:

```text
private_rhf/executed = true
private_rhf/materialized = false
private_rhf/converged = false
physics/endpoint_blocker = :physical_gausslet_private_rhf_not_converged
```

Do not:
- add a new SCF/DIIS strategy;
- run HFDMRG, DMRG, or CR2;
- add GTO/MWG supplement;
- compare to supplemented WL/QW H2 references;
- mark solver/export/public/HamV6 ready;
- pass final-basis self-overlap as downstream working data;
- mutate or reuse the H2 221 diagnostic route;
- add broad/default-runner tests.

Artifact/test expectations:
- Extend/update compact private-RHF fields only. Do not add matrices.
- The H2 endpoint test should assert the compact execution/convergence or
  nonconvergence status and a finite total energy only if materialized.
- Preserve the active He and H2 endpoint tests.

Line-count rule:
The active source/test/bin line-count rule still applies. For edits under
`src`, `test`, `bin`, and the CR2 generator script, final tracked diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs and blurb logs do not count. If source lines are added, pay for them by
deleting stale source/test pressure in the same pass.

Line-budget deletion candidates to inspect first:
- `test/nested/cartesian_ham_builder_one_center_config_smoke_runtests.jl`
  may be superseded by the He driver endpoint tests. Delete/shrink only if that
  is true.
- `test/nested/cartesian_pair_block_route_global_matrix_set_smoke_runtests.jl`
  and `test/nested/cartesian_pair_block_route_global_one_body_adapter_runtests.jl`
  are likely old global-overlap/one-body adapter pressure after the private
  global-overlap driver hook retirement. Delete/shrink only if live source
  callers or active tests do not require them.
- Do not delete the active H2 physical endpoint test or the He endpoint tests.

If no safe deletion is available, report exact live callers and stop rather
than adding net-positive source/test code.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 463 driver-artifact test with Julia-level timing:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

This test was about 79-81 seconds before actual SCF execution. If SCF is
attempted, report the SCF execution time separately if practical. If the test
becomes substantially slower, say why.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.217.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.217.md
```

Report:
- whether private RHF execution was attempted;
- whether it converged;
- statuses/blocker;
- electron count source and occupation policy;
- total energy and energy components if materialized;
- iteration count and residual diagnostics if available;
- final-density/one-step consistency status;
- what helper was generalized or wrapped;
- artifact fields changed;
- source/test/bin scoped line budget added/deleted/net;
- validation commands and elapsed times;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini
