Pass 216 - execute private physical H2 RHF diagnostic, or block precisely

Role:
You are `repo-doer@macmini` implementing one bounded private-RHF execution pass
for GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and
the unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `86c3b793 Add H2 physical RHF input contract`
- The H2 physical q5 gausslet-only driver path now has:
  - source plan available;
  - final basis available;
  - H1 available;
  - density/H1-J seam available;
  - private RHF input contract available;
  - explicit route/private-RHF electron count `2`;
  - final dimension `463`;
  - H1 lowest energy about `-0.7946609179724647`;
  - H1-J/self-Coulomb scalar about `0.45696639804337114`;
  - density gauge `:pre_final_localized_positive_weight`;
  - raw pair-factor convention `:raw_numerator`.
- The current endpoint blocker is:
  `:missing_physical_gausslet_private_rhf_execution`.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. This remains the
active PQS physics target. H2 221 remains source-box diagnostic only. He 419 is
the already validated atom endpoint.

Purpose:
Use the explicit physical-H2 RHF input contract to run the private RHF
diagnostic if the existing private RHF math can be adapted cleanly. This is a
Hamiltonian validator, not a product solver.

Default expected shape:
- Add a private physical-H2 RHF execution payload that consumes:

```text
_PQSDiatomicPhysicalGaussletRHFInputContractPayload
_PQSDiatomicPhysicalGaussletH1Payload
_PQSDiatomicPhysicalGaussletH1JPayload
```

- It may reuse existing private complete-core/shell RHF math if you add explicit
  accessors/predicates for the physical payloads. Do not fake old object kinds.
- It should produce either:
  - a converged private RHF diagnostic payload; or
  - a blocked/nonconverged payload with iteration/residual facts.

Required SCF/convergence contract if execution is attempted:
- closed-shell RHF only;
- electron count comes only from explicit `private_rhf_electron_count`;
- nocc = 1 for this H2 fixture;
- initial density may be H1-Aufbau;
- use existing private SCF controls only; do not invent a new SCF strategy;
- no DIIS tuning pass here;
- final returned density and final one-step/Fock diagnostics must correspond to
  the same density;
- convergence must be gated on final recomputed diagnostics, not just
  iteration-input checks;
- report density trace, idempotency residual, commutator residual if available,
  energy delta, and iteration count.

Expected success transition if converged:

```text
private_rhf/executed = true
private_rhf/materialized = true
private_rhf/converged = true
private_rhf/electron_count = 2
private_rhf/occupation_nocc = 1
private_rhf/total_energy = finite
private_rhf/commutator_residual = small
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_h2_gausslet_only_reference_comparison
comparison/ready = false
```

If execution runs but does not converge, use:

```text
private_rhf/executed = true
private_rhf/materialized = false
private_rhf/converged = false
physics/endpoint_blocker = :physical_gausslet_private_rhf_not_converged
```

If execution cannot be attempted because an adapter/convention is missing, use a
precise blocker such as:

```text
:missing_physical_gausslet_rhf_execution_adapter
:physical_gausslet_rhf_fock_convention_unreviewed
:missing_physical_gausslet_final_density_initialization
:missing_physical_gausslet_final_recomputed_diagnostics
```

Do not:
- add or tune a new SCF/DIIS algorithm;
- run HFDMRG, DMRG, or CR2;
- add GTO/MWG supplement;
- compare to supplemented WL/QW H2 references;
- mark solver/export/public/HamV6 ready;
- pass along final-basis self-overlap as downstream working data;
- mutate or reuse the H2 221 diagnostic route;
- add broad/default-runner tests.

Artifact/test expectations:
- Extend the H2 physical driver artifact with compact private-RHF execution
  fields only. Do not add large matrices or a report-field cloud.
- Keep the final-basis overlap as an identity-error diagnostic only.
- Update
  `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
  to assert the compact execution/convergence or nonconvergence fields.
- Preserve the active He and H2 endpoint tests.

Line-count rule:
The active source/test/bin line-count rule still applies. For edits under
`src`, `test`, `bin`, and the CR2 generator script, final tracked diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs and blurb logs do not count. If source lines are added, pay for them by
deleting stale source/test pressure in the same pass. Do not delete the active
He or H2 driver endpoint tests.

Line-budget deletion candidates to inspect first:
- `test/nested/cartesian_ham_builder_one_center_config_smoke_runtests.jl`
  may now be superseded by the He driver endpoint tests. Delete or shrink only
  if you verify the He endpoint tests cover the live one-center driver behavior.
- `test/nested/pqs_explicit_core_spacing_parent_axis_probe_runtests.jl`
  and `test/nested/pqs_raw_product_box_plan_probe_runtests.jl` are old probe
  scaffolds. Shrink/delete only if their active contracts are covered by H2
  endpoint or smaller module tests.

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

This test was about 80 seconds without package precompile after pass 215; a
manager rerun that included precompile printed about 139 seconds. If this pass
adds SCF execution, report the SCF phase time separately if practical.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.216.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.216.md
```

Report:
- whether private RHF execution was attempted;
- whether it converged;
- statuses/blocker;
- electron count source and occupation policy;
- total energy and energy components if materialized;
- iteration count and residual diagnostics if available;
- final-density/one-step consistency status;
- whether existing RHF helpers were adapted or left untouched;
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
