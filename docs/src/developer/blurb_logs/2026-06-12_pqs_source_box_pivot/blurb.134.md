Pass 134 - audit Be2 PQS complete-core/shell seam

Purpose:

Do a no-edit seam audit for the Be2/diatomic PQS route after the pass 133 Ham
payload fingerprint. The goal is to determine the smallest correct
implementation seam that would let Be2 PQS move from route skeleton/report
metadata toward the same route-owned complete-core/shell source-plan,
final-basis, H1, density-input, and Ham payload path that one-center PQS now
uses.

Why now:

Pass 133 proved the current Be2 PQS route reaches `cartesian_assembly(...)` but
the private Ham payload blocks because the upstream complete-core/shell payloads
are missing:

```text
:pqs_multilayer_shell_region_plan
:pqs_multilayer_shell_source_plan
:pqs_multilayer_complete_core_shell_final_basis
:pqs_multilayer_complete_core_shell_h1_payload
:axis_weights
:raw_pair_factor_terms
:coulomb_expansion
```

Before coding, we need to know whether the existing one-center helper path can
be generalized, whether a new diatomic/PQS source-plan producer is needed, and
where the route-owned object boundary should sit.

Governing constraints:

- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- `docs/src/developer/pqs_source_box_fixture_policy.md`
- `docs/src/developer/manager_doer_collaboration_contract_2026-05-26.md`
- source-box-first PQS remains the algorithmic framing.
- shell/support-row contraction is oracle/debug, not production PQS.
- retained diagnostic/self-integral weights are not IDA/quadrature weights.
- RHF is frozen as a private Hamiltonian validator only; HFDMRG is the serious
  HF/DMRG package.

Audit surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`
  - one-center complete-core/shell source-plan helper
  - one-center final-basis/H1 helper
  - diagnostic route payload and Ham payload handoff
- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- existing Be2/PQS route-driver report tests for route construction parameters
  and retained-rule vocabulary.

Questions to answer:

1. Which exact helper currently creates the one-center
   `pqs_multilayer_shell_region_plan`, source plan, Coulomb expansion, final
   basis, H1 payload, and density-input facts?
2. Which of those helpers are structurally one-center only, and why?
3. Does the Be2 staged route already carry source boxes, retained rules, parent
   axis bundles, center/nuclear metadata, and Coulomb expansion inputs in a
   structured form that a complete-core/shell producer could consume?
4. Is the missing seam best described as:
   - generalize the existing complete-core/shell source-plan producer;
   - add a diatomic/PQS source-plan producer;
   - add a route-owned adapter that maps existing Be2 route objects into the
     complete-core/shell producer;
   - or something else?
5. What is the smallest next implementation pass that advances the Be2 Ham
   payload without adding scalar report-field clouds or promoting debug/oracle
   paths?
6. What test should validate that pass, and what existing test pressure can be
   avoided or eventually shrunk?

Trust boundary:

- No source edits.
- No test edits.
- No commits unless you only write the response handoff.
- No WL payload implementation.
- No public API, exports, artifacts, hfdmrg, or CR2 execution.
- No RHF/SCF work.
- No route-global materialization or physics endpoint promotion.
- Do not synthesize scalar report aliases as a substitute for route-owned
  structured objects.
- Do not request interactive command approval during unattended baton work. If
  an approval would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- Read-only inspection should be enough.
- You may run `rg`, `sed`, `git show`, and other read-only commands.
- Do not run tests unless a specific short focused read-only probe is necessary
  to answer the seam question. If any command is expected to exceed 60 seconds,
  explain why in the response instead of running it by default.

Report back:

- Exact helpers/files inspected.
- Current one-center complete-core/shell payload dataflow.
- Current Be2/PQS available structured objects, if any.
- The true missing seam and why.
- Recommended smallest next implementation pass.
- Recommended focused validation for that pass.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
