Purpose:
  Make the driver consume the compact PQS H1/J seam before making the seam
  larger. Pass 089 produced `pqs_multilayer_complete_core_shell_h1_j_payload(...)`.
  The next step is driver-spine integration/audit, not RHF or fixture tuning.

Task:
  Audit where the compact PQS multilayer H1/J diagnostic payload belongs in
  the existing driver spine:

  ```text
  cartesian_assembly
  cartesian_report
  cartesian_materialization
  ```

  If the change is mechanical and small, add an internal optional
  assembly/report payload that can carry compact H1/J diagnostic status. If it
  is not obviously small, stop with the exact missing route object and the
  proposed integration point.

Current payload to integrate:

  ```text
  pqs_multilayer_complete_core_shell_h1_j_payload(...)
  ```

Relevant files to inspect first:

  ```text
  src/pqs_source_box_route_driver_helpers.jl
  src/pqs_source_box_route_driver_skeletons.jl
  src/pqs_multilayer_complete_core_shell_h1.jl
  docs/src/developer/pqs_near_term_final_basis_realization_plan.md
  test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
  test/nested/cartesian_report_stage_low_order_policy_runtests.jl
  ```

Guardrails:
  - Design/audit first. Do not implement unless the insertion point is clear
    and the patch stays small.
  - No RHF, SCF, density-matrix iteration, Fock construction, GTO, exports,
    artifacts, fixture-rule policy, side-13 rerun, q ladder, or acceptance
    promotion.
  - Do not extend the H1/J helper into an SCF route.
  - Do not use WL or old fixed-block data as active authority.
  - Keep `driver_route_materialized = false` unless a real driver-owned route
    object is created by this pass. A report-stage diagnostic hook is not a
    production route claim.
  - Do not add broad tests. Add or update a compact module/driver contract
    check only if it replaces or shrinks older probe/report glue.
  - No UI escalation; write `.agent_handoffs/ATTENTION.md` and stop if blocked.

Expected report:
  - where H1/J should enter the driver spine;
  - whether you implemented an internal optional assembly/report payload;
  - if not implemented, the precise missing route object;
  - what remains private/diagnostic;
  - whether an RHF contract is the next design boundary;
  - validation run;
  - deletion/shrinkage report:
      - what old test/probe/report glue became less necessary;
      - what was deleted or simplified;
      - if nothing was deleted, why no existing surface was made obsolete;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage.

Validation:
  - If audit/docs only: `git diff --check`.
  - If code changes: focused driver/module test touched by the change, focused
    H1 gate, load check, and `git diff --check`.

-- repo-manager@macmini
