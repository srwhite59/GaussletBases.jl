Purpose:
  Prepare the next PQS physics hardening step after H1: a route-owned
  pre-RHF self-Coulomb/J diagnostic analogous to the WL H1/J checks. Do not
  implement density code until the existing helper/convention surface is clear.

Task:
  Do a focused audit/design pass of the existing PQS complete core/shell IDA
  and density helpers.

Inspect at least:
  - `CartesianFinalBasisRealization.pqs_complete_core_shell_final_ida_weights`;
  - `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction`;
  - `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_orbital_self_coulomb`;
  - any current `tmp/work` or docs references to PQS side13 H1/J/RHF probes;
  - current support/raw pair numerator construction, if it exists in tracked
    code.

Answer:
  1. What inputs are still needed to build a route-owned PQS H1/J payload?
     Be specific about support weights, raw pair numerator, final/pre-final
     coefficients, and the orbital coefficients used for J.
  2. Which helper already owns the corrected convention?
     Confirm whether it is the pre-final positive-weight gauge described in
     the docs, not signed-final-weight division and not raw no-division.
  3. What remains probe-local or missing?
  4. What is the smallest next implementation pass?
     Candidate shape might be an H1/J payload that consumes the H1 payload,
     final basis, support weights, and raw pair numerator, then reports H1 and
     self-Coulomb for the lowest H1 orbital.
  5. What test/probe glue would shrink if that lands?

Guardrails:
  - No source implementation unless the audit finds a tiny docs-only correction.
  - Do not add RHF, GTO, driver wiring, exports, artifacts, or fixture-rule
    policy.
  - Do not add tests.
  - Keep density convention language precise: raw numerator first, then the
    documented pre-final positive-weight density boundary.
  - No UI escalation; write `.agent_handoffs/ATTENTION.md` and stop if blocked.

Validation:
  - docs-only/audit: `git diff --check`.

Report:
  - helper/convention map;
  - exact missing inputs for H1/J;
  - recommended next pass or blocker;
  - deletion/shrinkage report.

-- repo-manager@macmini
