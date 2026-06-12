Purpose:
  Correct pass 082 before acceptance. The H1 payload helper is functionally
  useful, but its returned object is too broad and carries consumed inputs /
  support-space intermediates.

Issue:
  `pqs_multilayer_complete_core_shell_h1_payload(...)` currently returns:

  ```text
  source_plan = plan
  final_basis
  support_kinetic
  support_nuclear_by_center
  final_kinetic
  final_nuclear_by_center
  final_hamiltonian
  h1
  summary
  ```

  The requested payload should expose final H1 assembly products plus compact
  summary/nonclaims. Returning source/final-basis inputs and support-space
  intermediates recreates the broad report-shape problem we have been trying to
  avoid.

Task:
  Trim the H1 payload return shape.

Keep in the returned payload:
  - `final_kinetic`;
  - `final_nuclear_by_center`;
  - `final_hamiltonian`;
  - `h1`;
  - `summary`;
  - compact metadata/nonclaims.

Remove from the returned payload unless a concrete live caller needs it:
  - `source_plan`;
  - `final_basis`;
  - `support_kinetic`;
  - `support_nuclear_by_center`.

Test adjustment:
  Keep the small axis-layer/origin-factor nuclear convention comparison, but do
  not broaden the H1 payload to support it. If the comparison needs the explicit
  support nuclear matrix, build it in that convention-check section with
  `pqs_multilayer_support_electron_nuclear_by_center_matrices(...)`, separate
  from the active H1 payload.

Do not:
  - add new features or physics values;
  - add IDA, density-density, RHF, GTO, driver wiring, exports, artifacts, or
    fixture policy;
  - add tests or broad metadata assertions;
  - request UI escalation. In unattended baton mode, write
    `.agent_handoffs/ATTENTION.md` and stop if permission is genuinely needed.

Validation:
  - focused H1 gate;
  - load check;
  - `git diff --check`.

Report:
  - final payload fields after trimming;
  - how the convention comparison gets its support nuclear matrix;
  - validation run;
  - deletion/shrinkage report.

-- repo-manager@macmini
