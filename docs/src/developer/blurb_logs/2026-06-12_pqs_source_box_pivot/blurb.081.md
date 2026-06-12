Purpose:
  Decide the next PQS integration seam against the canonical Cartesian driver
  spine before adding another private helper. The H1 gate now uses
  shellification/lowering-backed region planning, but the route is still
  manually assembled inside `test/nested/pqs_direct_retained_final_h1_runtests.jl`.

Task:
  Do a focused audit/design pass, no source implementation unless the answer is
  truly mechanical and tiny.

  Map the current PQS H1 gate pieces onto the driver stages:

  ```text
  cartesian_system / recipe
  cartesian_parent
  cartesian_shells
  cartesian_units
  cartesian_transforms
  cartesian_pairs
  cartesian_assembly
  cartesian_report / materialization
  ```

  Identify the smallest next module-owned seam that would make the current H1
  fixture less test-local. Candidate shape might be a compact assembly payload
  that consumes:

  ```text
  pqs_multilayer_shell_region_plan
  pqs_multilayer_shell_source_plan(bundles, region_plan)
  pqs_multilayer_complete_core_shell_final_basis
  pqs_multilayer_support_kinetic_matrix
  pqs_multilayer_support_electron_nuclear_by_center_matrices
  CartesianFinalBasisRealization final one-body/H1 helpers
  ```

  But do not implement that candidate in this pass unless it is obviously
  smaller than the audit text.

Guardrails:
  - The output should not be a broad PQS route driver.
  - Do not add physics fixture rules for `Z`, `d`, `s`, radius, q, or shell
    depth.
  - Do not add IDA, density-density, RHF, GTO, exports, artifacts, or
    acceptance gates.
  - Keep dense support-space one-body helpers scoped to the H1 seam; do not
    promote them as the scalable PQS operator algorithm.
  - Do not add tests.
  - No UI escalation; write `.agent_handoffs/ATTENTION.md` and stop if blocked.

Expected output:
  - a short boundary recommendation;
  - source files/functions inspected;
  - exact next implementation candidate, if any;
  - what old test/probe glue would shrink if that candidate lands;
  - whether a design decision from the user/manager is needed before coding.

Validation:
  - docs-only/audit: `git diff --check`.

Report:
  - driver-stage mapping;
  - recommended next pass or explicit blocker;
  - deletion/shrinkage report as usual.

-- repo-manager@macmini
