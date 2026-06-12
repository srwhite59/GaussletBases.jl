Purpose:
  Continue the PQS file-boundary cleanup before adding more physics. Pass 085
  separated shellification/lowering-backed region planning; the remaining
  `src/pqs_multilayer_shell_source_plan.jl` still starts with dense
  support-space one-body helpers that are conceptually H1 seam machinery, not
  source-plan realization.

Task:
  Do one same-module mechanical extraction of the support one-body helper
  section into a new file, if it remains low risk.

Candidate new file:

  `src/pqs_multilayer_support_one_body.jl`

Candidate functions to move:

  ```text
  _pqs_multilayer_support_product_matrix
  _pqs_multilayer_axis_tuple
  _pqs_multilayer_center_records
  _pqs_multilayer_center_property
  _pqs_multilayer_center_summary
  _pqs_multilayer_explicit_factor_terms
  _pqs_multilayer_centered_factor_terms
  _pqs_multilayer_term_first_factor_array
  _pqs_multilayer_validate_factor_terms
  _pqs_multilayer_support_electron_nuclear_matrix
  pqs_multilayer_support_kinetic_matrix
  pqs_multilayer_support_electron_nuclear_by_center_matrices
  ```

Expected include order if implemented:

  ```julia
  include("pqs_multilayer_shell_region_plan.jl")
  include("pqs_multilayer_support_one_body.jl")
  include("pqs_multilayer_shell_source_plan.jl")
  ```

  Adjust only if live dependencies require a different order. Preserve function
  names and return shapes.

Guardrails:
  - Same module only; do not create a new module.
  - No new physics, H1/J, density, RHF, driver wiring, exports, artifacts, or
    fixture-rule policy.
  - Do not redesign result types.
  - Do not add tests.
  - Do not optimize dense support-space algorithms in this pass; this is a file
    boundary pass only.
  - Keep support-space one-body helpers explicitly scoped to H1 seam/reference
    machinery, not general PQS operator assembly.
  - No UI escalation; write `.agent_handoffs/ATTENTION.md` and stop if blocked.

Stop condition:
  If moving both support public helpers creates circular include pressure with
  the source-plan type/fields, stop after moving only the private helper section
  or after writing the dependency audit. Do not force a tangled extraction.

Validation:
  - If code is mechanically split: focused H1 gate, load check, and
    `git diff --check`.
  - If audit only: `git diff --check`.

Report:
  - exactly what moved and what stayed;
  - include order used;
  - validation run;
  - deletion/shrinkage report:
      - what conceptual responsibility was separated;
      - what remains in `pqs_multilayer_shell_source_plan.jl`;
      - whether this makes the source-plan file less of a private route island;
      - remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini
