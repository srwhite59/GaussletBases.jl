Purpose:
  Finish the immediate PQS multi-layer file-boundary cleanup. Passes 085 and
  086 separated region planning and dense support one-body helpers. The source
  realization file still owns complete core/shell final-basis adaptation and
  the H1 payload, which are later-stage H1 seam concepts.

Task:
  Do one same-module mechanical extraction of the final-basis/H1 seam from
  `src/pqs_multilayer_shell_source_plan.jl`, if it remains low risk.

Candidate new file:

  `src/pqs_multilayer_complete_core_shell_h1.jl`

Candidate functions to move:

  ```text
  pqs_multilayer_complete_core_shell_final_basis
  pqs_multilayer_complete_core_shell_h1_payload
  _blocked_pqs_multilayer_complete_core_shell_final_basis
  ```

Expected include order if implemented:

  ```julia
  include("pqs_multilayer_shell_region_plan.jl")
  include("pqs_multilayer_support_one_body.jl")
  include("pqs_multilayer_shell_source_plan.jl")
  include("pqs_multilayer_complete_core_shell_h1.jl")
  ```

  Adjust only if live dependencies require it. Preserve function names and
  return shapes.

Guardrails:
  - Same module only; do not create a new module.
  - No new physics, H1/J diagnostics, density, RHF, driver wiring, exports,
    artifacts, or fixture-rule policy.
  - Do not redesign result types.
  - Do not add tests.
  - Do not optimize or rewrite dense support-space algorithms.
  - Keep this as file-boundary cleanup only.
  - No UI escalation; write `.agent_handoffs/ATTENTION.md` and stop if blocked.

Stop condition:
  If this extraction reveals that final-basis/H1 helpers should instead move to
  an existing module such as `CartesianFinalBasisRealization`, do not make that
  conceptual move in this pass. Stop with the dependency/audit finding, because
  that would be a design decision rather than a mechanical split.

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
      - whether `pqs_multilayer_shell_source_plan.jl` now mostly owns source
        realization plus explicit-box bridge only;
      - whether any later module-boundary decision remains;
      - remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini
