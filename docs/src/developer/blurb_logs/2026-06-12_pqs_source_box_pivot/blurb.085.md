Purpose:
  Prevent `src/pqs_multilayer_shell_source_plan.jl` from becoming a new private
  route island. Recent PQS direction corrections succeeded, but too many seams
  now live in one 1,100-line file.

Task:
  Do a focused file-boundary cleanup audit, with optional mechanical extraction
  only if it is obviously low-risk.

Current mixed responsibilities in `src/pqs_multilayer_shell_source_plan.jl`:

  ```text
  shellification/lowering-backed region-plan records and constructor
  explicit-box bridge layer specs
  common PQS source realization
  dense support-space product/kinetic/nuclear helpers
  complete core/shell final-basis adapter
  H1 assembly payload
  ```

Expected audit output:
  1. Propose a small file split under `src/`, using existing module/include
     style. Example shape:
       - `pqs_multilayer_shell_region_plan.jl`
       - `pqs_multilayer_shell_source_realization.jl`
       - `pqs_multilayer_support_one_body.jl`
       - `pqs_multilayer_h1_payload.jl`
     Adjust names to match repo style.
  2. Identify dependencies/order of includes.
  3. Identify which extraction, if any, can be done mechanically now without
     changing behavior.
  4. If you do implement a mechanical split, do only one small coherent
     extraction and preserve function names/return shapes.

Guardrails:
  - No new physics, H1/J, density, RHF, GTO, driver wiring, exports, artifacts,
    or fixture-rule policy.
  - Do not redesign result types in this pass.
  - Do not add tests.
  - Do not move code across modules unless the dependency boundary is obvious;
    a same-module file split is preferred first.
  - Keep support-space dense one-body helpers explicitly scoped to H1 seam
    machinery.
  - No UI escalation; write `.agent_handoffs/ATTENTION.md` and stop if blocked.

Validation:
  - If docs/audit only: `git diff --check`.
  - If code is mechanically split: focused H1 gate, load check, and
    `git diff --check`.

Report:
  - recommended file split and include order;
  - whether you implemented any mechanical extraction;
  - validation run;
  - deletion/shrinkage report:
      - what conceptual responsibility was separated;
      - what remains in the original file;
      - what old/private route island pressure remains.

-- repo-manager@macmini
