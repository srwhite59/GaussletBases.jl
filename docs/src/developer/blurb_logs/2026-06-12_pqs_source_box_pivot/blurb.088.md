Purpose:
  Move the next PQS H1/J ingredient out of tmp/work probe-local code without
  promoting RHF, acceptance, or fixture tuning. The near-term plan says the
  remaining route-owned gap for the H1/J density seam is the producer for
  support weights, the support raw pair numerator, and H1 orbital coefficients.
  Start with the first two.

Task:
  Add narrow route-owned support-density input helpers for a complete
  core/shell multi-layer PQS source plan.

Suggested file:

  `src/pqs_multilayer_support_density.jl`

Suggested public/private surface:

  ```text
  pqs_multilayer_support_weights(plan; axis_weights)
  pqs_multilayer_support_pair_raw_numerator_matrix(plan; raw_pair_factor_terms, coulomb_expansion)
  ```

  Adjust names only if the repo style strongly suggests better ones. These
  helpers should use the same support ordering as the H1 seam:

  ```text
  core_support_states followed by shell_support_states
  ```

  They should support the current one-center/common-axis case cleanly. If
  accepting `x/y/z` axis-specific weights or raw pair-factor terms is cheap and
  obvious, do it; otherwise keep the accepted input shape explicit.

Convention:
  - `support_weights` are product support weights, e.g.
    `wx[ix] * wy[iy] * wz[iz]`.
  - `support_pair_raw_numerator_matrix` contracts
    `pair_factor_terms_raw` with positive Coulomb expansion coefficients.
  - Do not divide by weights in the raw numerator helper.
  - Do not use density-normalized `pair_factor_terms` as authority.
  - Do not materialize final density interaction inside these helpers.

Validation:
  - Update or copy the existing `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl`
    so it consumes the new route-owned helpers instead of its probe-local
    `_support_weights` and `_support_pair_raw_numerator_matrix` implementations.
  - Run that tmp/work probe and report the H1/J diagnostic values.
  - Run the focused H1 gate, load check, and `git diff --check`.

Guardrails:
  - No RHF, no acceptance gate, no permanent test unless you can delete or
    shrink older coverage in the same pass.
  - No fixture tuning, q ladder, side-13 rerun, GTO, exports, artifacts, or
    driver wiring.
  - No broad metadata/report expansion.
  - No signed-final-weight density rule and no raw-no-division final density
    rule.
  - Keep this as input production for the reviewed pre-final positive-weight
    density interaction.
  - No UI escalation; write `.agent_handoffs/ATTENTION.md` and stop if blocked.

Report:
  - functions/files added or edited;
  - whether tmp/work probe-local duplicate helpers were bypassed or deleted;
  - H1/J probe result and whether it still matches the previous summary;
  - validation run;
  - deletion/shrinkage report:
      - what tmp/work or route-local duplication became unnecessary;
      - what remains probe-local, especially H1 orbital coefficient extraction;
      - why no permanent test was added, or what old coverage was shrunk if one
        was added.

-- repo-manager@macmini
