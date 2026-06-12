Purpose:
  Complete the narrow route-owned H1/J diagnostic seam without promoting RHF,
  acceptance, or fixture tuning. Pass 088 moved support weights and raw support
  pair numerator production out of tmp/work. The remaining probe-local accepted
  input is lowest H1 orbital coefficient extraction.

Task:
  Add a compact route-owned H1/J payload helper that consumes the existing
  complete core/shell H1 payload plus final basis and support-density inputs.

Suggested location:

  `src/pqs_multilayer_complete_core_shell_h1.jl`

Suggested function:

  ```text
  pqs_multilayer_complete_core_shell_h1_j_payload(
      plan;
      final_basis,
      h1_payload,
      axis_weights,
      raw_pair_factor_terms,
      coulomb_expansion,
      metadata = (;),
  )
  ```

  Adjust the signature if a clearer repo-local style exists, but keep the
  concept narrow.

Expected behavior:
  - validate the source plan, final basis, and H1 payload;
  - build support weights with `pqs_multilayer_support_weights(...)`;
  - build support raw numerator with
    `pqs_multilayer_support_pair_raw_numerator_matrix(...)`;
  - call
    `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction(...)`;
  - compute the lowest final H1 orbital coefficients from the final Hamiltonian
    already carried by the H1 payload;
  - call
    `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)`;
  - return a compact result containing status, blocker, density interaction,
    self-Coulomb record, and a small summary.

Convention guardrails:
  - Use the reviewed pre-final positive-weight density gauge.
  - Do not implement signed-final-weight division.
  - Do not implement raw-no-division final density.
  - Do not use density-normalized `pair_factor_terms` as authority.
  - Do not add RHF, acceptance, driver wiring, exports, artifacts, GTO, side-13
    rerun, q ladder, or fixture tuning.
  - Do not add a permanent test unless you shrink older coverage in the same
    pass.
  - Keep rejected-convention comparisons in tmp/work only.
  - No UI escalation; write `.agent_handoffs/ATTENTION.md` and stop if blocked.

Validation:
  - Update or copy `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl`
    so the accepted H1/J path consumes the new payload helper.
  - Keep the rejected signed-final-weight and raw-no-division comparisons
    probe-local if they are still useful.
  - Run the tmp/work H1/J probe, focused H1 gate, load check, and
    `git diff --check`.

Report:
  - functions/files added or edited;
  - H1/J probe result and whether it matches the previous accepted summary;
  - validation run;
  - deletion/shrinkage report:
      - what probe-local duplication became unnecessary;
      - what remains probe-local;
      - why no permanent test was added, or what older coverage was shrunk if
        one was added;
      - what remains before any RHF or fixture-rule work.

-- repo-manager@macmini
