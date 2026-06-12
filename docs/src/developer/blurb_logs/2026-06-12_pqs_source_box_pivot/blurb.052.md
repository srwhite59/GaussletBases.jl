Purpose:
  Add the reviewed pre-final density-interaction object and final-orbital
  consumption rule for the complete PQS route, still before RHF.

Context:
  Pass 051 established the accepted density gauge:

  ```text
  pre-final localized weights: positive
  pre-final positive-weight J: 0.6397851751855723
  fixed-block oracle J:        0.6397857768997106
  final signed-weight J:       rejected
  raw no-division J:           rejected
  ```

  The remaining blocker before RHF is:

  ```text
  :missing_reviewed_pre_final_density_interaction_consumption_for_rhf
  ```

Exact task:
  Add a compact pre-final density-interaction seam for the complete core/shell
  route and validate how final-basis orbital coefficients consume it.

Implementation target:
  A narrow helper/object may be added, for example:

  ```text
  pqs_complete_core_shell_pre_final_density_interaction(...)
  pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)
  ```

  Naming can differ if the codebase suggests a better local name, but the
  concept must be clear:

  ```text
  density interaction lives in pre-final/localized positive-weight gauge
  final orbital coefficients are mapped to that gauge explicitly
  ```

Expected coefficient rule to audit:
  If:

  ```text
  final_coefficients = pre_final_coefficients * combined_lowdin_cleanup
  ```

  then a final-basis orbital coefficient vector `c_final` should be represented
  in pre-final density gauge by:

  ```text
  c_prefinal = combined_lowdin_cleanup * c_final
  ```

  Verify this reconstruction in the probe before using it.

Required checks:
  - pre-final weights finite and positive;
  - pre-final pair matrix finite/symmetric;
  - final-to-pre-final coefficient reconstruction error is small;
  - H1 self-Coulomb J via the explicit final-to-pre-final rule matches the
    pass 051 pre-final J and the same-geometry fixed-block oracle;
  - signed final-weight division remains rejected/nonclaim;
  - RHF not materialized.

Do not:
  - run RHF;
  - add GTO, driver wiring, exports, or artifacts;
  - use fixed-block pair_sum as active authority;
  - use final signed weights as density weights;
  - add broad tests. Use a focused probe first unless you are replacing old
    coverage.

Probe:
  Create ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl
  tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe_summary.txt
  ```

  Include:
  - final dimension;
  - pre-final weight diagnostics;
  - pre-final pair matrix diagnostics;
  - final-to-pre-final reconstruction error;
  - H1 self-Coulomb J and oracle delta;
  - nonclaim flags for RHF/GTO/driver/export/artifact.

Validation:
  - focused consumption probe;
  - if source changed, run `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact helper names;
  - J diagnostic and fixed-oracle delta;
  - whether the final-to-pre-final consumption rule is accepted or blocked;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.052.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
