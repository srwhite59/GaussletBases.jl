Purpose:
  Run one more complete core/shell PQS He RHF scaling probe at `q=11`, putting
  the final basis near 2000 functions.

Context:
  The q=5/q=7/q=9 ladder is improving cleanly:

  ```text
  q=5  final dim 223   H1 -1.8476619225061284  RHF -2.7213372828531668
  q=7  final dim 561   H1 -1.9599740171825744  RHF -2.810068050134403
  q=9  final dim 1115  H1 -1.9869007472289404  RHF -2.8499091618019303
  ```

  q=9 is within `0.011770833810308634 Ha` of the He HF reference and ran in
  about `7.82 s` for the measured probe phases.

Exact task:
  Run a probe-only larger complete core/shell PQS He RHF fixture:

  ```text
  current_box:    (1:13)^3
  inner_box:      (2:12)^3
  raw dims:       11 x 11 x 11
  expected support count:          2197
  expected core support count:     1331
  expected shell support count:    866
  expected shell retained count:   602
  expected final dimension:        1933
  ```

  Reuse the accepted final/pre-final RHF convention:

  ```text
  c_prefinal = combined_lowdin_cleanup * c_final
  G_pre      = 2 * diag(V_pre * n_pre) - rho_pre .* V_pre
  F_final    = H_final + L' * G_pre * L
  ```

Probe artifacts:
  Write ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_he_rhf_q11_probe.jl
  tmp/work/pqs_complete_core_shell_he_rhf_q11_probe_summary.txt
  ```

Report:
  Include:
  - actual dimensions and support counts;
  - final overlap identity error;
  - Z=1 H1 lowest energy if cheap;
  - Z=2 H1 lowest energy and error vs `-2.0`;
  - H1 self-Coulomb J for the Z=2 occupied orbital if available;
  - RHF convergence, iterations, one-electron energy, electron-electron energy,
    total energy, and error vs He HF reference `-2.8616799956122388788`;
  - timing split for final-basis/support build, one-electron build,
    density-interaction build, and RHF solve;
  - comparison against q=5, q=7, and q=9;
  - whether q=11 looks like a good permanent gate candidate, whether q=9 is a
    better gate because it is cheaper, or whether another fixture/design review
    is needed.

Do not:
  - add a permanent test in this pass;
  - add broad metadata/report assertions;
  - use signed final weights or raw no-division density;
  - use fixed-block pair data as active authority;
  - add GTO, driver wiring, exports, or artifacts.

Source changes:
  Prefer no source changes. If q=11 exposes hard-coded size assumptions, make
  only the smallest parameterization needed and report it. If it exposes a
  design blocker, stop and write `ATTENTION.md`.

Validation:
  - run the q=11 RHF probe;
  - if source changed, run `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - whether source changed or this stayed probe-only;
  - q=11 H1/J/RHF numbers and timing;
  - comparison to q=5/q=7/q=9;
  - whether the final/pre-final RHF convention still appears accepted;
  - recommendation for q=9 vs q=11 as a future compact gate candidate;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.056.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
