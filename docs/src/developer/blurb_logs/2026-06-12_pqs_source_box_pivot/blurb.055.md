Purpose:
  Push the complete core/shell PQS He RHF scaling ladder one more step, from
  the successful `7^3` direct core to a `9^3` direct core plus one surrounding
  shell.

Context:
  Pass 054 showed the q=7 route is coherent and improves the compact fixture:

  ```text
  q=7 current_box:    (1:9)^3
  q=7 inner_box:      (2:8)^3
  q=7 final dim:      561
  q=7 Z=2 H1:         -1.9599740171825744
  q=7 J:              1.24311548900175
  q=7 RHF total:      -2.810068050134403
  q=7 He HF error:    +0.051611945477835874
  q=7 total probe:    4.102156 s
  ```

  The next question is whether the same route remains practical at a larger
  final-basis dimension and whether H1/RHF continue improving.

Exact task:
  Run a probe-only larger complete core/shell PQS He RHF fixture:

  ```text
  current_box:    (1:11)^3
  inner_box:      (2:10)^3
  raw dims:       9 x 9 x 9
  expected support count:          1331
  expected core support count:     729
  expected shell support count:    602
  expected shell retained count:   386
  expected final dimension:        1115
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
  tmp/work/pqs_complete_core_shell_he_rhf_q9_probe.jl
  tmp/work/pqs_complete_core_shell_he_rhf_q9_probe_summary.txt
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
  - comparison against q=5/pass 053 and q=7/pass 054;
  - whether the next move should be a permanent smoke gate, a larger probe, or
    fixture/design review.

Do not:
  - add a permanent test in this pass;
  - add broad metadata/report assertions;
  - use signed final weights or raw no-division density;
  - use fixed-block pair data as active authority;
  - add GTO, driver wiring, exports, or artifacts.

Source changes:
  Prefer no source changes. If hard-coded size assumptions block q=9, make only
  the smallest parameterization needed and report it. If the needed change is a
  broader design refactor, stop and write `ATTENTION.md`.

Validation:
  - run the q=9 RHF probe;
  - if source changed, run `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - whether source changed or this stayed probe-only;
  - q=9 H1/J/RHF numbers and timing;
  - comparison to q=5 and q=7;
  - whether the final/pre-final RHF convention still appears accepted;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.055.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
