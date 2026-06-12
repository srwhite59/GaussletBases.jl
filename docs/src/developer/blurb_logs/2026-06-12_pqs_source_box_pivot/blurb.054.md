Purpose:
  Probe whether the complete core/shell PQS route scales from a `5^3` direct
  core plus one surrounding shell to a `7^3` direct core plus one surrounding
  shell, still without promoting the RHF result to acceptance.

Context:
  Pass 053 proved route coherence for the compact fixture:

  ```text
  current_box:    (1:7)^3
  inner_box:      (2:6)^3
  raw dims:       5 x 5 x 5
  final dim:      223
  Z=2 H1:         -1.8476619225061284
  RHF total:      -2.7213372828531668
  He HF error:    +0.14034271275907217
  ```

  The result is physically sane but too compact to be an acceptance baseline.

Exact task:
  Run a probe-only larger complete core/shell PQS He RHF fixture:

  ```text
  current_box:    (1:9)^3
  inner_box:      (2:8)^3
  raw dims:       7 x 7 x 7
  expected support count:          729
  expected core support count:     343
  expected shell support count:    386
  expected shell retained count:   218
  expected final dimension:        561
  ```

  Reuse the same final/pre-final RHF convention accepted in pass 053:

  ```text
  c_prefinal = combined_lowdin_cleanup * c_final
  G_pre      = 2 * diag(V_pre * n_pre) - rho_pre .* V_pre
  F_final    = H_final + L' * G_pre * L
  ```

Probe artifacts:
  Write ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_he_rhf_q7_probe.jl
  tmp/work/pqs_complete_core_shell_he_rhf_q7_probe_summary.txt
  ```

Report:
  Include:
  - actual dimensions and support counts;
  - final overlap identity error;
  - Z=1 H1 lowest energy if cheap, for comparison with the earlier H1 route;
  - Z=2 H1 lowest energy and error vs `-2.0`;
  - H1 self-Coulomb J for the Z=2 occupied orbital if available;
  - RHF convergence, iterations, one-electron energy, electron-electron energy,
    total energy, and error vs He HF reference `-2.8616799956122388788`;
  - timing split for final-basis/support build, one-electron build,
    density-interaction build, and RHF solve;
  - comparison against the pass 053 compact fixture values;
  - whether this looks ready for a permanent smoke gate, a larger probe, or a
    convention audit.

Do not:
  - add a permanent test in this pass;
  - add broad metadata/report assertions;
  - use signed final weights or raw no-division density;
  - use fixed-block pair data as active authority;
  - add GTO, driver wiring, exports, or artifacts.

Source changes:
  Prefer no source changes. If hard-coded `5^3` assumptions block the probe,
  make only the smallest parameterization needed for the complete core/shell
  route and explain exactly what was hard-coded. If the needed change becomes a
  design refactor, stop and write `ATTENTION.md`.

Validation:
  - run the q=7 RHF probe;
  - if source changed, run `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - whether source changed or this stayed probe-only;
  - q=7 H1/J/RHF numbers and timing;
  - comparison to pass 053 q=5 compact fixture;
  - whether the final/pre-final RHF convention still appears accepted;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.054.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
