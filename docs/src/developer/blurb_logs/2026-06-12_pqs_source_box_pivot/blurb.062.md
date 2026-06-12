Purpose:
  Run the side13 multi-layer PQS He RHF probe now that H1 and J are coherent.
  Keep it probe-only and non-acceptance.

Context:
  Pass 061 validated the pre-RHF diagnostics for the WL-aligned side13
  multi-layer PQS fixture:

  ```text
  parent count: 13
  mapping:      AsinhMapping(c = 0.1, s = 1.0, tail = 10)
  core:         (4:10)^3
  outer box:    (1:13)^3
  shell layers: 3
  final dim:    1549
  Z=2 H1:      -1.9755618232013417
  J:            1.2169264388860319
  ```

  H1/J are close to the WL side13 diagnostics, and the pre-final density
  interaction has positive weights and a finite symmetric pair matrix.

Exact task:
  Build an ignored `tmp/work` probe for restricted closed-shell He RHF on the
  side13 multi-layer PQS final basis.

  Use the same convention as the prior complete core/shell RHF probes:

  ```text
  final one-electron basis: orthonormal final basis
  density interaction:     pre-final positive-weight gauge
  c_prefinal = combined_lowdin_cleanup * c_final
  G_pre      = 2 * diag(V_pre * n_pre) - rho_pre .* V_pre
  F_final    = H_final + L' * G_pre * L
  ```

  Use `Z = 2` at the origin and compare to:

  ```text
  He HF reference:       -2.8616799956122388788
  WL side13 RHF probe:   -2.8364979997009137
  ```

Artifacts:
  Write ignored artifacts:

  ```text
  tmp/work/pqs_multilayer_shell_side13_rhf_probe.jl
  tmp/work/pqs_multilayer_shell_side13_rhf_probe_summary.txt
  ```

Report:
  Include:
  - final dimension and shell-layer counts;
  - final overlap identity error;
  - H1/J values;
  - RHF convergence, iterations, one-electron energy, electron-electron energy,
    total energy;
  - error vs He HF reference;
  - comparison to WL side13 RHF;
  - density trace/electron count in final and pre-final gauges;
  - timing split for plan/final basis, one-body/H1, density interaction, RHF;
  - explicit nonclaims: no GTO, no driver route, no export/artifact route, no
    signed-final-weight density, no raw no-division density, no fixed-block
    pair authority.

Decision rule:
  If RHF converges and is physically sane, report whether this should remain a
  reference probe, become a future compact gate candidate, or be followed by a
  fixture-rule review. Do not promote it in this pass.

Do not:
  - add permanent tests;
  - add GTO, driver wiring, exports, or artifacts;
  - use signed final weights or raw no-division density;
  - use fixed-block pair data as active authority;
  - start a new spacing/Z rule study.

Docs:
  Update the PQS near-term plan only if the RHF result is concise and useful.
  Keep it labeled as a probe/non-acceptance result.

Validation:
  - run the RHF probe;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited;
  - RHF result and timing;
  - comparison to WL side13 and He HF reference;
  - whether the final/pre-final RHF convention remains accepted;
  - whether this remains probe-only or suggests a future gate candidate;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.062.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
