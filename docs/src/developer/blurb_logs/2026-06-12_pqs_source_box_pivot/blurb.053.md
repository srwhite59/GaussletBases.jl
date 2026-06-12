Purpose:
  Run the first complete core/shell PQS He RHF probe, using the reviewed
  pre-final density-interaction seam. This is exploratory physics validation,
  not an acceptance test yet.

Context:
  Pass 052 accepted the density consumption rule:

  ```text
  density interaction gauge: pre-final localized positive-weight gauge
  final orbital consumption: c_prefinal = combined_lowdin_cleanup * c_final
  H1 self-Coulomb J:         0.6397851751855723
  fixed-block oracle J:      0.6397857768997106
  ```

  The remaining question is whether a closed-shell RHF iteration can consume
  the final one-electron Hamiltonian and the pre-final density interaction
  coherently.

Exact task:
  Build an ignored `tmp/work` probe for complete core/shell PQS He RHF on the
  current 223-function final basis.

  Use the existing complete core/shell final-basis and one-electron machinery.
  Use a one-center He atom with nuclear charge `Z = 2` at the origin. Use the
  pre-final density interaction from pass 052 for electron-electron terms.

Required RHF convention:
  The occupied spatial orbital is solved in the orthonormal final basis, but
  the electron-electron term is evaluated in the pre-final density gauge.

  For a final-basis orbital `c_final`:

  ```text
  c_prefinal = combined_lowdin_cleanup * c_final
  rho_final  = c_final * c_final'
  rho_pre    = c_prefinal * c_prefinal'
  n_pre      = diag(rho_pre)
  ```

  Use the two-index closed-shell density-density energy:

  ```text
  E_one = 2 * tr(rho_final * H_final)
  E_ee  = 2 * n_pre' * V_pre * n_pre - sum(rho_pre .* V_pre .* rho_pre)
  E     = E_one + E_ee
  ```

  The Fock-like matrix diagonalized in the final basis should include the
  chain-rule transformation of the pre-final electron-electron response:

  ```text
  G_pre   = 2 * diag(V_pre * n_pre) - rho_pre .* V_pre
  F_final = H_final + L' * G_pre * L
  L       = combined_lowdin_cleanup
  ```

  If this formula appears inconsistent with the existing He/WL RHF convention
  or with the matrix dimensions, stop and write `ATTENTION.md`; do not invent a
  shortcut.

Probe artifacts:
  Write ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_he_rhf_probe.jl
  tmp/work/pqs_complete_core_shell_he_rhf_probe_summary.txt
  ```

Report:
  Include:
  - final dimension and support counts;
  - final overlap identity error;
  - one-electron lowest orbital energy for the Z=2 Hamiltonian;
  - H1 self-Coulomb diagnostic for the Z=2 occupied orbital if available;
  - RHF convergence status and iteration count;
  - RHF one-electron energy, electron-electron energy, total energy;
  - comparison with the He HF reference `-2.8616799956122388788`;
  - whether the result is physically sane enough to consider promotion later;
  - timing split for final-basis build, density interaction build, and RHF
    solve;
  - explicit nonclaims: no GTO, no driver route, no export/artifact route, no
    signed-final-weight density, no raw no-division density.

Do not:
  - add a permanent test in this pass unless it replaces or shrinks older PQS
    oracle coverage;
  - add broad report fields or metadata checks;
  - use fixed-block pair data as active authority;
  - use signed final weights as density weights;
  - use raw projected no-division density;
  - add GTO, PQS driver wiring, exports, or artifacts.

Validation:
  - run the RHF probe;
  - if source changed, run `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - whether source changed or this stayed probe-only;
  - final RHF result and HF-reference error;
  - whether the final/pre-final Fock convention was accepted or blocked;
  - timing summary;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.053.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
