Pass 107 - cross-check RHF one-step convention against one-orbital diagnostic

Baseline:

- Current pushed HEAD should include `abfd0a6d Add PQS RHF one-step payload`.
- The new private helper is
  `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`.
- Existing convention authority is
  `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)`
  and the helper it uses in
  `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`.

Task:

Add one tiny focused convention cross-check, preferably in:

`test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`

Goal:

For a closed-shell one-orbital final density, verify the new one-step
density-matrix path agrees with the existing one-orbital self-Coulomb
diagnostic on the same density interaction and final orbital.

Suggested shape:

- Reuse the synthetic payload helper already in the one-step test if possible.
- Choose a normalized occupied final orbital in the tiny synthetic fixture,
  for example `c = [1.0, 0.0]`.
- Build spin-summed closed-shell density:
  `final_density = 2.0 .* (c * c')`.
- Call:
  `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`
- Call:
  `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_orbital_self_coulomb(
      density_interaction,
      c,
  )`
  or the module-qualified equivalent available through `GaussletBases`.
- Assert:
  - one-step status is materialized;
  - self-Coulomb status is materialized;
  - `one_step.two_body_energy ≈ self_coulomb.self_coulomb`;
  - one-step density convention remains
    `:spin_summed_closed_shell_final_density`;
  - contraction rule remains
    `:pre_final_restricted_direct_minus_exchange_from_orbital_density`;
  - `scf_materialized === false`;
  - `rhf_converged === false`.

Decision rule:

- This should be test-only.
- Change production code only if the test exposes an obvious typo in the
  private one-step helper.
- If the comparison exposes a genuine convention mismatch, stop and report
  without forcing the test.

Exclusions:

- Do not add SCF.
- Do not add a new Fock convention.
- Do not route-wire this helper.
- Do not add report aliases, driver options, exports, artifacts, GTO, IDA/MWG,
  or fixture promotion.
- Do not run the heavy source-box dry-run.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- Whether this was test-only or required a typo fix.
- The observed self-Coulomb/two-body-energy comparison.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
