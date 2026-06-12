Review 089:

Accepted. The pass adds a compact route-owned H1/J payload that replaces the
accepted probe-local wiring for lowest-orbital extraction, density interaction,
and self-Coulomb contraction. It does not add RHF, a permanent acceptance gate,
driver wiring, fixture tuning, exports, or artifacts.

What I checked:

- `pqs_multilayer_complete_core_shell_h1_j_payload(...)` validates the source
  plan, complete core/shell final basis, and existing H1 payload.
- It uses `pqs_multilayer_support_weights(...)` and
  `pqs_multilayer_support_pair_raw_numerator_matrix(...)`.
- It calls the reviewed
  `pqs_complete_core_shell_pre_final_density_interaction(...)` and
  `pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)` helpers.
- It extracts the lowest H1 orbital from the final Hamiltonian carried by the
  H1 payload, and the reconstructed H1 energy matches the H1 payload energy to
  roundoff in the probe.
- It does not implement signed-final-weight division, raw-no-division final
  density, density-normalized pair-factor authority, RHF, GTO, driver wiring,
  exports, artifacts, q ladders, side-13 reruns, or fixture tuning.

Validation rerun by manager:

- `julia --project=. tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl`
  passed.
- The probe reports H1 energy about `-0.48047934800387126`,
  self-Coulomb `0.6397851751855723`, and fixed-oracle self-Coulomb delta about
  `-6.02e-7`.
- The rejected signed-final-weight and raw-no-division comparisons remain
  probe-local and visibly wrong.
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  passed, 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion / shrinkage:

- The accepted tmp/work H1/J path no longer manually wires support weights,
  support raw numerator, density interaction, lowest-orbital extraction, and
  self-Coulomb as separate route steps.
- Remaining probe-local pieces are the rejected-convention diagnostics and
  their local comparisons.
- No permanent test was added because this is still a diagnostic seam and the
  existing focused H1 gate plus tmp/work H1/J probe are the appropriate
  validation surfaces.

Caveat:

- The helper summary uses `density_density_materialized = true` for the
  pre-final H1/J density interaction. Treat that as the H1/J diagnostic
  interaction, not as a full RHF density route or fixture acceptance claim.

-- repo-manager@macmini
