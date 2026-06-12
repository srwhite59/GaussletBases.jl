Review 088:

Accepted. The pass added the intended route-owned support-density input
helpers without promoting density interaction, RHF, acceptance, driver wiring,
or fixture tuning.

What I checked:

- `src/pqs_multilayer_support_density.jl` adds
  `pqs_multilayer_support_weights(...)` and
  `pqs_multilayer_support_pair_raw_numerator_matrix(...)`.
- Both helpers consume the complete core/shell support ordering:
  `core_support_states` followed by `shell_support_states`.
- The support-weight helper builds product support weights.
- The raw-numerator helper contracts raw pair-factor terms with positive
  Coulomb coefficients, does not divide by weights, and does not materialize
  final density interaction.
- `src/GaussletBases.jl` includes the new file after support one-body and
  before source-plan/H1 payload code.
- No permanent test, RHF path, driver wiring, export, artifact, fixture tuning,
  or broad report surface was added.

Validation rerun by manager:

- `julia --project=. tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl`
  passed and reported `support_density_input_source =
  pqs_multilayer_support_density_helpers`.
- The compact H1/J probe still reports H1 energy about
  `-0.48047934800387126` and self-Coulomb
  `0.6397851751855723`.
- The fixed-oracle self-Coulomb comparison remains within about
  `6.02e-7`.
- The rejected signed-final-weight and raw-no-division conventions still report
  visibly wrong values, so the convention guard remains useful.
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  passed, 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion / shrinkage:

- The tmp/work density-consumption probe no longer needs its probe-local
  support-weight or raw-support-pair-numerator builders for the accepted path.
- The normalized-pair and signed-final-weight comparisons remain probe-local
  because they are rejected-convention diagnostics.
- H1 orbital coefficient extraction remains probe-local and is the next
  obvious input-production gap.
- No permanent test was added because this pass only replaces probe-local input
  production; the existing H1 gate plus the tmp/work H1/J probe are the right
  validation surfaces for now.

-- repo-manager@macmini
