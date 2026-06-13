Pass 116 response - RHF Fock/residual convention audit

Audit scope:
- Read-only production audit of `src/pqs_multilayer_complete_core_shell_rhf.jl`
  and `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`.
- Used existing ignored local summaries under `tmp/work/pqs_compact_route_owned_rhf_*`.
- Added one ignored local algebraic check:
  `tmp/work/pqs_rhf_fock_convention_check.jl`.
- No production source, tracked tests, or tracked docs changed.

Commands/results:
- `julia --project=. tmp/work/pqs_rhf_fock_convention_check.jl`
  - `energy=1.7141550055673576`
  - `finite_difference=0.007269751023386561`
  - `fock_directional=0.007269750922258428`
  - `abs_error=1.0112813364493434e-10`
  - `relative_error=1.3910811157027015e-8`
- `git status --short --branch`
  - `## main...origin/main`

Density convention conclusion:
- `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)` expects
  `final_density` to be the spin-summed closed-shell density
  `P = 2 C_occ C_occ'`.
- The input contract currently allows only even-electron closed-shell RHF with
  `occupancy = 2`, `nocc = electron_count ÷ 2`, and no fractional/open-shell
  support.
- Dividing `final_density` by occupancy to form `orbital_density_final` is
  correct for the existing one-orbital convention, because
  `pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)` consumes the
  spatial/orbital density `D = C C'`, not the spin-summed `P`.
- The current
  `coulomb_pre = 2*Diagonal(V*n) - rho .* V`
  is consistent with the existing restricted direct-minus-exchange diagnostic
  when `rho` is the closed-shell spatial density matrix `D_pre = P_pre / 2`.

Energy/Fock consistency:
- `two_body_energy = 0.5 * tr(P * coulomb_final)` matches the existing
  one-orbital self-Coulomb diagnostic and is the natural multi-orbital
  closed-shell generalization for this private density-density interaction:
  `E2(D) = 2 n' V n - sum_ij D_ij V_ij D_ij`.
- Under the fixed `final_to_pre_final` map, the implemented
  `effective_fock_matrix = h1 + coulomb_final` is the derivative of the
  reported `tr(P*h1) + 0.5*tr(P*coulomb_final)` energy with respect to the
  spin-summed final density `P`.
- I did not find a factor-of-two formula bug in the one-step payload.
- Label caveat: this is still a private density-density diagnostic Fock for the
  pre-final interaction gauge, not production route Coulomb/Fock. The existing
  `private_diagnostic_only = true` label is important. Before more SCF work, I
  recommend adding explicit residual/metric labels rather than renaming the
  formula itself.

Residual definition:
- `density_change = norm(P_next - P_current, Inf)` is a useful fixed-point
  iteration diagnostic, but it is not enough by itself to judge RHF stationarity
  in this pre-final-interaction setup.
- A future focused probe should measure an ordinary final-basis commutator
  residual, for example a norm of `[F, P]` or `[F, D]`, alongside idempotency and
  trace checks. This should be reported as a diagnostic residual, not promoted
  to route behavior.
- The ordinary final-basis metric assumption is still the right assumption for
  this path: the SCF diagonalization consumes the final H1/final Fock in the
  orthonormal final retained basis. The pre-final history is in the interaction
  construction map, not a generalized-overlap SCF metric.

Likely nonconvergence cause:
- The trace/sweep evidence does not look like a factor-of-two/Fock convention
  issue:
  - energy is monotone decreasing;
  - no two-cycle signal was recorded;
  - density trace stays at 4;
  - simple damping changes residual scale but does not converge;
  - the algebraic finite-difference check confirms Fock/energy consistency.
- The immediate blocker is more likely insufficient residual characterization
  than a formula convention bug. The current pass can say the fixed-point
  density delta has plateaued, but it cannot say whether the final density is
  near stationary because no commutator residual is measured.
- If the commutator residual is also stuck, then the next likely cause becomes
  missing acceleration/mixing beyond simple scalar damping. Do not add that
  until the residual probe says the SCF problem is genuinely nonstationary.

Recommended next pass:
- Add a small private/local diagnostic residual probe for the compact RHF SCF
  path that records:
  - `norm(F*P - P*F, Inf)` or equivalent in the ordinary final basis;
  - trace error;
  - closed-shell idempotency error for `D = P/2`;
  - the existing fixed-point density delta and energy delta.
- Keep it local/ignored or private diagnostic-only. Do not wire the route
  driver, add damping/DIIS, or change production formulas in that pass.

Exact formula/label bug:
- Formula bug found: none.
- Label issue found: no incorrect current label, but future SCF work should add
  explicit `residual_metric` / `density_change_rule` labels so
  `density_change` is not mistaken for a full RHF stationarity criterion.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none.
- simplified: none in this no-edit audit.
- quarantined: the algebraic check is local ignored `tmp/work` only.
- not deleted because: the private RHF one-step/SCF helpers remain the active
  diagnostic seam for the compact complete core/shell path.
- exact remaining caller/blocker: `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`
  currently reports only fixed-point `density_change`/`energy_change`; the
  compact real probe remains blocked by `:scf_not_converged` without a measured
  ordinary-final-basis commutator residual.

-- repo-doer@macmini
