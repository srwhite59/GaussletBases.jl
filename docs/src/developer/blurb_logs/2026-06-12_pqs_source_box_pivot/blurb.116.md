Pass 116 - RHF Fock/residual convention audit, no production code

Baseline:

- Current pushed HEAD should include `89eb11d6 Record PQS RHF damping sweep`.
- Compact route-owned private RHF SCF still does not converge; simple damping
  did not fix it.

Task:

Audit the private RHF one-step Fock/residual convention before any production
damping/mixing or route wiring.

Scope:

- No production code changes.
- No tracked tests/docs unless the audit result itself is a small tracked
  design note and you are confident it is necessary. Prefer response-only.
- Local ignored `tmp/work` scripts are allowed if needed for small algebraic
  checks, but do not rerun the heavy compact route unless absolutely necessary.

Questions to answer:

1. Density convention:
   - In `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`, is the
     supplied `final_density` a spin-summed closed-shell density `P = 2 C C'`?
   - Is dividing by occupancy to get `orbital_density_final` correct for the
     existing `_pqs_complete_core_shell_restricted_one_orbital_interaction_energy`
     convention?
   - Is the direct/exchange factor in `coulomb_pre =
     2*Diagonal(V*n) - rho .* V` consistent with a Fock operator whose energy is
     `0.5 * tr(P * G(P))`?

2. Energy/Fock consistency:
   - Does `two_body_energy = 0.5 * tr(P * coulomb_final)` match the existing
     one-orbital self-Coulomb diagnostic and generalize sensibly to multi-orbital
     closed-shell density?
   - Is the Fock matrix the derivative of the reported energy under this
     convention, or only a diagnostic effective operator?
   - If only diagnostic, state that explicitly and recommend label changes
     before more SCF work.

3. Residual definition:
   - Is `density_change = norm(P_next - P_current, Inf)` enough to judge SCF
     convergence in this nonorthogonal-history/pre-final interaction setup?
   - Should a future probe measure a commutator/residual such as `[F, P]` in the
     ordinary final basis?
   - Is the ordinary final-basis metric assumption still valid here?

4. Likely cause of nonconvergence:
   - Based on code and local trace/sweep evidence, is the plateau more likely a
     missing damping/acceleration problem, a factor-of-two/Fock convention issue,
     or an insufficient residual criterion?

Relevant files:

- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- ignored local summaries under `tmp/work/pqs_compact_route_owned_rhf_*`

Decision rules:

- Do not add damping/mixing.
- Do not change formulas unless the audit exposes a clear bug; if a clear bug
  is found, stop and report the exact corrective task rather than patching in
  this pass.
- Do not wire the route driver.
- Do not add report aliases/options/public API.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- `git status --short --branch`
- If local algebraic checks are run, report commands and results.

Report back:

- Audit conclusions for density convention, energy/Fock consistency, residual
  definition, and likely nonconvergence cause.
- Recommended next pass.
- Any exact formula/label bug if found.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
