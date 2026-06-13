Pass 118 - add residual diagnostics to private RHF SCF payload

Baseline:

- Current pushed HEAD should include `139c0d33 Record PQS RHF residual probe`.
- Pass 117 showed the compact private RHF SCF state has a real ordinary
  final-basis commutator residual around `1.3e-5`.

Task:

Add residual diagnostics to the private RHF SCF payload. Do not add damping,
DIIS, route wiring, or public/report behavior.

Implementation surface:

- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- focused SCF test:
  `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`

Goal:

Make the private SCF payload report stationarity diagnostics directly:

- density trace error;
- closed-shell idempotency error for `D = P / occupancy`;
- ordinary final-basis commutator residual, e.g.
  `norm(F * P - P * F, Inf)` or equivalent with `D`;
- explicit labels:
  - `density_change_rule = :fixed_point_spin_summed_density_inf_norm`
  - `residual_metric = :ordinary_final_basis_commutator_inf_norm`
  - `idempotency_rule = :closed_shell_spatial_density_idempotency`
  - `orbital_metric = :ordinary_orthonormal_final_basis`

Recommended implementation:

- Add a small private helper to compute compact residual diagnostics from:
  - final/spin-summed density `P`;
  - effective Fock matrix `F`;
  - electron count;
  - occupancy.
- The helper should return a compact NamedTuple, not a matrix payload.
- Use it in the converged return after final one-step recomputation.
- For the `:scf_not_converged` path, if a final density and one-step/Fock are
  available, include the same compact residual diagnostics in the blocked
  summary. If doing so would complicate the path too much, at least add labels
  making residuals unavailable explicit.
- Keep summaries/metadata free of large matrices.

Tests:

- Update the focused SCF test.
- In the existing tiny self-consistent fixture, assert:
  - commutator residual is near zero;
  - trace error is near zero;
  - idempotency error is near zero;
  - labels match the expected symbols;
  - route/report/export/artifact/public nonclaims remain false.
- Add a tiny nonstationary synthetic case only if it is easy and does not
  broaden the test much.

Exclusions:

- Do not add damping, mixing, DIIS, or acceleration.
- Do not rerun the compact route probe in this pass.
- Do not route-wire RHF.
- Do not add report aliases/options/public API.
- Do not touch GTO, IDA/MWG, exports, artifacts, or production route behavior.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- Helper/field names added.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
