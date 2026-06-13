Pass 118 manager review

Accepted.

The implementation adds compact residual diagnostics to the private SCF payload
without changing the RHF formula, adding damping/DIIS, or wiring any route/report
surface. The diagnostic labels are now explicit:

- `density_change_rule = :fixed_point_spin_summed_density_inf_norm`
- `residual_metric = :ordinary_final_basis_commutator_inf_norm`
- `idempotency_rule = :closed_shell_spatial_density_idempotency`
- `orbital_metric = :ordinary_orthonormal_final_basis`

Converged returns compute residual diagnostics after the final one-step
recomputation. Terminal nonconverged paths attempt the same recomputation and
otherwise carry blocked residual diagnostics explicitly. Summaries do not carry
large matrices.

Validation repeated by manager:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  passed: 53/53.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed.
- `git diff --check` passed.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: residual state is now a compact nested summary instead of a flat
  field cloud or local-only convention.
- quarantined: SCF remains private diagnostic/prototype behavior; no route,
  report, public API, export, or artifact behavior was added.
- not deleted because: this is the active private diagnostic seam needed before
  any SCF-control experiment.
- exact remaining caller/blocker: compact real private SCF still needs a
  residual-aware local probe from the updated payload before choosing any
  acceleration/mixing implementation.

-- repo-manager@macmini
