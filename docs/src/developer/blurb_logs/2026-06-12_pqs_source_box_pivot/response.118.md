Pass 118 response - add residual diagnostics to private RHF SCF payload

Files changed:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`

Helper/field names added:
- Added private helper:
  `_pqs_multilayer_complete_core_shell_rhf_scf_residual_diagnostics(...)`
- The helper returns a compact `NamedTuple` with:
  - `status`
  - `blocker`
  - `density_trace`
  - `density_trace_error`
  - `closed_shell_idempotency_error`
  - `commutator_residual`
  - `spatial_commutator_residual`
  - `density_change_rule = :fixed_point_spin_summed_density_inf_norm`
  - `residual_metric = :ordinary_final_basis_commutator_inf_norm`
  - `idempotency_rule = :closed_shell_spatial_density_idempotency`
  - `orbital_metric = :ordinary_orthonormal_final_basis`
- SCF summaries now carry:
  - top-level label fields `density_change_rule`, `residual_metric`,
    `idempotency_rule`, `orbital_metric`;
  - compact nested `residual_diagnostics`.

Behavior/shape notes:
- Converged SCF return computes residual diagnostics after the final one-step
  recomputation.
- The `:scf_not_converged` terminal path now attempts one final one-step
  recomputation on the returned final density. If that succeeds, the blocked
  summary reports residual diagnostics against the final density and marks
  `final_one_step_recomputed = true`.
- Earlier blocked paths keep residual diagnostics explicit but unavailable via
  blocked `residual_diagnostics` status, for example `:missing_final_density`.
- No large matrices were added to summaries or metadata.
- No damping, DIIS, acceleration, route wiring, report aliases, public API,
  exports, artifacts, GTO, or IDA/MWG behavior was added.

Focused test update:
- The existing tiny self-consistent SCF fixture now asserts:
  - commutator residual near zero;
  - density trace error near zero;
  - closed-shell idempotency error near zero;
  - residual labels match expected symbols;
  - route/report/export/artifact/public nonclaims remain false.
- The missing-contract blocked case now asserts labels remain present and
  residual diagnostics are explicitly blocked as unavailable.

Validation:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  - Passed: `53 / 53`
  - Note: this run precompiled `GaussletBases` first.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`
- `git diff --check`
  - Passed.

Git status:
```text
## main...origin/main
 M src/pqs_multilayer_complete_core_shell_rhf.jl
 M test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

Deletion/shrinkage report:
- deleted: none.
- simplified: residual state is now a compact nested summary instead of a new
  flat field cloud.
- quarantined: no new tracked quarantine; existing ignored `tmp/work` probes
  remain local-only.
- not deleted because: the private RHF SCF path still needs residual-aware SCF
  control before any route-facing shrinkage or removal.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` now reports
  residual diagnostics, but the compact real fixture remains blocked by
  `:scf_not_converged`; next blocker is private SCF control/acceleration design,
  not formula convention or route wiring.

-- repo-doer@macmini
