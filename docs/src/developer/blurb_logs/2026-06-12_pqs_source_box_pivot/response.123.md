Pass 123 response - tighten RHF SCF final convergence gating

Files changed:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`

Exact final convergence gating behavior:
- `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` still performs
  iteration-input convergence checks during the loop.
- In the converged branch it then recomputes the final one-step payload on the
  returned `final_density`.
- It now requires the recomputed final diagnostics to pass before returning
  `:materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload` with
  `rhf_converged = true`:
  - `final_one_step_density_matches_final_density`;
  - final commutator residual `<= residual_atol`;
  - final trace error `<= trace_atol`;
  - final idempotency error `<= idempotency_atol`.
- If iteration-input convergence passed but final recomputed diagnostics fail,
  the helper now returns the blocked SCF payload with blocker
  `:scf_not_converged`, carrying:
  - returned/final density;
  - recomputed final-one-step payload;
  - residual diagnostics through the existing blocked-summary path;
  - iteration history and compact DIIS counters.

Additional status/blocker labels:
- None added.
- The correction reuses existing labels:
  - materialized SCF status only after final diagnostics pass;
  - blocked SCF status with `:scf_not_converged` if final diagnostics fail.

Test update:
- The focused synthetic converged case now explicitly asserts the recomputed
  final residual diagnostics are within `residual_atol`, `trace_atol`, and
  `idempotency_atol`.
- Existing fixed-point and Fock-DIIS synthetic coverage still passes.
- No tracked compact route-smoke test was added.

Validation:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  - Passed: `91 / 91`
  - Note: run precompiled `GaussletBases`.
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
- simplified: convergence authority is now the final recomputed residual
  diagnostics, not only the iteration-input checks.
- quarantined: no new tracked quarantine; ignored compact probes remain in
  `tmp/work`.
- not deleted because: pass-122 private control/DIIS implementation is still
  under manager review and needs compact route validation before cleanup.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` now gates
  `rhf_converged = true` on final recomputed diagnostics; the remaining blocker
  is validating the tracked Fock-DIIS path on the compact route fixture, with
  route-driver integration still intentionally blocked.

-- repo-doer@macmini
