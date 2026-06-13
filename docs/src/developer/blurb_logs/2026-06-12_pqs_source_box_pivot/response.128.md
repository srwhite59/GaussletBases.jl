Pass 128 response - default Fock-DIIS history 8 compact route-smoke confirmation

Scope:
- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No tolerance changes.
- No route wiring/report fields/public API/exports/artifacts/GTO/IDA/MWG.
- No fixture promotion.
- No `hfdmrg` or CR2 comparison.

Artifact paths:
- Probe script:
  `tmp/work/pqs_compact_route_owned_rhf_default_fock_diis_probe.jl`
- Summary:
  `tmp/work/pqs_compact_route_owned_rhf_default_fock_diis_probe_summary.txt`
- Table:
  `tmp/work/pqs_compact_route_owned_rhf_default_fock_diis_probe.tsv`

Command/results:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_default_fock_diis_probe.jl`
  - elapsed: `107.467666458` seconds
- `git status --short --branch`
  - `## main...origin/main`

Probe table:
```text
resolved_max_history  converged  iter  final_or_last_total_energy  density_change       update_density_change  energy_change          commutator_residual  spatial_commutator  trace_error            idempotency_error      diis_used  fallbacks  solve_failures  coeff_pathologies  fallback_iterations  final_gate_blocked
8                     true       34    -10.032119189804888         8.3320104149464669e-9 3.0626759672802706e-7 2.8421709430404007e-14 2.9455127052713248e-9 1.4727563526356624e-9 1.7763568394002505e-15 8.3266726846886741e-17 30         3          0               3                  8,9,10              false
```

Summary:
- The omitted-history default path resolved `max_history = 8`.
- The compact route-smoke fixture converged under the strict gates.
- The run reproduced the explicit history-8 behavior from passes 125 and 126:
  - converged at iteration 34;
  - density delta `8.3320104149464669e-09`;
  - commutator residual `2.9455127052713248e-09`;
  - spatial commutator residual `1.4727563526356624e-09`;
  - energy change `2.8421709430404007e-14`;
  - trace/idempotency at machine precision;
  - fallback iterations `8, 9, 10`;
  - no solve failures;
  - final recomputed diagnostics did not block convergence.

Recommended next pass:
- No-edit route-RHF adoption audit.
- Identify the exact private driver slot and compact request object needed to
  call the private RHF diagnostic helper, while keeping:
  - route/report/public behavior unchanged;
  - serious HF outside this PQS route-smoke seam;
  - no fixture promotion.
- Do not implement route wiring until that audit pins down the smallest private
  boundary and the fields that can stay diagnostic-only.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none.
- simplified: confirmed that callers can omit `max_history` and get the
  validated private default.
- quarantined: default-history probe script/table/summary remain ignored
  `tmp/work` artifacts.
- not deleted because: local probe artifacts still document the private-control
  evidence chain for manager review.
- exact remaining caller/blocker:
  tracked private Fock-DIIS now converges the compact route-smoke fixture through
  the default control path; the remaining blocker is a no-edit audit for the
  smallest private route-driver request/slot, not SCF control behavior.

-- repo-doer@macmini
