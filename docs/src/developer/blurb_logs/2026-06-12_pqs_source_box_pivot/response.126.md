Pass 126 response - confirm tracked Fock-DIIS history 8

Scope:
- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No private default changes.
- No tolerance changes.
- No route wiring/report fields/public API/exports/artifacts/GTO/IDA/MWG.
- No fixture promotion or serious-HF comparison.

Artifact paths:
- Probe script:
  `tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_history8_confirm_probe.jl`
- Summary:
  `tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_history8_confirm_probe_summary.txt`
- Table:
  `tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_history8_confirm_probe.tsv`

Command/results:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_history8_confirm_probe.jl`
  - elapsed: `105.134415291` seconds
- `git status --short --branch`
  - `## main...origin/main`

Probe table:
```text
max_history  converged  iter  final_or_last_total_energy  density_change       update_density_change  energy_change          commutator_residual  spatial_commutator  trace_error            idempotency_error      diis_used  fallbacks  solve_failures  coeff_pathologies  fallback_iterations  final_gate_blocked
8            true       34    -10.032119189804888         8.3320104149464669e-9 3.0626759672802706e-7 2.8421709430404007e-14 2.9455127052713248e-9 1.4727563526356624e-9 1.7763568394002505e-15 8.3266726846886741e-17 30         3          0               3                  8,9,10              false
```

Summary:
- The tracked private Fock-DIIS path with `max_history = 8` again converged
  under strict gates.
- It converged at iteration 34, matching pass 125.
- Final diagnostics:
  - density delta `8.3320104149464669e-09`;
  - commutator residual `2.9455127052713248e-09`;
  - spatial commutator residual `1.4727563526356624e-09`;
  - energy change `2.8421709430404007e-14`;
  - trace/idempotency at machine precision.
- DIIS behavior:
  - used count `30`;
  - fallback count `3`;
  - solve failures `0`;
  - coefficient pathologies `3`;
  - fallback iterations `8, 9, 10`.
- Final convergence was not blocked by recomputed diagnostics.

Recommended next pass:
- Make a small tracked private-control default change:
  - change the private Fock-DIIS default `max_history` from `6` to `8`;
  - keep `diis_regularization = 1.0e-12`;
  - keep `diis_coefficient_max_abs = 25.0`;
  - keep strict tolerances unchanged.
- Use focused synthetic validation only.
- Do not route-wire RHF, promote the fixture, or add public/report behavior.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none.
- simplified: this confirms history 8 as the only default-control change worth
  making next.
- quarantined: history-8 confirmation script/table/summary remain ignored
  `tmp/work` artifacts.
- not deleted because: prior ignored probes remain useful comparisons until the
  private default change is reviewed.
- exact remaining caller/blocker:
  tracked `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` converges
  the compact fixture with explicit `:fock_diis, max_history = 8`, while the
  private control default remains `max_history = 6`; the remaining blocker is a
  small tracked default update, not tolerance or route design.

-- repo-doer@macmini
