Pass 125 response - tracked Fock-DIIS asymptote probe

Scope:
- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No route wiring/report fields/public API/exports/artifacts/GTO/IDA/MWG.
- No fixture promotion.
- No production tolerance/default changes.

Artifact paths:
- Probe script:
  `tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_asymptote_probe.jl`
- Summary:
  `tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_asymptote_probe_summary.txt`
- Table:
  `tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_asymptote_probe.tsv`

Command/results:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_asymptote_probe.jl`
  - elapsed: `103.647220416` seconds
- `git status --short --branch`
  - `## main...origin/main`

Probe table:
```text
max_history  converged  iter  final_or_last_total_energy  density_change       update_density_change  energy_change          commutator_residual  spatial_commutator  trace_error           idempotency_error     diis_used  fallbacks  solve_failures  coeff_pathologies  final_gate_blocked
6            false      200   -10.032119190006792         1.2225907675134096e-7 3.4892496003635642e-8 2.6645352591003757e-13 3.120558440816712e-8 1.560279220408356e-8 0                    5.5511151231257827e-17 199        0          0               0                  false
8            true       34    -10.032119189804888         8.3320104149464669e-9 3.0626759672802706e-7 2.8421709430404007e-14 2.9455127052713248e-9 1.4727563526356624e-9 1.7763568394002505e-15 8.3266726846886741e-17 30         3          0               3                  false
```

Summary:
- `max_history = 6` did not converge by 200 iterations.
  - It remained near the same plateau as pass 124:
    - density delta `1.22e-7`;
    - commutator residual `3.12e-8`.
  - No fallbacks or DIIS failures occurred.
- `max_history = 8` converged by iteration 34 under the strict gates.
  - density delta `8.33e-9`;
  - commutator residual `2.95e-9`;
  - energy change `2.84e-14`;
  - trace/idempotency at machine precision.
  - It had 3 coefficient-pathology fallbacks, no solve failures, and still
    converged.
- Neither case was blocked by final recomputed diagnostics. The history-8 case
  materialized cleanly after final diagnostics passed.

Recommended next pass:
- Do one local confirmatory tracked Fock-DIIS probe before changing defaults:
  - `max_history = 8`;
  - `max_iterations = 100`;
  - same strict gates and regularization;
  - preserve coefficient-pathology guard;
  - optionally record the fallback iterations.
- If the confirmatory run again converges near iteration 34 with bounded
  fallbacks, then a small tracked follow-up can change the private Fock-DIIS
  default `max_history` from 6 to 8.
- Do not loosen tolerances and do not route-wire RHF yet.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none.
- simplified: this narrows the next decision to confirming `max_history = 8`,
  not changing tolerances or route behavior.
- quarantined: asymptote probe script/table/summary remain ignored `tmp/work`
  artifacts.
- not deleted because: prior ignored probes remain useful comparisons until the
  private control default is settled.
- exact remaining caller/blocker:
  tracked `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` converges
  the compact fixture with `:fock_diis, max_history = 8`, but default
  `max_history = 6` still blocks at the strict `1e-8` gates; a confirmatory
  local history-8 run should precede any default change.

-- repo-doer@macmini
