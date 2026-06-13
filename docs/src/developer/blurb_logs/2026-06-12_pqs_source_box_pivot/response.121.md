Pass 121 response - local compact PQS Fock-DIIS experiment

Scope:
- Local ignored `tmp/work` scripts/summaries only.
- No tracked source/test/doc edits.
- No production/private DIIS code.
- No route wiring, report fields, public API, exports, artifacts, GTO, IDA/MWG,
  or fixture promotion.

Artifact paths:
- Script:
  `tmp/work/pqs_compact_route_owned_rhf_fock_diis_experiment.jl`
- Summary:
  `tmp/work/pqs_compact_route_owned_rhf_fock_diis_experiment_summary.txt`
- Table:
  `tmp/work/pqs_compact_route_owned_rhf_fock_diis_experiment.tsv`

Commands/results:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_fock_diis_experiment.jl`
  - First attempt exposed a typo in the ignored local script regularization
    path (`UniformScaling` broadcast). I fixed only that local script.
  - Corrected run elapsed: `114.94605325000001` seconds.
- `git status --short --branch`
  - `## main...origin/main`

Sweep table:
```text
max_history  reg       converged  iter  final_total_energy   fixed_point_dP       update_dP            dE_last             comm_residual        spatial_comm       trace_err            idem_err             diis_failures  coeff_pathologies  coeff_max_abs       coeff_l1            last_update
4            1e-12     false      100   -10.032183104737362  4.6552219618439805e-5 4.6527879701524499e-5 1.8985341121435795e-7 1.3857439112453795e-5 6.9287195562268977e-6 4.4408920985006262e-16 5.5511151231257827e-17 0              92                 39.047867290411766 103.38636684853395 diis_coefficient_pathology_fallback
6            1e-12     false      100   -10.032119189981447  1.1441678879342554e-7 3.2634555308108659e-8 2.4158453015843406e-13 2.9183295714024782e-8 1.4591647857011391e-8 4.4408920985006262e-16 5.5511151231257827e-17 0              0                  0.16690620700826006 1.0               fock_diis
4            1e-10     false      100   -10.032173934363737  4.5073946112139662e-5 1.4818217207557272e-5 5.619017251490277e-8  1.3266824339955297e-5 6.6334121699776483e-6 8.8817841970012523e-16 4.8572257327350599e-17 0              0                  0.41209731630803104 1.0               fock_diis
6            1e-10     false      100   -10.03217206552965   4.4692302840348974e-5 9.4794662461039625e-6 3.5372897499996725e-8 1.3127585394331898e-5 6.5637926971659488e-6 1.7763568394002505e-15 4.163336342344337e-17 0              0                  0.3399675573396599  1.0131077475808445 fock_diis
```

Best run:
- `max_history = 6`, `diis_regularization = 1.0e-12`.
- It did not meet the strict `1.0e-8` density/residual thresholds within 100
  iterations, but it reduced the commutator residual from the previous
  `~1.3e-5` scale to `2.9183295714024782e-8`.
- It also reduced the fixed-point density delta to
  `1.1441678879342554e-7` and the DIIS update delta to
  `3.2634555308108659e-8`.
- Energy change was essentially stalled at `2.4158453015843406e-13`.
- Trace and idempotency remained at machine precision.
- No DIIS solve failures or coefficient pathologies occurred in this best run.

Interpretation:
- Fock DIIS substantially improves the actual commutator residual for this
  compact route-smoke fixture.
- The `max_history = 4`, `reg = 1e-12` case was unstable under the local
  coefficient guard and fell back 92 times due to coefficient pathologies
  (`max_abs ≈ 39`, `l1 ≈ 103`).
- The stronger regularization `1e-10` was stable but did not materially improve
  the residual over fixed-point behavior.
- The useful configuration is therefore not generic "turn on DIIS"; it is a
  bounded Fock-DIIS path with history 6, small regularization, and coefficient
  pathology fallback.

Recommended next pass:
- Add a bounded private SCF-control implementation pass, still private and
  diagnostic-only:
  - compact control payload;
  - `mixing_kind = :fock_diis`;
  - `max_history = 6`;
  - `diis_regularization = 1.0e-12`;
  - `diis_start_iteration = 2`;
  - coefficient pathology guard, e.g. fallback if max abs coefficient exceeds
    about `25`;
  - convergence requiring density delta, energy delta, commutator residual,
    trace error, and idempotency error.
- Keep route-driver integration blocked.
- Keep tracked tests synthetic; use ignored compact route probes for real
  fixture performance/residual validation.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none.
- simplified: no tracked code changed; the experiment narrows the next
  implementation to bounded Fock DIIS rather than generic mixing.
- quarantined: DIIS experiment script/table/summary remain ignored `tmp/work`
  artifacts.
- not deleted because: prior ignored residual/damping probes remain useful
  comparisons until the private SCF-control implementation is reviewed.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` still has no
  built-in control payload or DIIS history. Local Fock DIIS can reduce the
  compact fixture commutator residual to `~2.9e-8`, but no tracked private
  implementation exists yet.

-- repo-doer@macmini
