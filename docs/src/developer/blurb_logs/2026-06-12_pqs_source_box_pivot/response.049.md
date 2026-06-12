Pass 049 response: PQS complete core-shell final IDA weights

Summary:
- Added the narrow final IDA weight seam for the complete core/shell final basis.
- Created and ran an ignored probe for the 223-function fixture.
- Did not build density-density, RHF, GTO, driver wiring, exports, or artifacts.

Files edited:
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  - added `pqs_complete_core_shell_final_ida_weights(final_basis, support_weights; near_zero_atol, metadata)`
- Ignored artifacts created:
  - `tmp/work/pqs_complete_core_shell_final_ida_weights_probe.jl`
  - `tmp/work/pqs_complete_core_shell_final_ida_weights_probe_summary.txt`

Helper behavior:
- Validates:
  - `pqs_complete_core_shell_final_basis` status/materialization;
  - support-weight length equals `core_support_count + shell_support_count`;
  - support weights are finite;
  - final coefficient row count matches support count;
  - projected final weights are finite.
- Computes:
  - `final_ida_weights = transpose(final_coefficients) * support_weights`
  - compact min/max/sum/sign/near-zero diagnostics.
- Records nonclaim flags:
  - fixed-block weights are not active authority;
  - raw source weights are not final weights;
  - boundary-shell diagnostic weights are not final weights;
  - density-density/RHF/driver/export/artifact are not materialized.

Probe fixture:
- `current_box = (1:7, 1:7, 1:7)`
- `inner_box = (2:6, 2:6, 2:6)`
- `raw_source_dims = (5, 5, 5)`
- final dimension: `223`
- support row order: `core_then_shell`
- support weight source: product PGDG integral weights on the complete support row order.

Final IDA weight numbers:
- support weight count: `343`
- final IDA weight count: `223`
- support weight sum: `1.5847117713178820e+03`
- final weight min: `-1.8967313490013488e+01`
- final weight max: `5.8136426512094616e+01`
- final weight sum: `-6.2222631936233967e+01`
- positive count: `111`
- negative count: `112`
- near-zero count at `1.0e-12`: `0`

Fixed-block oracle comparison:
- Same-geometry fixed-block weights were reachable as oracle-only.
- Native fixed-block weights:
  - count: `223`
  - min: `3.9540709965461518e-01`
  - max: `1.3495587524797234e+01`
  - sum: `1.2745426664962772e+03`
- After current-gauge alignment through the fixed/current cross overlap:
  - min: `-1.8967313490013606e+01`
  - max: `5.8136426512094779e+01`
  - sum: `-6.2222631936235622e+01`
  - max delta versus the new final IDA weights: `6.8922645368729718e-13`
- The fixed/current final-subspace identity error remains `5.8875461168328425e-04`, so the fixed block remains an oracle/comparison surface, not active route data.

Validation:
- `julia --project=. tmp/work/pqs_complete_core_shell_final_ida_weights_probe.jl`
  - passed
  - note: this process first precompiled `GaussletBases` after the source edit.
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 28 checks in 2.8s
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
- `git diff --check`
  - passed

Deletion/shrinkage report:
- No old production code or tests became obsolete in this pass. This adds the missing final-weight seam needed before any density-density/RHF work can be interpreted.
- Nothing was deleted or simplified because there was no existing complete core/shell final IDA weight helper to replace.
- No permanent test was added. The validation is an ignored `tmp/work` probe plus the existing compact H1 gate after the source edit.
- Remaining stale or duplicate surfaces to retire next:
  - any PQS weight probe or note that treats raw source weights or boundary-shell diagnostic weights as final IDA weights;
  - future density-density work should consume `pqs_complete_core_shell_final_ida_weights` rather than deriving a separate final-weight convention;
  - fixed-block weights should remain oracle-only and should not be threaded into the active route.

-- repo-doer@macmini
