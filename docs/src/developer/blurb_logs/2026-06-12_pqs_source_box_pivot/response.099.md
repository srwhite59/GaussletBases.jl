Pass 099 response

No-edit cleanup audit completed.

Files changed:

- none before writing this response

Read-only inspection commands used:

- `rg --files tmp/work`
- `rg -n ... tmp/work test/nested/pqs_direct_retained_final_h1_runtests.jl src/... docs/... .agent_handoffs -g '!*.tmp'`
- targeted `sed -n ...` reads of the requested source/test/log surfaces
- `git check-ignore -v ...` for concrete ignored probe candidates
- `git ls-files ...` for tracked surfaces
- `git status --short --branch`

Cleanup inventory

delete-local-ignored-now:

- `tmp/work/pqs_direct_retained_final_h1_timing_probe.jl`
- `tmp/work/pqs_direct_retained_final_h1_probe_summary.txt`
  - reason: pass 096 timing question is answered; focused test is short enough, and the timing probe duplicates the tracked fixture construction without adding current contract coverage.
  - remaining caller/blocker: none; ignored local timing artifact only.

- `tmp/work/pqs_complete_core_shell_h1_probe.jl`
- `tmp/work/pqs_complete_core_shell_h1_probe_summary.txt`
  - reason: the tracked focused H1 test now covers the source-plan/final-basis/H1 path, and pass 095 moved the same final-basis/H1 construction into the driver-owned diagnostic payload.
  - remaining caller/blocker: none; ignored local development probe only.

- `tmp/work/pqs_complete_core_shell_final_density_j_probe.jl`
- `tmp/work/pqs_complete_core_shell_final_density_j_probe_summary.txt`
  - reason: this records the rejected final signed IDA-weight division convention. The accepted H1/J path now uses the pre-final positive-weight gauge, and the rejected convention should not continue as active probe pressure.
  - remaining caller/blocker: none; ignored local convention-rejection artifact only.

delete-tracked-next-pass:

- none as an immediate tracked deletion.
  - reason: there is not yet a compact tracked driver H1/J smoke, so deleting tracked focused H1 assertions now would remove the only small tracked source-plan/final-basis/H1 contract gate.

- candidate after a driver H1/J smoke exists: narrow parts of `test/nested/pqs_direct_retained_final_h1_runtests.jl`:
  - `explicit_box_plan` / `explicit_box_final_basis` construction and equivalence assertions.
  - support nuclear duplicate construction comparing `gaussian_factor_terms_by_center` and `axis_layers`.
  - reason: both are development-era seam checks. They helped validate the explicit-box bridge and nuclear-factor convention, but they are not the driver-owned H1/J route authority.
  - remaining caller/blocker: keep until a tracked driver H1/J smoke covers the accepted path and manager decides whether explicit-box/nuclear convention checks should be retained as oracle/debug or removed.

quarantine/debug-only:

- `tmp/work/pqs_complete_core_shell_density_gauge_probe.jl`
- `tmp/work/pqs_complete_core_shell_density_gauge_probe_summary.txt`
- `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl`
- `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe_summary.txt`
  - reason: these contain the empirical convention trail for the accepted pre-final positive-weight density gauge and fixed-oracle J comparison. Their construction is now duplicated by `pqs_multilayer_complete_core_shell_h1_j_payload(...)` plus the driver density-input helper, but the summaries remain useful debug history until one compact tracked driver H1/J smoke exists.
  - remaining caller/blocker: no code caller; blocked only by lack of a tracked driver H1/J smoke carrying final dimension, H1 energy, self-Coulomb, gauge, and non-promotion flags.

- `tmp/work/pqs_complete_core_shell_final_ida_weights_probe.jl`
- `tmp/work/pqs_complete_core_shell_final_ida_weights_probe_summary.txt`
  - reason: useful for future IDA-weight semantics and for remembering that signed final projected weights are not the H1/J diagnostic gauge. Not active H1/J route pressure.
  - remaining caller/blocker: keep debug-only until IDA/MWG work is explicitly in scope.

- `tmp/work/pqs_complete_core_shell_one_body_oracle_probe.jl`
- `tmp/work/pqs_complete_core_shell_one_body_oracle_probe_summary.txt`
  - reason: useful oracle/debug trail for the corrected PGDG Gaussian-factor nuclear convention and fixed-block comparison. The tracked H1 test already checks the accepted convention at smaller cost.
  - remaining caller/blocker: keep only if manager wants a local historical oracle; otherwise can be deleted after compact driver smoke and retained H1 test cleanup.

- `tmp/work/pqs_multilayer_shell_side13_h1_smoke.jl`
- `tmp/work/pqs_multilayer_shell_side13_h1_smoke_summary.txt`
- `tmp/work/pqs_multilayer_shell_side13_j_probe.jl`
- `tmp/work/pqs_multilayer_shell_side13_j_probe_summary.txt`
- `tmp/work/pqs_multilayer_shell_side13_rhf_probe.jl`
- `tmp/work/pqs_multilayer_shell_side13_rhf_probe_summary.txt`
- `tmp/work/he_side13_h1j_probe.jl`
- `tmp/work/he_side13_h1j_probe_result.txt`
  - reason: larger side-13/Hydrogenic/He/RHF exploratory probes. They are not active driver H1/J contract tests and should not pressure this cleanup path. RHF remains explicitly out of scope.
  - remaining caller/blocker: keep quarantined only for future physics/RHF review; delete if manager wants to purge ignored exploratory artifacts.

keep-active-contract:

- `src/pqs_source_box_route_driver_helpers.jl`
  - keep `_PQSCompleteCoreShellDiagnosticRoutePayload`, `_pqs_source_box_route_driver_complete_core_shell_density_inputs(...)`, and `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)`.
  - reason: this is now the driver-owned private H1/J diagnostic path.

- `src/pqs_multilayer_complete_core_shell_h1.jl`
  - keep `pqs_multilayer_complete_core_shell_final_basis(...)`, `pqs_multilayer_complete_core_shell_h1_payload(...)`, and `pqs_multilayer_complete_core_shell_h1_j_payload(...)`.
  - reason: these are the active module-level construction payloads consumed by the driver diagnostic path.

- `src/pqs_multilayer_support_density.jl`
  - keep `pqs_multilayer_support_weights(...)` and `pqs_multilayer_support_pair_raw_numerator_matrix(...)`.
  - reason: these are the accepted support-density input builders for the private H1/J diagnostic. They enforce raw numerator terms and do not promote density-normalized pair terms as authority.

- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - keep for now.
  - reason: it is the small tracked gate for shellification-backed source plan, complete core/shell final basis, H1 payload, nuclear convention, final dimension `223`, H1 energy, and fixed-block H1 oracle comparison.
  - exact remaining blocker: no tracked driver H1/J smoke exists yet; deleting or heavily shrinking this test now would leave no small tracked seam validation.

Probe-local constructions now duplicated by the driver-owned H1/J diagnostic:

- repeated `_pqs_h1_test_bundle(7)` / `_CartesianNestedAxisBundles3D(...)` fixture construction;
- shellification/lowering/region/source-plan construction for the same 7x7x7, 5x5x5 fixture;
- complete core/shell final-basis construction;
- H1 payload construction;
- support weights from PGDG axis weights;
- raw support pair numerator matrix from `pair_factor_terms_raw`;
- pre-final positive-weight density interaction;
- H1 lowest-orbital self-Coulomb diagnostic.

Comparisons that remain valuable as oracle/debug:

- fixed-block H1 energy comparison in `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
- PGDG Gaussian-factor nuclear convention check, but likely not in both tracked test and ignored one-body oracle probe forever;
- fixed-oracle J comparison from the density gauge/pre-final consumption probes, until one compact driver H1/J smoke records the accepted self-Coulomb and gauge;
- side13/Hydrogenic comparisons only as quarantined future physics probes, not routine route-driver cleanup gates.

Compact driver H1/J smoke recommendation:

- Yes, one compact driver H1/J smoke should be added later, after manager assigns it.
- It should exercise the accepted one-center source-box dry-run/assembly path that manager validated in pass 097, then assert only compact facts:
  - `complete_core_shell_h1_j_diagnostic_status == :materialized_pqs_multilayer_complete_core_shell_h1_j_payload`;
  - internal `complete_core_shell_diagnostic_route_payload.density_inputs.status == :available_complete_core_shell_density_inputs`;
  - final dimension `223`;
  - H1 energy and self-Coulomb snapshot for the accepted driver fixture;
  - density gauge `:pre_final_localized_positive_weight`;
  - RHF/GTO/export/artifact flags remain false.
- It should replace ignored local H1/J and density-consumption probe pressure, and it can justify a later shrink of explicit-box and duplicate nuclear-convention assertions in the tracked H1 fixture. It should not replace the scientific/oracle H1 comparison until manager explicitly accepts that contraction.

Git status:

```text
## main...origin/main
```

Deletion/shrinkage forecast:

- deleted: none in this no-edit pass.
- simplified: next cleanup can delete local ignored timing/H1/rejected-final-density probes and reduce duplicate tracked fixture assertions once a compact driver H1/J smoke exists.
- quarantined: pre-final density gauge history, final IDA-weight semantics, one-body oracle history, side13 J/RHF exploratory probes.
- not deleted because: no-edit audit; tracked focused H1 test remains the only small tracked seam validation; compact driver H1/J smoke is not present yet.
- exact remaining caller/blocker: the remaining blocker for tracked shrinkage is absence of a compact tracked driver H1/J smoke over the accepted one-center source-box diagnostic path.

-- repo-doer@macmini
