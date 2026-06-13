Pass 133 response - Be2 PQS Ham payload blocker fingerprint

Commit:
- `3192f810 Fingerprint Be2 PQS Ham payload blocker`

Observed Be2 PQS Ham payload fingerprint:
- Staged Be2 PQS route reaches `cartesian_assembly(...)`.
- `assembly.object_kind == :cartesian_assembly`.
- `assembly.route_skeleton.route_family === :pqs_source_box`.
- `payload = assembly.complete_core_shell_diagnostic_route_payload` exists.
- `ham = payload.complete_core_shell_ham_payload` exists.

Observed payload status:
```text
payload.status = :blocked_missing_complete_core_shell_h1_j_route_inputs
payload.blocker = :missing_complete_core_shell_h1_j_route_inputs
payload.missing_inputs = (
    :pqs_multilayer_shell_region_plan,
    :pqs_multilayer_shell_source_plan,
    :pqs_multilayer_complete_core_shell_final_basis,
    :pqs_multilayer_complete_core_shell_h1_payload,
    :axis_weights,
    :raw_pair_factor_terms,
    :coulomb_expansion,
)
payload.source_payload.status = :blocked_missing_complete_core_shell_source_plan_inputs
payload.h1_j_payload.status = :blocked_missing_complete_core_shell_h1_j_route_inputs
```

Observed Ham payload status:
```text
ham.status = :blocked_complete_core_shell_ham_payload
ham.blocker = :missing_complete_core_shell_ham_inputs
ham.missing_inputs = (
    :pqs_multilayer_complete_core_shell_final_basis,
    :pqs_multilayer_complete_core_shell_h1_payload,
    :pqs_complete_core_shell_final_one_electron_hamiltonian,
    :complete_core_shell_density_inputs,
    :complete_core_shell_h1_j_diagnostic_payload,
    :pqs_complete_core_shell_pre_final_density_interaction,
)
```

Files changed:
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Validation:
```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```
- Passed: 16/16.
- Reported test time: `10.8s`.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```
- Passed: `load ok`.

```text
git diff --check
```
- Passed.

Focused pre-test probe:
- A focused Julia-level `@elapsed` status probe reported `elapsed_s=8.558704917`.
- This did not exceed 60s and was not precompilation-dominated.

Git status after commit:
```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:
- deleted: none.
- simplified: the Be2 PQS blocker is now a compact tracked fingerprint instead of an inferred readiness gap from broad report aliases.
- quarantined: test is standalone and not added to the default nested runner; no WL payload, export/artifact, public API, RHF/SCF, hfdmrg, CR2 execution, fixture promotion, or IDA/MWG semantic change.
- not deleted because: existing report/materialization tests remain compatibility coverage and this pass only adds the missing focused fingerprint.
- exact remaining caller/blocker: Be2 PQS still lacks a route-owned complete core/shell source plan, final basis, H1 payload/final Hamiltonian, density inputs, H1/J diagnostic payload, and pre-final density interaction for the private Ham payload seam.

-- repo-doer@macmini
