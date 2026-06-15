# Pass 252 response - independent H2 PQS private RHF diagnostic

Implemented the narrow private RHF diagnostic seam repair. No new solver was
added.

## Initial blocker

The first focused RHF run reached the H1-J route but blocked before execution:

```text
private_rhf/input_contract_status = blocked_pqs_physical_gausslet_rhf_input_contract
private_rhf/input_contract_blocker = unsupported_physical_gausslet_fixture_role
private_rhf/execution_status = blocked_pqs_physical_gausslet_private_rhf_execution
```

After recognizing the new independent-H2 stage roles, the input contract became
available but the shared RHF helper still rejected the independent route kind:

```text
private_rhf/input_contract_status = available_pqs_physical_gausslet_rhf_input_contract
private_rhf/execution_status = blocked_pqs_physical_gausslet_private_rhf_execution
private_rhf/execution_blocker = missing_rhf_input_contract
```

## Repair

Narrow source changes:

- Added a private helper that classifies the independent H2 PQS artifact roles:
  readiness, final-basis diagnostic, H1 diagnostic, and H1-J density diagnostic.
- Reused that classifier in artifact reporting so stage-specific independent
  inputs still get `source_backed_fixed_source_oracle_used=false` and
  `retained_transform_authority=:pqs_source_box_construction`.
- Allowed the physical-gausslet RHF input contract and shared private RHF helper
  to accept the independent H2 PQS route kind when the artifact role is one of
  those independent diagnostic roles.

## Final RHF result

Focused driver/artifact check using
`h2_pqs_q5_independent_source_box_r4_h1_j.jl` plus
`run_private_rhf=true`, `private_rhf_electron_count=2` passed:

```text
route/artifact_role = independent_h2_pqs_h1_j_density_diagnostic
route/fake_pqs_enabled = false
route/source_backed_fixed_source_oracle_used = false
route/retained_transform_authority = pqs_source_box_construction

private_rhf/input_contract_status =
  available_pqs_physical_gausslet_rhf_input_contract
private_rhf/execution_status =
  materialized_pqs_physical_gausslet_private_rhf_execution
private_rhf/requested = true
private_rhf/executed = true
private_rhf/materialized = true
private_rhf/converged = true
private_rhf/electron_count = 2
private_rhf/occupation_policy = closed_shell_rhf
```

RHF diagnostics:

```text
total_energy = -1.1589735957658853
one_body_energy = -1.5609752182694798
two_body_energy = 0.4020016225035945
iteration_count = 8
density_trace = 1.9999999999999993
idempotency_residual = 3.122502256758253e-17
commutator_residual = 3.965093933744335e-9
energy_delta = 2.5224267119483557e-13
final_density_one_step_consistency_status = reviewed_recomputed
```

The route remains diagnostic-only:

```text
physics/endpoint_ready = false
physics/endpoint_blocker = missing_h2_gausslet_only_reference_comparison
route/exports_materialized = false
route/public_api = false
```

I did not compare this independent PQS RHF value to supplemented WL/QW
references and did not use the fake-PQS route as evidence.

## Files changed

- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`

## Line budget and deletion offset

Scoped `src + test + bin` line budget:

```text
29 insertions, 31 deletions
net -2
```

Deletion/shrinkage report:

```text
deleted:
- stale assembly-stage terminal RouteCore pair/readiness mirror assertions
- stale assembly-stage summary identity mirrors for terminal scaffold and unit
  inventory

simplified:
- independent-H2 artifact-role recognition is centralized in one private helper
  instead of remaining a single readiness-role check

quarantined:
- private RHF remains a diagnostic/private route seam; no supplements, CR2,
  export/HamV6, public API, or public solver readiness was added

not deleted because:
- compact terminal assembly route/scaffold/deferral smoke remains active
- private RHF artifact fields are existing diagnostic outputs needed for this
  seam validation

exact remaining caller/blocker:
- private RHF materializes for the independent H2 PQS diagnostic route; endpoint
  readiness remains blocked by `:missing_h2_gausslet_only_reference_comparison`
```

## Validation

Commands run:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. tmp/work/pass252_independent_h2_private_rhf_probe.jl
git diff --check
git status --short --branch
```

Results:

```text
package load: passed
focused independent H2 PQS private RHF driver/artifact check: passed, elapsed_s=83.717773666
git diff --check: passed
```

The focused private RHF run exceeded 60 seconds as expected because it rebuilds
the independent H2 final basis, H1, H1-J density interaction, and then runs the
private RHF diagnostic. I did not run broad tests.

Final git status:

```text
## main...origin/main
 M src/pqs_multilayer_complete_core_shell_rhf.jl
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
```

-- repo-doer@macmini
