# Pass 249 response - independent H2 PQS H1 one-body seam

Implemented the narrow H1 reporting/artifact seam. The existing physical
gausslet H1 computation already worked for the independent H2 PQS final basis;
the code change only exposes existing H1 status facts in the artifact and keeps
the endpoint blocker from being overwritten by an unrequested H1-J payload.

## Result

Focused driver/artifact check with `run_h1=true`,
`run_h1_j=false`, and `run_private_rhf=false` passed:

```text
h1_status = materialized_pqs_physical_gausslet_h1_solve
h1_materialized = true
final_dimension = 471
h1_lowest_energy = -0.7946037173365885
h1_hamiltonian_matrix_finite = true
h1_hamiltonian_symmetry_error = 1.0658141036401503e-14
physics_endpoint_ready = false
physics_endpoint_blocker = missing_physical_gausslet_h1_j_builder
```

Support/final one-body statuses now appear in the artifact target group:

```text
support_kinetic_status =
  materialized_pqs_physical_gausslet_support_kinetic_matrix
support_electron_nuclear_status =
  materialized_pqs_physical_gausslet_support_electron_nuclear_by_center_matrix_set
final_kinetic_status =
  materialized_pqs_physical_gausslet_final_one_body_matrix
final_electron_nuclear_status =
  materialized_pqs_physical_gausslet_final_electron_nuclear_by_center
```

## Guard fields

The focused artifact check asserted:

```text
source_backed_fixed_source_oracle_used = false
fake_pqs_enabled = false
retained_transform_authority = pqs_source_box_construction
h1_j_materialized = false
private_rhf_materialized = false
exports_materialized = false
public_api = false
```

I did not enable H1-J, density-density interaction, RHF, supplements, CR2,
exports, or public API.

## Files changed

- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`

## Line budget and deletion offset

Scoped `src + test + bin` line budget:

```text
17 insertions, 26 deletions
net -9
```

Deletion/shrinkage report:

```text
deleted:
- stale unit-stage selected-terminal-sidecar count/materialization assertion
  cloud

simplified:
- endpoint blocker update now preserves an existing stronger blocker when the
  H1-J payload is merely not requested

quarantined:
- H1 remains private/diagnostic route readiness; no H1-J, RHF, supplement,
  CR2, export, or public API readiness was added

not deleted because:
- compact selected-terminal-sidecar status/private-development smoke remains as
  the active staged-test contract
- H1 target artifact fields are needed to validate the current H1 seam

exact remaining caller/blocker:
- independent H2 PQS H1 is materialized; endpoint readiness remains blocked by
  `:missing_physical_gausslet_h1_j_builder`
```

## Validation

Commands run:

```text
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
julia --project=. tmp/work/pass249_independent_h2_h1_probe.jl
git diff --check
git status --short --branch
```

Results:

```text
package load: passed, elapsed_s=58.3981265
focused independent H2 PQS H1 driver/artifact check: passed, elapsed_s=79.511969542
git diff --check: passed
```

The focused H1 driver check exceeded 60 seconds as expected because it rebuilds
the independent H2 final basis and materializes support/final one-body kinetic
and electron-nuclear operators while keeping H1-J and RHF off.

Final git status:

```text
## main...origin/main
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/nested/cartesian_unit_stage_low_order_policy_runtests.jl
```

-- repo-doer@macmini
