# Pass 241 response - independent H2 PQS shared-shell realization payload

Implemented the narrow shared-shell realization payload for the independent H2
PQS route. The complete source plan remains blocked; no final basis or physics
payload was assembled.

## Files changed

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`

## Shared-Shell Realization

Added private helpers:

```text
_pqs_source_box_route_driver_independent_h2_shared_shell_region_descriptors(...)
_pqs_source_box_route_driver_independent_h2_shared_shell_realization_payload(...)
_pqs_source_box_route_driver_independent_h2_shared_shell_realization(...)
```

The payload materializes only the two shared-shell realization blocks:

```text
shared_shell_realization_status = :available_independent_pqs_shared_shell_realization_payload
shared_shell_realization_blocker = nothing
shared_shell_realization_counts = (98, 98)
source_plan_blocker = :missing_independent_pqs_complete_core_shell_source_plan_assembly
```

Per shared shell, the internal payload carries:

```text
raw_source_plan
retained_rule
support_indices
support_states
shell_projection
lowdin_cleanup
shell_final_coefficients
retained_count
coefficient_shape
realized_overlap_identity_error
```

`:atom_contact_core` remains descriptor/identity-like only. No dense core
identity/source matrix was materialized.

## Artifact / Report Fields

Added compact target fields:

```text
target/shared_shell_realization_status
target/shared_shell_realization_blocker
target/shared_shell_realization_counts
target/shared_shell_realization_identity_errors
```

Existing report-facing route shape remains blocked:

```text
target/source_plan_status = :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
target/source_plan_blocker = :missing_independent_pqs_complete_core_shell_source_plan_assembly
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_independent_pqs_complete_core_shell_source_plan_assembly
```

## Fake / Source-Backed Evidence

The independent route still gates out the source-backed candidate. The focused
artifact check captured driver stdout and asserted:

```text
!occursin("diatomic.fixed_source", driver_output)
fake_pqs/enabled == false
route/source_backed_fixed_source_oracle_used == false
```

Result:

```text
source_backed_candidate_output_seen=false
```

No fake-PQS/WL fixed-source retained transforms or coefficient matrices were
used.

## Validation

Package load:

```text
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
```

Final result after the last source edit:

```text
load ok
elapsed_s=59.917775084
```

Focused independent-input artifact/readiness check:

```text
julia --project=. -e 'using JLD2; ... include("bin/cartesian_ham_builder.jl") ...'
```

Result:

```text
shared_shell_realization_status=available_independent_pqs_shared_shell_realization_payload
source_plan_blocker=missing_independent_pqs_complete_core_shell_source_plan_assembly
shared_shell_realization_counts=(98, 98)
shared_shell_realization_identity_errors=(1.2667660635192445e-14, 5.2966890053049083e-14)
physics_endpoint_blocker=missing_independent_pqs_complete_core_shell_source_plan_assembly
source_backed_candidate_output_seen=false
elapsed_s=78.1625405
```

The focused artifact check exceeded 60 seconds because it ran the driver route
construction after source edits.

Whitespace:

```text
git diff --check
```

passed.

## Scoped Line Budget

Scoped `src + test + bin` line budget:

```text
251	2	src/pqs_source_box_diatomic_complete_core_shell.jl
21	0	src/pqs_source_box_route_driver_helpers.jl
8	0	src/pqs_source_box_route_driver_reporting.jl
```

Total:

```text
280 added / 2 deleted, net +278
```

Line-budget exception granted by user/manager in:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/blurb.241_exception.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.241_attention.md
```

Reason: this pass stayed on the narrow shared-shell realization seam, while
pass 239 already retired the large projected-shell integration pressure. Chasing
net-negative here would have required broader cleanup that the pass explicitly
did not assign.

## Forbidden Surfaces Avoided

No complete source-plan assembly, final basis, H1, H1-J, RHF, supplements, CR2,
export, public API, fake-PQS/WL coefficient paths, or fixed-source retained
transforms were added.

## Git Status

```text
## main...origin/main [ahead 1]
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
```

## Deletion/Shrinkage Accounting

deleted:
- none.

simplified:
- the independent route blocker advanced from broad numerical materialization
  to `:missing_independent_pqs_complete_core_shell_source_plan_assembly`;
- compact shared-shell realization status/count/error fields replace the need
  for a wider report-key cloud in this pass.

quarantined:
- fake-PQS H2 463 source-backed reproduction remains separate;
- old projected-shell machinery is used only as an internal mathematical
  adapter fed by route-owned independent H2 support/source boxes.

not deleted because:
- the line-budget exception explicitly avoided broader cleanup in this
  numerical seam;
- `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl` is
  classified for careful shrinkage, not wholesale deletion;
- fake-PQS endpoint coverage and current raw product/final-basis/lowering
  contract tests still protect live contracts.

exact remaining caller/blocker:
- `:missing_independent_pqs_complete_core_shell_source_plan_assembly`

-- repo-doer@macmini
