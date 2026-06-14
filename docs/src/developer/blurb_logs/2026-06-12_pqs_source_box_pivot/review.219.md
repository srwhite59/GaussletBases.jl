Review 219 - accepted

Pass 219 added a compact no-matrix WL/old-QW H2 gausslet-only reference
candidate for the physical H2 q5 PQS endpoint.

Accepted commit:

```text
5e864d69 Add H2 WL reference candidate readiness
```

Result:

```text
comparison/wl_reference_candidate_status =
  :available_wl_h2_gausslet_only_reference_candidate
comparison/wl_reference_final_dimension = 463
comparison/wl_reference_retained_transform_kind =
  :white_lindsey_old_qw_gausslet_retained_transform
comparison/wl_reference_supplement_policy = :none
comparison/wl_reference_label = "WL/QW H2 R=4 gausslet-only 463"
comparison/old_supplemented_wl_qw_scalar_references_blocked = true
comparison/ready = false
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_wl_h2_gausslet_only_reference_values
```

Review notes:

- The new candidate is readiness metadata only. It does not materialize WL
  matrices or scalar values.
- The candidate checks the common H2 physical support/retained inventory and
  keeps old supplemented WL/QW scalar references blocked for this no-supplement
  endpoint.
- The route-state global safe-term test shrink was acceptable: it collapsed
  repeated assertions into helpers while preserving overlap, kinetic, position,
  and x2 smoke coverage.

Validation:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. test/nested/cartesian_pair_block_route_state_global_safe_terms_runtests.jl
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
git diff --check
git diff --cached --check
```

The focused H2 endpoint passed 87/87 in about 87.2 seconds.

Line budget:

```text
source/test/bin added:   315
source/test/bin deleted: 316
net:                    -1
```

Deletion/shrinkage:

```text
deleted: none
simplified:
  test/nested/cartesian_pair_block_route_state_global_safe_terms_runtests.jl
  collapsed repeated route-state adapter assertions into shared helpers
quarantined:
  old supplemented WL/QW H2 scalar references remain blocked as direct
  no-supplement comparison values
not deleted because:
  one-center driver smoke still has a live integration caller
exact remaining caller/blocker:
  :missing_wl_h2_gausslet_only_reference_values
```

-- repo-manager@macmini
