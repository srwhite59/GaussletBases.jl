Review 227 - accepted

Pass 227 correctly materialized the route-owned, matrix-free H2 GTO supplement
representation.

Accepted behavior:

```text
supplement_policy = :none
  representation status: :not_requested
  orbital count: 0

supplement_policy = :mwg_residual_gto
  request status: :available_pqs_physical_gausslet_supplement_request
  representation status: :available_pqs_physical_gausslet_gto_supplement_representation
  representation object kind: :cartesian_gaussian_shell_supplement_representation
  center count: 2
  orbital count: 18
  matrices materialized: false
  provider blocks materialized: false
```

The preflight blocker advanced as intended:

```text
:missing_gto_supplement_representation
  -> :missing_provider_gto_supplement_blocks
```

No GTO/GTO, mixed gausslet/GTO, MWG residual, density-density, supplemented
scalar, public API, export, CR2, HFDMRG, DMRG, or production route behavior was
added.

The source/test/bin line budget stayed negative:

```text
added:   351
deleted: 354
net:      -3
```

Manager validation:

```text
git diff --check
  passed

julia --project=. -e 'using GaussletBases; println("load ok")'
  passed

julia --project=. -e 'Meta.parse("begin\n" * read("test/nested/pqs_source_box_route_driver_report_runtests.jl", String) * "\nend"); println("route report test parse ok")'
  passed

julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
  passed: 115/115
  elapsed_s = 143.38063675
```

I did not rerun the mixed-GTO provider test locally; doer reported it passed
1889 checks in 85.594604875 seconds. The H2 endpoint and package load cover the
new route wiring.

One cleanup note for the next pass: now that `supplement_representation` exists
as its own artifact group, the request payload should not remain the long-term
authority for representation status. The next pass should shrink that duplicate
request/report surface while auditing the provider-block seam.

-- repo-manager@macmini
