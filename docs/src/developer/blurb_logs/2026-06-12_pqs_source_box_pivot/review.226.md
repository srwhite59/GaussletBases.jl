Review 226 - accepted

Pass 226 added the right supplement boundary: a compact, private, matrix-free
request payload for the H2 R=4 physical gausslet endpoint.

Accepted behavior:

```text
supplement_policy = :none
  request/preflight status: :not_requested
  blocker: nothing

supplement_policy = :mwg_residual_gto
  request status: :blocked_pqs_physical_gausslet_mwg_residual_gto_request
  preflight status: :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
  blocker: :missing_gto_supplement_representation
```

No GTO/GTO, mixed gausslet/GTO, MWG residual, density-density, supplemented
scalar, CR2, HFDMRG, export, or public API behavior was added. The accepted
no-supplement H2 endpoint values and WL/PQS deltas remain the live physics
gate.

The deletion work is also accepted. The pass removed 751 lines of stale
route-report materialization/assertion pressure and kept the scoped
source/test/bin diff line-negative:

```text
added:   336
deleted: 781
net:    -445
```

Manager validation:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
  passed

julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
  passed: 121/121
  elapsed_s = 139.531038583

julia --project=. -e 'Meta.parse("begin\n" * read("test/nested/pqs_source_box_route_driver_report_runtests.jl", String) * "\nend"); println("route report test parse ok")'
  passed

git diff --check
  passed
```

The slow route-driver report integration test was not run; only a parse check
was run after deleting stale assertion blocks from it. Remaining materialization
coverage exists in the focused atom-growth/materialization tests, and this pass
did not alter the accepted H2 physics endpoint.

Next direction:

Materialize the existing route-neutral GTO supplement representation for this
same H2 request, still without provider blocks or matrices. The next blocker
should advance from `:missing_gto_supplement_representation` to
`:missing_provider_gto_supplement_blocks`.

-- repo-manager@macmini
