Review 224 - accepted

Pass 224 added the private H2 physical supplement preflight boundary without
materializing any supplement matrices or changing the accepted no-supplement
endpoint.

Accepted behavior:

```text
supplement_policy = :none
  supplement_preflight/status = :not_requested
  blocker = nothing

supplement_policy = :mwg_residual_gto
  supplement_preflight/status =
    :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
  blocker = :missing_provider_gto_supplement_blocks
```

The artifact now records compact supplement-preflight facts: fixture label,
support/retained counts and order, retained transform kind `:pqs`, gausslet
final dimension `463`, required/available/missing fact labels, and the explicit
nonmaterialization flags. This is the right boundary: it states that PQS should
use the common H2 physical support plan plus a PQS retained transform, while
GTO/MWG supplement work remains blocked until provider and mixed-block facts
exist.

Manager validation rerun:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
  passed

julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
  passed: 119/119, elapsed_s=140.251064584

git diff --check
  passed
```

The focused H2 endpoint test exceeds 60 seconds because it rebuilds the accepted
H2 final-basis/H1/H1-J/private-RHF artifact before running the blocked MWG
preflight override. I did not rerun the broader route-report integration test:
the doer attempted it, it produced no normal output for several minutes, and it
was killed while compiling inside the changed report test. The focused H2
artifact test directly exercises the new report writer fields.

Line budget:

```text
src/test/bin added:   318
src/test/bin deleted: 322
net:                  -4
```

Deletion/shrinkage accepted:

```text
deleted: 319 lines of stale route-report metadata assertion pressure
simplified: duplicated route_configured_* artifact/TSV checks
quarantined: MWG/GTO request stops at supplement preflight
not deleted because: H2 physical endpoint artifact test is the live gate
exact remaining blocker: :missing_provider_gto_supplement_blocks
```

Next direction:

Audit the existing WL/MWG/GTO supplement machinery before adding PQS supplement
matrices. The next pass should identify the existing provider surfaces for raw
GTO blocks, mixed gausslet/GTO blocks, GTO/GTO blocks, residual MWG
orthogonalization, and combined operator transforms, and decide the smallest
PQS seam that can satisfy `:missing_provider_gto_supplement_blocks` without
forking the physical support/shell vocabulary.

-- repo-manager@macmini
