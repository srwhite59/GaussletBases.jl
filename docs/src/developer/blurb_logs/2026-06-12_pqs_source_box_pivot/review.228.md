Review 228 - accepted

Pass 228 did the right cleanup before the next materialization step.

Accepted changes:

```text
removed:
  supplement_request_representation_status
  supplement_request_representation_object_kind
  supplement_request/representation_status artifact field

kept:
  supplement_request = request/preflight metadata
  supplement_representation = sole representation authority
```

Current H2 MWG supplement status:

```text
request status:        :available_pqs_physical_gausslet_supplement_request
representation status: :available_pqs_physical_gausslet_gto_supplement_representation
orbital count:         18
preflight status:      :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
preflight blocker:     :missing_provider_gto_supplement_blocks
```

No provider blocks, matrices, residual MWG representation, supplemented scalar
values, CR2/HFDMRG/export/public behavior, or final-basis self-overlap working
data were added.

The provider-block audit is useful and accepted. The next live missing facts
are:

```text
CPB/source-box coverage for H2 physical support units
local row/source ordering for each support unit
support-row to retained-463 transform placement
mixed gausslet-GTO placement/accumulation rule
GTO/GTO self-block rule
raw moment matrices for MWG residualization
```

The scoped source/test/bin diff was a pure shrink:

```text
added:    0
deleted: 32
net:    -32
```

Manager validation:

```text
git diff --check
  passed

julia --project=. -e 'using GaussletBases; println("load ok")'
  passed

julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
  passed: 115/115
  elapsed_s = 143.836901333
```

Next direction:

Before building provider blocks, inspect the actual H2 physical source-plan and
final-basis payload fields and identify whether CPB/source-box coverage,
support-row ordering, and support-row-to-retained transform placement already
exist under stable names. If request-side provider missing-fact fields are still
duplicating the preflight role, shrink them while doing the audit.

-- repo-manager@macmini
