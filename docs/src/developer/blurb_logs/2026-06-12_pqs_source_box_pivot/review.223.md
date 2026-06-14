Review 223 - accepted

Pass 223 retired the stale component-smoke/CR2 sidecar surface before new
supplement work begins.

Accepted deletion:

```text
removed include from src/CartesianContractedParentMetrics.jl
deleted src/cartesian_contracted_parent_metrics/component_smoke_sidecars.jl
deleted test/nested/pqs_component_route_report_adapter_runtests.jl
removed include from test/nested/integration_runtests.jl
removed the component-smoke helper/call block from
  test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
```

The deletion is conceptually correct. The removed surface mixed source-box
diagnostic route-shadow vocabulary with final-residual MWG/CR2 sidecar concepts.
It had no non-test source caller and no accepted endpoint contract. Keeping it
would have pulled the next supplement seam toward an obsolete report/sidecar
contract.

Validation rerun by manager:

```text
rg -n "component_route_smoke|component_smoke_sidecars|pqs_component_route_report_adapter" src test
  no live hits

julia --project=. -e 'for path in (...); Meta.parseall(read(path, String)); println("parse ok ", path); end'
  parsed:
    test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
    test/nested/integration_runtests.jl
    src/CartesianContractedParentMetrics.jl

julia --project=. -e 'using GaussletBases; println("load ok")'
  passed

git diff --check
git diff --cached --check
  passed
```

Line budget:

```text
0 added
4103 deleted
net -4103
```

Committed and pushed:

```text
23990c8c Retire component smoke sidecar surface
```

Next direction:

Add the first private H2 physical supplement preflight/payload boundary, but keep
it metadata/readiness only. The accepted H2 no-supplement endpoint must remain
unchanged. The new boundary should identify the common physical support plan,
the PQS retained transform kind, the requested `:mwg_residual_gto` supplement
policy, and the exact missing facts before any GTO/MWG matrix materialization.

-- repo-manager@macmini
