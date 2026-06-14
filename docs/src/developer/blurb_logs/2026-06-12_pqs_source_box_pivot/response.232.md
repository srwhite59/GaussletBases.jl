# Pass 232 response - independent H2 PQS support plan authority

Implemented the pass as an honest blocked support-plan authority layer for the
independent H2 PQS target.

## Files changed

- `src/pqs_source_box_route_driver_skeletons.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `src/pqs_source_box_route_driver_helpers.jl`

## Support-plan status

- `target/support_plan_status = :blocked_independent_pqs_support_region_plan`
- `target/support_plan_blocker = :missing_independent_pqs_support_region_materializer`
- `target/support_plan_authority = :pqs_source_box_route_geometry_pending_materializer`
- `target/support_counts_generated = false`
- `target/support_counts_source = :target_constants_pending_support_region_materializer`

The target support constants remain:

- `:atom_contact_core => 275`
- `:shared_shell_1 => 578`
- `:shared_shell_2 => 362`

Those counts are still target constants, not generated support regions.

## Authority/provenance

The independent route now carries a compact `support_plan` summary through the
target inventory/report path. It records the intended support-plan authority as
route geometry pending a real support-region materializer. It does not claim a
source plan, retained transform, or physics endpoint.

Guard fields remain intact:

- `fake_pqs/enabled = false`
- `route/source_backed_fixed_source_oracle_used = false`
- `physics/endpoint_ready = false`
- `target/retained_counts = ()`

## Forbidden surfaces avoided

No atom-contact retained rule, `shared_shell_2` retained rule, fake/WL
coefficient import, final basis, H1, H1-J, RHF, supplement, CR2, export,
public API, or WL/QW physics comparison was added. The fake-PQS route was not
changed.

## Validation

Focused independent-input artifact/readiness check:

```sh
julia --project=. -e 'using JLD2; dir = mktempdir(); outfile = joinpath(dir, "h2_independent_pqs_support_plan.jld2"); saved_args = copy(ARGS); empty!(ARGS); append!(ARGS, ["test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl", "outfile=$(repr(outfile))", "save_tsv=false"]); t = @elapsed try include("bin/cartesian_ham_builder.jl") finally empty!(ARGS); append!(ARGS, saved_args) end; jldopen(outfile, "r") do f; @assert f["fake_pqs/enabled"] == false; @assert f["route/source_backed_fixed_source_oracle_used"] == false; @assert f["physics/endpoint_ready"] == false; @assert Tuple(f["target/support_counts"]) == (275, 578, 362); @assert Tuple(f["target/retained_counts"]) == (); @assert f["target/support_plan_status"] === :blocked_independent_pqs_support_region_plan; @assert f["target/support_plan_blocker"] === :missing_independent_pqs_support_region_materializer; @assert f["target/support_plan_authority"] === :pqs_source_box_route_geometry_pending_materializer; @assert f["target/support_counts_generated"] == false; println("support_plan_status=", f["target/support_plan_status"]); println("support_plan_blocker=", f["target/support_plan_blocker"]); println("support_counts_source=", f["target/support_counts_source"]); end; println("elapsed_s=", t)'
```

Result:

```text
support_plan_status=blocked_independent_pqs_support_region_plan
support_plan_blocker=missing_independent_pqs_support_region_materializer
support_counts_source=target_constants_pending_support_region_materializer
elapsed_s=67.634487375
```

The command exceeded 60 seconds because the edited source required driver
assembly/artifact writing on the independent input. The driver report also
showed `source_backed_fixed_source_oracle_used=false`,
`final_basis_materialized=false`, `h1_materialized=false`,
`h1_j_materialized=false`, and `rhf_materialized=false`.

Package load:

```sh
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
```

Result:

```text
load ok
elapsed_s=0.675291625
```

Whitespace:

```sh
git diff --check
```

Result: passed.

## Line budget

```text
git diff --numstat -- src test bin
3       10      src/pqs_source_box_diatomic_complete_core_shell.jl
0       8       src/pqs_source_box_route_driver_helpers.jl
10      6       src/pqs_source_box_route_driver_reporting.jl
10      0       src/pqs_source_box_route_driver_skeletons.jl
```

Scoped total: 23 added / 24 deleted, net -1.

## Git status

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M src/pqs_source_box_route_driver_skeletons.jl
```

## Deletion/shrinkage report

deleted:
- stale source-plan candidate artifact/report aliases from the independent
  target readiness surface;
- unused source-plan summary aliases that no longer belong on the visible
  blocked target report.

simplified:
- independent target reporting now exposes a support-plan fingerprint instead
  of another source-plan-candidate field cloud.

quarantined:
- independent H2 PQS remains a blocked target/readiness route only;
- fake-PQS remains quarantined and unchanged.

not deleted because:
- the fake-PQS endpoint is still the accepted golden driver regression;
- the independent target input remains the live placeholder for the next real
  PQS H2 route.

exact remaining caller/blocker:
- `:missing_independent_pqs_support_region_materializer` blocks generated
  support regions;
- downstream retained-rule blockers remain
  `:missing_independent_pqs_atom_contact_core_retained_rule` and
  `:missing_independent_pqs_shared_shell_2_retained_rule`.

-- repo-doer@macmini
