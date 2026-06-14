# Pass 234 response - independent H2 PQS support-region plan

Implemented the independent H2 PQS support-region materializer/fingerprint
under the manager exception addendum to the line-budget rule.

## Files changed

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `src/pqs_source_box_route_driver_skeletons.jl`

No commits or pushes.

## Support-plan result

The independent H2 target now generates its support-region fingerprint from
route geometry via `CartesianShellification.raw_terminal_geometry(...)`.

Artifact-facing result:

```text
target/support_plan_status = :available_independent_pqs_support_region_plan
target/support_plan_authority = :cartesian_shellification_route_geometry
target/support_counts_generated = true
target/support_counts_source = :cartesian_shellification_route_geometry
target/support_counts = (275, 578, 362)
target/support_plan_coverage_complete = true
target/support_plan_duplicate_count = 0
target/support_plan_missing_count = 0
target/support_plan_outside_count = 0
```

Generated grouping:

- `:atom_contact_core`: two `:atom_local_core` regions plus one
  `:midpoint_slab`, `125 + 125 + 25 = 275`;
- `:shared_shell_1`: outer `:shared_molecular_shell`, `578`;
- `:shared_shell_2`: inner `:shared_molecular_shell`, `362`.

The shared-shell order is explicitly outside-in by outer-box volume, not raw
terminal emission order.

## Source APIs used

- `CartesianShellification.raw_terminal_geometry(...)`
- `_nested_axis_pgdg(...)` to read parent axis-bundle centers from the existing
  parent object
- existing `parent.atom_locations`, `parent.bond_axis`, and route skeleton
  parent-axis counts

No fake-PQS/WL coefficient matrix, fixed-source retained transform, or
`bond_aligned_diatomic_nested_fixed_source(...)` support authority was used.

## Guardrails

Still blocked/unchanged:

- `fake_pqs/enabled = false`
- `route/source_backed_fixed_source_oracle_used = false`
- `physics/endpoint_ready = false`
- `target/retained_counts = ()`
- `target/expected_final_dimension = nothing`

No atom-contact retained rule, `shared_shell_2` retained rule, source-plan
authority beyond support-region generation, final basis, H1, H1-J, RHF,
supplement, CR2, export, public API, or WL/QW physics comparison was added.

## Validation

Package load:

```sh
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
```

Result:

```text
load ok
elapsed_s=58.297125083
```

Focused independent-input artifact/readiness check:

```sh
julia --project=. -e 'using JLD2; dir = mktempdir(); outfile = joinpath(dir, "h2_independent_pqs_support_plan.jld2"); saved_args = copy(ARGS); empty!(ARGS); append!(ARGS, ["test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl", "outfile=$(repr(outfile))", "save_tsv=false"]); t = @elapsed try include("bin/cartesian_ham_builder.jl") finally empty!(ARGS); append!(ARGS, saved_args) end; jldopen(outfile, "r") do f; @assert f["fake_pqs/enabled"] == false; @assert f["route/source_backed_fixed_source_oracle_used"] == false; @assert f["physics/endpoint_ready"] == false; @assert Tuple(f["target/support_counts"]) == (275, 578, 362); @assert Tuple(f["target/retained_counts"]) == (); @assert f["target/support_plan_status"] === :available_independent_pqs_support_region_plan; @assert f["target/support_plan_authority"] === :cartesian_shellification_route_geometry; @assert f["target/support_counts_generated"] == true; @assert f["target/support_counts_source"] === :cartesian_shellification_route_geometry; @assert f["target/support_plan_coverage_complete"] == true; @assert f["target/support_plan_duplicate_count"] == 0; @assert f["target/support_plan_missing_count"] == 0; @assert f["target/support_plan_outside_count"] == 0; println("support_plan_status=", f["target/support_plan_status"]); println("support_plan_authority=", f["target/support_plan_authority"]); println("support_counts_generated=", f["target/support_counts_generated"]); println("support_counts_source=", f["target/support_counts_source"]); println("support_plan_coverage_complete=", f["target/support_plan_coverage_complete"]); end; println("elapsed_s=", t)'
```

Result:

```text
support_plan_status=available_independent_pqs_support_region_plan
support_plan_authority=cartesian_shellification_route_geometry
support_counts_generated=true
support_counts_source=cartesian_shellification_route_geometry
support_plan_coverage_complete=true
elapsed_s=66.569249209
```

The focused driver command exceeded 60 seconds because it included driver
assembly and JLD2 artifact writing after source edits. The driver report kept
`source_backed_fixed_source_oracle_used=false`,
`final_basis_materialized=false`, `h1_materialized=false`,
`h1_j_materialized=false`, and `rhf_materialized=false`.

Whitespace:

```sh
git diff --check
```

Result: passed.

## Line budget

This pass uses the manager exception recorded in `blurb.234.md`.

```text
git diff --numstat -- src test bin
110     0       src/pqs_source_box_diatomic_complete_core_shell.jl
8       0       src/pqs_source_box_route_driver_reporting.jl
0       10      src/pqs_source_box_route_driver_skeletons.jl
```

Scoped total: 118 added / 10 deleted, net +108.

## Git status

Before writing this response:

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M src/pqs_source_box_route_driver_skeletons.jl
```

## Carrying-cost accounting

deleted:
- pass-232 hard-coded blocked `support_plan` field cloud from the route
  skeleton.

simplified:
- independent target support counts now come from the generated support-plan
  fingerprint when available, instead of remaining target constants;
- report output now exposes compact coverage/disjointness counters from the
  support-plan summary.

quarantined:
- independent H2 PQS remains a support-plan/readiness route only;
- fake-PQS remains quarantined and unchanged.

not deleted because:
- fake-PQS source-backed WL/QW reproduction remains the accepted golden
  regression;
- complete-core-shell support-window helpers are still called by the diagnostic
  path and were not safe deletion fuel for this pass.

exact remaining caller/blocker:
- source plan and retained transform remain blocked by
  `:missing_independent_pqs_atom_contact_core_retained_rule`,
  `:missing_independent_pqs_shared_shell_2_retained_rule`, and
  `:missing_independent_pqs_physical_source_plan_materializer`.

-- repo-doer@macmini
