# Pass 231 response - independent H2 PQS target/readiness surface

## Files changed

- `src/pqs_source_box_route_driver_skeletons.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl`

No commits or pushes.

## New route/input

Added driver input:

```text
test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl
```

New route kind:

```julia
:bond_aligned_diatomic_independent_pqs_source_box_core_shell
```

This route is target/readiness-only. It records the common H2 physical support
vocabulary and blocks before source-plan, final-basis, H1, H1-J, RHF, supplement,
CR2, export, or public API work.

## Artifact/readiness fields

The independent target writes the intended guard fields:

```text
fake_pqs/enabled = false
route/source_backed_fixed_source_oracle_used = false
route/retained_transform_authority = :pqs_source_box_construction
physics/source_backed_fixed_source_oracle_used = false
physics/retained_transform_authority = :pqs_source_box_construction
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_independent_pqs_atom_contact_core_retained_rule
```

The target group records support metadata only:

```text
target/support_counts = (275, 578, 362)
target/retained_counts = ()
target/expected_final_dimension = nothing
```

It does not claim the fake-PQS/WL retained counts `(251, 98, 114)` as
independent PQS retained counts.

## Blockers

Primary blocker:

```julia
:missing_independent_pqs_atom_contact_core_retained_rule
```

Secondary blocker:

```julia
:missing_independent_pqs_shared_shell_2_retained_rule
```

Broader source-plan blocker:

```julia
:missing_independent_pqs_physical_source_plan_materializer
```

The source-plan authority status is blocked as
`:blocked_pqs_source_box_construction_authority`.

## Validation

Package load:

```text
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
```

Result after final edit:

```text
load ok
elapsed_s=0.676889291
```

Focused readiness artifact check:

```text
julia --project=. -e 'using JLD2; ... include("bin/cartesian_ham_builder.jl") ... assert independent-PQS guard fields ...'
```

Final result:

```text
route_kind=bond_aligned_diatomic_independent_pqs_source_box_core_shell
target_status=blocked_independent_pqs_source_box_target_readiness
source_plan_status=blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
target_retained_counts=()
elapsed_s=128.704891958
```

The final focused driver command exceeded 60s because the edited package
precompiled again and the driver printed the full report. The timed assembly
inside the driver was 22.35s, with final-basis/H1/H1-J/RHF disabled.

Two earlier focused runs failed with `UndefVarError: independent_target not
defined`; those caught a predicate scoped to the wrong helper. I fixed that
before the final passing run.

Other validation:

```text
git diff --check
```

passed.

## Line budget

`git diff --numstat -- src test bin`:

```text
56    6    src/pqs_source_box_diatomic_complete_core_shell.jl
0     31   src/pqs_source_box_route_driver_helpers.jl
27    106  src/pqs_source_box_route_driver_reporting.jl
71    23   src/pqs_source_box_route_driver_skeletons.jl
```

Tracked diff subtotal: 154 added, 166 deleted. The new untracked driver input is
11 lines, so actual `src + test + bin` total is 165 added, 166 deleted, net -1.

## Git status

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M src/pqs_source_box_route_driver_skeletons.jl
?? test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl
```

## Carrying-cost accounting

deleted: redundant no-payload default field clouds in the target report helpers.

simplified: missing-payload summaries now rely on existing writer-side
`get(..., default)` handling instead of duplicating many defaults.

quarantined: the fake-PQS route remains quarantined separately; this pass adds a
distinct fake-free target/readiness route.

not deleted because: the fake-PQS source-backed WL/QW reproduction remains a
useful golden regression and was not mutated.

exact remaining caller/blocker: independent H2 PQS cannot materialize a source
plan until `atom_contact_core` has an independent retained rule; `shared_shell_2`
also lacks an independent retained rule for the old fake-route count `114`.

-- repo-doer@macmini
