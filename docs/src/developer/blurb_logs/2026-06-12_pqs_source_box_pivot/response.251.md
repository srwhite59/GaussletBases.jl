# Pass 251 response - clarify independent H2 PQS driver input taxonomy

Implemented the taxonomy/input cleanup only. No source or route behavior
changed, and I did not run slow H2 driver validations.

## New input files

Added three tiny include/override driver input variants:

```text
test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl
test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl
test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl
```

The existing base input remains the no-physics readiness input:

```text
test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl
```

Variant flags:

```text
final_basis: run_final_basis=true, run_h1=false, run_h1_j=false, run_private_rhf=false
H1:          run_final_basis=true, run_h1=true,  run_h1_j=false, run_private_rhf=false
H1-J:        run_final_basis=true, run_h1=true,  run_h1_j=true,  run_private_rhf=false
```

## Manifest changes

Updated `docs/src/developer/cartesian_driver_endpoint_manifest.md` to
distinguish:

```text
independent H2 PQS target/readiness input, endpoint not ready
independent H2 PQS final-basis diagnostic, dimension 471, not solver-ready
independent H2 PQS H1 diagnostic, dimension 471, not solver-ready
independent H2 PQS H1-J density diagnostic, dimension 471, not solver-ready
fake-PQS source-backed WL/QW 463 reproduction, still fake and separate
```

Removed stale manifest references to the deleted
`h2_pqs_q5_source_box_diagnostic_r4` input/test row and replaced the generic
planned independent-H2 row with the explicit stage taxonomy.

## Guardrails

No source files changed. No route/source-plan/final-basis/H1/H1-J behavior
changed. Fake-PQS remains a separate source-backed WL/QW reproduction row.
RHF/private RHF, supplements, CR2, exports, and public API readiness were not
touched.

## Line budget

Scoped `src + test + bin`, counting the new untracked driver input files:

```text
24 insertions, 27 deletions
net -3
```

Manifest/input taxonomy cost before deletion offset:

```text
new input files: 24 insertions
manifest: 4 insertions, 2 deletions
taxonomy subtotal: 28 insertions, 2 deletions, net +26
```

Including the stale staged-test deletion offset:

```text
overall touched taxonomy/test cleanup: 28 insertions, 29 deletions, net -1
```

Deletion/shrinkage report:

```text
deleted:
- stale manifest row for deleted `h2_pqs_q5_source_box_diagnostic_r4` input/test
- stale report-stage terminal route-core mirror assertion block

simplified:
- independent H2 PQS manifest taxonomy now names readiness, final-basis, H1,
  and H1-J stage inputs explicitly instead of one generic planned row

quarantined:
- fake-PQS remains a separate fake/source-backed WL/QW reproduction row

not deleted because:
- the base independent input remains useful as the no-physics readiness input
- compact report-stage terminal route/scaffold/materialization smoke remains

exact remaining caller/blocker:
- no caller needed the stale deleted H2 diagnostic manifest row; the actual
  deleted input/test files were already absent
```

## Validation

Commands run:

```text
julia --project=. -e 'variants = [("test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl", (true, false, false, false), :independent_h2_pqs_final_basis_diagnostic), ("test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl", (true, true, false, false), :independent_h2_pqs_h1_diagnostic), ("test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl", (true, true, true, false), :independent_h2_pqs_h1_j_density_diagnostic)]; for (path, flags, role) in variants; include(path); @assert (run_final_basis, run_h1, run_h1_j, run_private_rhf) == flags; @assert artifact_role === role; @assert route_family === :pqs_source_box; @assert route_kind === :bond_aligned_diatomic_independent_pqs_source_box_core_shell; println(path, " flags=", flags, " artifact_role=", role); end'
git diff --check
rg -n "[ \t]$" test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl docs/src/developer/cartesian_driver_endpoint_manifest.md test/nested/cartesian_report_stage_low_order_policy_runtests.jl
git status --short --branch
```

Results:

```text
include smoke: passed for all three new variants
git diff --check: passed
direct trailing-whitespace search: no matches
```

I did not run package load because no Julia source changed and the include
smoke did not require loading `GaussletBases`.

Final git status:

```text
## main...origin/main
 M docs/src/developer/cartesian_driver_endpoint_manifest.md
 M test/nested/cartesian_report_stage_low_order_policy_runtests.jl
?? test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl
?? test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl
?? test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl
```

-- repo-doer@macmini
