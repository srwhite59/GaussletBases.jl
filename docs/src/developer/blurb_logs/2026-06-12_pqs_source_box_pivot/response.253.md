# Pass 253 response - independent H2 PQS private RHF input taxonomy

Implemented the taxonomy-only private RHF input pass. I did not change RHF
numerical behavior and did not run the slow H2 RHF route.

## New input

Added:

```text
test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl
```

It includes the H1-J diagnostic input and sets:

```text
artifact_role = :independent_h2_pqs_private_rhf_diagnostic
physics_endpoint_ready = false
physics_endpoint_blocker = :private_rhf_diagnostic_not_public_solver_contract
run_final_basis = true
run_h1 = true
run_h1_j = true
run_private_rhf = true
private_rhf_electron_count = 2
```

## Manifest and Classifier

Updated `docs/src/developer/cartesian_driver_endpoint_manifest.md` with a row
for the independent H2 PQS private RHF diagnostic stage. The row marks it as
not endpoint-ready and not public solver/export-ready.

Updated the compact independent-H2 artifact-role classifier so
`:independent_h2_pqs_private_rhf_diagnostic` receives the same independent
guard behavior:

```text
source_backed_fixed_source_oracle_used = false
retained_transform_authority = :pqs_source_box_construction
```

Fake-PQS remains separate; the fake source-backed WL/QW reproduction role was
not changed.

## Files changed

- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`

## Line budget

Scoped `src + test + bin`, counting the new untracked input file:

```text
11 insertions, 17 deletions
net -6
```

Tracked diff plus manifest:

```text
docs/src/developer/cartesian_driver_endpoint_manifest.md: 1 insertion
src/pqs_source_box_route_driver_helpers.jl: 1 insertion
test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl: 17 deletions
```

Deletion/shrinkage report:

```text
deleted:
- stale assembly-stage terminal RouteCore summary/status mirror assertions

simplified:
- private RHF can now be requested through a tiny explicit driver input instead
  of command-line overrides

quarantined:
- private RHF remains diagnostic-only; no supplements, CR2/export, HamV6,
  public API, public solver readiness, or provider-block work was added

not deleted because:
- compact terminal assembly route/scaffold/deferral smoke remains active
- fake-PQS taxonomy remains separate and still has its manifest row

exact remaining caller/blocker:
- independent H2 PQS private RHF now has an explicit input; endpoint readiness
  remains blocked by `:private_rhf_diagnostic_not_public_solver_contract`
```

## Validation

Commands run:

```text
julia --project=. -e 'include("test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl"); @assert route_family === :pqs_source_box; @assert route_kind === :bond_aligned_diatomic_independent_pqs_source_box_core_shell; @assert artifact_role === :independent_h2_pqs_private_rhf_diagnostic; @assert (run_final_basis, run_h1, run_h1_j, run_private_rhf) == (true, true, true, true); @assert private_rhf_electron_count == 2; @assert physics_endpoint_ready == false; @assert physics_endpoint_blocker === :private_rhf_diagnostic_not_public_solver_contract; println("private_rhf_input_smoke ok")'
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
rg -n "[ \t]$" test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl docs/src/developer/cartesian_driver_endpoint_manifest.md src/pqs_source_box_route_driver_helpers.jl test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
git status --short --branch
```

Results:

```text
include/flag smoke: passed
package load: passed
git diff --check: passed
direct trailing-whitespace search: no matches
```

Final git status:

```text
## main...origin/main
 M docs/src/developer/cartesian_driver_endpoint_manifest.md
 M src/pqs_source_box_route_driver_helpers.jl
 M test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
?? test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl
```

-- repo-doer@macmini
