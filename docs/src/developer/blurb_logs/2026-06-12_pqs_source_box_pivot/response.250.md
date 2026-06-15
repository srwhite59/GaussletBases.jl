# Pass 250 response - independent H2 PQS H1-J density diagnostic

The existing physical-gausslet H1-J diagnostic path works for the independent
H2 PQS final basis. I did not add source code. This pass is a focused
validation plus staged-test cleanup pass.

## Result

Focused driver/artifact check with `run_h1=true`, `run_h1_j=true`, and
`run_private_rhf=false` passed:

```text
h1_status = materialized_pqs_physical_gausslet_h1_solve
h1_materialized = true
h1_j_status = materialized_pqs_physical_gausslet_h1_j_payload
h1_j_materialized = true
final_dimension = 471
h1_lowest_energy = -0.7946037173365885
h1_hamiltonian_matrix_finite = true
h1_hamiltonian_symmetry_error = 1.0658141036401503e-14
physics_endpoint_ready = false
physics_endpoint_blocker = missing_physical_gausslet_rhf_or_solver_contract
```

Density/H1-J diagnostics:

```text
density_interaction_status =
  materialized_pqs_physical_gausslet_pre_final_density_interaction
density_gauge = pre_final_localized_positive_weight
raw_pair_factor_convention = raw_numerator
support_weight_count = 1215
support_weights_all_positive = true
support_raw_pair_shape = (1215, 1215)
support_raw_pair_finite = true
pre_final_pair_matrix_shape = (471, 471)
pre_final_pair_matrix_finite = true
pre_final_pair_matrix_symmetry_error = 2.55351295663786e-15
h1_j_self_coulomb = 0.4569117646737236
```

This is a diagnostic density/H1-J result, not solver readiness.

## Guard fields

The focused artifact check asserted:

```text
source_backed_fixed_source_oracle_used = false
fake_pqs_enabled = false
retained_transform_authority = pqs_source_box_construction
private_rhf_materialized = false
exports_materialized = false
public_api = false
```

I did not enable RHF/private RHF, supplements, CR2, exports, public API, or
fake-PQS comparisons.

## Files changed

- `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`

No source files changed.

## Line budget and deletion offset

Scoped `src + test + bin` line budget:

```text
0 insertions, 12 deletions
net -12
```

Deletion/shrinkage report:

```text
deleted:
- stale transform-stage selected-terminal-sidecar materialization mirror
  assertions

simplified:
- no source simplification was needed; the H1-J path already materialized

quarantined:
- H1-J remains a diagnostic density interaction path, not RHF/solver/export
  readiness

not deleted because:
- compact transform selected-sidecar count smoke remains as the active staged
  contract
- H1-J artifact facts are existing live diagnostic outputs and were needed for
  the focused seam validation

exact remaining caller/blocker:
- independent H2 PQS H1-J is materialized; endpoint readiness remains blocked
  by `:missing_physical_gausslet_rhf_or_solver_contract`
```

## Validation

Commands run:

```text
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
julia --project=. tmp/work/pass250_independent_h2_h1_j_probe.jl
git diff --check
git status --short --branch
```

Results:

```text
package load: passed, elapsed_s=0.650674166
focused independent H2 PQS H1-J driver/artifact check: passed, elapsed_s=78.692066458
git diff --check: passed
```

The focused H1-J driver check exceeded 60 seconds as expected because it
materializes the independent final basis, H1, support weights, raw pair
numerator, and pre-final density interaction diagnostics while keeping RHF and
exports off.

Final git status:

```text
## main...origin/main
 M test/nested/cartesian_transform_stage_low_order_policy_runtests.jl
```

-- repo-doer@macmini
