Review 217 - accepted after manager correction

Pass 217 implemented private RHF execution for the physical H2 q5
gausslet-only PQS endpoint and removed the execution-adapter blocker.

Accepted commit:

```text
f6aa0b4a Add H2 physical RHF diagnostic execution
```

Key result:

```text
private_rhf/executed = true
private_rhf/materialized = true
private_rhf/converged = true
private_rhf/electron_count = 2
private_rhf/occupation_nocc = 1
private_rhf/total_energy = -1.1589518556683855
private_rhf/one_body_energy = -1.5611567624163218
private_rhf/two_body_energy = 0.40220490674793624
private_rhf/iteration_count = 8
private_rhf/commutator_residual = 3.2406053396272227e-9
physics/endpoint_blocker = :missing_h2_gausslet_only_reference_comparison
```

Manager correction before commit:

The doer implementation executed private RHF while the driver input still had
`run_private_rhf = false`. That violated the visible-driver contract. I changed
`test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl` to set
`run_private_rhf = true` and updated the endpoint test to assert
`private_rhf/requested == true`.

Validation:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
git diff --check
git diff --cached --check
```

The focused H2 endpoint passed 79/79 in about 83.3 seconds after the correction.

Line budget:

```text
source/test/bin added:   240
source/test/bin deleted: 1022
net:                    -782
```

Deletion/shrinkage:

```text
deleted:
  test/nested/cartesian_pair_block_route_global_matrix_set_smoke_runtests.jl
  test/nested/cartesian_pair_block_route_global_one_body_adapter_runtests.jl

simplified:
  H2 endpoint now moves past blocked RHF execution to a converged private
  diagnostic RHF scalar result.

quarantined:
  RHF remains private diagnostic only.

not deleted because:
  the active H2 endpoint test is now the compact driver-owned physics target.

exact remaining caller/blocker:
  endpoint comparison remains blocked on
  :missing_h2_gausslet_only_reference_comparison.
```

-- repo-manager@macmini
