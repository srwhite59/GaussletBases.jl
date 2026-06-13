Manager review for pass 171

Accepted.

The pass adds the missing consumer-facing PQS dense final-basis
density-density interaction matrix while preserving the existing pre-final
provenance fields. The formula is the expected congruence:

```julia
final_interaction_matrix =
    transpose(final_to_pre_final_coefficients) *
    pre_final_pair_matrix *
    final_to_pre_final_coefficients
```

The focused test verifies the orientation with deterministic final-density
vectors by checking that the final-basis quadratic form matches the pre-final
quadratic form after `d_pre = coefficients * d_final`.

Line budget is satisfied across `src`, `test`, and the tracked generator
surface: 47 added, 581 deleted, net -534. The deleted
`white_lindsey_materialized_seed_runtests.jl` was old one-center/materialized
seed scaffolding, and the live route-configured diatomic WL artifact path is
now validated elsewhere.

Solver/export readiness correctly remains false. This makes the private
artifact HF-input-readable for CR2 experimentation, but it does not promote a
GaussletBases solver contract.

-- repo-manager@macmini
