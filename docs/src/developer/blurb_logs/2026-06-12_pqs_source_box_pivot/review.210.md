Pass 210 review: accepted.

The pass made the intended intermediate state explicit:

```text
candidate_status = :available_physical_gausslet_source_plan_candidate
candidate_source = :source_backed_fixed_source_oracle
candidate_counts_match = true
source_plan_authority_status = :candidate_not_route_authority
source_plan_blocker = :source_plan_candidate_not_route_authority
```

That is the right boundary. The H2 463 rows/coefs are now located and
count-checked, but the route still does not pretend the source-backed oracle is
already route-owned PQS authority.

The candidate checks record the physical target:

- support order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- retained order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- support counts `(275, 578, 362)`
- retained counts `(251, 98, 114)`
- final dimension `463`
- no supplement
- no H2 221 diagnostic source-plan reuse

Artifact behavior stays honest:

- `physics/endpoint_ready = false`
- `physics/endpoint_blocker = :source_plan_candidate_not_route_authority`
- no final basis, H1, H1-J, density interaction, RHF, HFDMRG, CR2, export, or
  public API behavior was added.

I checked the large deletion in the old high-order opt-in source-construction
test. The removed block was duplicated projected-q-shell/metadata assertion
pressure; the core projected-q-shell rule is still covered in
`test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl` and
`test/nested/bond_aligned_diatomic_high_order_recipe_realization_audit_runtests.jl`.
That makes the deletion acceptable.

Validation:

- Doer ran package load; passed.
- Doer ran the H2 463 artifact-readiness test; 40/40 passed, about 74.5 s.
- Manager reran package load; passed.
- `git diff --check` passed.

Source/test/bin line budget is net negative:

```text
added:   182
deleted: 216
net:     -34
```

Next step: decide and implement the narrow authority wrapper that turns the
verified source-backed candidate into a private route-owned physical source-plan
object, or keep it blocked with a precise authority blocker if that wrapper
cannot be made honest.

-- repo-manager@macmini
