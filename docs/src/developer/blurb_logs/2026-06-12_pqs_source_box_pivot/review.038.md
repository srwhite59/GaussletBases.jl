Review 038: accepted.

The ignored H1 probe gives the right promotion candidate. It assembled the
explicit PQS final-basis H1 seam from direct retained overlap, kinetic, and
by-center nuclear blocks; transformed those into the shell-realized final
basis; assembled the final one-electron Hamiltonian; and matched the
shell-support oracle to roundoff.

Key probe facts:

```text
source dims/count: 5 x 5 x 5 / 125
boundary/final retained count: 98 / 98
final overlap identity error: 2.24e-14
direct retained flags: overlap/kinetic/nuclear all true
raw source one-body blocks materialized on active path: false
Hamiltonian max error vs shell-support oracle: 1.33e-15
H1 lowest eigenvalue delta vs oracle: 3.75e-16
ordinary symmetric eigensolve used: true
generalized overlap solve used: false
```

This is enough to promote one compact H1 gate, provided it replaces/shrinks old
oracle pressure. The identified shrink target is the slow integration section
around `_pqs_current_route_safe_term_matrices(...)` and
`_pqs_current_route_safe_term_authority_comparison(...)` in
`test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`.

Next pass should add the compact direct-retained final-H1 gate and remove or
shrink that old helper-vocabulary integration coverage. Keep the private CCPM
helpers as oracle/debug code for now; this is a test/coverage replacement pass,
not a source deletion pass.

-- repo-manager@macmini
