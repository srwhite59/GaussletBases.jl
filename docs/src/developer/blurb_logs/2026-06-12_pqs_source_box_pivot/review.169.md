Manager review for pass 169

Accepted.

The audit made the right distinction: the old White-Lindsey materialized seed
and one-center route-configured bundle are not honest Be2 comparators for the
CR2 artifact. The next usable source is the route-configured diatomic
atom-growth White-Lindsey path with:

- `low_order_shellization_policy = :atom_growth_complete_rectangular`
- `route_configured_diatomic_ham_interaction_treatment = :ggt_nearest`
- `save_ham_artifact = true`
- explicit `white_lindsey_expansion`

The populated WL schema must be labeled as a final-basis
ordinary-Cartesian/Qiu-White Hamiltonian bundle. Its two-body object is a
final-basis density-density interaction matrix, not PQS pre-final/source-box
pair data. Supplement/residual-GTO, correction/EGOI, Qiu-White atom-local HF,
MWG/IDA, solver/export, and HFDMRG remain unavailable placeholders.

Next pass should implement only the generator-side WL population. It should not
edit tracked `src/` or `test/`, and it should not try to delete the old
White-Lindsey seed path in the same pass. Keep that deletion candidate in view
for a later cleanup once the route-configured diatomic WL population is stable.

-- repo-manager@macmini
