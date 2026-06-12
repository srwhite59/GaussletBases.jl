Review result:

Accepted as the focused inventory-shape cleanup requested by `blurb.019`.
The pass removed the most concrete side-dependent tuple specialization from
the hot decomposed WL inventory result without adding compatibility shims or
new tests.

What changed:

- `unit_keys` is now `Vector{Symbol}` instead of `NTuple{N,Symbol}`.
- `unit_summaries` is now a vector of compact unit summaries instead of an
  `N`-element tuple.
- `pair_summaries` is now a compact count/status `NamedTuple`:
  `status`, `pair_count`, `detailed_pair_summaries_materialized`,
  `detailed_pair_summaries_source`, and `retained_ranges_available`.
- Detailed pair range checks now use live `unit_pairs` iteration.

Before/after shape result:

```text
field                 before side-7       before side-15      after both
unit_keys             NTuple{27,Symbol}   NTuple{131,Symbol}  Vector{Symbol}
unit_summaries        27-element tuple    131-element tuple   Vector{NamedTuple{...}}
pair_summaries        NTuple{378,...}     Symbol              compact NamedTuple
```

Physics/timing validation reviewed:

- Be S+P route still matches the old nested/QW oracle:
  `-14.574514244574639` versus `-14.574514244574694`, delta
  `5.5067062021407764e-14 Ha`.
- final dimension remains `636`.
- retained gausslet dimension remains `615`.
- units / pairs remain `131 / 8646`.
- fallback flags remain false for full-parent CPB, direct Cartesian, ordinary
  Cartesian IDA, raw GTO final density-density, and generalized final solve.

Test/coverage review:

The focused inventory test was updated rather than expanded. It no longer
treats `pair_summaries` as detailed range authority. It checks vector-backed
`unit_keys` / `unit_summaries`, compact `pair_summaries` status/count, and
uses `inventory.unit_pairs` for detailed pair range facts. No new test was
added.

Remaining specialization candidates:

- `pair_keys` still has a small/large shape split: vector-backed for small
  inventories and an omitted sentinel for large inventories.
- parent axis bundle objects still ride through metadata instead of a
  `ParentAxisContext3D`.
- atom+GTO route results still mix compute payloads with audit/report fields.
- factorized sidecar and residual moment bundles remain `NamedTuple`-shaped
  compute concepts.

Validation run by doer:

- `julia --project=. test/nested/cartesian_white_lindsey_decomposed_unit_pair_inventory_runtests.jl`
- `julia --project=. tmp/work/atom_gto_specialization_shape_audit.jl`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Additional validation run by manager:

- `julia --project=. test/nested/cartesian_white_lindsey_decomposed_unit_pair_inventory_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage review:

The pass deleted hot-route materialization of detailed
`pair_summaries::NTuple{N,...}` and removed tuple construction for hot
`unit_keys` and `unit_summaries`. It shrank the focused inventory test away
from pair-summary tuple vocabulary. No new tests or compatibility adapters
were added.

Next target:

Rerun warm/cold Be S+P timing and the specialization-shape audit after this
inventory cleanup, then choose the next compile-pressure cleanup from evidence.
Do not assume this change materially reduced cold route time until measured.
