Purpose:

Remove type-level inventory-size specialization from the hot decomposed WL
inventory result. Do this before adding larger precompile workloads or broader
route structs.

Why now:

The specialization audit found that the production data path is already
vector/table based:

```text
retained_units storage  Vector{RetainedUnitRecord}
unit_pairs storage      UnitPairIndexTable
```

but the hot inventory result still creates tuple-shaped report fields:

```julia
unit_keys = Tuple(unit.unit_key for unit in units)
unit_summaries = Tuple(...)
pair_summaries = Tuple(...)  # for inventories below the compact limit
```

This makes the side-7 synthetic precompile workload compile
`NTuple{27,Symbol}` and 27-element unit summary tuples, while side-15 Be needs
`NTuple{131,Symbol}` and 131-element summaries. The side-7 workload also
compiles a 378-pair `pair_summaries` tuple that production side-15 omits.

Exact task:

1. Update `_white_lindsey_decomposed_unit_pair_inventory_result(...)` and its
   helpers in:

   ```text
   src/cartesian_pair_block_materialization/white_lindsey_decomposed_unit_pair_inventory.jl
   ```

   so the hot result no longer stores inventory-size-specialized tuple
   summaries.

2. Preferred shape:
   - `unit_keys` should be vector-backed or count/status only, not
     `NTuple{N,Symbol}`.
   - `unit_summaries` should be vector-backed or an explicit compact summary
     object, not `Tuple{Vararg{...,N}}`.
   - `pair_summaries` should not materialize `NTuple{N,...}` in the hot route,
     even for small precompile inventories. Use count/status summary by
     default, or move detailed pair summaries behind an explicit audit helper.
   - Keep `retained_units::Vector{RetainedUnitRecord}` and
     `UnitPairIndexTable` as production data.

3. Preserve useful probe/test access without preserving the old type shape.
   Existing callers that iterate `unit_keys` or `unit_summaries` can usually
   work with vectors. If a test asserts tuple equality only to preserve helper
   vocabulary, shrink it instead of maintaining a tuple just for the test.

4. Do not create a broad struct framework in this pass. If you introduce a
   struct, it should be a small inventory-summary concept only, such as
   `DecomposedWLInventorySummary`, and it must replace an existing field cloud.

5. Rerun the specialization-shape audit after the change. The goal is that
   side-7 and side-15 no longer differ by `unit_keys::NTuple{N,...}` or
   `unit_summaries::Tuple{Vararg{...,N}}`. Report remaining differences.

6. Rerun enough physics/timing validation to show behavior stayed fixed:
   - Be S+P route or warm/cold probe if practical;
   - `using GaussletBases`;
   - `git diff --check`.

Trust boundary:

No public APIs, exports, route defaults, new acceptance tests by default, PQS,
ECP, Be2, high-l Be, H2, Cr, full-parent CPB fallback, direct Cartesian
fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or
generalized final-basis solve.

Test policy:

Do not add tests by default. If existing focused tests break because they
assert tuple-shaped summaries, update or shrink them to the live contract:
counts, statuses, sample first/last units, ranges, and iteration behavior.
Avoid adding broad metadata tests.

Decision rules:

- If the inventory summary shape can be changed locally with focused test
  updates, do it.
- If many unrelated callers depend on tuple-shaped summaries as a public
  contract, stop and report the smallest compatibility strategy.
- If changing `pair_summaries` exposes a genuinely live oracle need, move the
  detailed pair summaries to an explicit audit helper rather than keeping them
  on the hot result.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- exact new shape/types for `unit_keys`, `unit_summaries`, and
  `pair_summaries`;
- before/after specialization-shape audit table;
- whether any tests were shrunk because they only preserved tuple vocabulary;
- Be route/timing or physics validation run;
- `using GaussletBases` and `git diff --check`;
- deletion/shrinkage report.
