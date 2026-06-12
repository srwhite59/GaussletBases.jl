Completed the atom+GTO final-basis cold compile specialization audit. No production code was changed.

Files inspected:

- `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
- `src/cartesian_pair_block_materialization/route_global_mixed_gto_blocks.jl`
- `src/cartesian_pair_block_materialization/white_lindsey_decomposed_unit_pair_inventory.jl`
- `src/cartesian_pair_block_materialization/route_global_one_body_adapter.jl`
- `src/cartesian_retained_units/records.jl`
- `src/cartesian_unit_pairs/records.jl`
- `src/precompile_workloads.jl`

Generated ignored developer artifact:

- `tmp/work/atom_gto_specialization_shape_audit.jl`
- `tmp/work/atom_gto_specialization_shape_audit_summary.txt`

The probe compares the side-7 synthetic precompile shape with the side-15 Be route shape without loading the Be BasisSets supplement, running final density assembly, or solving RHF.

Shape evidence:

| field | side-7 synthetic | side-15 Be shape |
| --- | ---: | ---: |
| parent axis counts | `(7, 7, 7)` | `(15, 15, 15)` |
| retained dimension | `223` | `615` |
| unit count | `27` | `131` |
| pair count | `378` | `8646` |
| coefficient matrix size | `(343, 223)` | `(3375, 615)` |
| unit pairs storage | `UnitPairIndexTable` | `UnitPairIndexTable` |
| retained units storage | `Vector{RetainedUnitRecord}` | `Vector{RetainedUnitRecord}` |
| `unit_keys` type | `NTuple{27, Symbol}` | `NTuple{131, Symbol}` |
| `unit_summaries` type | tuple with 27 elements encoded | tuple with 131 elements encoded |
| `pair_summaries` type | `NTuple{378, ...}` | `Symbol` |
| factorized sidecar result type | stable `NamedTuple` | stable `NamedTuple` |

Top suspected specialization sources, ranked:

1. Inventory summaries encode unit count in the type.

   `_white_lindsey_decomposed_unit_pair_inventory_result(...)` stores `unit_keys = Tuple(unit.unit_key for unit in units)` and `unit_summaries = Tuple(...)`. The side-7 workload compiles `NTuple{27, Symbol}` and a 27-element `unit_summaries` tuple, while Be side-15 needs `NTuple{131, Symbol}` and a 131-element summary tuple. This directly explains why the small precompile workload does not cover the larger route shape.

2. Small inventories still materialize pair summaries as a giant value-shaped tuple.

   The side-7 synthetic workload has `pair_summaries_type = NTuple{378, ...}`. The side-15 route avoids this only because `_WHITE_LINDSEY_COMPACT_PAIR_SUMMARY_LIMIT = 1024` makes `pair_summaries = :too_many_pairs_for_compact_summary`. This means the precompile workload itself may be compiling a tuple-heavy diagnostic shape that production side-15 does not use. That is bad precompile coverage and carrying cost.

3. Route result objects carry both compute payloads and audit/report fields in one large `NamedTuple`.

   `_white_lindsey_atom_gto_route_result(...)` returns the full staged route payloads plus counts, statuses, fallback flags, and metadata. The property set is stable, but the value types vary across success/blocking states and carry inventory/mixed/final-density result shapes through helper signatures. This is a likely source of broad specialization because downstream helpers often receive full staged objects even when they only need matrices or compact status fields.

4. Metadata carries full parent axis bundle objects inside `NamedTuple` metadata.

   The inventory metadata type includes `parent_axis_bundle_object = (; x = bundle, y = bundle, z = bundle)`. That is a large typed object embedded in report metadata. It is useful for construction, but conceptually this is parent-axis context, not summary metadata. It should move out of audit/report metadata into a compact route-owned context object.

5. Repeated axis triples are represented as fresh `NamedTuple`s everywhere.

   `axis_inputs` carries `parent_axis_counts`, `parent_axis_bundle_object`, `overlap_1d`, `kinetic_1d`, `position_1d`, and `x2_1d` as separate `x/y/z` `NamedTuple`s. Fixed axis triples are fine as tuples by policy, but these are repeated stable route concepts and should be represented once as an axis context rather than copied through keyword surfaces.

Stable concepts that should become structs:

- `ParentAxisContext3D`
  - Replace the `axis_inputs` field cloud and the repeated keyword groups `parent_axis_counts`, `parent_axis_bundle_object`, `overlap_1d`, `kinetic_1d`, `position_1d`, and `x2_1d`.
  - This should be the first struct candidate because the route repeatedly passes these fields into overlap, kinetic, nuclear, moment, and density-density calls.

- `DecomposedWLInventorySummary`
  - Replace `unit_keys::NTuple{N, Symbol}` and `unit_summaries::Tuple{... Vararg ..., N}` in hot inventory result objects with vector-backed or count-only summary data.
  - Keep detailed unit summaries available only in explicit audit/report paths.

- `DecomposedWLFactorizedBasis3D`
  - The factorized result type is already stable, but it is still a `NamedTuple` carrying a real compute sidecar. Turning it into a small struct would clarify that it is production compute data, not report metadata.

- `ResidualMomentMatrices`
  - Replace the six-field position/x2 moment bundle passed into residual MWG construction. This would prevent the residual path from treating axis moment matrix groups as a loose field cloud.

- `AtomGTOFinalBasisMatrices`
  - A lean compute product for final overlap, final Hamiltonian, and final density-density matrix. This should be separate from the giant route report object consumed by probes.

Precompile recommendation:

- Do not broaden the precompile workload to side-15 or Be-like acceptance scale yet.
- The side-7 workload misses side-15 because the hot route still has type-level unit-count summaries. A larger workload would compile one more shape rather than solve the shape-specialization problem.
- First fix the inventory/report split:
  - keep `retained_units::Vector{RetainedUnitRecord}` and `UnitPairIndexTable` as production data;
  - replace `unit_keys` and `unit_summaries` tuple fields in the hot inventory result with vector-backed summaries or count/status summaries;
  - force pair summaries to remain compact/count-only in precompile and production route paths.
- After that, rerun the side-7 precompile workload and Be warm/cold probe. If cold time still remains high, then consider a small additional workload targeting the remaining shape-stable compute objects.

Tiny cleanup made:

- None. The audit found plausible local cleanup targets, but changing inventory result shape is a real contract change across route consumers and should be a separate implementation blurb.

Validation run:

- `julia --project=. tmp/work/atom_gto_specialization_shape_audit.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old production code, tests, metadata, or compatibility path became unnecessary in this audit-only pass.
- Nothing was deleted or simplified because the identified changes affect live route result contracts and should be handled as an explicit cleanup pass.
- No tests were added. The only new artifact is an ignored `tmp/work` developer probe.
- Remaining stale/duplicate surfaces to retire next:
  - tuple-valued `unit_keys` and `unit_summaries` in hot inventory results;
  - small-inventory `pair_summaries::NTuple{N,...}` in precompile/workload paths;
  - full parent axis bundle objects carried as metadata instead of parent-axis context;
  - giant route result objects used as both compute payloads and audit reports.

-- repo-doer@macmini
