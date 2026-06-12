Removed inventory-size-specialized tuple summaries from the hot decomposed WL inventory result.

Files changed:

- `src/cartesian_pair_block_materialization/white_lindsey_decomposed_unit_pair_inventory.jl`
- `test/nested/cartesian_white_lindsey_decomposed_unit_pair_inventory_runtests.jl`
- `docs/src/developer/numerical_contracts.md`
- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.019.md`

Generated ignored developer artifacts:

- `tmp/work/atom_gto_specialization_shape_audit.jl`
- `tmp/work/atom_gto_specialization_shape_audit_summary.txt`
- updated Be S+P probe artifacts under `tmp/work`

New hot-result shape:

- `unit_keys`: `Vector{Symbol}`
- `unit_summaries`: `Vector{NamedTuple{(:unit_key, :unit_index, :unit_kind, :stratum_kind, :source_cpb_index, :source_cpb_count, :dimension_status, :dimension, :column_range_status, :column_range)}}`
- `pair_summaries`: compact `NamedTuple`:
  - `status::Symbol`
  - `pair_count::Int64`
  - `detailed_pair_summaries_materialized::Bool`
  - `detailed_pair_summaries_source::Symbol`
  - `retained_ranges_available::Bool`
- `pair_keys`: now vector-backed for inventories below the compact limit.
- Production data remains `retained_units::Vector{RetainedUnitRecord}` and `unit_pairs::UnitPairIndexTable` or the existing compatible pair source.

Before/after specialization-shape audit:

| field | before side-7 | before side-15 | after side-7 | after side-15 |
| --- | ---: | ---: | ---: | ---: |
| parent axis counts | `(7, 7, 7)` | `(15, 15, 15)` | `(7, 7, 7)` | `(15, 15, 15)` |
| retained dimension | `223` | `615` | `223` | `615` |
| unit count | `27` | `131` | `27` | `131` |
| pair count | `378` | `8646` | `378` | `8646` |
| `unit_keys` type | `NTuple{27, Symbol}` | `NTuple{131, Symbol}` | `Vector{Symbol}` | `Vector{Symbol}` |
| `unit_summaries` type | 27-element tuple | 131-element tuple | `Vector{NamedTuple{...}}` | `Vector{NamedTuple{...}}` |
| `pair_summaries` type | `NTuple{378, ...}` | `Symbol` | compact `NamedTuple` | compact `NamedTuple` |
| coefficient matrix size | `(343, 223)` | `(3375, 615)` | `(343, 223)` | `(3375, 615)` |

Remaining expected differences:

- retained dimension, unit count, pair count, parent product size, and coefficient matrix size still differ because the physical route shapes differ.
- metadata still carries typed parent-axis bundle objects.
- the factorized sidecar result is still a `NamedTuple`, though its type did not vary between side-7 and side-15 in the audit.
- the giant atom+GTO route result still carries compute payloads and report fields together.

Test update/shrinkage:

- Updated `test/nested/cartesian_white_lindsey_decomposed_unit_pair_inventory_runtests.jl`.
- Removed assertions that treated `pair_summaries` as the detailed pair-range authority.
- Replaced them with:
  - vector-backed `unit_keys` / `unit_summaries` checks;
  - compact `pair_summaries` count/status checks;
  - detailed range checks by iterating `inventory.unit_pairs`, which is the live pair inventory.
- No new test was added.

Physics/timing validation:

- Be S+P route probe completed:
  - `rhf_total_energy = -14.574514244574639`
  - old nested/QW oracle total = `-14.574514244574694`
  - delta = `5.5067062021407764e-14 Ha`
  - final dimension = `636`
  - retained gausslet dimension = `615`
  - units/pairs = `131 / 8646`
  - no full-parent CPB, direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or generalized final solve.
- The cached-process route timing from this validation run was `31.609283375s` route elapsed and `34.203769083s` total elapsed. This was not a fresh cold comparison; it was a behavior/physics validation after the source edit.

Validation run:

- `julia --project=. test/nested/cartesian_white_lindsey_decomposed_unit_pair_inventory_runtests.jl`
- `julia --project=. tmp/work/atom_gto_specialization_shape_audit.jl`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- Deleted the hot-route materialization of detailed `pair_summaries::NTuple{N,...}`.
- Removed tuple construction for hot `unit_keys` and `unit_summaries`.
- Shrunk the focused inventory contract test so detailed pair ranges are validated from `unit_pairs` iteration instead of preserving pair-summary tuple vocabulary.
- No compatibility adapter was added.
- No new test was added.
- Remaining stale/duplicate surfaces to retire next:
  - `pair_keys` still has a small/large shape split because large inventories return an omitted sentinel; if it shows up in compile traces, make it a compact count/sample summary too.
  - parent axis bundles still ride through metadata instead of a `ParentAxisContext3D`.
  - atom+GTO route results still mix compute payloads with audit/report fields.

-- repo-doer@macmini
