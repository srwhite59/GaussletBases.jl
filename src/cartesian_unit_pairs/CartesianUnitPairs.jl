"""
    CartesianUnitPairs

Internal metadata layer for Cartesian retained-unit pair planning.

This module owns one step:

    retained-unit plan -> upper-triangular retained-unit pair metadata

Pair inventories start from `CartesianRetainedUnits.RetainedUnitRecord` objects.
They do not start from shellification regions, CPBs, or lowering contracts.

This module does not build source operators, pair blocks, Hamiltonian data,
exports, or artifacts.
"""
module CartesianUnitPairs

using ..CartesianRetainedUnits
using ..CartesianRouteCore

const CRU = CartesianRetainedUnits
const CRC = CartesianRouteCore

export UnitPairPolicy,
       MetadataOnlyUnitPairs,
       UnitPairRecord,
       UnitPairPlan,
       unit_pair_plan,
       unit_pairs,
       unavailable_summary,
       summary,
       route_core_pair_inventory

# File organization:
#
# records.jl
#     Unit-pair policy, record, and plan types.
#
# pair_inventory.jl
#     Retained-unit to upper-triangular pair inventory rules.
#
# summaries.jl
#     Compact pair-plan summaries for tests and reports.
include("records.jl")
include("pair_inventory.jl")
include("summaries.jl")

end # module CartesianUnitPairs
