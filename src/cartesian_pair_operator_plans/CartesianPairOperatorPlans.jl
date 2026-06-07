"""
    CartesianPairOperatorPlans

Internal metadata layer for Cartesian pair-operator planning.

This module owns one step:

    retained-unit pair metadata -> pair-operator construction-plan metadata

The input is `CartesianUnitPairs.UnitPairPlan`. Pair-operator planning starts
from retained-unit pairs; it does not start from shellification regions, raw
CPBs, or terminal lowering contracts.

This module does not build source operator blocks, final pair blocks,
Hamiltonian data, exports, or artifacts.
"""
module CartesianPairOperatorPlans

using ..CartesianRouteCore
using ..CartesianUnitPairs

const CRC = CartesianRouteCore
const CUP = CartesianUnitPairs

export PairOperatorPlanPolicy,
       MetadataOnlyPairOperatorPlans,
       PairOperatorPlanRecord,
       PairOperatorPlan,
       pair_operator_plan,
       pair_operator_records,
       unavailable_summary,
       summary,
       route_core_pair_operator_plan_inventory

# File organization:
#
# records.jl
#     Pair-operator policy, record, and plan types.
#
# plan_inventory.jl
#     Unit-pair to metadata-only pair-operator planning rules.
#
# summaries.jl
#     Compact pair-operator plan summaries for tests and reports.
include("records.jl")
include("plan_inventory.jl")
include("summaries.jl")

end # module CartesianPairOperatorPlans
