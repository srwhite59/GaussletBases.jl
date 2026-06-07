"""
    CartesianPairBlockMaterialization

Internal metadata/preflight layer for Cartesian pair-block materialization.

This module owns the future transition:

    pair-operator plan metadata -> numerical pair blocks

For now it creates metadata-only readiness records. It does not materialize
source operator blocks, final pair blocks, Hamiltonian data, exports, or
artifacts.
"""
module CartesianPairBlockMaterialization

using ..CartesianCPB
using ..CartesianUnitPairs
using ..CartesianPairOperatorPlans

const CPB = CartesianCPB
const CUP = CartesianUnitPairs
const CPOP = CartesianPairOperatorPlans

export PairBlockMaterializationPolicy,
       MetadataOnlyPairBlockMaterialization,
       PairBlockMaterializationRecord,
       PairBlockMaterializationResult,
       PairBlockMaterializationPlan,
       pair_block_materialization_plan,
       pair_block_materialization_records,
       direct_direct_overlap_block,
       unavailable_summary,
       summary

# File organization:
#
# records.jl
#     Pair-block materialization policy, record, and plan types.
#
# preflight.jl
#     Metadata-only pair-operator to pair-block readiness rules.
#
# summaries.jl
#     Compact summaries for tests and reports.
#
# direct_overlap.jl
#     First tiny direct/direct overlap pair-block pilot.
include("records.jl")
include("preflight.jl")
include("summaries.jl")
include("direct_overlap.jl")

end # module CartesianPairBlockMaterialization
