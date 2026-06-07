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

using ..CartesianPairOperatorPlans

const CPOP = CartesianPairOperatorPlans

export PairBlockMaterializationPolicy,
       MetadataOnlyPairBlockMaterialization,
       PairBlockMaterializationRecord,
       PairBlockMaterializationPlan,
       pair_block_materialization_plan,
       pair_block_materialization_records,
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
include("records.jl")
include("preflight.jl")
include("summaries.jl")

end # module CartesianPairBlockMaterialization
