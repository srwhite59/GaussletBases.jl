"""
    CartesianPairBlockMaterialization

Internal metadata/preflight layer for Cartesian pair-block materialization.

This module owns the transition:

    pair-operator plan metadata -> numerical pair blocks

It currently provides metadata/preflight records for pair-block
materialization, plus local direct/direct one-body pair-block pilots for:

    :overlap
    :position_x, :position_y, :position_z
    :x2_x, :x2_y, :x2_z
    :kinetic

The direct/direct one-body selector is local pair-block materialization only.
Signs and prefactors for supplied 1D factors are owned by the caller-provided
parent/axis blocks.

PQS/PQS raw-source pair support currently has metadata preflight plus narrow
raw source-space overlap, position, x2, and kinetic pilots from caller-supplied
1D factors. It does not build shell projection, Lowdin realization, or final
retained PQS pair blocks.

This module still does not assemble full operators, Hamiltonians, Ham bundles,
exports, artifacts, final retained PQS pair blocks, or White--Lindsey blocks.
"""
module CartesianPairBlockMaterialization

using ..CartesianCPB
using ..CartesianUnitPairs
using ..CartesianPairOperatorPlans
using ..CartesianRetainedUnitTransformContracts
using ..CartesianRawProductSources

const CPB = CartesianCPB
const CUP = CartesianUnitPairs
const CPOP = CartesianPairOperatorPlans
const CRTC = CartesianRetainedUnitTransformContracts
const CRPS = CartesianRawProductSources

export PairBlockMaterializationPolicy,
       MetadataOnlyPairBlockMaterialization,
       PairBlockMaterializationRecord,
       PairBlockMaterializationResult,
       PairBlockMaterializationBatchResult,
       PairBlockMaterializationPlan,
       pair_block_materialization_plan,
       pair_block_materialization_records,
       pqs_source_pair_overlap_block,
       pqs_source_pair_position_block,
       pqs_source_pair_x2_block,
       pqs_source_pair_kinetic_block,
       pqs_source_pair_overlap_blocks,
       pqs_source_pair_position_blocks,
       pqs_source_pair_x2_blocks,
       pqs_source_pair_kinetic_blocks,
       pqs_source_pair_one_body_block,
       pqs_source_pair_one_body_blocks,
       direct_direct_overlap_block,
       direct_direct_overlap_blocks,
       direct_direct_position_block,
       direct_direct_position_blocks,
       direct_direct_x2_block,
       direct_direct_x2_blocks,
       direct_direct_kinetic_block,
       direct_direct_kinetic_blocks,
       direct_direct_one_body_block,
       direct_direct_one_body_blocks,
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
# pqs_source_overlap.jl
#     Tiny PQS/PQS raw source-space safe-term pilots.
#
# direct_overlap.jl
#     First tiny direct/direct overlap pair-block pilot.
#
# direct_position.jl
#     Direct/direct position pair-block pilot.
#
# direct_x2.jl
#     Direct/direct x2 pair-block pilot.
#
# direct_kinetic.jl
#     Direct/direct kinetic pair-block pilot.
#
# direct_one_body.jl
#     Direct/direct one-body term selector.
#
# pqs_source_one_body.jl
#     PQS/PQS raw source-space one-body term selector.
include("records.jl")
include("preflight.jl")
include("summaries.jl")
include("direct_overlap.jl")
include("pqs_source_overlap.jl")
include("direct_position.jl")
include("direct_x2.jl")
include("direct_kinetic.jl")
include("direct_one_body.jl")
include("pqs_source_one_body.jl")

end # module CartesianPairBlockMaterialization
