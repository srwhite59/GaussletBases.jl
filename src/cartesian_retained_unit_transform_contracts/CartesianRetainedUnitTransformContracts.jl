"""
    CartesianRetainedUnitTransformContracts

Internal metadata layer for retained-unit transform planning.

This module owns one step:

    retained-unit records -> planned transform-contract metadata

The contracts produced here describe how each retained unit will eventually get
from source rows/modes to final retained columns. They do not materialize
matrices, Lowdin data, coefficient maps, pair inventories, operator blocks,
Hamiltonians, artifacts, or reports.
"""
module CartesianRetainedUnitTransformContracts

using ..CartesianRetainedUnits

const CRU = CartesianRetainedUnits

export RetainedUnitTransformContractPolicy,
       MetadataOnlyRetainedUnitTransformContracts,
       RetainedUnitTransformContract,
       RetainedUnitTransformContractPlan,
       retained_unit_transform_contract_plan,
       transform_contracts,
       unavailable_summary,
       summary

# File organization:
#
# records.jl
#     Transform-contract policy, record, and plan types.
#
# unit_contracts.jl
#     Retained-unit to transform-contract metadata rules.
#
# summaries.jl
#     Compact transform-contract plan summaries for tests and reports.
include("records.jl")
include("unit_contracts.jl")
include("summaries.jl")

end # module CartesianRetainedUnitTransformContracts
