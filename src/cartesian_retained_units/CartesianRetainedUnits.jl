"""
    CartesianRetainedUnits

Internal metadata layer for Cartesian retained-unit planning.

This module owns one step:

    selected terminal-lowering contracts -> planned final retained-unit records

The records produced here are column-owning metadata for later construction and
pair planning. This module chooses retained-unit granularity, for example one
direct unit for a direct slab, one PQS unit for a complete shell, or many
White--Lindsey boundary-stratum units for one complete shell.

This module does not build numerical transforms, coefficient maps, Lowdin
matrices, pair inventories, operator blocks, Hamiltonian data, artifacts, or
reports.
"""
module CartesianRetainedUnits

using ..CartesianCPB
using ..CartesianTerminalLowering
using ..CartesianRouteCore

const CPB = CartesianCPB
const CTL = CartesianTerminalLowering
const CRC = CartesianRouteCore

export RetainedUnitPolicy,
       MetadataOnlyRetainedUnits,
       RetainedUnitRecord,
       RetainedUnitPlan,
       retained_unit_plan,
       units,
       unavailable_summary,
       summary,
       route_core_final_units

# File organization:
#
# records.jl
#     Retained-unit policy, record, and plan types.
#
# lower_contract_units.jl
#     Selected lowering-contract to retained-unit metadata rules.
#
# summaries.jl
#     Compact retained-unit plan summaries for tests and reports.
include("records.jl")
include("lower_contract_units.jl")
include("summaries.jl")

end # module CartesianRetainedUnits
