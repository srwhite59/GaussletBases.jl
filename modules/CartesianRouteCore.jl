"""
    CartesianRouteCore

Standalone load-path prototype for the Cartesian route-core typed contract
layer.

This bare module mirrors the nested `GaussletBases.CartesianRouteCore` module
for direct load/precompile experiments. It depends only on the standalone
`CartesianCPB` load-path module and reuses the existing route-core source files.
"""
module CartesianRouteCore

using CartesianCPB

import CartesianCPB: role, support_count

const CPB = CartesianCPB
const CoordinateProductBox = CartesianCPB.CoordinateProductBox
const cpb = CartesianCPB.cpb
const filled_cpb = CartesianCPB.filled_cpb
const slab_cpb = CartesianCPB.slab_cpb
const complete_shell_boundary_strata = CartesianCPB.complete_shell_boundary_strata
const intervals = CartesianCPB.intervals
const shape = CartesianCPB.shape
const codimension = CartesianCPB.codimension

export CoordinateProductBox,
       OwnedSupport,
       ShellificationRegion,
       LoweringSource,
       IntermediateRetainedSpace,
       ShellRealization,
       FinalRetainedUnit,
       UnitPair,
       UnitPairInventory,
       cpb,
       filled_cpb,
       slab_cpb,
       complete_shell_support,
       complete_shell_boundary_strata,
       shellification_region,
       lowering_source,
       white_lindsey_boundary_strata_lowering,
       pqs_filled_source_lowering,
       intermediate_retained_space,
       boundary_product_mode_count,
       shell_realization,
       trivial_shell_realization,
       pqs_shell_realization,
       final_retained_unit,
       unit_pair_inventory,
       intervals,
       shape,
       codimension,
       role,
       support_count,
       owned_support,
       source_cpbs,
       lowering_recipe,
       unit_keys,
       pair_keys,
       pair_entries,
       final_units

const _ROUTE_CORE_SOURCE_DIR =
    joinpath(@__DIR__, "..", "src", "cartesian_route_core")

include(joinpath(_ROUTE_CORE_SOURCE_DIR, "shellification_regions.jl"))
include(joinpath(_ROUTE_CORE_SOURCE_DIR, "lowering_sources.jl"))
include(joinpath(_ROUTE_CORE_SOURCE_DIR, "retained_spaces.jl"))
include(joinpath(_ROUTE_CORE_SOURCE_DIR, "unit_pairs.jl"))
include(joinpath(_ROUTE_CORE_SOURCE_DIR, "pair_operator_plans.jl"))

end # module CartesianRouteCore
