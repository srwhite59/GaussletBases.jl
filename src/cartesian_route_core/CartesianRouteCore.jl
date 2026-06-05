module CartesianRouteCore

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

include("coordinate_product_boxes.jl")
include("shellification_regions.jl")
include("lowering_sources.jl")
include("retained_spaces.jl")
include("unit_pairs.jl")

end # module CartesianRouteCore
