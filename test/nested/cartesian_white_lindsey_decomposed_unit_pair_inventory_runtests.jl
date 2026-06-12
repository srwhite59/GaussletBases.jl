using Test
using GaussletBases

const WLDInvCPBM = GaussletBases.CartesianPairBlockMaterialization
const WLDInvCPB = GaussletBases.CartesianCPB
const WLDInvCRU = GaussletBases.CartesianRetainedUnits
const WLDInvCRC = GaussletBases.CartesianRouteCore
const WLDInvCUP = GaussletBases.CartesianUnitPairs

function _wld_inventory_unit(
    unit_key::Symbol,
    unit_index::Int,
    source_cpb,
    stratum_kind::Symbol,
    column_range;
    dimension = isnothing(column_range) ? nothing : length(column_range),
    dimension_status = isnothing(dimension) ? :not_materialized : :available,
    column_range_status = isnothing(column_range) ? :not_materialized : :available,
)
    return WLDInvCRU.RetainedUnitRecord(
        unit_key,
        unit_index,
        :white_lindsey_boundary_stratum_retained_unit,
        Symbol(unit_key, "_contract"),
        Symbol(unit_key, "_terminal_region"),
        :synthetic_terminal_region,
        :synthetic_terminal_region,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        WLDInvCRC.owned_cpb(source_cpb),
        (source_cpb,),
        unit_index,
        dimension_status,
        dimension,
        column_range_status,
        column_range,
        nothing,
        false,
        (; stratum_kind, source_cpb_index = unit_index),
    )
end

function _wld_inventory_pair(left_unit, right_unit)
    return WLDInvCUP.UnitPairRecord(
        (left_unit.unit_key, right_unit.unit_key),
        1,
        :white_lindsey_boundary_stratum_pair,
        left_unit,
        right_unit,
        left_unit.unit_index,
        right_unit.unit_index,
        left_unit.unit_key,
        right_unit.unit_key,
        left_unit.unit_kind,
        right_unit.unit_kind,
        nothing,
        false,
        (;),
    )
end

@testset "White-Lindsey decomposed unit-pair inventory" begin
    left_cpb = WLDInvCPB.slab_cpb(
        1:1,
        2:6,
        2:6;
        role = :wld_inventory_left_facet,
        metadata = (; stratum_kind = :facet_cpb),
    )
    right_cpb = WLDInvCPB.cpb(
        7:7,
        1:1,
        2:6;
        role = :wld_inventory_right_edge,
        metadata = (; stratum_kind = :edge_cpb),
    )
    left_unit =
        _wld_inventory_unit(:wld_inventory_left, 1, left_cpb, :facet_cpb, 1:3)
    right_unit =
        _wld_inventory_unit(:wld_inventory_right, 2, right_cpb, :edge_cpb, 4:6)
    pair = _wld_inventory_pair(left_unit, right_unit)

    inventory =
        WLDInvCPBM.white_lindsey_decomposed_unit_pair_inventory((pair,))
    @test inventory.status ==
          :available_white_lindsey_decomposed_unit_pair_inventory
    @test isnothing(inventory.blocker)
    @test inventory.unit_count == 2
    @test inventory.pair_count == 1
    @test inventory.retained_dimension == 6
    @test inventory.retained_dimension_status ==
          :available_from_decomposed_wl_unit_column_ranges
    @test inventory.retained_unit_column_ranges_materialized
    @test inventory.decomposed_unit_pair_column_ranges_available
    @test inventory.term_compatibility.overlap
    @test inventory.term_compatibility.kinetic
    @test inventory.term_compatibility.electron_nuclear_by_center
    @test inventory.unit_keys isa Vector{Symbol}
    @test inventory.unit_keys == [:wld_inventory_left, :wld_inventory_right]
    @test inventory.unit_summaries isa Vector
    @test inventory.pair_summaries.status ==
          :omitted_from_hot_decomposed_wl_inventory
    @test inventory.pair_summaries.pair_count == 1
    @test !inventory.pair_summaries.detailed_pair_summaries_materialized
    first_pair = first(inventory.unit_pairs)
    @test first_pair.left_unit.column_range == 1:3
    @test first_pair.right_unit.column_range == 4:6
    @test !inventory.full_parent_window_cpb_used
    @test !inventory.direct_cartesian_product_assembly_used
    @test !inventory.ordinary_cartesian_ida_operators_used
    @test !inventory.local_operator_blocks_materialized
    @test !inventory.global_matrices_materialized
    @test !inventory.hamiltonian_data_materialized
    @test !hasproperty(inventory, :dense_block)
    @test !hasproperty(inventory, :global_matrix)

    missing_range_unit = _wld_inventory_unit(
        :wld_inventory_missing_range,
        3,
        right_cpb,
        :edge_cpb,
        nothing;
        dimension = 3,
        dimension_status = :available,
    )
    blocked_range = WLDInvCPBM.white_lindsey_decomposed_unit_pair_inventory(
        (_wld_inventory_pair(left_unit, missing_range_unit),),
    )
    @test blocked_range.status ==
          :blocked_white_lindsey_decomposed_unit_pair_inventory
    @test blocked_range.blocker == :missing_retained_unit_column_ranges
    @test !blocked_range.retained_unit_column_ranges_materialized
    @test !blocked_range.decomposed_unit_pair_column_ranges_available
    @test isnothing(blocked_range.retained_dimension)

    missing_source = WLDInvCPBM.white_lindsey_decomposed_unit_pair_inventory(nothing)
    @test missing_source.status ==
          :blocked_white_lindsey_decomposed_unit_pair_inventory
    @test missing_source.blocker ==
          :missing_decomposed_wl_unit_pair_inventory_source
    @test missing_source.unit_count == 0
    @test missing_source.pair_count == 0
    @test !missing_source.decomposed_unit_pair_column_ranges_available

    seed_report = GaussletBases._white_lindsey_low_order_materialized_seed_report()
    seed_inventory =
        WLDInvCPBM.white_lindsey_decomposed_unit_pair_inventory(seed_report)
    @test seed_inventory.status ==
          :available_white_lindsey_decomposed_unit_pair_inventory
    @test isnothing(seed_inventory.blocker)
    @test seed_inventory.source_kind ==
          :white_lindsey_low_order_materialized_seed_ranges
    @test seed_inventory.unit_count == 27
    @test seed_inventory.pair_count == 378
    @test seed_inventory.retained_dimension == 223
    @test first(seed_inventory.unit_summaries).stratum_kind == :direct_core
    @test first(seed_inventory.unit_summaries).column_range == 1:125
    @test first(seed_inventory.pair_keys) ==
          (:white_lindsey_seed_direct_core, :white_lindsey_seed_direct_core)
    @test seed_inventory.retained_unit_column_ranges_materialized
    @test seed_inventory.decomposed_unit_pair_column_ranges_available
    @test seed_inventory.term_compatibility.overlap
    @test seed_inventory.term_compatibility.kinetic
    @test seed_inventory.term_compatibility.electron_nuclear_by_center
    @test all(summary -> summary.source_cpb_count == 1, seed_inventory.unit_summaries)
    @test seed_inventory.unit_keys isa Vector{Symbol}
    @test seed_inventory.unit_summaries isa Vector
    @test seed_inventory.pair_summaries.status ==
          :omitted_from_hot_decomposed_wl_inventory
    @test seed_inventory.pair_summaries.pair_count == seed_inventory.pair_count
    @test !seed_inventory.pair_summaries.detailed_pair_summaries_materialized
    @test all(
        pair -> !isnothing(pair.left_unit.column_range) &&
                !isnothing(pair.right_unit.column_range),
        seed_inventory.unit_pairs,
    )
    @test seed_inventory.metadata.seed_inventory_source ==
          :white_lindsey_low_order_materialized_seed_inventory
    @test !seed_inventory.metadata.fixed_block_operator_matrices_used
    @test !seed_inventory.full_parent_window_cpb_used
    @test !seed_inventory.direct_cartesian_product_assembly_used
    @test !seed_inventory.ordinary_cartesian_ida_operators_used
    @test !seed_inventory.global_matrices_materialized

    shellification_source =
        WLDInvCPBM.white_lindsey_shellification_decomposed_unit_pair_inventory(
            ntuple(_ -> collect(-3.0:3.0), 3),
            ((0.0, 0.0, 0.0),);
            metadata = (; q = 5),
            parent_axis_counts = (7, 7, 7),
            parent_axis_bundle_object = nothing,
        )
    shellification_inventory = shellification_source.inventory
    @test shellification_source.status ==
          :available_white_lindsey_decomposed_unit_pair_inventory
    @test isnothing(shellification_source.blocker)
    @test shellification_source.shellification_backed_decomposed_wl_inventory
    @test !shellification_source.low_order_materialized_seed_inventory_used
    @test shellification_inventory.source_kind ==
          :cartesian_shellification_retained_unit_pair_plan
    @test shellification_inventory.metadata.shellification_unit_pair_plan_source ==
          :upper_triangular_unit_index_table
    @test shellification_inventory.unit_pairs isa WLDInvCUP.UnitPairIndexTable
    @test !shellification_inventory.metadata.rich_unit_pair_records_stored
    @test !shellification_inventory.metadata.route_core_pair_sidecars_duplicated
    @test shellification_inventory.unit_count == 27
    @test shellification_inventory.pair_count == 378
    @test shellification_inventory.retained_dimension == 223
    @test first(shellification_inventory.unit_pairs).pair_key ==
          (first(shellification_inventory.unit_summaries).unit_key,
              first(shellification_inventory.unit_summaries).unit_key)
    @test first(shellification_inventory.unit_summaries).stratum_kind == :direct_core
    @test first(shellification_inventory.unit_summaries).column_range == 1:125
    @test !shellification_source.full_parent_window_cpb_used
    @test !shellification_source.direct_cartesian_product_assembly_used
    @test !shellification_source.ordinary_cartesian_ida_operators_used
end
