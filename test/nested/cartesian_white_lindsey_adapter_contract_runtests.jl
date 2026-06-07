using Test
using GaussletBases

const CPBMForLWAdapter = GaussletBases.CartesianPairBlockMaterialization
const CRUForLWAdapter = GaussletBases.CartesianRetainedUnits
const CRCForLWAdapter = GaussletBases.CartesianRouteCore
const CPBForLWAdapter = GaussletBases.CartesianCPB

function _lw_adapter_retained_unit(
    unit_key::Symbol,
    unit_index::Int,
    source_cpb,
    stratum_kind::Symbol,
    source_cpb_index::Int,
)
    return CRUForLWAdapter.RetainedUnitRecord(
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
        CRCForLWAdapter.owned_cpb(source_cpb),
        (source_cpb,),
        source_cpb_index,
        :not_materialized,
        nothing,
        :not_materialized,
        nothing,
        nothing,
        false,
        (; stratum_kind, source_cpb_index),
    )
end

function _lw_adapter_descriptor_units()
    facet_source = CPBForLWAdapter.slab_cpb(
        1:1,
        1:3,
        1:3;
        role = :lw_adapter_test_facet_source_cpb,
        metadata = (; stratum_kind = :facet_cpb, source_cpb_index = 1),
    )
    edge_source = CPBForLWAdapter.cpb(
        4:4,
        2:2,
        1:3;
        role = :lw_adapter_test_edge_source_cpb,
        metadata = (; stratum_kind = :edge_cpb, source_cpb_index = 2),
    )
    corner_source = CPBForLWAdapter.cpb(
        4:4,
        3:3,
        3:3;
        role = :lw_adapter_test_corner_source_cpb,
        metadata = (; stratum_kind = :corner_cpb, source_cpb_index = 3),
    )
    return (
        _lw_adapter_retained_unit(
            :lw_adapter_test_facet_unit,
            1,
            facet_source,
            :facet_cpb,
            1,
        ),
        _lw_adapter_retained_unit(
            :lw_adapter_test_edge_unit,
            2,
            edge_source,
            :edge_cpb,
            2,
        ),
        _lw_adapter_retained_unit(
            :lw_adapter_test_corner_unit,
            3,
            corner_source,
            :corner_cpb,
            3,
        ),
    )
end

@testset "CartesianPairBlockMaterialization White-Lindsey corner coefficients" begin
    facet_unit, edge_unit, corner_unit = _lw_adapter_descriptor_units()
    facet_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            facet_unit,
        )
    edge_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            edge_unit,
        )
    corner_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            corner_unit,
        )

    corner_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            corner_descriptor,
        )
    @test corner_coefficients.object_kind ==
          :white_lindsey_boundary_stratum_unit_coefficients
    @test corner_coefficients.status ==
          :materialized_white_lindsey_corner_unit_coefficients
    @test isnothing(corner_coefficients.blocker)
    @test corner_coefficients.unit_key == :lw_adapter_test_corner_unit
    @test corner_coefficients.stratum_kind == :corner_cpb
    @test corner_coefficients.source_cpb_role ==
          :lw_adapter_test_corner_source_cpb
    @test corner_coefficients.source_cpb_shape == (1, 1, 1)
    @test corner_coefficients.fixed_axis_coordinates == (
        (; axis = :x, coordinate = 4),
        (; axis = :y, coordinate = 3),
        (; axis = :z, coordinate = 3),
    )
    @test corner_coefficients.planned_old_kernel == :_nested_corner_piece
    @test corner_coefficients.coefficient_space == :source_cpb_support_local
    @test !corner_coefficients.parent_row_indices_available
    @test size(corner_coefficients.coefficient_matrix) == (1, 1)
    @test corner_coefficients.coefficient_matrix[1, 1] == 1.0
    @test corner_coefficients.source_support_row_count == 1
    @test corner_coefficients.retained_column_count == 1
    @test corner_coefficients.source_support_row_index == 1
    @test corner_coefficients.retained_column_index == 1
    @test corner_coefficients.nonzero_count == 1
    @test corner_coefficients.nonzero_values == (1.0,)
    @test corner_coefficients.coefficient_maps_materialized
    @test !corner_coefficients.source_operator_blocks_materialized
    @test !corner_coefficients.final_pair_blocks_materialized
    @test !corner_coefficients.operator_blocks_materialized
    @test !corner_coefficients.hamiltonian_data_materialized
    @test !corner_coefficients.artifacts_materialized

    corner_unit_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            corner_unit,
        )
    @test corner_unit_coefficients.status ==
          :materialized_white_lindsey_corner_unit_coefficients
    @test corner_unit_coefficients.coefficient_matrix == [1.0;;]

    facet_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            facet_descriptor,
        )
    @test facet_coefficients.status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test facet_coefficients.blocker ==
          :white_lindsey_unit_coefficients_not_implemented_for_stratum
    @test facet_coefficients.stratum_kind == :facet_cpb
    @test isnothing(facet_coefficients.coefficient_matrix)
    @test !facet_coefficients.coefficient_maps_materialized

    edge_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            edge_descriptor,
        )
    @test edge_coefficients.status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test edge_coefficients.blocker ==
          :white_lindsey_unit_coefficients_not_implemented_for_stratum
    @test edge_coefficients.stratum_kind == :edge_cpb
    @test isnothing(edge_coefficients.coefficient_matrix)
    @test !edge_coefficients.coefficient_maps_materialized

    bad_descriptor = (;
        object_kind = :white_lindsey_boundary_stratum_unit_adapter_descriptor,
        status = :available_metadata_only_white_lindsey_unit_adapter_descriptor,
        unit_key = :lw_adapter_test_bad_corner_unit,
        stratum_kind = :corner_cpb,
        planned_old_kernel = :_nested_corner_piece,
    )
    bad_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            bad_descriptor,
        )
    @test bad_coefficients.status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test bad_coefficients.blocker ==
          :white_lindsey_corner_source_cpb_not_support_local
    @test !bad_coefficients.coefficient_maps_materialized
end

@testset "CartesianPairBlockMaterialization White-Lindsey seed oracle summary" begin
    oracle_summary =
        CPBMForLWAdapter.white_lindsey_materialized_seed_oracle_summary()

    @test oracle_summary.object_kind ==
          :white_lindsey_materialized_seed_oracle_summary
    @test oracle_summary.status ==
          :available_white_lindsey_materialized_seed_oracle_summary
    @test oracle_summary.oracle_role == :validation_oracle_only
    @test !oracle_summary.route_authority
    @test !oracle_summary.adapter_authority
    @test oracle_summary.private_development_only
    @test oracle_summary.seed_report_kind ==
          :white_lindsey_low_order_materialized_seed_report
    @test oracle_summary.seed_report_status == :private_development_seed
    @test oracle_summary.route_family == :white_lindsey_low_order
    @test oracle_summary.shellization_source == :white_lindsey_one_center_seed
    @test oracle_summary.packet_kernel == :factorized_direct
    @test oracle_summary.packet_inventory_available

    @test oracle_summary.retained_dimension == 223
    @test oracle_summary.retained_unit_count == 4
    @test oracle_summary.unit_keys == (
        :low_order_core_direct,
        :low_order_face_interiors,
        :low_order_edges,
        :low_order_corners,
    )
    @test oracle_summary.unit_roles ==
          (:direct_core, :face_interiors, :edges, :corners)
    @test oracle_summary.retained_unit_kinds == (
        :white_lindsey_direct_core,
        :white_lindsey_face_interior_2d_products,
        :white_lindsey_edge_1d_side_functions,
        :white_lindsey_corner_direct_single_sites,
    )
    @test oracle_summary.retained_counts == (
        low_order_core_direct = 125,
        low_order_face_interiors = 54,
        low_order_edges = 36,
        low_order_corners = 8,
    )
    @test oracle_summary.retained_ranges == (
        low_order_core_direct = 1:125,
        low_order_face_interiors = 126:179,
        low_order_edges = 180:215,
        low_order_corners = 216:223,
    )
    @test oracle_summary.piece_counts ==
          (core = 1, faces = 6, edges = 12, corners = 8)
    @test oracle_summary.support_counts ==
          (core = 125, shell = 218, total_source = 343)
    @test oracle_summary.seed_retained_counts == (
        core = 125,
        faces = 54,
        edges = 36,
        corners = 8,
        shell = 98,
        total = 223,
    )
    @test oracle_summary.fixed_block_ready
    @test oracle_summary.overlap_ready
    @test oracle_summary.retained_basis_integral_weights_ready
    @test oracle_summary.weight_semantics == :retained_basis_integral_weights

    @test oracle_summary.operator_inventory_available
    @test oracle_summary.operator_source == :nested_fixed_block
    @test oracle_summary.operator_terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test all(values(oracle_summary.one_body_operator_matrix_available))
    @test oracle_summary.fixed_block_operator_matrix_sizes.overlap == (223, 223)
    @test oracle_summary.fixed_block_operator_matrix_sizes.kinetic == (223, 223)
    @test oracle_summary.all_operator_matrices_finite
    @test oracle_summary.operator_symmetric_ready.overlap
    @test oracle_summary.operator_symmetric_ready.kinetic
    @test oracle_summary.overlap_identity_ready
    @test !oracle_summary.operator_pairs_materialized
    @test !oracle_summary.electron_electron_materialized
    @test !oracle_summary.new_adapter_coefficient_maps_materialized
    @test !oracle_summary.new_adapter_pair_blocks_materialized
    @test !oracle_summary.hamiltonian_data_materialized
    @test !oracle_summary.artifacts_materialized

    @test !haskey(oracle_summary, :fixture)
    @test !haskey(oracle_summary, :fixed_block)
    @test !haskey(oracle_summary, :operator_matrices)
end
