using Test
using GaussletBases

const CPBMForLWSeedOracle = GaussletBases.CartesianPairBlockMaterialization

@testset "CartesianPairBlockMaterialization White-Lindsey seed oracle summary" begin
    oracle_summary =
        CPBMForLWSeedOracle.white_lindsey_materialized_seed_oracle_summary()

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
