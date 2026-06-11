# Runtime role: compact combined gausslet+GTO matrix assembly check.
#
# This test protects the numerical placement convention for one tiny combined
# route-global overlap/Hamiltonian assembly. The H + GTO readiness test owns the
# real decomposed WL blocker.

using Test
using LinearAlgebra
using GaussletBases

const CombinedGTOLayoutCPBM = GaussletBases.CartesianPairBlockMaterialization

@testset "route-global combined GTO matrix assembly" begin
    layout = CombinedGTOLayoutCPBM.route_global_combined_gto_basis_layout(
        (;
            status = :available_white_lindsey_decomposed_unit_pair_inventory,
            source_kind = :test_decomposed_wl_unit_pair_inventory,
            unit_count = 1,
            pair_count = 1,
            retained_dimension = 2,
        ),
        (;
            supplement_kind = :test_cartesian_shell,
            orbitals = ((; label = "g1"),),
        ),
    )
    bundle = (;
        mixed_blocks = (;
            overlap = (; dense_block = [0.1; 0.2;;]),
            kinetic = (; dense_block = [0.3; 0.4;;]),
            position_x = (; dense_block = [0.5; 0.6;;]),
            position_y = (; dense_block = [0.7; 0.8;;]),
            position_z = (; dense_block = [0.9; 1.0;;]),
            x2_x = (; dense_block = [1.1; 1.2;;]),
            x2_y = (; dense_block = [1.3; 1.4;;]),
            x2_z = (; dense_block = [1.5; 1.6;;]),
        ),
        gto_blocks = (;
            overlap = (; dense_block = [1.5;;]),
            kinetic = (; dense_block = [3.0;;]),
            position_x = (; dense_block = [2.0;;]),
            position_y = (; dense_block = [2.1;;]),
            position_z = (; dense_block = [2.2;;]),
            x2_x = (; dense_block = [2.3;;]),
            x2_y = (; dense_block = [2.4;;]),
            x2_z = (; dense_block = [2.5;;]),
        ),
        mixed_nuclear_by_center_blocks = ((; dense_block = [-0.2; -0.3;;]),),
        gto_nuclear_by_center_blocks = ((; dense_block = [-0.7;;]),),
    )
    nuclear_result = (;
        matrix_results = (
            (;
                matrix = [-1.0 0.0; 0.0 -0.5],
                center_index = 1,
                center_key = :proton_a,
                center_location = (0.0, 0.0, 0.0),
                nuclear_charge = 1.0,
                nuclear_charge_applied = false,
            ),
        ),
    )

    assembled =
        CombinedGTOLayoutCPBM.route_global_combined_gto_one_electron_matrices(
            layout;
            overlap_result = (; matrix = Matrix{Float64}(I, 2, 2)),
            kinetic_result = (; matrix = [1.0 0.1; 0.1 2.0]),
            electron_nuclear_by_center_results = nuclear_result,
            gto_bundle = bundle,
            mixed_gausslet_row_range = 1:2,
        )

    @test assembled.status ==
          :materialized_route_global_combined_gto_one_electron_matrices
    @test isnothing(assembled.blocker)
    @test assembled.gausslet_retained_range == 1:2
    @test assembled.gto_supplement_range == 3:3
    @test assembled.total_combined_dimension == 3
    @test assembled.block_layout_keys ==
          (:gausslet_gausslet, :gausslet_gto, :gto_gausslet, :gto_gto)
    @test assembled.overlap_matrix ≈ [
        1.0 0.0 0.1
        0.0 1.0 0.2
        0.1 0.2 1.5
    ]
    @test assembled.hamiltonian_matrix ≈ [
        0.0 0.1 0.1
        0.1 1.5 0.1
        0.1 0.1 2.3
    ]
    @test assembled.overlap_matrix_shape == (3, 3)
    @test assembled.hamiltonian_matrix_shape == (3, 3)
    @test assembled.overlap_symmetry_error == 0.0
    @test assembled.overlap_positive_definite
    @test assembled.nuclear_charge_application_stage ==
          :route_global_combined_gto_hamiltonian_assembly
    @test assembled.nuclear_charge_applied_at_hamiltonian_assembly
    @test assembled.centers_summed_at_hamiltonian_assembly
    @test !assembled.direct_cartesian_product_assembly_used
    @test !assembled.ordinary_cartesian_ida_operators_used

    moments =
        CombinedGTOLayoutCPBM.route_global_combined_gto_residual_moment_matrices(
            layout;
            position_x_result = (; matrix = [0.0 0.0; 0.0 1.0]),
            position_y_result = (; matrix = [0.0 0.1; 0.1 0.0]),
            position_z_result = (; matrix = [1.0 0.0; 0.0 0.0]),
            x2_x_result = (; matrix = [1.0 0.2; 0.2 1.5]),
            x2_y_result = (; matrix = [1.1 0.0; 0.0 1.6]),
            x2_z_result = (; matrix = [1.2 0.1; 0.1 1.7]),
            gto_bundle = bundle,
            mixed_gausslet_row_range = 1:2,
        )

    @test moments.status ==
          :materialized_route_global_combined_gto_residual_moment_matrices
    @test isnothing(moments.blocker)
    @test moments.raw_moment_matrix_fields ==
          (:position_x, :position_y, :position_z, :x2_x, :x2_y, :x2_z)
    @test all(shape -> shape == (3, 3), values(moments.moment_matrix_shapes))
    @test moments.position_x ≈ [
        0.0 0.0 0.5
        0.0 1.0 0.6
        0.5 0.6 2.0
    ]
    @test moments.x2_z ≈ [
        1.2 0.1 1.5
        0.1 1.7 1.6
        1.5 1.6 2.5
    ]
    @test !moments.raw_gto_density_density_used_as_final_operator
    @test !moments.direct_cartesian_product_assembly_used
    @test !moments.ordinary_cartesian_ida_operators_used
end
