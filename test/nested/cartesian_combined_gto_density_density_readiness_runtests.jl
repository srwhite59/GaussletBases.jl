using Test
using LinearAlgebra
using GaussletBases

const CombinedGTODensityCPBM = GaussletBases.CartesianPairBlockMaterialization

@testset "combined GTO final-basis density-density blocks on missing MWG residual representation" begin
    gausslet_dimension = 2
    raw_supplement_count = 1
    retained_supplement_count = 1
    raw_dimension = gausslet_dimension + raw_supplement_count
    final_dimension = gausslet_dimension + retained_supplement_count

    combined_matrices = (;
        status = :materialized_route_global_combined_gto_one_electron_matrices,
        gausslet_retained_dimension = gausslet_dimension,
        gausslet_retained_range = 1:gausslet_dimension,
        gto_supplement_range = (gausslet_dimension + 1):raw_dimension,
        gto_supplement_orbital_count = raw_supplement_count,
        total_combined_dimension = raw_dimension,
    )
    final_basis_projection = (;
        status = :materialized_route_global_combined_gto_final_basis_projection,
        raw_supplement_count,
        retained_supplement_count,
        final_dimension,
        transformation_matrix = Matrix{Float64}(I, raw_dimension, final_dimension),
    )
    gausslet_density_density = (;
        status = :materialized_route_global_density_density_interaction_matrix,
        matrix = [1.0 0.2; 0.2 0.9],
    )

    result =
        CombinedGTODensityCPBM.route_global_combined_gto_final_basis_density_density_matrix(
            combined_matrices,
            final_basis_projection;
            gausslet_density_density_result = gausslet_density_density,
        )

    @test result.status ===
          :blocked_route_global_combined_gto_final_basis_density_density_matrix
    @test result.blocker ===
          :missing_residual_gto_to_mwg_effective_gaussian_representation
    @test isnothing(result.gausslet_density_density_matrix)
    @test result.gausslet_density_density_matrix_shape ==
          (gausslet_dimension, gausslet_dimension)
    @test result.gausslet_density_density_matrix_matches_existing_wl_block
    @test !result.gausslet_gausslet_dense_payload_retained_in_blocked_summary
    @test result.final_dimension == final_dimension
    @test result.raw_supplement_count == raw_supplement_count
    @test result.retained_supplement_count == retained_supplement_count
    @test result.missing_density_density_blocks ==
          (
              :mixed_gausslet_residual_density_density,
              :residual_residual_density_density,
          )
    @test !result.final_density_density_matrix_materialized
    @test !result.raw_gto_density_density_used_as_final_operator
    @test !result.generalized_overlap_final_solve
    @test !result.full_parent_window_cpb_used
    @test !result.direct_cartesian_product_assembly_used
    @test !result.ordinary_cartesian_ida_operators_used
end
