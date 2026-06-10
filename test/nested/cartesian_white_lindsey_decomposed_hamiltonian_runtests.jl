using Test
using GaussletBases

const WLHamiltonianCPBM = GaussletBases.CartesianPairBlockMaterialization

@testset "decomposed WL one-electron Hamiltonian assembly" begin
    kinetic = (;
        status = :materialized_route_global_kinetic_matrix,
        global_kinetic_matrix_materialized = true,
        matrix = [1.0 0.1; 0.1 2.0],
    )
    by_center = (;
        status = :materialized_route_global_electron_nuclear_by_center_matrix_set,
        global_electron_nuclear_by_center_matrices_materialized = true,
        nuclear_charge_applied = false,
        centers_summed = false,
        by_center_matrix_count = 1,
        matrix_results = (
            (;
                center_index = 1,
                nuclear_charge = 2.0,
                matrix = [-0.4 -0.05; -0.05 -0.25],
            ),
        ),
    )

    assembled =
        WLHamiltonianCPBM.white_lindsey_decomposed_one_electron_hamiltonian(
            kinetic,
            by_center,
        )

    @test assembled.status ==
          :materialized_white_lindsey_decomposed_one_electron_hamiltonian
    @test isnothing(assembled.blocker)
    @test assembled.matrix ≈ [0.2 0.0; 0.0 1.5]
    @test assembled.nuclear_charge_application_convention ==
          :multiply_unit_charge_nuclear_attraction_by_recorded_charge
    @test assembled.nuclear_charge_applied_at_hamiltonian_assembly
    @test assembled.centers_summed_at_hamiltonian_assembly
    @test assembled.by_center_matrices_remain_separated
    @test assembled.hamiltonian_data_materialized
    @test assembled.global_hamiltonian_data_materialized
    @test !assembled.full_parent_window_cpb_used
    @test !assembled.direct_cartesian_product_assembly_used
    @test !assembled.ordinary_cartesian_ida_operators_used

    blocked =
        WLHamiltonianCPBM.white_lindsey_decomposed_one_electron_hamiltonian(
            merge(kinetic, (; global_kinetic_matrix_materialized = false)),
            by_center,
        )
    @test blocked.status ==
          :blocked_white_lindsey_decomposed_one_electron_hamiltonian
    @test blocked.blocker == :missing_route_global_kinetic_matrix
    @test isnothing(blocked.matrix)
    @test !blocked.hamiltonian_data_materialized
end
