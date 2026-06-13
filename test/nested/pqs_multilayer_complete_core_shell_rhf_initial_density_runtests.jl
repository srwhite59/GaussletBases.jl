using Test
using LinearAlgebra
using GaussletBases

function _pqs_rhf_initial_density_synthetic_payloads(;
    dimension = 3,
    electron_count = 2,
)
    h1_diagonal = [-1.0, 0.5, 2.0][1:dimension]
    h1_payload = (;
        object_kind = :pqs_multilayer_complete_core_shell_h1_payload,
        status = :materialized_pqs_multilayer_complete_core_shell_h1_payload,
        blocker = nothing,
        final_hamiltonian = (;
            hamiltonian_matrix = Matrix(Diagonal(h1_diagonal)),
        ),
        summary = (final_dimension = dimension,),
    )
    input_contract =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_input_contract(
            ;
            source_plan = (;
                object_kind = :pqs_multilayer_shell_source_plan,
                status = :available_pqs_multilayer_shell_source_plan,
            ),
            final_basis = (;
                object_kind = :pqs_complete_core_shell_final_basis,
                status = :available_pqs_complete_core_shell_final_basis,
                final_retained_count = dimension,
            ),
            h1_payload,
            density_inputs = (;
                status = :available_complete_core_shell_density_inputs,
                axis_weights = (x = [1.0], y = [1.0], z = [1.0]),
                raw_pair_factor_terms = (
                    x = ones(2, 1, 1),
                    y = ones(2, 1, 1),
                    z = ones(2, 1, 1),
                ),
                summary = (term_count = 2,),
            ),
            coulomb_expansion = (coefficients = [1.0, 0.5],),
            electron_count,
            fixture_role = :route_smoke,
        )
    return (; input_contract, h1_payload)
end

@testset "PQS complete core-shell RHF initial density payload" begin
    payloads = _pqs_rhf_initial_density_synthetic_payloads()
    initial_density =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_initial_density_payload(
            ;
            payloads...,
            metadata = (; fixture = :synthetic_rhf_initial_density),
        )

    @test initial_density.object_kind ==
          :pqs_multilayer_complete_core_shell_rhf_initial_density_payload
    @test initial_density.status ==
          :materialized_pqs_multilayer_complete_core_shell_rhf_initial_density_payload
    @test initial_density.blocker === nothing
    @test isempty(initial_density.missing_inputs)
    @test initial_density.occupied_orbital_coefficients ≈
          reshape([1.0, 0.0, 0.0], 3, 1)
    @test initial_density.final_density ≈ Diagonal([2.0, 0.0, 0.0])
    @test initial_density.eigenvalues ≈ [-1.0, 0.5, 2.0]
    @test initial_density.occupied_eigenvalues ≈ [-1.0]
    @test initial_density.electron_trace ≈ 2.0
    @test tr(initial_density.final_density) ≈ payloads.input_contract.electron_count
    @test initial_density.final_density ≈ transpose(initial_density.final_density)
    @test initial_density.final_density^2 ≈
          2.0 .* initial_density.final_density
    @test initial_density.summary.initial_density_source === :h1_aufbau
    @test initial_density.summary.initial_density_materialized
    @test !initial_density.summary.scf_materialized
    @test !initial_density.summary.rhf_converged
    @test !initial_density.summary.rhf_energy_materialized
    @test !initial_density.summary.driver_route_materialized
    @test !initial_density.summary.exports_materialized
    @test !initial_density.summary.artifacts_materialized

    compact_payload =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_initial_density_payload(
            ;
            input_contract = payloads.input_contract,
            h1_j_payload = (; h1_payload = payloads.h1_payload),
        )
    @test compact_payload.status == initial_density.status
    @test compact_payload.final_density ≈ initial_density.final_density

    missing_contract =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_initial_density_payload(
            ;
            h1_payload = payloads.h1_payload,
        )
    @test missing_contract.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_initial_density_payload
    @test missing_contract.blocker == :missing_rhf_input_contract
    @test missing_contract.missing_inputs == (:rhf_input_contract,)

    insufficient_dimension_payloads =
        _pqs_rhf_initial_density_synthetic_payloads(
            ;
            dimension = 1,
            electron_count = 4,
        )
    insufficient_dimension =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_initial_density_payload(
            ;
            insufficient_dimension_payloads...,
        )
    @test insufficient_dimension.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_initial_density_payload
    @test insufficient_dimension.blocker ==
          :insufficient_final_dimension_for_occupation
end
