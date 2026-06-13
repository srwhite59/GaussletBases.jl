using Test
using LinearAlgebra
using GaussletBases

function _pqs_rhf_one_step_synthetic_payloads()
    dimension = 2
    h1_payload = (;
        object_kind = :pqs_multilayer_complete_core_shell_h1_payload,
        status = :materialized_pqs_multilayer_complete_core_shell_h1_payload,
        blocker = nothing,
        final_hamiltonian = (;
            hamiltonian_matrix = [-1.0 0.0; 0.0 0.5],
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
            electron_count = 2,
            fixture_role = :route_smoke,
        )
    density_interaction = (;
        object_kind = :pqs_complete_core_shell_pre_final_density_interaction,
        status =
            :materialized_pqs_complete_core_shell_pre_final_density_interaction,
        blocker = nothing,
        density_gauge = :pre_final_localized_positive_weight,
        pre_final_weight_count = dimension,
        final_to_pre_final_coefficients =
            Matrix{Float64}(I, dimension, dimension),
        final_to_pre_final_reconstruction_error = 0.0,
        pre_final_pair_matrix = [2.0 0.5; 0.5 1.0],
    )
    return (; input_contract, h1_payload, density_interaction)
end

@testset "PQS complete core-shell RHF one-step payload" begin
    payloads = _pqs_rhf_one_step_synthetic_payloads()
    final_density = [2.0 0.0; 0.0 0.0]

    available =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_one_step_payload(
            ;
            payloads...,
            final_density,
            metadata = (; fixture = :synthetic_rhf_one_step),
        )

    @test available.object_kind ==
          :pqs_multilayer_complete_core_shell_rhf_one_step_payload
    @test available.status ==
          :materialized_pqs_multilayer_complete_core_shell_rhf_one_step_payload
    @test available.blocker === nothing
    @test isempty(available.missing_inputs)
    @test available.final_density == final_density
    @test available.effective_fock_matrix ≈ [1.0 0.0; 0.0 1.5]
    @test available.fock_matrix == available.effective_fock_matrix
    @test available.one_body_energy ≈ -2.0
    @test available.two_body_energy ≈ 2.0
    @test available.total_energy ≈ 0.0
    @test available.summary.fixture_role === :route_smoke
    @test available.summary.density_convention ===
          :spin_summed_closed_shell_final_density
    @test available.summary.contraction_rule ===
          :pre_final_restricted_direct_minus_exchange_from_orbital_density
    @test available.summary.fock_materialized
    @test available.summary.one_step_energy_materialized
    @test !available.summary.scf_materialized
    @test !available.summary.rhf_converged
    @test !available.summary.rhf_energy_materialized
    @test !available.summary.driver_route_materialized
    @test !available.summary.exports_materialized
    @test !available.summary.artifacts_materialized

    occupied_orbital = [1.0, 0.0]
    one_orbital_density =
        2.0 .* (occupied_orbital * transpose(occupied_orbital))
    one_orbital_payload =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_one_step_payload(
            ;
            payloads...,
            final_density = one_orbital_density,
        )
    final_basis_module = GaussletBases.CartesianFinalBasisRealization
    self_coulomb =
        final_basis_module.pqs_complete_core_shell_pre_final_orbital_self_coulomb(
            payloads.density_interaction,
            occupied_orbital,
        )
    @test one_orbital_payload.status ==
          :materialized_pqs_multilayer_complete_core_shell_rhf_one_step_payload
    @test self_coulomb.status ==
          :materialized_pqs_complete_core_shell_pre_final_orbital_self_coulomb
    @test one_orbital_payload.two_body_energy ≈ self_coulomb.self_coulomb
    @test one_orbital_payload.summary.density_convention ===
          :spin_summed_closed_shell_final_density
    @test one_orbital_payload.summary.contraction_rule ===
          :pre_final_restricted_direct_minus_exchange_from_orbital_density
    @test one_orbital_payload.summary.scf_materialized === false
    @test one_orbital_payload.summary.rhf_converged === false

    compact_payload =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_one_step_payload(
            ;
            input_contract = payloads.input_contract,
            h1_j_payload = (;
                h1_payload = payloads.h1_payload,
                density_interaction = payloads.density_interaction,
            ),
            final_density,
        )
    @test compact_payload.status == available.status
    @test compact_payload.effective_fock_matrix ≈ available.effective_fock_matrix

    missing_contract =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_one_step_payload(
            ;
            h1_payload = payloads.h1_payload,
            density_interaction = payloads.density_interaction,
            final_density,
        )
    @test missing_contract.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_one_step_payload
    @test missing_contract.blocker == :missing_rhf_input_contract
    @test missing_contract.missing_inputs == (:rhf_input_contract,)

    nonsymmetric =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_one_step_payload(
            ;
            payloads...,
            final_density = [2.0 0.1; 0.0 0.0],
        )
    @test nonsymmetric.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_one_step_payload
    @test nonsymmetric.blocker == :nonsymmetric_final_density

    trace_mismatch =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_one_step_payload(
            ;
            payloads...,
            final_density = [1.0 0.0; 0.0 0.0],
        )
    @test trace_mismatch.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_one_step_payload
    @test trace_mismatch.blocker == :electron_trace_mismatch
    @test trace_mismatch.summary.final_density_trace ≈ 1.0
end
