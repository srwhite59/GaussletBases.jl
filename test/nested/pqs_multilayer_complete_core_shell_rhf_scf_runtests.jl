using Test
using LinearAlgebra
using GaussletBases

function _pqs_rhf_scf_synthetic_payloads()
    dimension = 2
    h1_payload = (;
        object_kind = :pqs_multilayer_complete_core_shell_h1_payload,
        status = :materialized_pqs_multilayer_complete_core_shell_h1_payload,
        blocker = nothing,
        final_hamiltonian = (;
            hamiltonian_matrix = Diagonal([-1.0, 0.5]),
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
        pre_final_pair_matrix = zeros(dimension, dimension),
    )
    return (; input_contract, h1_payload, density_interaction)
end

@testset "PQS complete core-shell RHF SCF payload" begin
    payloads = _pqs_rhf_scf_synthetic_payloads()
    scf = GaussletBases._pqs_multilayer_complete_core_shell_rhf_scf_payload(
        ;
        payloads...,
        max_iterations = 3,
        density_atol = 1.0e-10,
        energy_atol = 1.0e-12,
        metadata = (; fixture = :synthetic_rhf_scf),
    )

    @test scf.object_kind ==
          :pqs_multilayer_complete_core_shell_rhf_scf_payload
    @test scf.status ==
          :materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload
    @test scf.blocker === nothing
    @test isempty(scf.missing_inputs)
    @test scf.final_density ≈ Diagonal([2.0, 0.0])
    @test tr(scf.final_density) ≈ payloads.input_contract.electron_count
    @test abs.(scf.occupied_orbital_coefficients) ≈ reshape([1.0, 0.0], 2, 1)
    @test scf.final_one_step_payload.total_energy ≈ -2.0
    @test scf.final_one_step_payload.final_density ≈ scf.final_density
    @test length(scf.iteration_records) == 1
    @test scf.iteration_records[1].iteration == 1
    @test scf.iteration_records[1].density_change ≈ 0.0
    @test scf.iteration_records[1].energy_change === nothing
    @test scf.iteration_records[1].converged
    @test scf.summary.fixture_role === :route_smoke
    @test scf.summary.iteration_count == 1
    @test scf.summary.converged_iteration == 1
    @test scf.summary.first_iteration_energy_change_rule ===
          :energy_converged_when_previous_energy_missing
    @test scf.summary.density_change_rule ===
          :fixed_point_spin_summed_density_inf_norm
    @test scf.summary.residual_metric ===
          :ordinary_final_basis_commutator_inf_norm
    @test scf.summary.idempotency_rule ===
          :closed_shell_spatial_density_idempotency
    @test scf.summary.orbital_metric === :ordinary_orthonormal_final_basis
    @test scf.summary.residual_diagnostics.status ===
          :materialized_pqs_multilayer_complete_core_shell_rhf_scf_residual_diagnostics
    @test scf.summary.residual_diagnostics.density_trace_error ≈ 0.0
    @test scf.summary.residual_diagnostics.closed_shell_idempotency_error ≈ 0.0
    @test scf.summary.residual_diagnostics.commutator_residual ≈ 0.0
    @test scf.summary.residual_diagnostics.spatial_commutator_residual ≈ 0.0
    @test scf.summary.residual_diagnostics.density_change_rule ===
          scf.summary.density_change_rule
    @test scf.summary.residual_diagnostics.residual_metric ===
          scf.summary.residual_metric
    @test scf.summary.residual_diagnostics.idempotency_rule ===
          scf.summary.idempotency_rule
    @test scf.summary.residual_diagnostics.orbital_metric ===
          scf.summary.orbital_metric
    @test scf.summary.final_one_step_recomputed === true
    @test scf.summary.final_one_step_density_matches_final_density === true
    @test scf.summary.final_total_energy ≈
          scf.final_one_step_payload.total_energy
    @test scf.summary.rhf_materialized
    @test scf.summary.rhf_converged
    @test !scf.summary.driver_route_materialized
    @test !scf.summary.route_report_materialized
    @test !scf.summary.exports_materialized
    @test !scf.summary.artifacts_materialized
    @test scf.summary.public_api === false

    missing_contract =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_scf_payload(
            ;
            h1_payload = payloads.h1_payload,
            density_interaction = payloads.density_interaction,
        )
    @test missing_contract.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload
    @test missing_contract.blocker == :missing_rhf_input_contract
    @test missing_contract.missing_inputs == (:rhf_input_contract,)
    @test missing_contract.summary.density_change_rule ===
          :fixed_point_spin_summed_density_inf_norm
    @test missing_contract.summary.residual_metric ===
          :ordinary_final_basis_commutator_inf_norm
    @test missing_contract.summary.idempotency_rule ===
          :closed_shell_spatial_density_idempotency
    @test missing_contract.summary.orbital_metric ===
          :ordinary_orthonormal_final_basis
    @test missing_contract.summary.residual_diagnostics.status ===
          :blocked_pqs_multilayer_complete_core_shell_rhf_scf_residual_diagnostics
    @test missing_contract.summary.residual_diagnostics.blocker ===
          :missing_final_density

    missing_density_interaction =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_scf_payload(
            ;
            input_contract = payloads.input_contract,
            h1_payload = payloads.h1_payload,
        )
    @test missing_density_interaction.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload
    @test missing_density_interaction.blocker == :missing_density_interaction
    @test missing_density_interaction.missing_inputs == (:density_interaction,)
end
