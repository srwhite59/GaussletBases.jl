using Test
using LinearAlgebra
using GaussletBases

function _pqs_rhf_input_contract_synthetic_inputs()
    dimension = 3
    source_plan = (;
        object_kind = :pqs_multilayer_shell_source_plan,
        status = :available_pqs_multilayer_shell_source_plan,
        blocker = nothing,
        layer_count = 1,
    )
    final_basis = (;
        object_kind = :pqs_complete_core_shell_final_basis,
        status = :available_pqs_complete_core_shell_final_basis,
        blocker = nothing,
        final_retained_count = dimension,
    )
    h1_payload = (;
        object_kind = :pqs_multilayer_complete_core_shell_h1_payload,
        status = :materialized_pqs_multilayer_complete_core_shell_h1_payload,
        blocker = nothing,
        final_hamiltonian = (;
            hamiltonian_matrix = Matrix{Float64}(I, dimension, dimension),
        ),
        h1 = (final_dimension = dimension, lowest_energy = -1.0),
        summary = (final_dimension = dimension, h1_solve_materialized = true),
    )
    density_inputs = (;
        status = :available_complete_core_shell_density_inputs,
        blocker = nothing,
        axis_weights = (x = [1.0], y = [1.0], z = [1.0]),
        raw_pair_factor_terms = (
            x = ones(2, 1, 1),
            y = ones(2, 1, 1),
            z = ones(2, 1, 1),
        ),
        summary = (term_count = 2,),
    )
    coulomb_expansion = (coefficients = [0.5, 0.25],)
    return (;
        source_plan,
        final_basis,
        h1_payload,
        density_inputs,
        coulomb_expansion,
    )
end

@testset "PQS complete core-shell RHF input contract" begin
    inputs = _pqs_rhf_input_contract_synthetic_inputs()
    available =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_input_contract(
            ;
            inputs...,
            electron_count = 2,
            fixture_role = :route_smoke,
            metadata = (; fixture = :synthetic_rhf_input_contract),
        )
    @test available.object_kind ==
          :pqs_multilayer_complete_core_shell_rhf_input_contract
    @test available.status ==
          :available_pqs_multilayer_complete_core_shell_rhf_input_contract
    @test available.blocker === nothing
    @test isempty(available.missing_inputs)
    @test available.electron_count == 2
    @test available.occupation.nocc == 1
    @test available.occupation.occupancy == 2
    @test available.occupation.occupation_vector == (2,)
    @test available.fixture_role === :route_smoke
    @test available.summary.coulomb_term_count_compatible
    @test !available.summary.fock_materialized
    @test !available.summary.scf_materialized
    @test !available.summary.rhf_energy_materialized
    @test !available.summary.driver_route_materialized
    @test !available.summary.exports_materialized
    @test !available.summary.artifacts_materialized

    missing_electron_count =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_input_contract(
            ;
            inputs...,
            fixture_role = :route_smoke,
        )
    @test missing_electron_count.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_input_contract
    @test missing_electron_count.blocker == :missing_electron_count
    @test missing_electron_count.missing_inputs == (:electron_count,)
    @test missing_electron_count.occupation === nothing

    odd_electron_count =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_input_contract(
            ;
            inputs...,
            electron_count = 3,
            fixture_role = :route_smoke,
        )
    @test odd_electron_count.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_input_contract
    @test odd_electron_count.blocker == :unsupported_open_shell_rhf_input
    @test isempty(odd_electron_count.missing_inputs)
    @test odd_electron_count.occupation === nothing

    missing_fixture_role =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_input_contract(
            ;
            inputs...,
            electron_count = 2,
        )
    @test missing_fixture_role.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_input_contract
    @test missing_fixture_role.blocker == :missing_fixture_role
    @test missing_fixture_role.missing_inputs == (:fixture_role,)

    missing_h1_payload =
        GaussletBases._pqs_multilayer_complete_core_shell_rhf_input_contract(
            ;
            source_plan = inputs.source_plan,
            final_basis = inputs.final_basis,
            h1_payload = nothing,
            density_inputs = inputs.density_inputs,
            coulomb_expansion = inputs.coulomb_expansion,
            electron_count = 2,
            fixture_role = :route_smoke,
        )
    @test missing_h1_payload.status ==
          :blocked_pqs_multilayer_complete_core_shell_rhf_input_contract
    @test missing_h1_payload.blocker == :missing_h1_payload
    @test missing_h1_payload.missing_inputs == (:h1_payload,)
end
