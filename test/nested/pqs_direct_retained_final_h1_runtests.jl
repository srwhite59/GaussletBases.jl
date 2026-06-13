using Test
using GaussletBases

const PQSH1CSH = GaussletBases.CartesianShellification
const PQSH1CTL = GaussletBases.CartesianTerminalLowering

function _pqs_h1_test_bundle(count::Int)
    xmax = 8.0
    tail = 10.0
    endpoint = (count - 1) / 2
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count,
            mapping = AsinhMapping(
                a = 0.25,
                s = asinh(xmax / 0.25) / (endpoint - xmax / tail),
                tail_spacing = tail,
            ),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
end

function _pqs_h1_complete_fixture()
    expansion = coulomb_gaussian_expansion(doacc = false)
    current_box = (1:11, 1:11, 1:11)
    inner_box = (4:8, 4:8, 4:8)
    bundle11 = _pqs_h1_test_bundle(11)
    bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle11, bundle11, bundle11)
    parent_index_axes = ntuple(_ -> collect(1:11), 3)
    shellification = PQSH1CSH.shellify(
        parent_index_axes,
        (6.0, 6.0, 6.0),
        PQSH1CSH.OneCenterShellification(core_side = 5, q = 5),
        metadata = (; fixture = :pqs_direct_retained_final_h1_complete_core_shell),
    )
    lowering = PQSH1CTL.lower_terminal_regions(
        shellification,
        PQSH1CTL.PQSLowering(q = 5),
        metadata = (; fixture = :pqs_direct_retained_final_h1_complete_core_shell),
    )
    region_plan = GaussletBases.pqs_multilayer_shell_region_plan(
        shellification,
        lowering;
        metadata = (; fixture = :pqs_direct_retained_final_h1_complete_core_shell),
    )
    plan = GaussletBases.pqs_multilayer_shell_source_plan(
        bundles,
        region_plan;
        bond_axis = :z,
        term_coefficients = Float64.(expansion.coefficients),
        metadata = (;
            fixture = :pqs_direct_retained_final_h1_complete_core_shell,
            route_authority = :shellification_lowering_region_plan,
        ),
    )
    final_basis = GaussletBases.pqs_multilayer_complete_core_shell_final_basis(
        plan;
        metadata = (; fixture = :pqs_direct_retained_final_h1_complete_core_shell),
    )
    return (;
        expansion,
        bundle = bundle11,
        current_box,
        inner_box,
        raw_source_dims = (5, 5, 5),
        plan,
        region_plan,
        final_basis,
    )
end

@testset "PQS fixed-q complete core-shell inventory gate" begin
    fixture = _pqs_h1_complete_fixture()
    @test fixture.plan.status == :available_pqs_multilayer_shell_source_plan
    @test fixture.region_plan.status == :available_pqs_multilayer_shell_region_plan
    @test fixture.plan.source_kind ===
          :shellification_backed_repeated_one_cell_projected_q_shell_layers
    @test fixture.region_plan.core_box == fixture.inner_box
    @test fixture.region_plan.outer_box == fixture.current_box
    @test fixture.region_plan.summary.outer_support_coverage
    @test fixture.current_box == (1:11, 1:11, 1:11)
    @test fixture.inner_box == (4:8, 4:8, 4:8)
    @test fixture.raw_source_dims == (5, 5, 5)
    @test fixture.region_plan.summary.shell_layer_count == 3
    @test fixture.plan.layer_count == 3
    @test all(record -> record.raw_source_dims == (5, 5, 5), fixture.plan.shell_records)
    @test all(record -> record.retained_count == 98, fixture.plan.shell_records)
    @test all(
        record -> record.source_mode_shape_source === :terminal_lowering_contract,
        fixture.plan.shell_records,
    )
    @test fixture.plan.summary.fixed_source_mode_shape_used
    @test fixture.plan.summary.source_mode_shape_sources ==
          (:terminal_lowering_contract,)
    @test fixture.final_basis.core_support_count == 125
    @test fixture.final_basis.shell_support_count == 1206
    @test fixture.final_basis.shell_final_retained_count == 294
    @test fixture.final_basis.final_retained_count == 419
    @test fixture.final_basis.final_overlap_identity_error < 1.0e-10
    @test fixture.final_basis.source_plan_status ==
          :available_pqs_multilayer_shell_source_plan

    center = (;
        center_key = :origin,
        center_index = 1,
        location = (0.0, 0.0, 0.0),
        charge = 2.0,
    )
    h1_payload = GaussletBases.pqs_multilayer_complete_core_shell_h1_payload(
        fixture.plan;
        final_basis = fixture.final_basis,
        coulomb_expansion = fixture.expansion,
        center_records = (center,),
        gaussian_factor_terms_by_center =
            fixture.bundle.pgdg_intermediate.gaussian_factor_terms,
        metadata = (; fixture = :pqs_fixed_q_he_h1_gate),
    )
    hamiltonian = h1_payload.final_hamiltonian.hamiltonian_matrix
    h1 = h1_payload.h1

    @test h1_payload.status ==
          :materialized_pqs_multilayer_complete_core_shell_h1_payload
    @test all(isfinite, hamiltonian)
    @test hamiltonian ≈ transpose(hamiltonian) atol = 1.0e-10 rtol = 0.0
    @test h1.final_dimension == 419
    @test h1.solve_kind === :ordinary_symmetric
    @test isfinite(h1.lowest_energy)
    @test h1.lowest_energy < -1.0
    @test h1.generalized_overlap_solve_materialized == false
    @test h1_payload.summary.h1_solve_materialized
    @test h1_payload.summary.density_density_materialized == false
    @test h1_payload.summary.rhf_materialized == false
    @test h1.exports_materialized == false
    @test h1.artifacts_materialized == false

    provenance =
        GaussletBases.CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(
            fixture.plan.bundles;
            expected_term_count = length(fixture.expansion.coefficients),
        )
    h1_j_payload = GaussletBases.pqs_multilayer_complete_core_shell_h1_j_payload(
        fixture.plan;
        final_basis = fixture.final_basis,
        h1_payload,
        axis_weights = provenance.axis_weights,
        raw_pair_factor_terms = provenance.raw_axis_pair_factor_terms,
        coulomb_expansion = fixture.expansion,
        metadata = (; fixture = :pqs_fixed_q_he_h1j_gate),
    )
    density_interaction = h1_j_payload.density_interaction
    h1_j_summary = h1_j_payload.summary

    @test provenance.diagnostics.interaction_path === :ida_gausslet_source_box
    @test h1_j_payload.status ==
          :materialized_pqs_multilayer_complete_core_shell_h1_j_payload
    @test h1_j_summary.final_dimension == 419
    @test h1_j_summary.h1_energy ≈ h1.lowest_energy atol = 1.0e-12 rtol = 0.0
    @test h1_j_summary.h1_energy_reconstruction_error < 1.0e-10
    @test density_interaction.status ==
          :materialized_pqs_complete_core_shell_pre_final_density_interaction
    @test density_interaction.density_gauge === :pre_final_localized_positive_weight
    @test density_interaction.pre_final_pair_matrix_finite
    @test density_interaction.pre_final_weights_all_positive
    @test isfinite(h1_j_summary.self_coulomb)
    @test h1_j_summary.self_coulomb > 0.0
    @test h1_j_summary.rhf_materialized == false
    @test h1_j_summary.gto_materialized == false
    @test h1_j_summary.driver_route_materialized == false
    @test h1_j_summary.exports_materialized == false
    @test h1_j_summary.artifacts_materialized == false
end
