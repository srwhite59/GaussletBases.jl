using Test
using LinearAlgebra
using GaussletBases

const PQSH1CFBR = GaussletBases.CartesianFinalBasisRealization
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
    current_box = (1:7, 1:7, 1:7)
    inner_box = (2:6, 2:6, 2:6)
    bundle7 = _pqs_h1_test_bundle(7)
    bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle7, bundle7, bundle7)
    explicit_box_plan = GaussletBases.pqs_multilayer_shell_source_plan(
        bundles,
        inner_box,
        current_box;
        bond_axis = :z,
        term_coefficients = Float64.(expansion.coefficients),
        metadata = (;
            fixture = :pqs_direct_retained_final_h1_complete_core_shell,
        ),
    )
    parent_index_axes = ntuple(_ -> collect(1:7), 3)
    shellification = PQSH1CSH.shellify(
        parent_index_axes,
        (4.0, 4.0, 4.0),
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
    explicit_box_final_basis = GaussletBases.pqs_multilayer_complete_core_shell_final_basis(
        explicit_box_plan;
        metadata = (; fixture = :pqs_direct_retained_final_h1_complete_core_shell),
    )
    return (;
        expansion,
        bundle7,
        metrics = plan.metrics,
        current_box,
        inner_box,
        raw_source_dims = (5, 5, 5),
        plan,
        explicit_box_plan,
        shellification,
        lowering,
        region_plan,
        core_states = plan.core_support_states,
        shell_states = plan.shell_support_states,
        all_states = vcat(plan.core_support_states, plan.shell_support_states),
        final_basis,
        explicit_box_final_basis,
    )
end

function _pqs_h1_fixed_block_oracle_energy(bundle, expansion)
    fixed_block = one_center_atomic_full_parent_fixed_block(
        bundle;
        expansion,
        nside = 5,
    )
    hamiltonian = Matrix{Float64}(fixed_block.kinetic) - Matrix{Float64}(fixed_block.gaussian_sum)
    return first(eigvals(Symmetric((hamiltonian + transpose(hamiltonian)) ./ 2)))
end

@testset "PQS complete core-shell final H1 gate" begin
    fixture = _pqs_h1_complete_fixture()
    @test fixture.plan.status == :available_pqs_multilayer_shell_source_plan
    @test fixture.region_plan.status == :available_pqs_multilayer_shell_region_plan
    @test fixture.plan.source_kind ===
          :shellification_backed_repeated_one_cell_projected_q_shell_layers
    @test fixture.explicit_box_plan.status ==
          :available_pqs_multilayer_shell_source_plan
    @test fixture.region_plan.core_box == fixture.inner_box
    @test fixture.region_plan.outer_box == fixture.current_box
    @test fixture.region_plan.summary.shell_layer_count ==
          fixture.explicit_box_plan.layer_count
    @test fixture.region_plan.summary.outer_support_coverage
    @test fixture.plan.summary.core_support_count ==
          fixture.explicit_box_plan.summary.core_support_count
    @test fixture.plan.summary.shell_support_count ==
          fixture.explicit_box_plan.summary.shell_support_count
    @test fixture.plan.summary.shell_final_retained_count ==
          fixture.explicit_box_plan.summary.shell_final_retained_count
    @test fixture.final_basis.final_retained_count ==
          fixture.explicit_box_final_basis.final_retained_count
    @test fixture.final_basis.final_overlap_identity_error < 1.0e-10
    @test fixture.final_basis.source_plan_status ==
          :available_pqs_multilayer_shell_source_plan
    @test fixture.final_basis.source_plan_final_basis_helper ==
          :pqs_multilayer_complete_core_shell_final_basis

    center = (;
        center_key = :origin,
        center_index = 1,
        location = (0.0, 0.0, 0.0),
        charge = 1.0,
    )
    support_kinetic =
        GaussletBases.pqs_multilayer_support_kinetic_matrix(fixture.plan)
    support_nuclear_by_center =
        GaussletBases.pqs_multilayer_support_electron_nuclear_by_center_matrices(
            fixture.plan;
            coulomb_expansion = fixture.expansion,
            center_records = (center,),
            gaussian_factor_terms_by_center =
                fixture.bundle7.pgdg_intermediate.gaussian_factor_terms,
        )
    support_nuclear = only(support_nuclear_by_center.records)
    support_nuclear_axis_layer =
        only(
            GaussletBases.pqs_multilayer_support_electron_nuclear_by_center_matrices(
                fixture.plan;
                coulomb_expansion = fixture.expansion,
                center_records = (center,),
                axis_layers = (
                    fixture.bundle7.pgdg_intermediate.auxiliary_layer,
                    fixture.bundle7.pgdg_intermediate.auxiliary_layer,
                    fixture.bundle7.pgdg_intermediate.auxiliary_layer,
                ),
            ).records,
        )
    @test support_nuclear_axis_layer.support_operator ≈
          support_nuclear.support_operator atol = 1.0e-12 rtol = 0.0

    final_nuclear = PQSH1CFBR.pqs_complete_core_shell_final_one_body_matrix(
        fixture.final_basis,
        support_nuclear.support_operator;
        term = :electron_nuclear_by_center,
        center_record = center,
        metadata = merge(
            support_nuclear.metadata,
            (;
                nuclear_factor_source = :pgdg_intermediate_gaussian_factor_terms,
                raw_base_layer_gaussian_factor_matrices_used = false,
            ),
        ),
    )

    final_kinetic = PQSH1CFBR.pqs_complete_core_shell_final_one_body_matrix(
        fixture.final_basis,
        support_kinetic;
        term = :kinetic,
    )
    final_hamiltonian = PQSH1CFBR.pqs_complete_core_shell_final_one_electron_hamiltonian(
        final_kinetic,
        (final_nuclear,),
    )
    h1 = PQSH1CFBR.pqs_complete_core_shell_final_h1_solve(final_hamiltonian)
    fixed_h1 =
        _pqs_h1_fixed_block_oracle_energy(fixture.bundle7, fixture.expansion)

    @test fixture.current_box == (1:7, 1:7, 1:7)
    @test fixture.inner_box == (2:6, 2:6, 2:6)
    @test fixture.raw_source_dims == (5, 5, 5)
    @test fixture.final_basis.core_support_count == 125
    @test fixture.final_basis.shell_support_count == 218
    @test fixture.final_basis.shell_final_retained_count == 98
    @test fixture.final_basis.final_retained_count == 223
    @test fixture.final_basis.final_overlap_identity_error < 1.0e-10
    @test final_nuclear.metadata.nuclear_factor_source ===
          :pgdg_intermediate_gaussian_factor_terms
    @test final_nuclear.metadata.raw_base_layer_gaussian_factor_matrices_used == false
    @test final_hamiltonian.hamiltonian_matrix_finite
    @test final_hamiltonian.hamiltonian_matrix ≈
          transpose(final_hamiltonian.hamiltonian_matrix) atol = 1.0e-12 rtol = 0.0
    @test h1.solve_kind === :ordinary_symmetric
    @test h1.lowest_energy ≈ -0.48047934800387226 atol = 1.0e-12 rtol = 0.0
    @test abs(h1.lowest_energy - fixed_h1) < 1.0e-6
    @test final_kinetic.old_fixed_block_matrix_authority_used == false
    @test final_nuclear.old_fixed_block_matrix_authority_used == false
    @test final_hamiltonian.metadata.old_fixed_block_matrix_authority_used == false
    @test final_hamiltonian.metadata.current_route_safe_term_matrices_used == false
    @test h1.generalized_overlap_solve_materialized == false
    @test h1.ida_data_materialized == false
    @test h1.density_density_materialized == false
    @test h1.rhf_materialized == false
    @test h1.driver_route_materialized == false
    @test h1.exports_materialized == false
    @test h1.artifacts_materialized == false
end
