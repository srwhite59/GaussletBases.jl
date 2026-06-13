using Test
using LinearAlgebra
using GaussletBases

const _PQS_BE2_HAM_PAYLOAD_FINGERPRINT_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

function _pqs_be2_ham_payload_fingerprint_components(;
    probe_parent_axis_construction = false,
)
    system_inputs = (;
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0,
        parent_axis_counts = (x = 9, y = 7, z = 9),
        map_backend = :pgdg_localized_experimental,
    )
    spacing_inputs = (;
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = nothing,
    )
    probe_inputs = (;
        probe_parent_axis_construction,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
        probe_raw_product_box_plans = false,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
    )
    route_inputs = (;
        route_family = :pqs_source_box,
        route_kind = :be2_cartesian_nesting_route_driver_spine,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = _PQS_BE2_HAM_PAYLOAD_FINGERPRINT_TERMS,
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities = (:support_row_oracle, :dense_parent_projection),
        white_lindsey_route_shape =
            (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks,
        white_lindsey_benchmark_role =
            :published_cartesian_baseline_for_pqs_comparison,
    )

    system = GaussletBases.cartesian_system(system_inputs)
    recipe = GaussletBases.cartesian_recipe(route_inputs)
    parent = GaussletBases.cartesian_parent(
        system,
        spacing_inputs,
        probe_inputs,
        recipe,
    )
    shells = GaussletBases.cartesian_shells(parent, spacing_inputs, recipe)
    units = GaussletBases.cartesian_units(parent, shells, probe_inputs, recipe)
    transforms = GaussletBases.cartesian_transforms(units, recipe)
    pairs = GaussletBases.cartesian_pair_terms(units, transforms, recipe)
    assembly = GaussletBases.cartesian_assembly(
        parent,
        shells,
        units,
        transforms,
        pairs,
        recipe,
    )
    return (;
        system,
        recipe,
        parent,
        shells,
        units,
        transforms,
        pairs,
        assembly,
    )
end

function _pqs_be2_ham_payload_fingerprint_assembly(;
    probe_parent_axis_construction = false,
)
    return _pqs_be2_ham_payload_fingerprint_components(;
        probe_parent_axis_construction,
    ).assembly
end

@testset "Be2 PQS probe-enabled Ham readiness fingerprint" begin
    assembly = _pqs_be2_ham_payload_fingerprint_assembly(
        probe_parent_axis_construction = :auto,
    )
    @test assembly.object_kind == :cartesian_assembly
    @test assembly.route_skeleton.route_family === :pqs_source_box

    source_plan_payload =
        assembly.diatomic_complete_core_shell_source_plan_payload
    final_basis_payload =
        assembly.diatomic_complete_core_shell_final_basis_payload
    h1_payload = assembly.diatomic_complete_core_shell_h1_payload
    ham_input_payload = assembly.diatomic_complete_core_shell_ham_input_payload
    hamiltonian_handoff_payload =
        assembly.diatomic_complete_core_shell_hamiltonian_handoff_payload
    consumer_contract_payload =
        assembly.diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
    readiness = assembly.diatomic_complete_core_shell_ham_readiness_payload

    source_plan = source_plan_payload.source_plan
    @test source_plan !== nothing
    @test source_plan.object_kind ==
          :pqs_diatomic_complete_core_shell_source_plan
    @test source_plan.object_kind !== :pqs_multilayer_shell_source_plan
    @test length(source_plan.core_support_indices) == 25
    @test length(source_plan.shell_support_indices) == 250
    @test source_plan.summary.core_support_count == 25
    @test source_plan.summary.shell_support_count == 250
    @test source_plan.summary.shell_retained_count == 196
    @test source_plan.summary.precleanup_retained_dimension == 221
    @test size(source_plan.shell_final_coefficients) == (250, 196)
    @test all(isfinite, source_plan.shell_final_coefficients)
    @test source_plan.support_order == (:product, :pqs_left, :pqs_right)
    @test source_plan.route_retained_order ==
          (:pqs_left, :pqs_right, :product)
    @test isempty(
        intersect(
            source_plan.core_support_indices,
            source_plan.shell_support_indices,
        ),
    )

    final_basis = final_basis_payload.final_basis
    @test final_basis_payload.status ==
          :available_diatomic_complete_core_shell_final_basis_payload
    @test final_basis !== nothing
    @test final_basis.final_retained_count == 221
    @test final_basis.support_row_order == :core_then_shell
    @test final_basis_payload.summary.final_dimension == 221
    @test final_basis_payload.summary.precleanup_retained_dimension == 221
    @test !final_basis_payload.summary.old_source_plan_object_kind
    @test final_basis_payload.summary.final_basis_materialized

    @test h1_payload.status == :available_diatomic_complete_core_shell_h1_payload
    @test h1_payload.summary.final_dimension == 221
    h1_matrix = h1_payload.final_hamiltonian.hamiltonian_matrix
    @test size(h1_matrix) == (221, 221)
    @test all(isfinite, h1_matrix)
    @test norm(h1_matrix - h1_matrix') <= 1.0e-10
    @test isfinite(h1_payload.summary.lowest_energy)
    @test isapprox(
        h1_payload.summary.lowest_energy,
        -0.27746109235228694;
        atol = 1.0e-12,
        rtol = 0.0,
    )

    @test ham_input_payload.status ==
          :available_diatomic_complete_core_shell_ham_input_payload
    @test ham_input_payload.summary.density_gauge ==
          :pre_final_localized_positive_weight
    @test ham_input_payload.summary.raw_pair_factor_convention == :raw_numerator
    @test ham_input_payload.summary.support_weight_count == 275
    @test ham_input_payload.summary.pre_final_pair_matrix_shape == (221, 221)
    @test ham_input_payload.summary.ham_input_materialized

    @test hamiltonian_handoff_payload.status ==
          :available_diatomic_complete_core_shell_hamiltonian_handoff_payload
    @test hamiltonian_handoff_payload.summary.private_inspect_only
    @test hamiltonian_handoff_payload.summary.hamiltonian_handoff_materialized
    @test hamiltonian_handoff_payload.one_body_hamiltonian === h1_matrix
    @test hamiltonian_handoff_payload.density_interaction ===
          ham_input_payload.density_interaction
    @test hamiltonian_handoff_payload.summary.nuclear_charges == (4.0, 4.0)
    @test hamiltonian_handoff_payload.summary.nuclear_coordinates ==
          ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))
    @test hamiltonian_handoff_payload.summary.nuclear_repulsion == 4.0
    @test hamiltonian_handoff_payload.summary.electron_count == 8
    @test hamiltonian_handoff_payload.summary.spin_sector ==
          :closed_shell_singlet

    @test consumer_contract_payload.status ==
          :available_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
    @test consumer_contract_payload.summary.private_inspector_ready
    @test consumer_contract_payload.summary.final_dimension == 221
    @test consumer_contract_payload.source_handoff ===
          hamiltonian_handoff_payload
    @test consumer_contract_payload.one_body_hamiltonian === h1_matrix
    @test consumer_contract_payload.summary.two_body_representation_kind ==
          :pre_final_density_interaction
    @test consumer_contract_payload.summary.density_gauge ==
          hamiltonian_handoff_payload.summary.density_gauge
    @test consumer_contract_payload.summary.raw_pair_factor_convention ==
          hamiltonian_handoff_payload.summary.raw_pair_factor_convention
    @test !consumer_contract_payload.summary.hfdmrg_density_density_ready
    @test !consumer_contract_payload.summary.hfdmrg_sliced_ready
    @test !consumer_contract_payload.summary.hamv6_export_ready
    @test !consumer_contract_payload.summary.cr2_ready

    @test readiness.status == :blocked_diatomic_complete_core_shell_ham_readiness
    @test readiness.blocker == :missing_hfdmrg_density_density_contract
    @test :diatomic_hamiltonian_consumer_contract in readiness.available_objects
    @test !(
        :diatomic_hamiltonian_consumer_contract in readiness.missing_objects
    )
end
