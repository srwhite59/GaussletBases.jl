using Test
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

@testset "Be2 PQS Ham payload readiness fingerprint" begin
    assembly = _pqs_be2_ham_payload_fingerprint_assembly()
    @test assembly.object_kind == :cartesian_assembly
    @test assembly.route_skeleton.route_family === :pqs_source_box

    payload = assembly.complete_core_shell_diagnostic_route_payload
    support_window_payload =
        assembly.diatomic_complete_core_shell_support_window_payload
    raw_box_route_payload = assembly.diatomic_raw_box_route_payload
    source_realization_payload =
        assembly.diatomic_complete_core_shell_source_realization_payload
    source_plan_payload =
        assembly.diatomic_complete_core_shell_source_plan_payload
    final_basis_payload =
        assembly.diatomic_complete_core_shell_final_basis_payload
    h1_payload = assembly.diatomic_complete_core_shell_h1_payload
    ham_input_payload = assembly.diatomic_complete_core_shell_ham_input_payload
    readiness = assembly.diatomic_complete_core_shell_ham_readiness_payload
    ham = payload.complete_core_shell_ham_payload

    @test support_window_payload.status ==
          :available_diatomic_complete_core_shell_support_windows
    @test support_window_payload.blocker === nothing
    @test support_window_payload.parent_dims == (9, 7, 9)
    @test !support_window_payload.parent_axis_bundle_object_available
    @test support_window_payload.source_box_windows.pqs_left ==
          (1:5, 1:5, 1:5)
    @test support_window_payload.source_box_windows.product ==
          (1:5, 1:5, 5:5)
    @test support_window_payload.source_box_windows.pqs_right ==
          (1:5, 1:5, 5:9)
    @test support_window_payload.source_mode_dims.pqs_left == (5, 5, 5)
    @test support_window_payload.source_mode_dims.product == (5, 5, 1)
    @test support_window_payload.source_mode_dims.pqs_right == (5, 5, 5)
    @test support_window_payload.retained_order ==
          (:pqs_left, :pqs_right, :product)
    @test support_window_payload.candidate_core_then_shell_support_order ==
          (:product, :pqs_left, :pqs_right)
    @test support_window_payload.retained_to_support_order_permutation_required
    @test support_window_payload.support_counts ==
          (pqs_left = 125, product = 25, pqs_right = 125)
    @test :raw_product_box_plan_objects in support_window_payload.missing_objects
    @test :pqs_axis_local_coefficients in support_window_payload.missing_objects
    @test :diatomic_complete_core_shell_source_plan_materializer in
          support_window_payload.missing_objects
    @test !support_window_payload.summary.support_states_materialized
    @test !support_window_payload.summary.raw_product_box_plans_materialized
    @test !support_window_payload.summary.source_coefficients_materialized
    @test !support_window_payload.summary.source_plan_materialized
    @test !support_window_payload.summary.final_basis_materialized
    @test !support_window_payload.summary.h1_materialized
    @test !support_window_payload.summary.h1_j_materialized
    @test !support_window_payload.summary.ham_payload_materialized
    @test !support_window_payload.summary.route_driver_public_surface
    @test !support_window_payload.summary.exports_materialized
    @test !support_window_payload.summary.artifacts_materialized

    @test raw_box_route_payload.status == :blocked_diatomic_raw_box_route_payload
    @test raw_box_route_payload.blocker == :missing_parent_axis_bundle_object
    @test raw_box_route_payload.producer === nothing
    @test raw_box_route_payload.producer_status == :not_available
    @test raw_box_route_payload.support_window_payload_status ==
          :available_diatomic_complete_core_shell_support_windows
    @test :diatomic_complete_core_shell_support_windows in
          raw_box_route_payload.available_objects
    @test :parent_axis_bundle_object in raw_box_route_payload.missing_objects
    @test :raw_product_box_plan_objects in raw_box_route_payload.missing_objects
    @test :pqs_axis_local_coefficients in raw_box_route_payload.missing_objects
    @test :product_doside_unit in raw_box_route_payload.missing_objects
    @test :pair_inventory in raw_box_route_payload.missing_objects
    @test raw_box_route_payload.summary.private_candidate_only
    @test !raw_box_route_payload.summary.raw_product_box_probe_authority
    @test !raw_box_route_payload.summary.source_plan_materialized
    @test !raw_box_route_payload.summary.final_basis_materialized
    @test !raw_box_route_payload.summary.h1_materialized
    @test !raw_box_route_payload.summary.h1_j_materialized
    @test !raw_box_route_payload.summary.ham_payload_materialized
    @test !raw_box_route_payload.summary.route_driver_public_surface
    @test !raw_box_route_payload.summary.exports_materialized
    @test !raw_box_route_payload.summary.artifacts_materialized

    @test source_realization_payload.status ==
          :blocked_diatomic_complete_core_shell_source_realization
    @test source_realization_payload.blocker == :missing_parent_axis_bundle_object
    @test source_realization_payload.route_family === :pqs_source_box
    @test source_realization_payload.system_classification === :bond_aligned_diatomic
    @test source_realization_payload.bond_axis === :x
    @test !source_realization_payload.parent_axis_bundle_object_available
    @test source_realization_payload.support_window_payload_status ==
          :available_diatomic_complete_core_shell_support_windows
    @test source_realization_payload.raw_box_route_payload_status ==
          :blocked_diatomic_raw_box_route_payload
    @test source_realization_payload.core_unit_key == :product
    @test source_realization_payload.shell_unit_keys == (:pqs_left, :pqs_right)
    @test source_realization_payload.retained_order ==
          (:pqs_left, :pqs_right, :product)
    @test source_realization_payload.support_order ==
          (:product, :pqs_left, :pqs_right)
    @test source_realization_payload.retained_to_support_order_permutation_required
    @test source_realization_payload.core_support_count == 25
    @test source_realization_payload.shell_support_counts ==
          (pqs_left = 125, pqs_right = 125)
    @test source_realization_payload.shell_support_count == 250
    @test source_realization_payload.shell_retained_counts ==
          (pqs_left = nothing, pqs_right = nothing)
    @test source_realization_payload.shell_retained_count === nothing
    @test source_realization_payload.precleanup_retained_dimension === nothing
    @test source_realization_payload.shell_final_coefficients_shape === nothing
    @test source_realization_payload.shell_coefficient_block_structure ==
          :block_diagonal_left_right_pqs
    @test source_realization_payload.bundles_role == :parent_axis_bundle
    @test source_realization_payload.object_kind_claim ==
          :not_pqs_multilayer_shell_source_plan
    @test :parent_axis_bundle_object in source_realization_payload.missing_objects
    @test !source_realization_payload.summary.source_plan_materialized
    @test !source_realization_payload.summary.returns_pqs_multilayer_shell_source_plan
    @test !source_realization_payload.summary.final_basis_materialized
    @test !source_realization_payload.summary.h1_materialized
    @test !source_realization_payload.summary.h1_j_materialized
    @test !source_realization_payload.summary.ham_payload_materialized
    @test !source_realization_payload.summary.route_driver_public_surface
    @test !source_realization_payload.summary.exports_materialized
    @test !source_realization_payload.summary.artifacts_materialized

    @test source_plan_payload.status ==
          :blocked_diatomic_complete_core_shell_source_plan
    @test source_plan_payload.blocker == :missing_parent_axis_bundle_object
    @test source_plan_payload.route_family === :pqs_source_box
    @test source_plan_payload.system_classification === :bond_aligned_diatomic
    @test source_plan_payload.bond_axis === :x
    @test !source_plan_payload.parent_axis_bundle_object_available
    @test source_plan_payload.source_plan === nothing
    @test source_plan_payload.source_plan_status ==
          :not_materialized_diatomic_complete_core_shell_source_plan
    @test source_plan_payload.support_window_payload === support_window_payload
    @test source_plan_payload.source_realization_payload ===
          source_realization_payload
    @test source_plan_payload.summary.support_window_payload_status ==
          :available_diatomic_complete_core_shell_support_windows
    @test source_plan_payload.summary.raw_box_route_payload_status ==
          :blocked_diatomic_raw_box_route_payload
    @test source_plan_payload.summary.source_realization_payload_status ==
          :blocked_diatomic_complete_core_shell_source_realization
    @test :parent_axis_bundle_object in source_plan_payload.missing_objects
    @test :diatomic_complete_core_shell_source_realization_payload in
          source_plan_payload.missing_objects
    @test !source_plan_payload.summary.source_plan_materialized
    @test !source_plan_payload.summary.final_basis_materialized
    @test !source_plan_payload.summary.h1_materialized
    @test !source_plan_payload.summary.h1_j_materialized
    @test !source_plan_payload.summary.ham_payload_materialized
    @test !source_plan_payload.summary.route_driver_public_surface
    @test !source_plan_payload.summary.exports_materialized
    @test !source_plan_payload.summary.artifacts_materialized

    @test final_basis_payload.status ==
          :blocked_diatomic_complete_core_shell_final_basis_payload
    @test final_basis_payload.blocker ==
          :missing_pqs_diatomic_complete_core_shell_source_plan
    @test final_basis_payload.source_plan === nothing
    @test final_basis_payload.source_plan_status == :not_available
    @test final_basis_payload.final_basis === nothing
    @test final_basis_payload.final_basis_status ==
          :not_materialized_diatomic_complete_core_shell_final_basis
    @test :pqs_diatomic_complete_core_shell_source_plan in
          final_basis_payload.missing_objects
    @test !final_basis_payload.summary.final_basis_materialized
    @test !final_basis_payload.summary.h1_materialized
    @test !final_basis_payload.summary.h1_j_materialized
    @test !final_basis_payload.summary.ham_payload_materialized
    @test !final_basis_payload.summary.route_driver_public_surface
    @test !final_basis_payload.summary.exports_materialized
    @test !final_basis_payload.summary.artifacts_materialized

    @test h1_payload.status == :blocked_diatomic_complete_core_shell_h1_payload
    @test h1_payload.blocker ==
          :missing_pqs_diatomic_complete_core_shell_source_plan
    @test h1_payload.source_plan === nothing
    @test h1_payload.source_plan_status == :not_available
    @test h1_payload.final_basis === nothing
    @test h1_payload.final_basis_status == :not_available
    @test h1_payload.support_kinetic_status ==
          :not_materialized_diatomic_support_kinetic_matrix
    @test h1_payload.support_electron_nuclear_status ==
          :not_materialized_diatomic_support_electron_nuclear_by_center
    @test h1_payload.h1_status ==
          :not_materialized_diatomic_complete_core_shell_h1
    @test !h1_payload.summary.h1_materialized
    @test !h1_payload.summary.h1_j_materialized
    @test !h1_payload.summary.density_density_materialized
    @test !h1_payload.summary.ham_payload_materialized
    @test !h1_payload.summary.rhf_materialized
    @test !h1_payload.summary.public_api
    @test !h1_payload.summary.exports_materialized
    @test !h1_payload.summary.artifacts_materialized

    @test ham_input_payload.status ==
          :blocked_diatomic_complete_core_shell_ham_input_payload
    @test ham_input_payload.blocker ==
          :missing_diatomic_complete_core_shell_source_plan
    @test ham_input_payload.source_plan === nothing
    @test ham_input_payload.source_plan_status == :not_available
    @test ham_input_payload.final_basis === nothing
    @test ham_input_payload.final_basis_status == :not_available
    @test ham_input_payload.h1_payload === h1_payload
    @test ham_input_payload.h1_payload_status ==
          :blocked_diatomic_complete_core_shell_h1_payload
    @test ham_input_payload.density_provenance_status ==
          :not_materialized_diatomic_complete_core_shell_density_provenance
    @test ham_input_payload.support_weights_status ==
          :not_materialized_diatomic_complete_core_shell_support_weights
    @test ham_input_payload.raw_pair_factor_status ==
          :not_materialized_diatomic_raw_pair_factor_terms
    @test ham_input_payload.support_pair_raw_numerator_status ==
          :not_materialized_diatomic_support_raw_pair_numerator_matrix
    @test ham_input_payload.density_interaction_status ==
          :not_materialized_diatomic_pre_final_density_interaction
    @test !ham_input_payload.summary.ham_input_materialized
    @test !ham_input_payload.summary.h1_j_materialized
    @test !ham_input_payload.summary.ham_payload_materialized
    @test !ham_input_payload.summary.rhf_materialized
    @test !ham_input_payload.summary.public_api
    @test !ham_input_payload.summary.exports_materialized
    @test !ham_input_payload.summary.artifacts_materialized

    @test readiness.status == :blocked_diatomic_complete_core_shell_ham_readiness
    @test readiness.blocker ==
          :missing_diatomic_complete_core_shell_source_plan_producer
    @test readiness.route_family === :pqs_source_box
    @test readiness.system_classification === :bond_aligned_diatomic
    @test readiness.bond_axis === :x
    @test !readiness.parent_axis_bundle_object_available

    @test readiness.center_summary.center_count == 2
    @test readiness.center_summary.nuclear_charges == (4, 4)
    @test readiness.source_box_summary.source_box_count == 3
    @test readiness.source_box_summary.source_box_keys ==
          (:pqs_left, :product, :pqs_right)
    @test readiness.retained_unit_summary.retained_unit_count == 3
    @test readiness.retained_unit_summary.unit_keys ==
          (:pqs_left, :pqs_right, :product)
    @test readiness.retained_unit_summary.retained_unit_kinds ==
          (:pqs, :pqs, :product_doside)
    @test readiness.pair_inventory_summary.pair_count == 6
    @test readiness.pair_inventory_summary.pair_family_counts ==
          (pqs_pqs = 3, pqs_product = 2, product_pqs = 0, product_product = 1)
    @test readiness.pair_inventory_summary.source_box_algorithmic_path_for_all_pairs

    @test :route_skeleton in readiness.available_objects
    @test :source_boxes in readiness.available_objects
    @test :retained_units in readiness.available_objects
    @test :pair_inventory in readiness.available_objects
    @test :diatomic_complete_core_shell_source_plan_producer in
          readiness.missing_objects
    @test :parent_axis_bundle_object in readiness.missing_objects
    @test !readiness.summary.public_api
    @test !readiness.summary.final_basis_materialized
    @test !readiness.summary.h1_materialized
    @test !readiness.summary.h1_j_materialized
    @test !readiness.summary.ham_payload_materialized
    @test !readiness.summary.rhf_materialized
    @test !readiness.summary.exports_materialized
    @test !readiness.summary.artifacts_materialized

    @test payload.status == :blocked_missing_complete_core_shell_h1_j_route_inputs
    @test payload.blocker == :missing_complete_core_shell_h1_j_route_inputs
    @test payload.missing_inputs == (
        :pqs_multilayer_shell_region_plan,
        :pqs_multilayer_shell_source_plan,
        :pqs_multilayer_complete_core_shell_final_basis,
        :pqs_multilayer_complete_core_shell_h1_payload,
        :axis_weights,
        :raw_pair_factor_terms,
        :coulomb_expansion,
    )
    @test payload.source_payload.status ==
          :blocked_missing_complete_core_shell_source_plan_inputs
    @test payload.h1_j_payload.status ==
          :blocked_missing_complete_core_shell_h1_j_route_inputs

    @test ham.object_kind == :pqs_source_box_complete_core_shell_ham_payload
    @test ham.status == :blocked_complete_core_shell_ham_payload
    @test ham.blocker == :missing_complete_core_shell_ham_inputs
    @test ham.missing_inputs == (
        :pqs_multilayer_complete_core_shell_final_basis,
        :pqs_multilayer_complete_core_shell_h1_payload,
        :pqs_complete_core_shell_final_one_electron_hamiltonian,
        :complete_core_shell_density_inputs,
        :complete_core_shell_h1_j_diagnostic_payload,
        :pqs_complete_core_shell_pre_final_density_interaction,
    )
    @test !ham.summary.public_api
    @test !ham.summary.exports_materialized
    @test !ham.summary.artifacts_materialized
    @test !ham.summary.rhf_product_surface
    @test !ham.summary.serious_hf_claim
end

@testset "Be2 PQS probe-enabled Ham readiness fingerprint" begin
    assembly = _pqs_be2_ham_payload_fingerprint_assembly(
        probe_parent_axis_construction = :auto,
    )
    @test assembly.object_kind == :cartesian_assembly
    @test assembly.route_skeleton.route_family === :pqs_source_box

    readiness = assembly.diatomic_complete_core_shell_ham_readiness_payload
    support_window_payload =
        assembly.diatomic_complete_core_shell_support_window_payload
    raw_box_route_payload = assembly.diatomic_raw_box_route_payload
    source_realization_payload =
        assembly.diatomic_complete_core_shell_source_realization_payload
    source_plan_payload =
        assembly.diatomic_complete_core_shell_source_plan_payload
    final_basis_payload =
        assembly.diatomic_complete_core_shell_final_basis_payload
    h1_payload = assembly.diatomic_complete_core_shell_h1_payload
    ham_input_payload = assembly.diatomic_complete_core_shell_ham_input_payload
    payload = assembly.complete_core_shell_diagnostic_route_payload
    ham = payload.complete_core_shell_ham_payload

    @test support_window_payload.status ==
          :available_diatomic_complete_core_shell_support_windows
    @test support_window_payload.blocker === nothing
    @test support_window_payload.parent_axis_bundle_object_available
    @test support_window_payload.parent_dims isa NTuple{3,Int}
    @test support_window_payload.source_box_windows.pqs_left isa
          NTuple{3,UnitRange{Int}}
    @test support_window_payload.source_box_windows.product isa
          NTuple{3,UnitRange{Int}}
    @test support_window_payload.source_box_windows.pqs_right isa
          NTuple{3,UnitRange{Int}}
    @test length.(support_window_payload.source_box_windows.pqs_left) ==
          support_window_payload.source_mode_dims.pqs_left
    @test length.(support_window_payload.source_box_windows.product) ==
          support_window_payload.source_mode_dims.product
    @test length.(support_window_payload.source_box_windows.pqs_right) ==
          support_window_payload.source_mode_dims.pqs_right
    @test support_window_payload.retained_order ==
          (:pqs_left, :pqs_right, :product)
    @test support_window_payload.candidate_core_then_shell_support_order ==
          (:product, :pqs_left, :pqs_right)
    @test support_window_payload.retained_to_support_order_permutation_required
    @test support_window_payload.support_counts ==
          (pqs_left = 125, product = 25, pqs_right = 125)
    @test :raw_product_box_plan_objects in support_window_payload.missing_objects
    @test :pqs_axis_local_coefficients in support_window_payload.missing_objects
    @test :diatomic_complete_core_shell_source_plan_materializer in
          support_window_payload.missing_objects
    @test !support_window_payload.summary.support_states_materialized
    @test !support_window_payload.summary.raw_product_box_plans_materialized
    @test !support_window_payload.summary.source_coefficients_materialized
    @test !support_window_payload.summary.source_plan_materialized
    @test !support_window_payload.summary.final_basis_materialized
    @test !support_window_payload.summary.h1_materialized
    @test !support_window_payload.summary.h1_j_materialized
    @test !support_window_payload.summary.ham_payload_materialized
    @test !support_window_payload.summary.route_driver_public_surface
    @test !support_window_payload.summary.exports_materialized
    @test !support_window_payload.summary.artifacts_materialized

    @test raw_box_route_payload.status == :available_diatomic_raw_box_route_payload
    @test raw_box_route_payload.blocker === nothing
    @test raw_box_route_payload.producer_status == :private_shadow_only
    @test raw_box_route_payload.producer !== nothing
    @test raw_box_route_payload.support_window_payload_status ==
          :available_diatomic_complete_core_shell_support_windows
    @test raw_box_route_payload.descriptor_summary.object_kind ==
          :pqs_pqs_product_safe_term_route_descriptor
    @test raw_box_route_payload.descriptor_summary.roles ==
          (:pqs_left, :pqs_right, :product)
    @test raw_box_route_payload.descriptor_summary.retained_dimension == 221
    @test raw_box_route_payload.descriptor_summary.expected_pair_count == 6
    left_raw_plan_summary =
        raw_box_route_payload.raw_product_box_plan_summary.pqs_left
    right_raw_plan_summary =
        raw_box_route_payload.raw_product_box_plan_summary.pqs_right
    @test left_raw_plan_summary.object_kind ==
          :cartesian_raw_product_box_plan_3d
    @test right_raw_plan_summary.object_kind ==
          :cartesian_raw_product_box_plan_3d
    @test left_raw_plan_summary.source_mode_dims == (5, 5, 5)
    @test right_raw_plan_summary.source_mode_dims == (5, 5, 5)
    @test left_raw_plan_summary.axis_local_coefficient_shapes ==
          ((5, 5), (5, 5), (5, 5))
    @test right_raw_plan_summary.axis_local_coefficient_shapes ==
          ((5, 5), (5, 5), (5, 5))
    left_raw_pqs_summary = raw_box_route_payload.raw_pqs_plan_summary.pqs_left
    right_raw_pqs_summary = raw_box_route_payload.raw_pqs_plan_summary.pqs_right
    @test left_raw_pqs_summary.representation == :orthogonal_raw_product_box
    @test right_raw_pqs_summary.representation == :orthogonal_raw_product_box
    @test left_raw_pqs_summary.boundary_selected_count == 98
    @test right_raw_pqs_summary.boundary_selected_count == 98
    @test raw_box_route_payload.product_unit_summary.kind == :product_doside
    @test raw_box_route_payload.product_unit_summary.support_count == 25
    @test raw_box_route_payload.product_unit_summary.support_state_count == 25
    @test raw_box_route_payload.product_unit_summary.coefficient_matrix_shape ==
          (25, 25)
    raw_pair_inventory_summary = raw_box_route_payload.pair_inventory_summary
    @test raw_pair_inventory_summary.every_pair_uses_source_box_algorithmic_policy
    @test raw_pair_inventory_summary.source_box_algorithmic_pair_count == 6
    @test :raw_product_box_plan_objects in raw_box_route_payload.available_objects
    @test :pqs_axis_local_coefficients in raw_box_route_payload.available_objects
    @test :product_doside_unit in raw_box_route_payload.available_objects
    @test :pair_inventory in raw_box_route_payload.available_objects
    @test :diatomic_complete_core_shell_source_plan_materializer in
          raw_box_route_payload.missing_objects
    @test raw_box_route_payload.summary.private_candidate_only
    @test !raw_box_route_payload.summary.raw_product_box_probe_authority
    @test !raw_box_route_payload.summary.source_plan_materialized
    @test !raw_box_route_payload.summary.final_basis_materialized
    @test !raw_box_route_payload.summary.h1_materialized
    @test !raw_box_route_payload.summary.h1_j_materialized
    @test !raw_box_route_payload.summary.ham_payload_materialized
    @test !raw_box_route_payload.summary.route_driver_public_surface
    @test !raw_box_route_payload.summary.exports_materialized
    @test !raw_box_route_payload.summary.artifacts_materialized

    @test source_realization_payload.status ==
          :available_diatomic_complete_core_shell_source_realization
    @test source_realization_payload.blocker === nothing
    @test source_realization_payload.route_family === :pqs_source_box
    @test source_realization_payload.system_classification === :bond_aligned_diatomic
    @test source_realization_payload.bond_axis === :x
    @test source_realization_payload.parent_axis_bundle_object_available
    @test source_realization_payload.support_window_payload_status ==
          :available_diatomic_complete_core_shell_support_windows
    @test source_realization_payload.raw_box_route_payload_status ==
          :available_diatomic_raw_box_route_payload
    @test source_realization_payload.core_unit_key == :product
    @test source_realization_payload.shell_unit_keys == (:pqs_left, :pqs_right)
    @test source_realization_payload.retained_order ==
          (:pqs_left, :pqs_right, :product)
    @test source_realization_payload.support_order ==
          (:product, :pqs_left, :pqs_right)
    @test source_realization_payload.retained_to_support_order_permutation_required
    @test source_realization_payload.route_retained_ranges ==
          raw_box_route_payload.producer.descriptor.expected_ranges
    @test source_realization_payload.source_plan_precleanup_ranges ==
          (product = 1:25, pqs_left = 26:123, pqs_right = 124:221)
    @test source_realization_payload.core_support_count == 25
    @test source_realization_payload.shell_support_counts ==
          (pqs_left = 125, pqs_right = 125)
    @test source_realization_payload.shell_support_count == 250
    @test source_realization_payload.shell_retained_counts ==
          (pqs_left = 98, pqs_right = 98)
    @test source_realization_payload.shell_retained_count == 196
    @test source_realization_payload.precleanup_retained_dimension == 221
    @test source_realization_payload.shell_final_coefficients_shape == (250, 196)
    @test source_realization_payload.shell_coefficient_block_structure ==
          :block_diagonal_left_right_pqs
    @test source_realization_payload.bundles_role == :parent_axis_bundle
    @test source_realization_payload.object_kind_claim ==
          :not_pqs_multilayer_shell_source_plan
    @test :diatomic_complete_core_shell_source_realization in
          source_realization_payload.available_objects
    @test :pqs_multilayer_shell_source_plan_adapter_contract in
          source_realization_payload.missing_objects
    @test !source_realization_payload.summary.source_plan_materialized
    @test !source_realization_payload.summary.returns_pqs_multilayer_shell_source_plan
    @test !source_realization_payload.summary.final_basis_materialized
    @test !source_realization_payload.summary.h1_materialized
    @test !source_realization_payload.summary.h1_j_materialized
    @test !source_realization_payload.summary.ham_payload_materialized
    @test !source_realization_payload.summary.route_driver_public_surface
    @test !source_realization_payload.summary.exports_materialized
    @test !source_realization_payload.summary.artifacts_materialized

    @test source_plan_payload.status ==
          :blocked_diatomic_complete_core_shell_source_plan
    @test source_plan_payload.blocker ==
          :missing_diatomic_complete_core_shell_final_basis_consumer
    @test source_plan_payload.route_family === :pqs_source_box
    @test source_plan_payload.system_classification === :bond_aligned_diatomic
    @test source_plan_payload.bond_axis === :x
    @test source_plan_payload.parent_axis_bundle_object_available
    source_plan = source_plan_payload.source_plan
    @test source_plan !== nothing
    @test source_plan.object_kind ==
          :pqs_diatomic_complete_core_shell_source_plan
    @test source_plan.object_kind !== :pqs_multilayer_shell_source_plan
    @test source_plan.status ==
          :available_pqs_diatomic_complete_core_shell_source_plan
    @test source_plan.blocker === nothing
    @test source_plan.bundles !== nothing
    @test source_plan.metrics !== nothing
    @test source_plan.core_unit_key == :product
    @test source_plan.shell_unit_keys == (:pqs_left, :pqs_right)
    @test length(source_plan.core_support_indices) == 25
    @test length(source_plan.core_support_states) == 25
    @test length(source_plan.shell_support_indices) == 250
    @test length(source_plan.shell_support_states) == 250
    @test isempty(
        intersect(
            source_plan.core_support_indices,
            source_plan.shell_support_indices,
        ),
    )
    @test length(unique(source_plan.shell_support_indices)) == 250
    @test size(source_plan.shell_final_coefficients) == (250, 196)
    @test all(isfinite, source_plan.shell_final_coefficients)
    @test source_plan.support_order == (:product, :pqs_left, :pqs_right)
    @test source_plan.route_retained_order ==
          (:pqs_left, :pqs_right, :product)
    @test source_plan.retained_pre_final_map.precleanup_ranges ==
          (product = 1:25, pqs_left = 26:123, pqs_right = 124:221)
    @test source_plan.convention_labels.old_source_plan_object_kind == false
    @test source_plan.summary.old_source_plan_object_kind == false
    @test source_plan.summary.core_support_count == 25
    @test source_plan.summary.shell_support_count == 250
    @test source_plan.summary.shell_retained_count == 196
    @test source_plan.summary.precleanup_retained_dimension == 221
    @test source_plan.summary.shell_final_coefficients_shape == (250, 196)
    @test source_plan.summary.shell_coefficient_block_structure ==
          :block_diagonal_left_right_pqs
    @test source_plan.summary.source_plan_materialized
    @test !source_plan.summary.final_basis_materialized
    @test !source_plan.summary.h1_materialized
    @test !source_plan.summary.h1_j_materialized
    @test !source_plan.summary.ham_payload_materialized
    @test !source_plan.summary.route_driver_public_surface
    @test !source_plan.summary.exports_materialized
    @test !source_plan.summary.artifacts_materialized
    @test source_plan_payload.source_plan_status ==
          :available_pqs_diatomic_complete_core_shell_source_plan
    @test source_plan_payload.support_window_payload === support_window_payload
    @test source_plan_payload.source_realization_payload ===
          source_realization_payload
    @test source_plan_payload.summary.support_window_payload_status ==
          :available_diatomic_complete_core_shell_support_windows
    @test source_plan_payload.summary.raw_box_route_payload_status ==
          :available_diatomic_raw_box_route_payload
    @test source_plan_payload.summary.source_realization_payload_status ==
          :available_diatomic_complete_core_shell_source_realization
    @test :diatomic_complete_core_shell_support_windows in
          source_plan_payload.available_objects
    @test :diatomic_raw_box_route_payload in source_plan_payload.available_objects
    @test :diatomic_complete_core_shell_source_realization in
          source_plan_payload.available_objects
    @test :pqs_diatomic_complete_core_shell_source_plan in
          source_plan_payload.available_objects
    @test :parent_axis_bundle_object in source_plan_payload.available_objects
    @test !in(:parent_axis_bundle_object, source_plan_payload.missing_objects)
    @test :diatomic_complete_core_shell_final_basis_consumer in
          source_plan_payload.missing_objects
    @test source_plan_payload.summary.source_plan_materialized
    @test !source_plan_payload.summary.final_basis_materialized
    @test !source_plan_payload.summary.h1_materialized
    @test !source_plan_payload.summary.h1_j_materialized
    @test !source_plan_payload.summary.ham_payload_materialized
    @test !source_plan_payload.summary.route_driver_public_surface
    @test !source_plan_payload.summary.exports_materialized
    @test !source_plan_payload.summary.artifacts_materialized

    @test final_basis_payload.status ==
          :available_diatomic_complete_core_shell_final_basis_payload
    @test final_basis_payload.blocker === nothing
    @test final_basis_payload.source_plan === source_plan
    @test final_basis_payload.source_plan_status ==
          :available_pqs_diatomic_complete_core_shell_source_plan
    @test final_basis_payload.final_basis !== nothing
    @test final_basis_payload.final_basis_status ==
          :available_pqs_complete_core_shell_final_basis
    final_basis = final_basis_payload.final_basis
    @test final_basis.status == :available_pqs_complete_core_shell_final_basis
    @test final_basis.support_row_order == :core_then_shell
    @test final_basis.core_support_count == 25
    @test final_basis.shell_support_count == 250
    @test final_basis.shell_final_retained_count == 196
    @test final_basis.final_retained_count == 221
    @test final_basis_payload.summary.final_dimension == 221
    @test final_basis_payload.summary.precleanup_retained_dimension == 221
    @test final_basis_payload.summary.core_support_count == 25
    @test final_basis_payload.summary.shell_support_count == 250
    @test final_basis_payload.summary.shell_retained_count == 196
    @test final_basis_payload.summary.support_row_order == :core_then_shell
    @test !final_basis_payload.summary.old_source_plan_object_kind
    @test final_basis_payload.summary.final_basis_materialized
    @test !final_basis_payload.summary.h1_materialized
    @test !final_basis_payload.summary.h1_j_materialized
    @test !final_basis_payload.summary.ham_payload_materialized
    @test !final_basis_payload.summary.route_driver_public_surface
    @test !final_basis_payload.summary.exports_materialized
    @test !final_basis_payload.summary.artifacts_materialized
    @test :diatomic_complete_core_shell_final_basis in
          final_basis_payload.available_objects
    @test :diatomic_complete_core_shell_h1_consumer in
          final_basis_payload.missing_objects

    @test h1_payload.status == :available_diatomic_complete_core_shell_h1_payload
    @test h1_payload.blocker === nothing
    @test h1_payload.source_plan === source_plan
    @test h1_payload.source_plan_status ==
          :available_pqs_diatomic_complete_core_shell_source_plan
    @test h1_payload.final_basis === final_basis
    @test h1_payload.final_basis_status ==
          :available_pqs_complete_core_shell_final_basis
    @test h1_payload.support_kinetic_status ==
          :materialized_diatomic_complete_core_shell_support_kinetic_matrix
    @test h1_payload.support_electron_nuclear_status ==
          :materialized_diatomic_complete_core_shell_support_electron_nuclear_by_center_matrix_set
    @test h1_payload.final_kinetic_status ==
          :materialized_pqs_complete_core_shell_final_one_body_matrix
    @test h1_payload.final_electron_nuclear_status ==
          :materialized_diatomic_final_electron_nuclear_by_center
    @test h1_payload.final_hamiltonian_status ==
          :materialized_pqs_complete_core_shell_final_one_electron_hamiltonian
    @test h1_payload.h1_status ==
          :materialized_pqs_complete_core_shell_final_h1_solve
    @test h1_payload.summary.final_dimension == 221
    @test isfinite(h1_payload.summary.lowest_energy)
    @test h1_payload.h1.final_dimension == 221
    @test h1_payload.summary.center_count == 2
    @test h1_payload.summary.support_row_order == :core_then_shell
    @test h1_payload.summary.source_plan_support_row_order ==
          :core_product_then_shell_left_right_pqs
    @test h1_payload.summary.h1_materialized
    @test !h1_payload.summary.h1_j_materialized
    @test !h1_payload.summary.density_density_materialized
    @test !h1_payload.summary.ham_payload_materialized
    @test !h1_payload.summary.rhf_materialized
    @test !h1_payload.summary.public_api
    @test !h1_payload.summary.exports_materialized
    @test !h1_payload.summary.artifacts_materialized
    @test :diatomic_complete_core_shell_h1_payload in
          h1_payload.available_objects
    @test isempty(h1_payload.missing_objects)

    @test ham_input_payload.status ==
          :available_diatomic_complete_core_shell_ham_input_payload
    @test ham_input_payload.blocker === nothing
    @test ham_input_payload.source_plan === source_plan
    @test ham_input_payload.source_plan_status ==
          :available_pqs_diatomic_complete_core_shell_source_plan
    @test ham_input_payload.final_basis === final_basis
    @test ham_input_payload.final_basis_status ==
          :available_pqs_complete_core_shell_final_basis
    @test ham_input_payload.h1_payload === h1_payload
    @test ham_input_payload.h1_payload_status ==
          :available_diatomic_complete_core_shell_h1_payload
    @test ham_input_payload.density_provenance_status ==
          :available_diatomic_complete_core_shell_density_provenance
    @test ham_input_payload.support_weights_status ==
          :materialized_diatomic_complete_core_shell_support_weights
    @test ham_input_payload.raw_pair_factor_status ==
          :available_diatomic_raw_pair_factor_terms
    @test ham_input_payload.support_pair_raw_numerator_status ==
          :materialized_diatomic_complete_core_shell_support_raw_pair_numerator_matrix
    @test ham_input_payload.density_interaction_status ==
          :materialized_pqs_complete_core_shell_pre_final_density_interaction
    @test ham_input_payload.summary.final_dimension == 221
    @test ham_input_payload.summary.pre_final_dimension == 221
    @test ham_input_payload.summary.density_gauge ==
          :pre_final_localized_positive_weight
    @test ham_input_payload.summary.raw_pair_factor_convention == :raw_numerator
    @test ham_input_payload.summary.support_row_order == :core_then_shell
    @test ham_input_payload.summary.source_plan_support_row_order ==
          :core_product_then_shell_left_right_pqs
    @test ham_input_payload.summary.support_weight_count == 275
    @test ham_input_payload.summary.support_raw_pair_shape == (275, 275)
    @test ham_input_payload.summary.pre_final_pair_matrix_shape == (221, 221)
    @test ham_input_payload.summary.ham_input_materialized
    @test !ham_input_payload.summary.h1_j_materialized
    @test !ham_input_payload.summary.ham_payload_materialized
    @test !ham_input_payload.summary.rhf_materialized
    @test !ham_input_payload.summary.public_api
    @test !ham_input_payload.summary.exports_materialized
    @test !ham_input_payload.summary.artifacts_materialized
    @test !ham_input_payload.summary.density_normalized_pair_terms_used_as_authority
    @test !ham_input_payload.summary.retained_diagnostic_weights_are_ida_weights
    @test :diatomic_complete_core_shell_ham_input_payload in
          ham_input_payload.available_objects
    @test :diatomic_complete_core_shell_pre_final_density_interaction in
          ham_input_payload.available_objects
    @test isempty(ham_input_payload.missing_objects)

    @test readiness.status == :blocked_diatomic_complete_core_shell_ham_readiness
    @test readiness.blocker == :missing_diatomic_ham_consumer_contract
    @test readiness.route_family === :pqs_source_box
    @test readiness.system_classification === :bond_aligned_diatomic
    @test readiness.bond_axis === :x
    @test readiness.parent_axis_bundle_object_available
    @test :parent_axis_bundle_object in readiness.available_objects
    @test :pqs_diatomic_complete_core_shell_source_plan in
          readiness.available_objects
    @test :diatomic_complete_core_shell_final_basis in readiness.available_objects
    @test :diatomic_complete_core_shell_h1_payload in readiness.available_objects
    @test :diatomic_complete_core_shell_ham_input_payload in
          readiness.available_objects
    @test !in(:parent_axis_bundle_object, readiness.missing_objects)
    @test !in(:diatomic_complete_core_shell_h1_consumer, readiness.missing_objects)
    @test !in(:diatomic_complete_core_shell_h1_j_consumer, readiness.missing_objects)
    @test :diatomic_ham_consumer_contract in readiness.missing_objects
    @test !in(
        :diatomic_complete_core_shell_final_basis_consumer,
        readiness.missing_objects,
    )
    @test readiness.summary.final_basis_materialized
    @test readiness.summary.h1_materialized
    @test readiness.summary.ham_input_materialized
    @test !readiness.summary.h1_j_materialized
    @test !readiness.summary.ham_payload_materialized
    @test !readiness.summary.rhf_materialized
    @test !readiness.summary.public_api
    @test !readiness.summary.exports_materialized
    @test !readiness.summary.artifacts_materialized

    @test readiness.source_box_summary.source_box_keys ==
          (:pqs_left, :product, :pqs_right)
    @test readiness.retained_unit_summary.unit_keys ==
          (:pqs_left, :pqs_right, :product)
    @test readiness.pair_inventory_summary.pair_family_counts ==
          (pqs_pqs = 3, pqs_product = 2, product_pqs = 0, product_product = 1)

    @test ham.object_kind == :pqs_source_box_complete_core_shell_ham_payload
    @test ham.status == :blocked_complete_core_shell_ham_payload
    @test ham.blocker == :missing_complete_core_shell_ham_inputs
    @test ham.missing_inputs == (
        :pqs_multilayer_complete_core_shell_final_basis,
        :pqs_multilayer_complete_core_shell_h1_payload,
        :pqs_complete_core_shell_final_one_electron_hamiltonian,
        :complete_core_shell_density_inputs,
        :complete_core_shell_h1_j_diagnostic_payload,
        :pqs_complete_core_shell_pre_final_density_interaction,
    )
    @test !ham.summary.public_api
    @test !ham.summary.exports_materialized
    @test !ham.summary.artifacts_materialized
    @test !ham.summary.rhf_product_surface
    @test !ham.summary.serious_hf_claim
end
