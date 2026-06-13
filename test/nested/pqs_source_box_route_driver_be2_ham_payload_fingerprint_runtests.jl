using Test
using GaussletBases

const _PQS_BE2_CCPM = GaussletBases.CartesianContractedParentMetrics

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

function _pqs_be2_axis_metrics(bundles)
    pgdg_x = GaussletBases._nested_axis_pgdg(bundles, :x)
    pgdg_y = GaussletBases._nested_axis_pgdg(bundles, :y)
    pgdg_z = GaussletBases._nested_axis_pgdg(bundles, :z)
    return (
        x = (
            overlap = pgdg_x.overlap,
            position = pgdg_x.position,
            x2 = pgdg_x.x2,
            weights = pgdg_x.weights,
            centers = pgdg_x.centers,
            kinetic = pgdg_x.kinetic,
            source = :nested_pgdg_axis,
        ),
        y = (
            overlap = pgdg_y.overlap,
            position = pgdg_y.position,
            x2 = pgdg_y.x2,
            weights = pgdg_y.weights,
            centers = pgdg_y.centers,
            kinetic = pgdg_y.kinetic,
            source = :nested_pgdg_axis,
        ),
        z = (
            overlap = pgdg_z.overlap,
            position = pgdg_z.position,
            x2 = pgdg_z.x2,
            weights = pgdg_z.weights,
            centers = pgdg_z.centers,
            kinetic = pgdg_z.kinetic,
            source = :nested_pgdg_axis,
        ),
    )
end

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
    source_plan_payload =
        assembly.diatomic_complete_core_shell_source_plan_payload
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
    @test source_plan_payload.summary.support_window_payload_status ==
          :available_diatomic_complete_core_shell_support_windows
    @test :parent_axis_bundle_object in source_plan_payload.missing_objects
    @test :diatomic_complete_core_shell_source_realization_contract in
          source_plan_payload.missing_objects
    @test !source_plan_payload.summary.source_plan_materialized
    @test !source_plan_payload.summary.final_basis_materialized
    @test !source_plan_payload.summary.h1_materialized
    @test !source_plan_payload.summary.h1_j_materialized
    @test !source_plan_payload.summary.ham_payload_materialized
    @test !source_plan_payload.summary.route_driver_public_surface
    @test !source_plan_payload.summary.exports_materialized
    @test !source_plan_payload.summary.artifacts_materialized

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
    source_plan_payload =
        assembly.diatomic_complete_core_shell_source_plan_payload
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

    @test source_plan_payload.status ==
          :blocked_diatomic_complete_core_shell_source_plan
    @test source_plan_payload.blocker ==
          :missing_diatomic_complete_core_shell_source_realization_contract
    @test source_plan_payload.route_family === :pqs_source_box
    @test source_plan_payload.system_classification === :bond_aligned_diatomic
    @test source_plan_payload.bond_axis === :x
    @test source_plan_payload.parent_axis_bundle_object_available
    @test source_plan_payload.source_plan === nothing
    @test source_plan_payload.source_plan_status ==
          :not_materialized_diatomic_complete_core_shell_source_plan
    @test source_plan_payload.support_window_payload === support_window_payload
    @test source_plan_payload.summary.support_window_payload_status ==
          :available_diatomic_complete_core_shell_support_windows
    @test :diatomic_complete_core_shell_support_windows in
          source_plan_payload.available_objects
    @test :parent_axis_bundle_object in source_plan_payload.available_objects
    @test !in(:parent_axis_bundle_object, source_plan_payload.missing_objects)
    @test :diatomic_complete_core_shell_source_realization_contract in
          source_plan_payload.missing_objects
    @test !source_plan_payload.summary.source_plan_materialized
    @test !source_plan_payload.summary.final_basis_materialized
    @test !source_plan_payload.summary.h1_materialized
    @test !source_plan_payload.summary.h1_j_materialized
    @test !source_plan_payload.summary.ham_payload_materialized
    @test !source_plan_payload.summary.route_driver_public_surface
    @test !source_plan_payload.summary.exports_materialized
    @test !source_plan_payload.summary.artifacts_materialized

    @test readiness.status == :blocked_diatomic_complete_core_shell_ham_readiness
    @test readiness.blocker ==
          :missing_diatomic_complete_core_shell_source_plan_producer
    @test readiness.route_family === :pqs_source_box
    @test readiness.system_classification === :bond_aligned_diatomic
    @test readiness.bond_axis === :x
    @test readiness.parent_axis_bundle_object_available
    @test :parent_axis_bundle_object in readiness.available_objects
    @test !in(:parent_axis_bundle_object, readiness.missing_objects)
    @test :diatomic_complete_core_shell_source_plan_producer in
          readiness.missing_objects

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

@testset "Be2 PQS raw-box producer fingerprint" begin
    components = _pqs_be2_ham_payload_fingerprint_components(
        probe_parent_axis_construction = :auto,
    )
    parent = components.parent
    assembly = components.assembly
    support_window_payload =
        assembly.diatomic_complete_core_shell_support_window_payload

    @test support_window_payload.status ==
          :available_diatomic_complete_core_shell_support_windows
    @test parent.parent_axis_bundle_object_available
    bundles = parent.parent_axis_bundle_object
    metrics = _pqs_be2_axis_metrics(bundles)

    producer = _PQS_BE2_CCPM._pqs_pqs_product_raw_box_route_producer(
        bundles,
        support_window_payload.source_box_windows.pqs_left,
        support_window_payload.source_box_windows.pqs_right,
        support_window_payload.source_box_windows.product,
        metrics;
        source_mode_dims = support_window_payload.source_mode_dims.pqs_left,
        route_name = :be2_pqs_raw_box_route_producer_fingerprint,
        parent_dims = support_window_payload.parent_dims,
        bond_axis = parent.bond_axis,
        metadata = (source = :be2_pqs_raw_box_route_producer_fingerprint,),
        provenance = (source = :be2_pqs_raw_box_route_producer_fingerprint,),
    )

    @test producer.object_kind == :pqs_pqs_product_raw_box_route_producer
    @test producer.status == :private_shadow_only
    @test producer.descriptor.object_kind ==
          :pqs_pqs_product_safe_term_route_descriptor
    @test producer.descriptor.roles == (:pqs_left, :pqs_right, :product)
    @test producer.descriptor.retained_dimension == 221
    @test producer.descriptor.expected_pair_count == 6

    left_raw = producer.raw_product_box_plans.pqs_left
    right_raw = producer.raw_product_box_plans.pqs_right
    @test left_raw.object_kind == :cartesian_raw_product_box_plan_3d
    @test right_raw.object_kind == :cartesian_raw_product_box_plan_3d
    @test left_raw.source_mode_dims == (5, 5, 5)
    @test right_raw.source_mode_dims == (5, 5, 5)
    @test length(left_raw.axis_local_coefficients) == 3
    @test length(right_raw.axis_local_coefficients) == 3
    @test all(axis -> size(left_raw.axis_local_coefficients[axis]) == (5, 5), 1:3)
    @test all(axis -> size(right_raw.axis_local_coefficients[axis]) == (5, 5), 1:3)

    left_pqs = _PQS_BE2_CCPM._pqs_raw_product_box_plan_view(
        producer.raw_pqs_plans.pqs_left,
    )
    right_pqs = _PQS_BE2_CCPM._pqs_raw_product_box_plan_view(
        producer.raw_pqs_plans.pqs_right,
    )
    @test left_pqs.representation == :orthogonal_raw_product_box
    @test right_pqs.representation == :orthogonal_raw_product_box
    @test left_pqs.boundary_selector.selected_count == 98
    @test right_pqs.boundary_selector.selected_count == 98
    @test length(left_pqs.axis_local_coefficients) == 3
    @test length(right_pqs.axis_local_coefficients) == 3

    product_unit = producer.product_unit
    @test product_unit.kind == :product_doside
    @test length(product_unit.support_indices) == 25
    @test length(product_unit.support_states) == 25
    @test size(product_unit.coefficient_matrix) == (25, 25)

    inventory_diagnostics = producer.all_pairs_inventory.diagnostics
    producer_diagnostics = producer.diagnostics
    @test inventory_diagnostics.every_pair_uses_source_box_algorithmic_policy
    @test inventory_diagnostics.source_box_algorithmic_pair_count == 6
    @test producer_diagnostics.private_shadow_only
    @test producer_diagnostics.raw_product_box_plan_built
    @test producer_diagnostics.retained_rule_built
    @test producer_diagnostics.route_descriptor_emitted
    @test !producer_diagnostics.dense_raw_source_box_pair_matrix_materialized
    dense_raw_pair_matrix_built_by_producer =
        producer_diagnostics.dense_raw_source_box_pair_matrix_materialized_by_producer
    @test !dense_raw_pair_matrix_built_by_producer
    @test producer_diagnostics.dense_raw_source_box_pair_matrices_validation_only
    @test !producer_diagnostics.shell_projection_used
    @test !producer_diagnostics.lowdin_cleanup_used
    @test !producer_diagnostics.support_local_pqs_oracle_used
    @test !producer_diagnostics.support_coefficient_matrix_used
    @test !producer_diagnostics.retained_pqs_weights_used
    @test producer_diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !producer_diagnostics.ida_weight_division_allowed
    @test !producer_diagnostics.packet_adoption
    @test !producer_diagnostics.fixed_block_routing
    @test !producer_diagnostics.qwhamiltonian_consumes
    @test !producer_diagnostics.public_default_consumes
end
