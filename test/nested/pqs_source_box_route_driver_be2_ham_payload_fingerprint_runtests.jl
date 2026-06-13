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

function _pqs_be2_ham_payload_fingerprint_assembly()
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
        probe_parent_axis_construction = false,
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
    return GaussletBases.cartesian_assembly(
        parent,
        shells,
        units,
        transforms,
        pairs,
        recipe,
    )
end

@testset "Be2 PQS Ham payload readiness fingerprint" begin
    assembly = _pqs_be2_ham_payload_fingerprint_assembly()
    @test assembly.object_kind == :cartesian_assembly
    @test assembly.route_skeleton.route_family === :pqs_source_box

    payload = assembly.complete_core_shell_diagnostic_route_payload
    readiness = assembly.diatomic_complete_core_shell_ham_readiness_payload
    ham = payload.complete_core_shell_ham_payload

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
