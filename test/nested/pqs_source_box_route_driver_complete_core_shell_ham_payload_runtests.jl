using Test
using LinearAlgebra
using GaussletBases

const _PQS_COMPLETE_CORE_SHELL_HAM_PAYLOAD_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

function _pqs_complete_core_shell_ham_payload_assembly()
    system_inputs = (;
        atom_symbols = ("Be",),
        nuclear_charges = (4,),
        atom_locations = ((0.0, 0.0, 0.0),),
        radius = 15.0,
        parent_axis_counts = (x = 7, y = 7, z = 7),
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
        route_kind = :one_center_pqs_source_box_complete_core_shell_ham_payload,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = _PQS_COMPLETE_CORE_SHELL_HAM_PAYLOAD_TERMS,
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

@testset "PQS source-box complete core/shell Ham payload" begin
    assembly = _pqs_complete_core_shell_ham_payload_assembly()
    route_payload = assembly.complete_core_shell_diagnostic_route_payload
    ham_payload = route_payload.complete_core_shell_ham_payload

    @test ham_payload.object_kind == :pqs_source_box_complete_core_shell_ham_payload
    @test ham_payload.status ==
          :materialized_pqs_source_box_complete_core_shell_ham_payload
    @test ham_payload.blocker === nothing
    @test ham_payload.route_family === :pqs_source_box
    @test ham_payload.dimension_summary.final_dimension == 223
    @test ham_payload.one_body_hamiltonian.status ==
          :materialized_pqs_complete_core_shell_final_one_electron_hamiltonian
    @test ham_payload.one_body_hamiltonian.hamiltonian_matrix_finite
    @test ham_payload.one_body_hamiltonian.hamiltonian_matrix_symmetry_error <= 1.0e-8
    @test ham_payload.density_interaction.status ==
          :materialized_pqs_complete_core_shell_pre_final_density_interaction
    @test ham_payload.electron_electron_representation ===
          :pre_final_density_interaction
    @test ham_payload.convention_labels.density_gauge ===
          :pre_final_localized_positive_weight
    @test ham_payload.convention_labels.raw_pair_factor_convention ===
          :raw_numerator
    @test ham_payload.convention_labels.support_row_order === :core_then_shell
    @test !ham_payload.convention_labels.signed_final_weight_division_used
    @test !ham_payload.convention_labels.raw_no_division_used
    @test !ham_payload.convention_labels.density_normalized_pair_terms_used_as_authority
    @test !ham_payload.convention_labels.public_api
    @test !ham_payload.convention_labels.exports_materialized
    @test !ham_payload.convention_labels.artifacts_materialized
    @test !ham_payload.convention_labels.rhf_product_surface
    @test !ham_payload.convention_labels.serious_hf_claim
end
