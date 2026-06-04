using Test
using GaussletBases

function _cartesian_route_diatomic_white_lindsey_report_for_test()
    return GaussletBases._pqs_source_box_route_driver_dry_run(
        route_family = :white_lindsey_low_order,
        route_kind = :be2_cartesian_diatomic_materializer_probe,
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0,
        parent_axis_counts = (x = 9, y = 7, z = 9),
        map_backend = :pgdg_localized_experimental,
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = nothing,
        probe_parent_axis_construction = :auto,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
        probe_raw_product_box_plans = :auto,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (:overlap, :position_x, :position_y, :position_z,
            :x2_x, :x2_y, :x2_z, :kinetic),
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities = (:support_row_oracle, :dense_parent_projection),
    )
end

@testset "route-configured diatomic materializer probe blocker" begin
    expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false)
    materialization = GaussletBases._pqs_source_box_route_driver_materialization(
        _cartesian_route_diatomic_white_lindsey_report_for_test();
        materialize_route = true,
        save_basis_artifact = false,
        save_ham_artifact = false,
        basisfile = "unused_basis.jld2",
        hamfile = "unused_ham.jld2",
        white_lindsey_expansion = expansion,
    )

    @test materialization.route_configured_system_classification ==
          :bond_aligned_diatomic
    @test materialization.route_configured_diatomic_materializer_probe_requested
    @test materialization.route_configured_diatomic_materializer_probe_status ==
          :blocked_missing_diatomic_materializer_contract
    @test !materialization.route_configured_diatomic_materializer_probe_materialized
    @test !materialization.route_configured_diatomic_materializer_probe_consumed
    @test materialization.route_configured_diatomic_materializer_probe_blocker ==
          :pending_route_configured_bond_aligned_diatomic_materializer_contract
    @test materialization.route_configured_diatomic_seed_fallback

    missing = materialization.route_configured_diatomic_materializer_missing_contract
    @test :parent_qw_basis_object_handoff in missing
    @test :parent_axis_bundle_object_handoff in missing
    @test :axis_bundle_backend_provenance in missing
    @test :shared_shell_layer_policy in missing
    @test :packet_kernel in missing
    @test :pgdg_requires_endcap_panel_owned_shared_shell_policy in missing
    @test !in(:coulomb_expansion_or_term_coefficients, missing)

    @test materialization.route_configured_materializer_backend_requested ==
          :pgdg_localized_experimental
    @test materialization.route_configured_materializer_backend_consumed === nothing
    @test materialization.route_configured_materializer_d_requested == 0.15
    @test materialization.route_configured_materializer_d_consumed === nothing
    @test materialization.route_configured_materializer_nside_requested == 5
    @test materialization.route_configured_materializer_nside_consumed === nothing

    @test materialization.status == :materialized_seed_report_available
    @test materialization.shellization_source == :white_lindsey_one_center_seed
    @test !materialization.route_configured_shellization_consumed
end
