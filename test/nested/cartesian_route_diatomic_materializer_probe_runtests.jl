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

@testset "route-configured diatomic materializer probe consumes parent handoff" begin
    expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false)
    report = _cartesian_route_diatomic_white_lindsey_report_for_test()
    @test report.route_materializer_payload.private_development_only
    @test report.route_materializer_payload.transient_only
    @test report.route_materializer_payload.parent_qw_basis_object_available
    @test report.route_materializer_payload.parent_axis_bundle_object_available
    @test report.route_materializer_payload.axis_bundle_backend ==
          :pgdg_localized_experimental

    materialization = GaussletBases._pqs_source_box_route_driver_materialization(
        report;
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
          :materialized_route_configured_bond_aligned_diatomic_shellization
    @test materialization.route_configured_diatomic_materializer_probe_materialized
    @test materialization.route_configured_diatomic_materializer_probe_consumed
    @test materialization.route_configured_diatomic_materializer_probe_blocker ===
          nothing
    @test !materialization.route_configured_diatomic_seed_fallback
    @test materialization.route_configured_diatomic_materializer_payload_available
    @test materialization.route_configured_diatomic_parent_qw_basis_object_handoff_available
    @test materialization.route_configured_diatomic_parent_axis_bundle_object_handoff_available
    @test materialization.route_configured_diatomic_axis_bundle_backend_handoff_available
    @test materialization.route_configured_diatomic_axis_bundle_backend_handoff ==
          :pgdg_localized_experimental
    @test materialization.route_configured_diatomic_shared_shell_layer_policy ==
          :endcap_panel_owned
    @test materialization.route_configured_diatomic_packet_kernel ==
          :factorized_direct
    @test materialization.route_configured_diatomic_policy_source ==
          :existing_endcap_panel_owned_pgdg_route

    missing = materialization.route_configured_diatomic_materializer_missing_contract
    @test isempty(missing)
    @test !in(:parent_qw_basis_object_handoff, missing)
    @test !in(:parent_axis_bundle_object_handoff, missing)
    @test !in(:axis_bundle_backend_provenance, missing)
    @test !in(:shared_shell_layer_policy, missing)
    @test !in(:packet_kernel, missing)
    @test !in(:pgdg_requires_endcap_panel_owned_shared_shell_policy, missing)
    @test !in(:coulomb_expansion_or_term_coefficients, missing)

    @test materialization.route_configured_materializer_backend_requested ==
          :pgdg_localized_experimental
    @test materialization.route_configured_materializer_backend_consumed ==
          :pgdg_localized_experimental
    @test materialization.route_configured_materializer_d_requested == 0.15
    @test materialization.route_configured_materializer_d_consumed == 0.15
    @test materialization.route_configured_materializer_nside_requested == 5
    @test materialization.route_configured_materializer_nside_consumed == 5

    probe_materialization =
        materialization.route_configured_diatomic_materializer_probe.materialization
    probe_options = probe_materialization.materializer_options
    @test probe_options.axis_bundle_backend == :pgdg_localized_experimental
    @test probe_options.shared_shell_layer_policy == :endcap_panel_owned
    @test probe_options.packet_kernel == :factorized_direct
    @test probe_options.term_coefficients_source == :coulomb_expansion_coefficients

    @test materialization.status ==
          :materialized_route_configured_diatomic_shellization_available
    @test materialization.materialized_report === nothing
    @test materialization.materialized_report_kind ==
          :cartesian_shellization_route_bond_aligned_diatomic_materialization
    @test materialization.shellization_source ==
          :route_configured_bond_aligned_diatomic_source
    @test materialization.route_configured_shellization_consumed
    @test materialization.seed_materialization_status ==
          :not_seed_route_configured_diatomic_shellization
    @test materialization.retained_dimension > 0
    @test materialization.basis_bundle_export_status ==
          :pending_route_configured_diatomic_basis_export
    @test materialization.basis_artifact_status == :not_requested
    @test !materialization.basis_artifact_written
    @test materialization.ham_artifact_status == :not_requested
    @test !materialization.ham_artifact_written
end
