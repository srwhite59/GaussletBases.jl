using Test
using GaussletBases

const _PQS_ROUTE_DRIVER_ONE_CENTER_PROBE_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

function _pqs_route_driver_be2_probe_report(; route_family = :pqs_source_box)
    return GaussletBases._pqs_source_box_route_driver_dry_run(
        ;
        route_family,
        route_kind = :be2_cartesian_nesting_route_driver_spine,
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
        probe_parent_axis_construction = false,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
        probe_raw_product_box_plans = false,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = _PQS_ROUTE_DRIVER_ONE_CENTER_PROBE_TERMS,
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities = (:support_row_oracle, :dense_parent_projection),
    )
end

function _pqs_route_driver_one_center_probe_report()
    return (;
        route_family = :white_lindsey_low_order,
        retained_dimension = 223,
        system_metadata = (;
            atom_symbols = ("Be",),
            nuclear_charges = (4,),
            atom_locations = ((0.0, 0.0, 0.0),),
            parent_axis_counts = (x = 7, y = 7, z = 7),
            parent_axis_counts_source = :manual_fixture,
            parent_box = (x = -3.0:3.0, y = -3.0:3.0, z = -3.0:3.0),
            map_backend = :pgdg_localized_experimental,
        ),
        recipe_metadata = (;
            route_kind = :one_center_low_order_probe,
            route_shape = (:standard_cartesian_units, :low_order_comx_coarsening),
            benchmark_role = :published_cartesian_baseline_for_pqs_comparison,
            n_s = 5,
            core_spacing = 0.15,
            reference_spacing = 1.0,
            tail_spacing = 10.0,
            parent_axis_probe_backend = :pgdg_localized_experimental,
        ),
    )
end

@testset "Route-driver one-center materializer probe" begin
    density_expansion = coulomb_gaussian_expansion(doacc = false)

    default_status = GaussletBases._pqs_source_box_route_driver_materialization(
        _pqs_route_driver_be2_probe_report();
        materialize_route = false,
        save_basis_artifact = false,
        save_ham_artifact = false,
        basisfile = "unused_basis.jld2",
        hamfile = "unused_ham.jld2",
    )
    @test default_status.route_configured_one_center_materializer_probe_status ==
          :not_requested
    @test !default_status.route_configured_one_center_materializer_probe_requested
    @test !default_status.route_configured_one_center_materializer_probe_materialized
    @test !default_status.route_configured_one_center_materializer_probe_consumed
    @test default_status.route_configured_one_center_materializer_probe_blocker === nothing

    be2_white_lindsey_status =
        GaussletBases._pqs_source_box_route_driver_materialization(
            _pqs_route_driver_be2_probe_report(route_family = :white_lindsey_low_order);
            materialize_route = false,
            probe_route_configured_one_center_materializer = true,
            save_basis_artifact = false,
            save_ham_artifact = false,
            basisfile = "unused_basis.jld2",
            hamfile = "unused_ham.jld2",
            white_lindsey_expansion = density_expansion,
        )
    @test be2_white_lindsey_status.route_configured_system_classification ==
          :bond_aligned_diatomic
    @test be2_white_lindsey_status.route_configured_one_center_materializer_probe_status ==
          :blocked_not_one_center
    @test be2_white_lindsey_status.route_configured_one_center_materializer_probe_blocker ==
          :route_config_not_one_center
    @test !be2_white_lindsey_status.route_configured_one_center_materializer_probe_materialized
    @test !be2_white_lindsey_status.route_configured_one_center_materializer_probe_consumed

    one_center_status = GaussletBases._pqs_source_box_route_driver_materialization(
        _pqs_route_driver_one_center_probe_report();
        materialize_route = false,
        probe_route_configured_one_center_materializer = true,
        save_basis_artifact = false,
        save_ham_artifact = false,
        basisfile = "unused_basis.jld2",
        hamfile = "unused_ham.jld2",
        white_lindsey_expansion = density_expansion,
    )
    @test one_center_status.route_configured_system_classification == :one_center
    @test one_center_status.route_configured_one_center_materializer_probe_status ==
          :materialized_route_configured_one_center_low_order
    @test one_center_status.route_configured_one_center_materializer_probe_materialized
    @test one_center_status.route_configured_one_center_materializer_probe_consumed
    @test one_center_status.route_configured_one_center_materializer_probe_blocker === nothing

    materialization =
        one_center_status.route_configured_one_center_materializer_probe.materialization
    @test materialization.object_kind ==
          :cartesian_shellization_route_one_center_materialization
    @test materialization.retained_dimension == 223
    @test materialization.route_configured_shellization_consumed
    @test materialization.shellification_plan_path_used
    @test materialization.calls_shellification_plan_materializer
    @test !materialization.calls_build_one_center_atomic_full_parent_shell_sequence
    @test !materialization.calls_white_lindsey_seed_fixture
    @test materialization.materializer_options.gausslet_backend ==
          :pgdg_localized_experimental
    @test materialization.materializer_options.d == 0.15
    @test materialization.materializer_options.nside == 5
    plan_summary = materialization.shellification_plan_summary
    @test plan_summary.object_kind == :cartesian_shellification_plan_private_summary
    @test plan_summary.source_kind == :one_center_atomic_full_parent_shellification_plan
    @test plan_summary.route_family == :white_lindsey_low_order
    @test plan_summary.system_classification == :one_center
    @test plan_summary.region_count == 2
    @test plan_summary.shell_region_count == 1
    @test plan_summary.direct_core_region_count == 1
    @test plan_summary.retained_dimension == materialization.retained_dimension
    @test plan_summary.coverage_complete
    @test one_center_status.route_configured_materializer_backend_requested ==
          :pgdg_localized_experimental
    @test one_center_status.route_configured_materializer_backend_consumed ==
          :pgdg_localized_experimental
    @test one_center_status.route_configured_materializer_d_requested == 0.15
    @test one_center_status.route_configured_materializer_d_consumed == 0.15
    @test one_center_status.route_configured_materializer_nside_requested == 5
    @test one_center_status.route_configured_materializer_nside_consumed == 5
    @test !materialization.public_default_behavior_changed
end
