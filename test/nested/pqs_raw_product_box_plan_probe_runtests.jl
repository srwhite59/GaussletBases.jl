using Test
using GaussletBases

@testset "PQS explicit-core-spacing raw product-box plan probe" begin
    metrics_module = GaussletBases.CartesianContractedParentMetrics

    manual_counts = (x = 9, y = 7, z = 9)
    unavailable_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            q_to_core_spacing_rule = :explicit_core_spacing_only,
        )
    default_skeleton =
        metrics_module._pqs_pqs_product_source_box_route_skeleton(
            q = 5,
            parent_axis_counts = manual_counts,
        )
    default_probe =
        metrics_module._pqs_explicit_core_spacing_route_raw_product_box_plan_probe(
            unavailable_setup,
            default_skeleton,
        )
    @test default_probe.status == :not_constructed_pending_facts
    @test default_probe.raw_product_box_plan_count == 0
    @test :explicit_core_spacing in default_probe.pending_facts

    explicit_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            core_spacing = 0.15,
        )
    explicit_readiness =
        metrics_module._pqs_standard_parent_axis_construction_readiness(
            explicit_setup;
            parent_axis_counts = manual_counts,
        )
    parent_axis_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(explicit_setup)
    selected_counts =
        metrics_module._pqs_source_box_route_parent_axis_counts_for_skeleton(
            explicit_setup,
            explicit_readiness,
            parent_axis_probe;
            manual_parent_axis_counts = manual_counts,
        )
    constructed_skeleton =
        metrics_module._pqs_pqs_product_source_box_route_skeleton(
            q = 5,
            parent_axis_counts = selected_counts.parent_axis_counts,
        )
    raw_probe =
        metrics_module._pqs_explicit_core_spacing_route_raw_product_box_plan_probe(
            explicit_setup,
            constructed_skeleton,
        )

    @test raw_probe.status == :constructed_raw_product_box_plan_metadata
    @test raw_probe.raw_product_box_plan_count == 3
    @test raw_probe.gausslet_backend == :pgdg_localized_experimental
    @test raw_probe.gausslet_backend_role == :default_development_pgdg_localized
    @test raw_probe.diagnostics.gausslet_backend_role ==
          :default_development_pgdg_localized
    @test raw_probe.all_pgdg_exact
    @test !raw_probe.any_numerical_reference_fallback
    @test isfinite(raw_probe.max_axis_overlap_error)
    @test raw_probe.max_axis_overlap_error < 1.0e-8
    @test raw_probe.pending_facts == ()
    @test Tuple(metadata.unit_key for metadata in raw_probe.unit_plan_metadata) ==
          (:pqs_left, :pqs_right, :product)
    @test all(
        metadata -> metadata.integration_contract == :pgdg_exact,
        raw_probe.unit_plan_metadata,
    )
    @test all(
        metadata -> !metadata.numerical_reference_fallback,
        raw_probe.unit_plan_metadata,
    )
    @test all(
        metadata -> !metadata.retained_rule_attached,
        raw_probe.unit_plan_metadata,
    )
    @test all(metadata -> !metadata.packet_adoption, raw_probe.unit_plan_metadata)
    @test raw_probe.unit_plan_metadata[1].source_box ==
          constructed_skeleton.source_boxes.pqs_left
    @test raw_probe.unit_plan_metadata[2].source_box ==
          constructed_skeleton.source_boxes.pqs_right
    @test raw_probe.unit_plan_metadata[3].source_box ==
          constructed_skeleton.source_boxes.product
    @test raw_probe.unit_plan_metadata[1].source_mode_count == 125
    @test raw_probe.unit_plan_metadata[2].source_mode_count == 125
    @test raw_probe.unit_plan_metadata[3].source_mode_count == 25

    manual_route_probe =
        metrics_module._pqs_explicit_core_spacing_route_raw_product_box_plan_probe(
            explicit_setup,
            default_skeleton,
        )
    @test manual_route_probe.status == :not_constructed_pending_facts
    @test :route_skeleton_parent_axis_counts_from_constructed_probe in
          manual_route_probe.pending_facts

    for diagnostics in (
        default_probe.diagnostics,
        raw_probe.diagnostics,
        manual_route_probe.diagnostics,
    )
        @test diagnostics.private_development_only
        @test !diagnostics.production_route
        @test !diagnostics.public_default_consumes
        @test !diagnostics.packet_adoption
        @test !diagnostics.fixed_block_routing
        @test !diagnostics.qwhamiltonian_consumes
        @test !diagnostics.hamiltonian_matrix_built
        @test !diagnostics.shell_projection_used
        @test !diagnostics.lowdin_cleanup_used
        @test !diagnostics.support_local_shell_row_algorithm
        @test !diagnostics.support_coefficient_matrix_used
        @test !diagnostics.retained_pqs_weights_used
        @test !diagnostics.retained_weight_division_allowed
        @test !diagnostics.repo_side_ray_id
        @test !diagnostics.mwg_ida_semantics_changed
        @test !diagnostics.ecp_terms_implemented
        @test !diagnostics.cr2_science_status_changed
    end
end
