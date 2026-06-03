using Test
using GaussletBases

@testset "PQS route skeleton parent-axis count selection" begin
    metrics_module = GaussletBases.CartesianContractedParentMetrics

    manual_counts = (x = 9, y = 7, z = 9)
    default_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
        )
    default_readiness =
        metrics_module._pqs_standard_parent_axis_construction_readiness(
            default_setup;
            parent_axis_counts = manual_counts,
        )
    manual_selection =
        metrics_module._pqs_source_box_route_parent_axis_counts_for_skeleton(
            default_setup,
            default_readiness,
            nothing;
            manual_parent_axis_counts = manual_counts,
        )
    @test manual_selection.object_kind ==
          :pqs_source_box_route_parent_axis_counts_for_skeleton
    @test manual_selection.status == :available
    @test manual_selection.parent_axis_counts == manual_counts
    @test manual_selection.parent_axis_counts_source == :manual_fixture
    @test manual_selection.parent_axis_counts_manual_fixture
    @test !manual_selection.parent_axis_counts_derived
    @test manual_selection.parent_axis_probe_status == :not_requested
    @test manual_selection.q_minimum_satisfied
    @test manual_selection.pending_facts == ()

    manual_skeleton =
        metrics_module._pqs_pqs_product_source_box_route_skeleton(
            q = 5,
            parent_axis_counts = manual_selection.parent_axis_counts,
        )
    @test manual_skeleton.parent_axis_counts == manual_counts
    @test manual_skeleton.source_boxes.pqs_right == (x = 1:5, y = 1:5, z = 5:9)
    @test manual_skeleton.source_boxes.product == (x = 1:5, y = 1:5, z = 5:5)

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
    explicit_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(explicit_setup)
    constructed_selection =
        metrics_module._pqs_source_box_route_parent_axis_counts_for_skeleton(
            explicit_setup,
            explicit_readiness,
            explicit_probe;
            manual_parent_axis_counts = manual_counts,
        )
    @test constructed_selection.status == :available
    @test constructed_selection.parent_axis_counts == (x = 21, y = 11, z = 11)
    @test constructed_selection.parent_axis_counts_source ==
          :constructed_parent_axis_probe
    @test !constructed_selection.parent_axis_counts_manual_fixture
    @test constructed_selection.parent_axis_counts_derived
    @test constructed_selection.parent_axis_probe_status ==
          :constructed_explicit_core_spacing_parent_axis_metadata
    @test constructed_selection.q_minimum_satisfied
    @test constructed_selection.pending_facts == ()

    constructed_skeleton =
        metrics_module._pqs_pqs_product_source_box_route_skeleton(
            q = 5,
            parent_axis_counts = constructed_selection.parent_axis_counts,
        )
    @test constructed_skeleton.parent_axis_counts == (x = 21, y = 11, z = 11)
    @test constructed_skeleton.source_boxes.pqs_right ==
          (x = 1:5, y = 1:5, z = 7:11)
    @test constructed_skeleton.source_boxes.product == (x = 1:5, y = 1:5, z = 6:6)
    @test constructed_skeleton.retained_counts == manual_skeleton.retained_counts
    @test constructed_skeleton.ranges == manual_skeleton.ranges

    unavailable_selection =
        metrics_module._pqs_source_box_route_parent_axis_counts_for_skeleton(
            default_setup,
            default_readiness,
            nothing;
            manual_parent_axis_counts = nothing,
        )
    @test unavailable_selection.status == :not_available_pending_facts
    @test unavailable_selection.parent_axis_counts === nothing
    @test unavailable_selection.parent_axis_counts_source == :unavailable
    @test :manual_parent_axis_counts_or_constructed_parent_axis_probe in
          unavailable_selection.pending_facts

    for diagnostics in (
        manual_selection.diagnostics,
        constructed_selection.diagnostics,
        unavailable_selection.diagnostics,
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
