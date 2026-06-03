using Test
using GaussletBases

@testset "PQS standard parent-axis construction readiness helper" begin
    metrics_module = GaussletBases.CartesianContractedParentMetrics

    default_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
        )
    default_readiness =
        metrics_module._pqs_standard_parent_axis_construction_readiness(default_setup)

    @test default_readiness.object_kind ==
          :pqs_standard_parent_axis_construction_readiness
    @test default_readiness.status == :not_ready_pending_facts
    @test !default_readiness.core_spacing_available
    @test !default_readiness.d_available
    @test !default_readiness.white_lindsey_spacing_facts_available
    @test default_readiness.charge_family == :homonuclear
    @test default_readiness.homonuclear
    @test !default_readiness.heteronuclear
    @test default_readiness.geometry.status == :axis_aligned_diatomic
    @test default_readiness.geometry.bond_axis == :x
    @test default_readiness.geometry.bond_length == 4.0
    @test default_readiness.geometry.existing_bond_aligned_api_geometry_ready
    @test default_readiness.extent_candidates.available
    @test default_readiness.extent_candidates.xmax_parallel == 5.0
    @test default_readiness.extent_candidates.xmax_transverse == 3.0
    @test default_readiness.parent_axis_counts === nothing
    @test default_readiness.parent_axis_counts_status ==
          :pending_helper_or_documented_rule
    @test !default_readiness.parent_axis_counts_manual_fixture
    @test !default_readiness.parent_axis_counts_derived
    @test !default_readiness.existing_parent_api_appears_applicable
    @test !default_readiness.standard_parent_axis_rule_ready
    @test !default_readiness.parent_axis_construction_ready
    @test !default_readiness.parent_axis_metadata_constructed
    @test :explicit_core_spacing_or_documented_q_to_core_spacing_rule in
          default_readiness.pending_facts
    @test :parent_axis_counts_or_documented_axis_count_rule in
          default_readiness.pending_facts
    @test :reviewed_parent_axis_constructor_call in default_readiness.pending_facts

    explicit_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            core_spacing = 0.15,
        )
    manual_readiness =
        metrics_module._pqs_standard_parent_axis_construction_readiness(
            explicit_setup;
            parent_axis_counts = (x = 9, y = 7, z = 9),
        )

    @test manual_readiness.core_spacing_available
    @test manual_readiness.d_available
    @test manual_readiness.white_lindsey_spacing_facts_available
    @test manual_readiness.core_spacing == 0.15
    @test manual_readiness.d == 0.15
    @test manual_readiness.mapping_s ≈ sqrt(0.6)
    @test manual_readiness.mapping_s_by_atom[1] ≈ sqrt(0.6)
    @test manual_readiness.mapping_s_by_atom[2] ≈ sqrt(0.6)
    @test manual_readiness.parent_axis_counts == (x = 9, y = 7, z = 9)
    @test manual_readiness.parent_axis_counts_status == :manual_fixture
    @test manual_readiness.parent_axis_counts_manual_fixture
    @test !manual_readiness.parent_axis_counts_derived
    @test manual_readiness.parent_axis_counts_derivation ==
          :manual_driver_fixture_not_standard_setup_derivation
    @test manual_readiness.existing_parent_api_appears_applicable
    homonuclear_candidate =
        manual_readiness.existing_parent_api_candidates.bond_aligned_homonuclear_qw_basis
    @test homonuclear_candidate.appears_applicable
    @test !manual_readiness.standard_parent_axis_rule_ready
    @test !manual_readiness.parent_axis_construction_ready
    @test !manual_readiness.parent_axis_metadata_constructed
    @test manual_readiness.construction_decision ==
          :readiness_only_no_parent_axis_construction_added
    @test :standard_parent_axis_count_rule_replacing_manual_fixture in
          manual_readiness.pending_facts
    @test :reviewed_parent_axis_constructor_call in manual_readiness.pending_facts

    heteronuclear_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 3),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            core_spacing = 0.15,
        )
    heteronuclear_readiness =
        metrics_module._pqs_standard_parent_axis_construction_readiness(
            heteronuclear_setup;
            parent_axis_counts = (9, 7, 9),
        )
    @test heteronuclear_readiness.charge_family == :heteronuclear
    @test !heteronuclear_readiness.homonuclear
    @test heteronuclear_readiness.heteronuclear
    @test heteronuclear_readiness.existing_parent_api_appears_applicable
    heteronuclear_candidate =
        heteronuclear_readiness.existing_parent_api_candidates.bond_aligned_heteronuclear_qw_basis
    @test heteronuclear_candidate.appears_applicable
    @test :atom_symbol_labels_for_heteronuclear_parent_api in
          heteronuclear_readiness.pending_facts
    @test !heteronuclear_readiness.standard_parent_axis_rule_ready

    for diagnostics in (
        default_readiness.diagnostics,
        manual_readiness.diagnostics,
        heteronuclear_readiness.diagnostics,
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
