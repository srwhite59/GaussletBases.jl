using Test
using GaussletBases

@testset "PQS standard source-box route setup helper" begin
    metrics_module = GaussletBases.CartesianContractedParentMetrics

    setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
        )

    @test setup.object_kind == :pqs_standard_source_box_route_setup
    @test setup.status == :private_development_setup
    @test setup.nuclear_charges == (4.0, 4.0)
    @test setup.atom_locations == ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))
    @test setup.atom_count == 2
    @test setup.q == 5
    @test setup.n_s == 5
    @test setup.n_s_source == :q_default
    @test setup.core_cube_side == 5
    @test setup.core_cube_side_rule == :q_for_odd_q_q_plus_one_for_even_q
    @test setup.parent_box == (x = (-5.0, 5.0), y = (-3.0, 3.0), z = (-3.0, 3.0))
    @test setup.parent_box_lengths == (x = 10.0, y = 6.0, z = 6.0)
    @test setup.parent_box_rule ==
          :minimal_axis_aligned_box_enclosing_radius_pads_around_atoms
    @test setup.core_spacing == 0.15
    @test setup.d == 0.15
    @test setup.mapping_s ≈ sqrt(0.6)
    @test setup.mapping_s_by_atom == (sqrt(0.6), sqrt(0.6))
    @test setup.spacing.q_to_core_spacing_rule_status ==
          :standard_n_s_core_spacing_default
    @test setup.spacing.provenance ==
          :white_lindsey_shared_shell_policy_core_spacing_1p2_over_4_ns_minus_3
    @test setup.spacing.core_spacing_source == :standard_n_s_default
    @test setup.spacing.core_spacing_default_formula ==
          :core_spacing_equals_1p2_over_4_times_n_s_minus_3
    @test setup.spacing.white_lindsey_formula_available_when_d_is_explicit
    @test setup.diagnostics.private_development_only
    @test !setup.diagnostics.production_route
    @test setup.diagnostics.n_s_equals_q
    @test setup.diagnostics.n_s_source == :q_default
    @test setup.diagnostics.physical_parent_box_minimal_radius_pad
    @test setup.diagnostics.core_spacing_source == :standard_n_s_default
    @test setup.diagnostics.standard_n_s_default_core_spacing_used
    @test setup.diagnostics.q_to_core_spacing_non_optimality_claim == :not_claimed
    @test setup.diagnostics.q_to_core_spacing_replaceable
    @test !setup.diagnostics.explicit_core_spacing_override_used
    @test !setup.diagnostics.public_default_consumes
    @test !setup.diagnostics.packet_adoption
    @test !setup.diagnostics.fixed_block_routing
    @test !setup.diagnostics.qwhamiltonian_consumes
    @test !setup.diagnostics.hamiltonian_matrix_built
    @test !setup.diagnostics.shell_projection_used
    @test !setup.diagnostics.lowdin_cleanup_used
    @test !setup.diagnostics.support_local_shell_row_algorithm
    @test !setup.diagnostics.support_coefficient_matrix_used
    @test !setup.diagnostics.retained_pqs_weights_used
    @test !setup.diagnostics.retained_weight_division_allowed
    @test !setup.diagnostics.repo_side_ray_id
    @test !setup.diagnostics.mwg_ida_semantics_changed
    @test !setup.diagnostics.ecp_terms_implemented
    @test !setup.diagnostics.cr2_science_status_changed

    even_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 4,
            radius = 1.5,
    )
    @test even_setup.n_s == 4
    @test even_setup.core_spacing ≈ 0.3
    @test even_setup.core_cube_side == 5
    @test even_setup.parent_box ==
          (x = (-3.5, 3.5), y = (-1.5, 1.5), z = (-1.5, 1.5))

    ns_override =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            n_s = 6,
            radius = 3.0,
        )
    @test ns_override.q == 5
    @test ns_override.n_s == 6
    @test ns_override.n_s_source == :explicit_override
    @test !ns_override.diagnostics.n_s_equals_q
    @test ns_override.core_cube_side == 5
    @test ns_override.core_spacing ≈ 0.1
    @test ns_override.spacing.core_spacing_source == :standard_n_s_default

    explicit =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            core_spacing = 0.15,
        )
    @test explicit.core_spacing == 0.15
    @test explicit.d == 0.15
    @test explicit.mapping_s ≈ sqrt(0.6)
    @test explicit.mapping_s_by_atom[1] ≈ sqrt(0.6)
    @test explicit.mapping_s_by_atom[2] ≈ sqrt(0.6)
    @test explicit.spacing.core_range_by_atom[1] ≈ sqrt(0.15 / 4.0)
    @test explicit.spacing.core_range_by_atom[2] ≈ sqrt(0.15 / 4.0)
    @test explicit.spacing.q_to_core_spacing_rule_status ==
          :explicit_core_spacing_override
    @test explicit.spacing.provenance ==
          :explicit_core_spacing_with_white_lindsey_mapping_s_sqrt_dZ
    @test explicit.spacing.core_spacing_source == :explicit_core_spacing_override
    @test explicit.spacing.core_spacing_default_formula === nothing
    @test explicit.diagnostics.explicit_core_spacing_override_used
    @test !explicit.diagnostics.standard_n_s_default_core_spacing_used
    @test explicit.diagnostics.q_to_core_spacing_non_optimality_claim == :not_claimed

    explicit_only =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            q_to_core_spacing_rule = :explicit_core_spacing_only,
        )
    @test explicit_only.core_spacing === nothing
    @test explicit_only.spacing.q_to_core_spacing_rule_status ==
          :explicit_core_spacing_required
    @test explicit_only.spacing.core_spacing_source == :unavailable

    heteronuclear =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 3),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            core_spacing = 0.15,
        )
    @test heteronuclear.mapping_s === nothing
    @test heteronuclear.mapping_s_by_atom[1] ≈ sqrt(0.6)
    @test heteronuclear.mapping_s_by_atom[2] ≈ sqrt(0.45)

    @test_throws ArgumentError metrics_module._pqs_standard_source_box_route_setup(
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        q = 1,
        radius = 3.0,
    )
    @test_throws ArgumentError metrics_module._pqs_standard_source_box_route_setup(
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        q = 5,
        radius = 0.0,
    )
    @test_throws ArgumentError metrics_module._pqs_standard_source_box_route_setup(
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        q = 5,
        radius = 3.0,
        n_s = 3,
    )
    @test_throws ArgumentError metrics_module._pqs_standard_source_box_route_setup(
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        q = 5,
        radius = 3.0,
        q_to_core_spacing_rule = :not_a_rule,
    )
    @test_throws ArgumentError metrics_module._pqs_standard_source_box_route_setup(
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        q = 5,
        radius = 3.0,
        core_spacing = 0.0,
    )
    @test_throws DimensionMismatch metrics_module._pqs_standard_source_box_route_setup(
        nuclear_charges = (4,),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        q = 5,
        radius = 3.0,
    )
end
