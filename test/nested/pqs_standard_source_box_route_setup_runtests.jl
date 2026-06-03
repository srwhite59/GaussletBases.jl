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
    @test setup.core_cube_side == 5
    @test setup.core_cube_side_rule == :q_for_odd_q_q_plus_one_for_even_q
    @test setup.parent_box == (x = (-5.0, 5.0), y = (-3.0, 3.0), z = (-3.0, 3.0))
    @test setup.parent_box_lengths == (x = 10.0, y = 6.0, z = 6.0)
    @test setup.parent_box_rule ==
          :minimal_axis_aligned_box_enclosing_radius_pads_around_atoms
    @test setup.core_spacing === nothing
    @test setup.d === nothing
    @test setup.mapping_s === nothing
    @test setup.mapping_s_by_atom === nothing
    @test setup.spacing.q_to_core_spacing_rule_status ==
          :unavailable_no_documented_q_to_core_spacing_formula
    @test setup.spacing.provenance ==
          :pending_documented_standard_rule_or_explicit_override
    @test setup.spacing.white_lindsey_formula_available_when_d_is_explicit
    @test setup.diagnostics.private_development_only
    @test !setup.diagnostics.production_route
    @test setup.diagnostics.n_s_equals_q
    @test setup.diagnostics.physical_parent_box_minimal_radius_pad
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
    @test even_setup.core_cube_side == 5
    @test even_setup.parent_box ==
          (x = (-3.5, 3.5), y = (-1.5, 1.5), z = (-1.5, 1.5))

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
    @test explicit.diagnostics.explicit_core_spacing_override_used
    @test explicit.diagnostics.q_to_core_spacing_non_optimality_claim == :not_claimed

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
        core_spacing = 0.0,
    )
    @test_throws DimensionMismatch metrics_module._pqs_standard_source_box_route_setup(
        nuclear_charges = (4,),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        q = 5,
        radius = 3.0,
    )
end
