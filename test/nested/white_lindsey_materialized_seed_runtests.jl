using LinearAlgebra
using Test
using GaussletBases

@testset "White-Lindsey materialized complete-shell seed fixture" begin
    fixture = GaussletBases._white_lindsey_low_order_materialized_seed_fixture()
    sequence = fixture.sequence
    fixed_block = fixture.fixed_block
    structure = fixture.structure
    inventory = fixture.inventory
    route_units = GaussletBases._white_lindsey_low_order_materialized_seed_route_units(fixture)
    operator_inventory =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_inventory(fixture)

    shell_retained_count = 98
    core_retained_count = 5^3
    total_retained_dimension = core_retained_count + shell_retained_count

    @test fixture.packet_kernel == :factorized_direct
    @test sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(sequence.shell_layers) == 1
    shell = only(sequence.shell_layers)
    @test shell isa GaussletBases._CartesianNestedCompleteShell3D

    @test length(shell.faces) == 6
    @test length(shell.edges) == 12
    @test length(shell.corners) == 8
    @test all(face -> face.side_first.retained_count == 3, shell.faces)
    @test all(face -> face.side_second.retained_count == 3, shell.faces)
    @test all(edge -> edge.side.retained_count == 3, shell.edges)

    @test length(shell.support_indices) == 7^3 - 5^3
    @test isempty(intersect(sequence.core_indices, shell.support_indices))
    @test sequence.working_box == (1:7, 1:7, 1:7)
    @test shell.provenance.source_box == (1:7, 1:7, 1:7)
    @test shell.provenance.next_inner_box == (2:6, 2:6, 2:6)
    @test shell.provenance.source_point_count == 7^3 - 5^3
    @test shell.provenance.retained_fixed_count == shell_retained_count

    @test sum(length, shell.face_column_ranges) == 54
    @test sum(length, shell.edge_column_ranges) == 36
    @test sum(length, shell.corner_column_ranges) == 8
    @test length(sequence.core_column_range) == core_retained_count
    @test length(only(sequence.layer_column_ranges)) == shell_retained_count
    @test first(sequence.core_column_range) == 1
    @test last(sequence.core_column_range) == core_retained_count
    @test only(sequence.layer_column_ranges) ==
          (core_retained_count + 1):total_retained_dimension

    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test fixed_block.shell === sequence
    @test size(fixed_block.coefficient_matrix, 2) == total_retained_dimension
    @test size(fixed_block.overlap) == (total_retained_dimension, total_retained_dimension)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    # These are retained-basis integral weights carried by the nested fixed block.
    @test length(fixed_block.weights) == total_retained_dimension
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0

    @test structure.parent_side_count == 7
    @test structure.working_box_side_count == 7
    @test structure.nside == 5
    @test structure.core_side_count == 5
    @test structure.shell_layer_count == 1
    @test structure.expected_shell_increment == shell_retained_count
    @test structure.expected_face_retained_count == 54
    @test structure.expected_edge_retained_count == 36
    @test structure.expected_corner_retained_count == 8
    @test structure.total_face_retained_count == 54
    @test structure.total_edge_retained_count == 36
    @test structure.total_corner_retained_count == 8
    @test structure.total_expected_gausslet_count == total_retained_dimension
    @test structure.total_actual_gausslet_count == total_retained_dimension
    @test structure.layers_match_expected

    @test inventory.object_kind == :white_lindsey_low_order_materialized_seed_inventory
    @test inventory.route_family == :white_lindsey_low_order
    @test inventory.status == :private_development_seed
    @test inventory.private_development_only
    @test inventory.packet_kernel == :factorized_direct
    @test inventory.parent_side_count == 7
    @test inventory.source_side_count == 7
    @test inventory.nside == 5
    @test inventory.piece_counts == (core = 1, faces = 6, edges = 12, corners = 8)
    @test inventory.support_counts == (core = 125, shell = 218, total_source = 343)
    @test inventory.retained_counts == (
        core = 125,
        faces = 54,
        edges = 36,
        corners = 8,
        shell = 98,
        total = 223,
    )
    @test inventory.retained_ranges.core == 1:125
    @test inventory.retained_ranges.shell == 126:223
    @test inventory.retained_ranges.faces ==
          Tuple((first(range) + 125):(last(range) + 125) for range in shell.face_column_ranges)
    @test inventory.retained_ranges.edges ==
          Tuple((first(range) + 125):(last(range) + 125) for range in shell.edge_column_ranges)
    @test inventory.retained_ranges.corners ==
          Tuple((first(range) + 125):(last(range) + 125) for range in shell.corner_column_ranges)
    @test inventory.materialized_shell_local_ranges.faces == Tuple(shell.face_column_ranges)
    @test inventory.materialized_shell_local_ranges.edges == Tuple(shell.edge_column_ranges)
    @test inventory.materialized_shell_local_ranges.corners == Tuple(shell.corner_column_ranges)
    @test inventory.fixed_block_ready
    @test inventory.overlap_ready
    @test inventory.retained_basis_integral_weights_ready
    @test inventory.weight_semantics == :retained_basis_integral_weights

    @test route_units.object_kind == :white_lindsey_low_order_materialized_seed_route_units
    @test route_units.route_family == :white_lindsey_low_order
    @test route_units.status == :private_development_seed
    @test route_units.packet_kernel == :factorized_direct
    @test route_units.retained_dimension == 223
    @test !route_units.operator_pairs_materialized
    @test isempty(route_units.pair_entries)
    @test route_units.pair_family_counts == (white_lindsey_low_order = 0,)
    @test route_units.weight_semantics == :retained_basis_integral_weights

    @test Tuple(unit.unit_key for unit in route_units.retained_units) == (
        :low_order_core_direct,
        :low_order_face_interiors,
        :low_order_edges,
        :low_order_corners,
    )
    @test Tuple(unit.retained_unit_kind for unit in route_units.retained_units) == (
        :white_lindsey_direct_core,
        :white_lindsey_face_interior_2d_products,
        :white_lindsey_edge_1d_side_functions,
        :white_lindsey_corner_direct_single_sites,
    )
    @test Tuple(unit.retained_rule_kind for unit in route_units.retained_units) == (
        :direct_parent_sites,
        :face_interior_2d_products_of_1d_retained_side_functions,
        :edge_1d_retained_side_functions,
        :corner_direct_single_site_pieces,
    )
    @test all(
        unit -> unit.weight_semantics == :retained_basis_integral_weights,
        route_units.retained_units,
    )
    @test route_units.unit_inventory.retained_counts == (
        low_order_core_direct = 125,
        low_order_face_interiors = 54,
        low_order_edges = 36,
        low_order_corners = 8,
    )
    @test route_units.unit_inventory.ranges == (
        low_order_core_direct = 1:125,
        low_order_face_interiors = 126:179,
        low_order_edges = 180:215,
        low_order_corners = 216:223,
    )
    @test route_units.standard_unit_inventory.retained_counts_materialized
    @test route_units.standard_unit_inventory.retained_ranges_materialized
    @test route_units.standard_unit_inventory.retained_dimension == 223
    @test route_units.standard_unit_inventory.pair_count == 0
    @test route_units.standard_unit_inventory.pair_family_counts ==
          (white_lindsey_low_order = 0,)

    @test operator_inventory.object_kind ==
          :white_lindsey_low_order_materialized_seed_operator_inventory
    @test operator_inventory.route_family == :white_lindsey_low_order
    @test operator_inventory.status == :private_development_seed
    @test operator_inventory.packet_kernel == :factorized_direct
    @test operator_inventory.operator_source == :nested_fixed_block
    @test operator_inventory.retained_dimension == 223
    @test operator_inventory.terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test all(==((223, 223)), values(operator_inventory.matrix_sizes))
    @test all(values(operator_inventory.finite_ready))
    @test operator_inventory.all_finite
    @test operator_inventory.overlap_identity_error < 1.0e-10
    @test operator_inventory.overlap_identity_ready
    @test operator_inventory.symmetric_ready.overlap
    @test operator_inventory.symmetric_ready.x2_x
    @test operator_inventory.symmetric_ready.x2_y
    @test operator_inventory.symmetric_ready.x2_z
    @test operator_inventory.symmetric_ready.kinetic
    @test !operator_inventory.operator_pairs_materialized
    @test !operator_inventory.electron_electron_materialized

    direct_inventory = GaussletBases._white_lindsey_low_order_materialized_seed_inventory(
        sequence,
        fixed_block,
        structure,
    )
    direct_route_units =
        GaussletBases._white_lindsey_low_order_materialized_seed_route_units(direct_inventory)
    direct_operator_inventory =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_inventory(fixed_block)
    @test direct_inventory.packet_kernel === nothing
    @test direct_route_units.packet_kernel === nothing
    @test direct_operator_inventory.packet_kernel === nothing

    support_fixture =
        GaussletBases._white_lindsey_low_order_materialized_seed_fixture(
            packet_kernel = :support_reference,
        )
    support_inventory = support_fixture.inventory
    support_route_units =
        GaussletBases._white_lindsey_low_order_materialized_seed_route_units(support_fixture)
    support_operator_inventory =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_inventory(
            support_fixture,
        )

    @test support_fixture.packet_kernel == :support_reference
    @test support_inventory.packet_kernel == :support_reference
    @test support_route_units.packet_kernel == :support_reference
    @test support_operator_inventory.packet_kernel == :support_reference
    @test support_inventory.retained_counts == inventory.retained_counts
    @test support_inventory.retained_ranges == inventory.retained_ranges
    @test support_inventory.retained_basis_integral_weights_ready
    @test support_route_units.retained_dimension == route_units.retained_dimension
    @test support_route_units.unit_inventory.retained_counts ==
          route_units.unit_inventory.retained_counts
    @test support_route_units.unit_inventory.ranges == route_units.unit_inventory.ranges
    @test Tuple(unit.unit_key for unit in support_route_units.retained_units) ==
          Tuple(unit.unit_key for unit in route_units.retained_units)
    @test Tuple(unit.retained_unit_kind for unit in support_route_units.retained_units) ==
          Tuple(unit.retained_unit_kind for unit in route_units.retained_units)
    @test support_operator_inventory.terms == operator_inventory.terms

    direct_matrices =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_matrices(fixed_block)
    support_matrices =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_matrices(
            support_fixture.fixed_block,
        )
    matrix_comparison_tolerance = 1.0e-10
    for term in operator_inventory.terms
        @test norm(
            getproperty(direct_matrices, term) - getproperty(support_matrices, term),
            Inf,
        ) <= matrix_comparison_tolerance
    end
end
