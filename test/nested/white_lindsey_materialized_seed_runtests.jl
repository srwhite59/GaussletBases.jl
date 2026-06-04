using LinearAlgebra
using Test
using GaussletBases

function _white_lindsey_materialized_seed_fixture()
    parent_side_count = 7
    nside = 5
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = parent_side_count,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    sequence = build_one_center_atomic_full_parent_shell_sequence(
        basis;
        expansion = expansion,
        nside = nside,
        gausslet_backend = :numerical_reference,
        refinement_levels = 0,
    )
    fixed_block = GaussletBases._nested_fixed_block(sequence, basis, :numerical_reference)
    structure = one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = parent_side_count,
        nside = nside,
    )
    return (; parent_side_count, nside, basis, sequence, fixed_block, structure)
end

@testset "White-Lindsey materialized complete-shell seed fixture" begin
    fixture = _white_lindsey_materialized_seed_fixture()
    sequence = fixture.sequence
    fixed_block = fixture.fixed_block
    structure = fixture.structure

    shell_retained_count = 98
    core_retained_count = 5^3
    total_retained_dimension = core_retained_count + shell_retained_count

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
end
