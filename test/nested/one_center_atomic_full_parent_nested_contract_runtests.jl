# Integration/slow test. Do not include in default nested runner.

@testset "One-center atomic full-parent nested contract" begin
    basis, sequence, audit = _one_center_atomic_full_parent_contract_fixture()
    count = length(basis)
    diagnostics = one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = count,
        nside = 7,
    )
    common_contract = GaussletBases._nested_glass_box_contract(diagnostics)
    count_only_27 = one_center_atomic_nested_structure_diagnostics(27; nside = 7)
    count_only_29 = one_center_atomic_nested_structure_diagnostics(29; nside = 7)
    count_only_5 = one_center_atomic_nested_structure_diagnostics(15; nside = 5)

    @test sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test sequence.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sequence.working_box == (1:count, 1:count, 1:count)
    @test audit.full_parent_working_box
    @test audit.support_count == count^3
    @test audit.expected_support_count == count^3
    @test audit.missing_row_count == 0
    @test audit.ownership_group_count_min == 1
    @test audit.ownership_group_count_max == 1
    @test audit.ownership_unowned_row_count == 0
    @test audit.ownership_multi_owned_row_count == 0

    @test diagnostics.parent_side_count == count
    @test diagnostics.working_box_side_count == count
    @test diagnostics.nside == 7
    @test diagnostics.core_side_count == 7
    @test diagnostics.shell_layer_count == 6
    @test diagnostics.expected_shell_increment == 7^3 - 5^3
    @test diagnostics.expected_shell_increment == 218
    @test diagnostics.expected_face_retained_count == 6 * 5^2
    @test diagnostics.expected_edge_retained_count == 12 * 5
    @test diagnostics.expected_corner_retained_count == 8
    @test diagnostics.total_face_retained_count == 6 * (6 * 5^2)
    @test diagnostics.total_edge_retained_count == 6 * (12 * 5)
    @test diagnostics.total_corner_retained_count == 6 * 8
    @test diagnostics.total_expected_gausslet_count == 7^3 + 6 * 218
    @test diagnostics.total_actual_gausslet_count == 7^3 + 6 * 218
    @test diagnostics.total_actual_gausslet_count == size(sequence.coefficient_matrix, 2)
    @test diagnostics.layers_match_expected
    @test length(diagnostics.layer_structures) == 6
    @test all(layer.face_retained_count == 150 for layer in diagnostics.layer_structures)
    @test all(layer.edge_retained_count == 60 for layer in diagnostics.layer_structures)
    @test all(layer.corner_retained_count == 8 for layer in diagnostics.layer_structures)
    @test all(layer.retained_dimension == 218 for layer in diagnostics.layer_structures)
    @test common_contract.fixed_dimension == diagnostics.total_actual_gausslet_count
    @test common_contract.contract_audit === diagnostics
    @test common_contract.layer_dimensions ==
          Int[layer.retained_dimension for layer in diagnostics.layer_structures]
    @test common_contract.layer_provenance ==
          GaussletBases._CartesianNestedShellLayerProvenance3D[
              shell.provenance for shell in sequence.shell_layers
          ]
    @test common_contract.leaf_count === nothing
    @test common_contract.layer_provenance[1].source_box == sequence.working_box
    @test common_contract.layer_provenance[1].source_point_count == 19^3 - 17^3
    @test common_contract.layer_provenance[end].next_inner_box == ((7:13), (7:13), (7:13))
    report_text = one_center_atomic_nested_structure_report(diagnostics)
    @test occursin("layer_1_source_box = ", report_text)
    @test occursin("layer_1_next_inner_box = ", report_text)
    @test occursin("layer_1_source_point_count = ", report_text)
    @test !occursin("leaf_count", report_text)

    @test count_only_5.expected_shell_increment == 5^3 - 3^3
    @test count_only_5.expected_shell_increment == 98
    @test count_only_5.layer_structures[1].provenance.source_point_count == 15^3 - 13^3

    @test count_only_27.parent_side_count == 27
    @test count_only_27.working_box_side_count == 27
    @test count_only_27.nside == 7
    @test count_only_27.core_side_count == 7
    @test count_only_27.shell_layer_count == 10
    @test count_only_27.expected_shell_increment == 218
    @test count_only_27.total_expected_gausslet_count == 343 + 10 * 218
    @test count_only_27.total_actual_gausslet_count == 343 + 10 * 218
    @test count_only_27.total_actual_gausslet_count == 2523
    @test count_only_27.layers_match_expected

    @test count_only_29.parent_side_count == 29
    @test count_only_29.working_box_side_count == 29
    @test count_only_29.shell_layer_count == 11
    @test count_only_29.expected_shell_increment == 218
    @test count_only_29.total_actual_gausslet_count == 343 + 11 * 218
end
