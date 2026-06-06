@testset "One-center atomic legacy-profile nested contract" begin
    basis, sequence, diagnostics, ownership = _one_center_atomic_legacy_profile_contract_fixture()
    common_contract = GaussletBases._nested_glass_box_contract(diagnostics)
    range_groups = UnitRange{Int}[sequence.core_column_range]
    append!(range_groups, sequence.layer_column_ranges)
    support_group_counts = Int[]
    for row in sequence.support_indices
        nzcols = findall(!iszero, @view sequence.coefficient_matrix[row, :])
        touched_groups = 0
        for range in range_groups
            any(col -> col in range, nzcols) && (touched_groups += 1)
        end
        push!(support_group_counts, touched_groups)
    end

    @test sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test sequence.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sequence.working_box == (2:14, 2:14, 2:14)
    @test length(sequence.support_indices) == 13^3
    @test ownership.min_group_count == 0
    @test ownership.max_group_count == 1
    @test ownership.unowned_row_count == length(basis)^3 - 13^3
    @test ownership.multi_owned_row_count == 0
    @test minimum(support_group_counts) == 1
    @test maximum(support_group_counts) == 1
    @test diagnostics.parent_side_count == length(basis)
    @test diagnostics.working_box_side_count == 13
    @test diagnostics.nside == 5
    @test diagnostics.core_side_count == 5
    @test diagnostics.shell_layer_count == 4
    @test diagnostics.expected_shell_increment == 98
    @test diagnostics.total_actual_gausslet_count == 5^3 + 4 * 98
    @test diagnostics.layers_match_expected
    @test common_contract.fixed_dimension == diagnostics.total_actual_gausslet_count
    @test common_contract.contract_audit === diagnostics
    @test common_contract.layer_dimensions == [98, 98, 98, 98]
    @test common_contract.layer_provenance ==
          GaussletBases._CartesianNestedShellLayerProvenance3D[
              shell.provenance for shell in sequence.shell_layers
          ]
    @test common_contract.leaf_count === nothing
    @test common_contract.layer_provenance[1].source_box == ((2:14), (2:14), (2:14))
    @test common_contract.layer_provenance[end].next_inner_box == ((6:10), (6:10), (6:10))

    count_only_legacy_ne = one_center_atomic_nested_structure_diagnostics(
        29;
        working_box_side_count = 27,
        nside = 7,
    )
    count_only_modern_ne = one_center_atomic_nested_structure_diagnostics(29; nside = 7)
    @test count_only_legacy_ne.parent_side_count == 29
    @test count_only_legacy_ne.working_box_side_count == 27
    @test count_only_legacy_ne.shell_layer_count == 10
    @test count_only_legacy_ne.expected_shell_increment == 218
    @test count_only_legacy_ne.total_actual_gausslet_count == 2523
    @test count_only_modern_ne.working_box_side_count == 29
    @test count_only_modern_ne.shell_layer_count == 11
    @test count_only_modern_ne.total_actual_gausslet_count == 2741
end
