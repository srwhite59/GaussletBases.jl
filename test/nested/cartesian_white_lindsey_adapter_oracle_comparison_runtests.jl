using Test

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

function _lw_adapter_oracle_pair_family(adapter_block)
    metadata = adapter_block.metadata
    return Symbol(
        String(metadata.left_stratum_kind),
        "__",
        String(metadata.right_stratum_kind),
    )
end

function _lw_adapter_oracle_overlap_comparison(adapter_block, oracle_summary)
    expected_adapter_shape = (
        adapter_block.metadata.left_retained_column_count,
        adapter_block.metadata.right_retained_column_count,
    )
    adapter_shape = size(adapter_block.block)
    oracle_shape_available =
        hasproperty(oracle_summary.fixed_block_operator_matrix_sizes, :overlap)
    oracle_shape = oracle_shape_available ?
                   oracle_summary.fixed_block_operator_matrix_sizes.overlap :
                   nothing
    oracle_shape_status =
        oracle_summary.overlap_ready && oracle_shape_available ?
        :available_global_retained_overlap_shape :
        :blocked_missing_global_retained_overlap_shape

    shape_status =
        adapter_shape == expected_adapter_shape &&
        oracle_shape_status === :available_global_retained_overlap_shape ?
        :local_pair_shape_matches_adapter_metadata_global_oracle_shape_available :
        :blocked_shape_metadata_mismatch
    status =
        shape_status ===
        :local_pair_shape_matches_adapter_metadata_global_oracle_shape_available ?
        :metadata_shape_only :
        :blocked_oracle_shape_comparison
    blocker =
        status === :metadata_shape_only ? nothing : :shape_metadata_mismatch

    return (;
        object_kind = :white_lindsey_adapter_oracle_comparison,
        term = adapter_block.term,
        pair_key = adapter_block.pair_key,
        pair_family = _lw_adapter_oracle_pair_family(adapter_block),
        adapter_materialization_path = adapter_block.metadata.materialization_path,
        adapter_shape,
        expected_adapter_shape,
        oracle_shape,
        oracle_shape_status,
        shape_status,
        max_abs_error = nothing,
        max_abs_error_status =
            :not_compared_local_seed_pair_block_not_available,
        symmetry_error = nothing,
        symmetry_error_status =
            :not_applicable_rectangular_local_pair_block,
        status,
        blocker,
        oracle_role = oracle_summary.oracle_role,
        route_authority = false,
        adapter_authority = false,
        local_pair_block_materialized =
            adapter_block.metadata.local_pair_block_materialized,
        source_operator_blocks_materialized =
            adapter_block.metadata.source_operator_blocks_materialized,
        final_pair_blocks_materialized =
            adapter_block.metadata.final_pair_blocks_materialized,
        operator_blocks_materialized =
            adapter_block.metadata.operator_blocks_materialized,
        hamiltonian_data_materialized =
            adapter_block.metadata.hamiltonian_data_materialized,
        artifacts_materialized = adapter_block.metadata.artifacts_materialized,
        dense_parent_parent_overlap_materialized =
            adapter_block.metadata.dense_parent_parent_overlap_materialized,
    )
end

function _lw_adapter_oracle_comparison_summary(comparisons)
    comparison_tuple = Tuple(comparisons)
    return (;
        object_kind = :white_lindsey_adapter_oracle_comparison_summary,
        comparison_count = length(comparison_tuple),
        terms = Tuple(comparison.term for comparison in comparison_tuple),
        pair_families =
            Tuple(comparison.pair_family for comparison in comparison_tuple),
        statuses = Tuple(comparison.status for comparison in comparison_tuple),
        metadata_shape_only_count =
            count(comparison -> comparison.status === :metadata_shape_only,
                comparison_tuple),
        value_comparison_count =
            count(comparison -> comparison.max_abs_error !== nothing,
                comparison_tuple),
        blocked_count =
            count(comparison -> !isnothing(comparison.blocker),
                comparison_tuple),
    )
end

@testset "CartesianPairBlockMaterialization White-Lindsey oracle comparison scaffold" begin
    fixture = _lw_adapter_prepared_facet_edge_fixture(;
        prefix = "lw_oracle_comparison",
    )
    overlap_result =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_overlap_block(
            fixture.real_pair_coefficients;
            parent_axis_counts = (7, 7, 7),
            overlap_1d = fixture.factors.overlap_1d,
        )
    oracle_summary =
        CPBMForLWAdapter.white_lindsey_materialized_seed_oracle_summary()

    comparison =
        _lw_adapter_oracle_overlap_comparison(overlap_result, oracle_summary)
    @test comparison.object_kind ==
          :white_lindsey_adapter_oracle_comparison
    @test comparison.term == :overlap
    @test comparison.pair_key == (
        :lw_oracle_comparison_real_facet_unit,
        :lw_oracle_comparison_real_edge_unit,
    )
    @test comparison.pair_family == :facet_cpb__edge_cpb
    @test comparison.adapter_materialization_path ==
          :white_lindsey_boundary_stratum_overlap_adapter
    @test comparison.adapter_shape == (9, 3)
    @test comparison.expected_adapter_shape == (9, 3)
    @test comparison.oracle_shape == (223, 223)
    @test comparison.oracle_shape_status ==
          :available_global_retained_overlap_shape
    @test comparison.shape_status ==
          :local_pair_shape_matches_adapter_metadata_global_oracle_shape_available
    @test comparison.status == :metadata_shape_only
    @test isnothing(comparison.blocker)
    @test isnothing(comparison.max_abs_error)
    @test comparison.max_abs_error_status ==
          :not_compared_local_seed_pair_block_not_available
    @test isnothing(comparison.symmetry_error)
    @test comparison.symmetry_error_status ==
          :not_applicable_rectangular_local_pair_block
    @test comparison.oracle_role == :validation_oracle_only
    @test !comparison.route_authority
    @test !comparison.adapter_authority
    @test comparison.local_pair_block_materialized
    @test comparison.source_operator_blocks_materialized
    @test comparison.final_pair_blocks_materialized
    @test !comparison.operator_blocks_materialized
    @test !comparison.hamiltonian_data_materialized
    @test !comparison.artifacts_materialized
    @test !comparison.dense_parent_parent_overlap_materialized

    summary = _lw_adapter_oracle_comparison_summary((comparison,))
    @test summary.object_kind ==
          :white_lindsey_adapter_oracle_comparison_summary
    @test summary.comparison_count == 1
    @test summary.terms == (:overlap,)
    @test summary.pair_families == (:facet_cpb__edge_cpb,)
    @test summary.statuses == (:metadata_shape_only,)
    @test summary.metadata_shape_only_count == 1
    @test summary.value_comparison_count == 0
    @test summary.blocked_count == 0
end
