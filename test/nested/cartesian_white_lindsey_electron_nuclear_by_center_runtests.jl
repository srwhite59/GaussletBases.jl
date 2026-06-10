using Test
using GaussletBases

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const WLENCPBM = GaussletBases.CartesianPairBlockMaterialization

function _wlen_collection_entry(result, left_column_range, right_column_range)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection_entry,
        entry_kind = :materialized_result,
        term = result.term,
        block_set_term = :electron_nuclear_by_center,
        result_term = result.term,
        source_space_term = nothing,
        pair_key = result.pair_key,
        pair_index = 1,
        selector_family = :white_lindsey_boundary_stratum,
        materialization_path = result.metadata.materialization_path,
        block_space = :final_local_space,
        block_shape = size(result.block),
        block_size = length(result.block),
        status = :materialized_local_one_body_block_collection_entry,
        blocker = nothing,
        left_column_range,
        right_column_range,
        center_index = result.metadata.center_index,
        center_key = result.metadata.center_key,
        center_location = result.metadata.center_location,
        nuclear_charge_recorded = result.metadata.nuclear_charge_recorded,
        nuclear_charge_applied = result.metadata.nuclear_charge_applied,
        materialized = result.materialized,
        result_available = true,
        skipped_record_available = false,
        result,
        skipped_record = nothing,
    )
end

function _wlen_collection(entry)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        status = :materialized_local_one_body_block_collection,
        blocker = nothing,
        terms = (:electron_nuclear_by_center,),
        requested_terms = (:electron_nuclear_by_center,),
        requested_materialize_terms = (:electron_nuclear_by_center,),
        materialized_terms = (:electron_nuclear_by_center,),
        deferred_terms = (),
        entries = (entry,),
        materialized_entries = (entry,),
        skipped_entries = (),
        entry_count = 1,
        materialized_entry_count = 1,
        skipped_entry_count = 0,
    )
end

@testset "White-Lindsey decomposed electron-nuclear by-center one-body" begin
    fixture = _lw_adapter_prepared_facet_edge_fixture(
        prefix = "lw_electron_nuclear_test",
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    parent_axis_bundle_object = (;
        x = fixture.doside_source_1d,
        y = fixture.doside_source_1d,
        z = fixture.doside_source_1d,
    )
    center_record = (;
        center_key = :proton_a,
        center_index = 1,
        nuclear_charge = 1.0,
        location = (0.0, 0.0, 0.0),
    )

    block = WLENCPBM.white_lindsey_boundary_stratum_one_body_block(
        fixture.real_pair_coefficients,
        :electron_nuclear_by_center;
        parent_axis_counts = (7, 7, 7),
        parent_axis_bundle_object,
        coulomb_expansion = expansion,
        center_record,
    )

    @test block.term == :electron_nuclear_by_center
    @test block.materialized
    @test block.source_operator_blocks_materialized
    @test block.final_pair_blocks_materialized
    @test !block.operator_blocks_materialized
    @test !block.hamiltonian_data_materialized
    @test block.metadata.materialization_path ===
          :white_lindsey_boundary_stratum_electron_nuclear_by_center_adapter
    @test block.metadata.by_center
    @test !block.metadata.centers_summed
    @test block.metadata.nuclear_charge_recorded
    @test !block.metadata.nuclear_charge_applied
    @test block.metadata.charge_application_stage ===
          :acceptance_or_hamiltonian_assembly
    @test !block.metadata.full_parent_window_cpb_used
    @test block.metadata.gaussian_expansion_loop === :inner_support_contraction
    @test size(block.block) == (
        fixture.real_pair_coefficients.left_retained_column_count,
        fixture.real_pair_coefficients.right_retained_column_count,
    )
    @test all(isfinite, block.block)

    collection = _wlen_collection(
        _wlen_collection_entry(
            block,
            1:size(block.block, 1),
            (size(block.block, 1) + 1):(size(block.block, 1) + size(block.block, 2)),
        ),
    )
    plan = WLENCPBM.one_body_electron_nuclear_by_center_placement_plan(
        collection;
        global_dimension = sum(size(block.block)),
    )
    @test plan.status == :placeable_local_one_body_placement_plan
    @test plan.record_count == 1
    @test plan.placeable_count == 1
    @test plan.placeable_records[1].center_index == 1
    @test plan.placeable_records[1].nuclear_charge_recorded
    @test !plan.placeable_records[1].nuclear_charge_applied

    global_result =
        WLENCPBM.one_body_global_electron_nuclear_by_center_matrix(plan)
    @test global_result.status ==
          :materialized_global_electron_nuclear_by_center_matrix
    @test global_result.by_center
    @test !global_result.centers_summed
    @test global_result.center_index == 1
    @test !global_result.nuclear_charge_applied
    @test global_result.charge_application_stage ===
          :acceptance_or_hamiltonian_assembly
    @test global_result.global_electron_nuclear_by_center_matrix_materialized
    @test !global_result.hamiltonian_data_materialized
    @test isapprox(
        global_result.matrix[
            1:size(block.block, 1),
            (size(block.block, 1) + 1):sum(size(block.block)),
        ],
        block.block;
        atol = 0.0,
        rtol = 0.0,
    )
    @test isapprox(
        global_result.matrix[
            (size(block.block, 1) + 1):sum(size(block.block)),
            1:size(block.block, 1),
        ],
        transpose(block.block);
        atol = 0.0,
        rtol = 0.0,
    )
end
