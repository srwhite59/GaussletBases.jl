# Runtime role: tiny smoke / local one-body placement-plan contract.
#
# This checks metadata-only placement records for a local overlap block
# collection. It does not assemble a global matrix, run route drivers, or touch
# White-Lindsey oracle fixtures.

using Test
using GaussletBases

const CPBMPlacement = GaussletBases.CartesianPairBlockMaterialization

function _placement_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return sum(entry -> hasproperty(entry, :count) ? entry.count : 1, matches)
end

function _placement_result(
    term::Symbol,
    pair_key::Tuple{Symbol,Symbol},
    block;
    source_operator_blocks_materialized::Bool,
    final_pair_blocks_materialized::Bool,
    metadata = (;),
)
    return CPBMPlacement.PairBlockMaterializationResult(
        term,
        pair_key,
        Matrix{Float64}(block),
        true,
        source_operator_blocks_materialized,
        final_pair_blocks_materialized,
        false,
        false,
        false,
        NamedTuple(metadata),
    )
end

function _placement_materialized_entry(result; block_set_term = :overlap)
    return CPBMPlacement._one_body_local_block_collection_entry(
        result;
        block_set_term,
    )
end

function _placement_skipped_entry(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    selector_family::Symbol,
    blocker::Symbol,
)
    return CPBMPlacement._one_body_local_block_collection_skipped_entry(
        (;
            requested_term = :overlap,
            pair_key,
            pair_index,
            selector_family,
            materialization_path = :synthetic_not_materialized,
            blocker,
        );
        block_set_term = :overlap,
    )
end

function _placement_collection(; include_direct_ranges::Bool = true)
    direct_metadata =
        include_direct_ranges ?
        (;
            pair_index = 1,
            selector_family = :direct_direct,
            materialization_path = :direct_direct_pair_block_materialization_pilot,
            left_column_range = 1:2,
            right_column_range = 3:4,
        ) :
        (;
            pair_index = 1,
            selector_family = :direct_direct,
            materialization_path = :direct_direct_pair_block_materialization_pilot,
        )
    direct = _placement_materialized_entry(
        _placement_result(
            :overlap,
            (:direct_left, :direct_right),
            [1.0 0.2; 0.2 1.1];
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = true,
            metadata = direct_metadata,
        ),
    )
    pqs = _placement_materialized_entry(
        _placement_result(
            :source_overlap,
            (:pqs_left, :pqs_right),
            [1.0 0.3; 0.4 1.2];
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = false,
            metadata = (;
                pair_index = 2,
                selector_family = :pqs_source_pair,
                materialization_path = :pqs_source_pair_preflight,
            ),
        ),
    )
    missing_ranges = _placement_materialized_entry(
        _placement_result(
            :overlap,
            (:missing_left, :missing_right),
            [2.0 0.1; 0.1 2.2];
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = true,
            metadata = (;
                pair_index = 3,
                selector_family = :direct_direct,
                materialization_path =
                    :direct_direct_pair_block_materialization_pilot,
            ),
        ),
    )
    lw_skip = _placement_skipped_entry(
        (:lw_left, :lw_right),
        4,
        :white_lindsey_boundary_stratum,
        :missing_white_lindsey_unit_pair,
    )
    unsupported_skip = _placement_skipped_entry(
        (:unsupported_left, :unsupported_right),
        5,
        :unsupported,
        :unsupported_pair_block_materialization_path,
    )
    materialized_entries = (direct, pqs, missing_ranges)
    skipped_entries = (lw_skip, unsupported_skip)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        terms = (:overlap,),
        requested_terms = (:overlap,),
        materialized_terms = (:overlap,),
        deferred_terms = (),
        entries = (materialized_entries..., skipped_entries...),
        materialized_entries,
        skipped_entries,
    )
end

@testset "CartesianPairBlockMaterialization local one-body placement plan" begin
    collection = _placement_collection()
    plan = CPBMPlacement.one_body_overlap_placement_plan(
        collection;
        global_dimension = 8,
    )

    @test plan.object_kind ==
          :cartesian_pair_block_local_one_body_placement_plan
    @test plan.status ==
          :partially_placeable_local_one_body_placement_plan
    @test plan.blocker === :blocked_placement_records_present
    @test plan.term === :overlap
    @test plan.global_dimension == 8
    @test plan.global_dimension_status === :available
    @test plan.record_count == 5
    @test plan.placeable_count == 1
    @test plan.blocked_count == 4
    @test _placement_count(
        plan.blocker_counts,
        :blocker,
        :source_space_block_requires_shell_realization,
    ) == 1
    @test _placement_count(plan.blocker_counts, :blocker, :missing_column_ranges) == 1
    @test _placement_count(
        plan.blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test _placement_count(
        plan.blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 1

    direct_record = only(plan.placeable_records)
    @test direct_record.term === :overlap
    @test direct_record.pair_key == (:direct_left, :direct_right)
    @test direct_record.pair_index == 1
    @test direct_record.selector_family === :direct_direct
    @test direct_record.block_space === :final_local_space
    @test direct_record.block_shape == (2, 2)
    @test direct_record.status === :placeable_local_one_body_placement_record
    @test isnothing(direct_record.blocker)
    @test direct_record.left_column_range == 1:2
    @test direct_record.right_column_range == 3:4
    @test direct_record.target_ranges.rows == 1:2
    @test direct_record.target_ranges.columns == 3:4
    @test direct_record.placeable_in_global_retained_operator
    @test direct_record.record_copies_block_matrix == false

    pqs_record = only(
        record for record in plan.blocked_records
        if record.selector_family === :pqs_source_pair
    )
    @test pqs_record.block_space === :source_space
    @test pqs_record.result_term === :source_overlap
    @test !pqs_record.placeable_in_global_retained_operator
    @test pqs_record.blocker === :source_space_block_requires_shell_realization
    @test pqs_record.not_placeable_reason ===
          :source_space_block_requires_shell_realization

    missing_record = only(
        record for record in plan.blocked_records
        if record.pair_key == (:missing_left, :missing_right)
    )
    @test missing_record.block_space === :final_local_space
    @test missing_record.blocker === :missing_column_ranges
    @test isnothing(missing_record.left_column_range)
    @test isnothing(missing_record.right_column_range)

    lw_record = only(
        record for record in plan.blocked_records
        if record.selector_family === :white_lindsey_boundary_stratum
    )
    unsupported_record = only(
        record for record in plan.blocked_records
        if record.selector_family === :unsupported
    )
    @test lw_record.blocker === :missing_white_lindsey_unit_pair
    @test lw_record.skipped_record.blocker === :missing_white_lindsey_unit_pair
    @test unsupported_record.blocker ===
          :unsupported_pair_block_materialization_path
    @test unsupported_record.skipped_record.blocker ===
          :unsupported_pair_block_materialization_path

    all_blocked = CPBMPlacement.one_body_placement_plan(
        _placement_collection(include_direct_ranges = false);
        term = :overlap,
    )
    @test all_blocked.status === :blocked_local_one_body_placement_plan
    @test all_blocked.placeable_count == 0
    @test all_blocked.blocked_count == 5
    @test all_blocked.global_dimension_status === :not_supplied
    @test _placement_count(
        all_blocked.blocker_counts,
        :blocker,
        :missing_column_ranges,
    ) == 2

    @test !plan.operator_matrix_materialized
    @test !plan.operator_blocks_materialized
    @test !plan.hamiltonian_data_materialized
    @test !plan.artifacts_materialized
    @test !plan.exports_materialized
    @test !plan.global_operator_assembled
    @test !plan.global_operator_blocks_materialized
    @test !plan.global_hamiltonian_data_materialized
    @test !plan.global_artifacts_materialized
    @test !plan.coulomb_materialized
    @test !plan.density_density_materialized
    @test !plan.ida_mwg_data_materialized
    @test !plan.pqs_lowdin_materialized
    @test !plan.pqs_shell_projection_materialized
    @test !plan.full_white_lindsey_route_assembled

    kinetic_plan = CPBMPlacement.one_body_kinetic_placement_plan(collection)
    @test kinetic_plan.term === :kinetic
    @test kinetic_plan.status === :empty_local_one_body_placement_plan
    @test kinetic_plan.record_count == 0

    position_x_plan = CPBMPlacement.one_body_position_x_placement_plan(collection)
    @test position_x_plan.term === :position_x
    @test position_x_plan.status === :empty_local_one_body_placement_plan
    @test position_x_plan.record_count == 0

    position_y_plan = CPBMPlacement.one_body_position_y_placement_plan(collection)
    @test position_y_plan.term === :position_y
    @test position_y_plan.status === :empty_local_one_body_placement_plan
    @test position_y_plan.record_count == 0

    position_z_plan = CPBMPlacement.one_body_position_z_placement_plan(collection)
    @test position_z_plan.term === :position_z
    @test position_z_plan.status === :empty_local_one_body_placement_plan
    @test position_z_plan.record_count == 0

    x2_x_plan = CPBMPlacement.one_body_x2_x_placement_plan(collection)
    @test x2_x_plan.term === :x2_x
    @test x2_x_plan.status === :empty_local_one_body_placement_plan
    @test x2_x_plan.record_count == 0

    x2_y_plan = CPBMPlacement.one_body_x2_y_placement_plan(collection)
    @test x2_y_plan.term === :x2_y
    @test x2_y_plan.status === :empty_local_one_body_placement_plan
    @test x2_y_plan.record_count == 0

    x2_z_plan = CPBMPlacement.one_body_x2_z_placement_plan(collection)
    @test x2_z_plan.term === :x2_z
    @test x2_z_plan.status === :empty_local_one_body_placement_plan
    @test x2_z_plan.record_count == 0

    @test_throws ArgumentError CPBMPlacement.one_body_placement_plan(
        collection;
        term = :coulomb,
    )
end
