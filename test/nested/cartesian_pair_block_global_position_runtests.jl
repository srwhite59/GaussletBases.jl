# Runtime role: tiny smoke / global position-axis assembly pilots.
#
# This validates position_x/y/z dense global assembly from already-placeable
# final-local placement records. It keeps position terms separated and does not
# cover Hamiltonians, overlap, kinetic/x2, Coulomb, IDA/MWG, PQS realization,
# exports, or White-Lindsey oracle paths.

using Test
using GaussletBases

const CPBMGlobalPosition = GaussletBases.CartesianPairBlockMaterialization

function _global_position_result(
    term::Symbol,
    pair_key::Tuple{Symbol,Symbol},
    block;
    source_operator_blocks_materialized::Bool,
    final_pair_blocks_materialized::Bool,
    metadata = (;),
)
    return CPBMGlobalPosition.PairBlockMaterializationResult(
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

function _global_position_axis_offset(term::Symbol)
    term === :position_x && return 0.0
    term === :position_y && return 10.0
    term === :position_z && return 20.0
    throw(ArgumentError("position test term must be position_x/y/z"))
end

function _global_position_entry(result; block_set_term)
    return CPBMGlobalPosition._one_body_local_block_collection_entry(
        result;
        block_set_term,
    )
end

function _global_position_blocks(term::Symbol)
    offset = _global_position_axis_offset(term)
    return (;
        diagonal = [
            1.5 + offset 0.4 + offset
            0.4 + offset 2.5 + offset
        ],
        off_diagonal = [
            3.0 + offset 4.0 + offset
            5.0 + offset 6.0 + offset
        ],
        source = [
            100.0 + offset 101.0 + offset
            102.0 + offset 103.0 + offset
        ],
        missing_ranges = [
            2.0 + offset 0.1 + offset
            0.1 + offset 2.2 + offset
        ],
    )
end

function _global_position_expected(term::Symbol)
    blocks = _global_position_blocks(term)
    diagonal = blocks.diagonal
    off_diagonal = blocks.off_diagonal
    return [
        diagonal[1, 1] diagonal[1, 2] off_diagonal[1, 1] off_diagonal[1, 2]
        diagonal[2, 1] diagonal[2, 2] off_diagonal[2, 1] off_diagonal[2, 2]
        off_diagonal[1, 1] off_diagonal[2, 1] 0.0 0.0
        off_diagonal[1, 2] off_diagonal[2, 2] 0.0 0.0
    ]
end

function _global_position_collection(
    term::Symbol;
    include_missing_ranges::Bool = true,
)
    blocks = _global_position_blocks(term)
    diagonal = _global_position_entry(
        _global_position_result(
            term,
            (:diag_left, :diag_left),
            blocks.diagonal;
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = true,
            metadata = (;
                pair_index = 1,
                selector_family = :direct_direct,
                materialization_path =
                    :direct_direct_pair_block_materialization_pilot,
                left_column_range = 1:2,
                right_column_range = 1:2,
            ),
        );
        block_set_term = term,
    )
    off_diagonal = _global_position_entry(
        _global_position_result(
            term,
            (:left, :right),
            blocks.off_diagonal;
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = true,
            metadata = (;
                pair_index = 2,
                selector_family = :direct_direct,
                materialization_path =
                    :direct_direct_pair_block_materialization_pilot,
                left_column_range = 1:2,
                right_column_range = 3:4,
            ),
        );
        block_set_term = term,
    )
    pqs_source = _global_position_entry(
        _global_position_result(
            Symbol("source_", String(term)),
            (:pqs_left, :pqs_right),
            blocks.source;
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = false,
            metadata = (;
                pair_index = 3,
                selector_family = :pqs_source_pair,
                materialization_path = :pqs_source_pair_preflight,
            ),
        );
        block_set_term = term,
    )
    materialized_entries = (diagonal, off_diagonal, pqs_source)
    if include_missing_ranges
        missing_ranges = _global_position_entry(
            _global_position_result(
                term,
                (:missing_left, :missing_right),
                blocks.missing_ranges;
                source_operator_blocks_materialized = true,
                final_pair_blocks_materialized = true,
                metadata = (;
                    pair_index = 4,
                    selector_family = :direct_direct,
                    materialization_path =
                        :direct_direct_pair_block_materialization_pilot,
                ),
            );
            block_set_term = term,
        )
        materialized_entries = (materialized_entries..., missing_ranges)
    end
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        terms = (term,),
        requested_terms = (term,),
        materialized_terms = (term,),
        deferred_terms = (),
        entries = materialized_entries,
        materialized_entries,
        skipped_entries = (),
    )
end

function _global_position_plan(
    term::Symbol;
    global_dimension = 4,
    include_missing_ranges::Bool = true,
)
    term === :position_x && return CPBMGlobalPosition.one_body_position_x_placement_plan(
        _global_position_collection(term; include_missing_ranges);
        global_dimension,
    )
    term === :position_y && return CPBMGlobalPosition.one_body_position_y_placement_plan(
        _global_position_collection(term; include_missing_ranges);
        global_dimension,
    )
    term === :position_z && return CPBMGlobalPosition.one_body_position_z_placement_plan(
        _global_position_collection(term; include_missing_ranges);
        global_dimension,
    )
    throw(ArgumentError("position test term must be position_x/y/z"))
end

function _global_position_matrix(term::Symbol, placement_plan)
    term === :position_x &&
        return CPBMGlobalPosition.one_body_global_position_x_matrix(placement_plan)
    term === :position_y &&
        return CPBMGlobalPosition.one_body_global_position_y_matrix(placement_plan)
    term === :position_z &&
        return CPBMGlobalPosition.one_body_global_position_z_matrix(placement_plan)
    throw(ArgumentError("position test term must be position_x/y/z"))
end

function _global_position_alias_matrix(term::Symbol, placement_plan)
    term === :position_x &&
        return CPBMGlobalPosition.global_position_x_matrix(placement_plan)
    term === :position_y &&
        return CPBMGlobalPosition.global_position_y_matrix(placement_plan)
    term === :position_z &&
        return CPBMGlobalPosition.global_position_z_matrix(placement_plan)
    throw(ArgumentError("position test term must be position_x/y/z"))
end

function _global_position_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return sum(entry -> hasproperty(entry, :count) ? entry.count : 1, matches)
end

@testset "CartesianPairBlockMaterialization global position matrices" begin
    for term in (:position_x, :position_y, :position_z)
        placement_plan = _global_position_plan(term)
        result = _global_position_matrix(term, placement_plan)
        alias_result = _global_position_alias_matrix(term, placement_plan)
        term_string = String(term)
        materialized_flag = Symbol("global_", term_string, "_matrix_materialized")

        @test result.object_kind ===
              Symbol("cartesian_pair_block_global_", term_string, "_matrix")
        @test result.status === Symbol("materialized_global_", term_string, "_matrix")
        @test isnothing(result.blocker)
        @test result.term === term
        @test result.global_dimension == 4
        @test result.matrix == _global_position_expected(term)
        @test alias_result.matrix == _global_position_expected(term)
        @test result.placeable_record_count == 2
        @test result.blocked_record_count == 2
        @test result.placed_block_count == 3
        @test result.skipped_block_count == 2
        @test result.symmetry === :symmetric

        @test _global_position_count(
            placement_plan.blocker_counts,
            :blocker,
            :source_space_block_requires_shell_realization,
        ) == 1
        @test _global_position_count(
            placement_plan.blocker_counts,
            :blocker,
            :missing_column_ranges,
        ) == 1
        @test !any(value -> value >= 100.0, result.matrix)

        missing_dimension_plan = _global_position_plan(
            term;
            global_dimension = nothing,
        )
        missing_dimension = _global_position_matrix(term, missing_dimension_plan)
        @test missing_dimension.status ===
              Symbol("blocked_global_", term_string, "_matrix")
        @test missing_dimension.blocker === :missing_global_dimension
        @test isnothing(missing_dimension.matrix)

        non_position = _global_position_matrix(
            term,
            merge(placement_plan, (; term = :overlap)),
        )
        @test non_position.status ===
              Symbol("blocked_global_", term_string, "_matrix")
        @test non_position.blocker ===
              Symbol("non_", term_string, "_placement_plan")

        placeable = only(
            record for record in placement_plan.placeable_records
            if record.pair_key == (:left, :right)
        )
        shape_mismatch_record = merge(
            placeable,
            (; target_ranges = (; rows = 1:3, columns = 3:4)),
        )
        shape_mismatch_plan = merge(
            placement_plan,
            (;
                placeable_records = (shape_mismatch_record,),
                placeable_count = 1,
                blocked_count = 0,
                blocked_records = (),
                record_count = 1,
            ),
        )
        shape_mismatch = _global_position_matrix(term, shape_mismatch_plan)
        @test shape_mismatch.status ===
              Symbol("blocked_global_", term_string, "_matrix")
        @test shape_mismatch.blocker === :local_block_shape_mismatch

        outside_record = merge(
            placeable,
            (; target_ranges = (; rows = 1:2, columns = 4:5)),
        )
        outside_plan = merge(
            placement_plan,
            (;
                placeable_records = (outside_record,),
                placeable_count = 1,
                blocked_count = 0,
                blocked_records = (),
                record_count = 1,
            ),
        )
        outside = _global_position_matrix(term, outside_plan)
        @test outside.status === Symbol("blocked_global_", term_string, "_matrix")
        @test outside.blocker === :target_ranges_outside_global_dimension

        empty_plan = merge(
            placement_plan,
            (;
                placeable_records = (),
                placeable_count = 0,
                blocked_count = placement_plan.record_count,
            ),
        )
        empty = _global_position_matrix(term, empty_plan)
        @test empty.status === Symbol("blocked_global_", term_string, "_matrix")
        @test empty.blocker === Symbol("no_placeable_", term_string, "_blocks")

        @test result.operator_matrix_materialized
        @test getproperty(result, materialized_flag)
        @test result.global_operator_assembled
        @test !result.operator_blocks_materialized
        @test !result.hamiltonian_data_materialized
        @test !result.artifacts_materialized
        @test !result.exports_materialized
        @test !result.global_operator_blocks_materialized
        @test !result.global_hamiltonian_data_materialized
        @test !result.global_artifacts_materialized
        @test !result.coulomb_materialized
        @test !result.density_density_materialized
        @test !result.ida_mwg_data_materialized
        @test !result.pqs_lowdin_materialized
        @test !result.pqs_shell_projection_materialized
        @test !result.full_white_lindsey_route_assembled
    end

    @test_throws ArgumentError CPBMGlobalPosition.one_body_placement_plan(
        _global_position_collection(:position_x);
        term = :coulomb,
    )
end
