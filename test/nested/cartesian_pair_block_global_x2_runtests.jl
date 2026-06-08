# Runtime role: tiny smoke / global x2-axis assembly pilots.
#
# This validates x2_x/y/z dense global assembly from already-placeable
# final-local placement records. It keeps x2 terms separated and does not cover
# Hamiltonians, overlap, kinetic/position, Coulomb, IDA/MWG, PQS realization,
# exports, or White-Lindsey oracle paths.

using Test
using GaussletBases

const CPBMGlobalX2 = GaussletBases.CartesianPairBlockMaterialization

function _global_x2_result(
    term::Symbol,
    pair_key::Tuple{Symbol,Symbol},
    block;
    source_operator_blocks_materialized::Bool,
    final_pair_blocks_materialized::Bool,
    metadata = (;),
)
    return CPBMGlobalX2.PairBlockMaterializationResult(
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

function _global_x2_axis_offset(term::Symbol)
    term === :x2_x && return 0.0
    term === :x2_y && return 20.0
    term === :x2_z && return 40.0
    throw(ArgumentError("x2 test term must be x2_x/y/z"))
end

function _global_x2_entry(result; block_set_term)
    return CPBMGlobalX2._one_body_local_block_collection_entry(
        result;
        block_set_term,
    )
end

function _global_x2_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return sum(entry -> hasproperty(entry, :count) ? entry.count : 1, matches)
end

function _global_x2_blocks(term::Symbol)
    offset = _global_x2_axis_offset(term)
    return (;
        diagonal = [
            7.5 + offset 0.8 + offset
            0.8 + offset 8.5 + offset
        ],
        off_diagonal = [
            9.0 + offset 10.0 + offset
            11.0 + offset 12.0 + offset
        ],
        source = [
            100.0 + offset 101.0 + offset
            102.0 + offset 103.0 + offset
        ],
        missing_ranges = [
            4.0 + offset 0.2 + offset
            0.2 + offset 4.4 + offset
        ],
    )
end

function _global_x2_expected(term::Symbol)
    blocks = _global_x2_blocks(term)
    diagonal = blocks.diagonal
    off_diagonal = blocks.off_diagonal
    return [
        diagonal[1, 1] diagonal[1, 2] off_diagonal[1, 1] off_diagonal[1, 2]
        diagonal[2, 1] diagonal[2, 2] off_diagonal[2, 1] off_diagonal[2, 2]
        off_diagonal[1, 1] off_diagonal[2, 1] 0.0 0.0
        off_diagonal[1, 2] off_diagonal[2, 2] 0.0 0.0
    ]
end

function _global_x2_collection(term::Symbol; include_missing_ranges::Bool = true)
    blocks = _global_x2_blocks(term)
    diagonal = _global_x2_entry(
        _global_x2_result(
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
    off_diagonal = _global_x2_entry(
        _global_x2_result(
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
    source_term = Symbol("source_", String(term))
    pqs_source = _global_x2_entry(
        _global_x2_result(
            source_term,
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
        missing_ranges = _global_x2_entry(
            _global_x2_result(
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

function _global_x2_plan(
    term::Symbol;
    global_dimension = 4,
    include_missing_ranges::Bool = true,
)
    term === :x2_x && return CPBMGlobalX2.one_body_x2_x_placement_plan(
        _global_x2_collection(term; include_missing_ranges);
        global_dimension,
    )
    term === :x2_y && return CPBMGlobalX2.one_body_x2_y_placement_plan(
        _global_x2_collection(term; include_missing_ranges);
        global_dimension,
    )
    term === :x2_z && return CPBMGlobalX2.one_body_x2_z_placement_plan(
        _global_x2_collection(term; include_missing_ranges);
        global_dimension,
    )
    throw(ArgumentError("x2 test term must be x2_x/y/z"))
end

function _global_x2_matrix(term::Symbol, placement_plan)
    term === :x2_x &&
        return CPBMGlobalX2.one_body_global_x2_x_matrix(placement_plan)
    term === :x2_y &&
        return CPBMGlobalX2.one_body_global_x2_y_matrix(placement_plan)
    term === :x2_z &&
        return CPBMGlobalX2.one_body_global_x2_z_matrix(placement_plan)
    throw(ArgumentError("x2 test term must be x2_x/y/z"))
end

function _global_x2_alias_matrix(term::Symbol, placement_plan)
    term === :x2_x && return CPBMGlobalX2.global_x2_x_matrix(placement_plan)
    term === :x2_y && return CPBMGlobalX2.global_x2_y_matrix(placement_plan)
    term === :x2_z && return CPBMGlobalX2.global_x2_z_matrix(placement_plan)
    throw(ArgumentError("x2 test term must be x2_x/y/z"))
end

@testset "CartesianPairBlockMaterialization global x2 matrices" begin
    for term in (:x2_x, :x2_y, :x2_z)
        placement_plan = _global_x2_plan(term)
        result = _global_x2_matrix(term, placement_plan)
        alias_result = _global_x2_alias_matrix(term, placement_plan)
        term_string = String(term)
        materialized_flag = Symbol("global_", term_string, "_matrix_materialized")

        @test result.object_kind ===
              Symbol("cartesian_pair_block_global_", term_string, "_matrix")
        @test result.status === Symbol("materialized_global_", term_string, "_matrix")
        @test isnothing(result.blocker)
        @test result.term === term
        @test result.global_dimension == 4
        @test result.matrix == _global_x2_expected(term)
        @test alias_result.matrix == _global_x2_expected(term)
        @test result.placeable_record_count == 2
        @test result.blocked_record_count == 2
        @test result.placed_block_count == 3
        @test result.skipped_block_count == 2
        @test result.symmetry === :symmetric

        @test _global_x2_count(
            placement_plan.blocker_counts,
            :blocker,
            :source_space_block_requires_shell_realization,
        ) == 1
        @test _global_x2_count(
            placement_plan.blocker_counts,
            :blocker,
            :missing_column_ranges,
        ) == 1
        @test !any(value -> value >= 100.0, result.matrix)

        missing_dimension_plan = _global_x2_plan(term; global_dimension = nothing)
        missing_dimension = _global_x2_matrix(term, missing_dimension_plan)
        @test missing_dimension.status ===
              Symbol("blocked_global_", term_string, "_matrix")
        @test missing_dimension.blocker === :missing_global_dimension
        @test isnothing(missing_dimension.matrix)

        non_x2 = _global_x2_matrix(
            term,
            merge(placement_plan, (; term = :overlap)),
        )
        @test non_x2.status === Symbol("blocked_global_", term_string, "_matrix")
        @test non_x2.blocker === Symbol("non_", term_string, "_placement_plan")

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
        shape_mismatch = _global_x2_matrix(term, shape_mismatch_plan)
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
        outside = _global_x2_matrix(term, outside_plan)
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
        empty = _global_x2_matrix(term, empty_plan)
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

    @test_throws ArgumentError CPBMGlobalX2.one_body_placement_plan(
        _global_x2_collection(:x2_x);
        term = :coulomb,
    )
end
