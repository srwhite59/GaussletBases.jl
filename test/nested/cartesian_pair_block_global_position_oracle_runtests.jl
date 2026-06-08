# Runtime role: focused oracle / global position_x equivalence.
#
# This compares selected White-Lindsey fixed-block position_x slices against
# the new local-block placement and global position_x assembly pilot. It is
# narrower than the broad White-Lindsey oracle comparison gate and covers
# position_x only.

using Test

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const CPBMGlobalPositionOracle = CPBMForLWAdapter

function _global_position_oracle_seed_source_1d(seed; expansion)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        seed.basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
end

function _global_position_oracle_ranges(seed_report)
    seed = seed_report.fixture
    shell = only(seed.sequence.shell_layers)
    face_index = findfirst(
        face -> (
            face.face_kind === :yz &&
            face.fixed_axis === :x &&
            face.fixed_side === :low
        ),
        shell.faces,
    )
    edge_index = findfirst(
        edge -> (
            edge.free_axis === :z &&
            edge.fixed_axes == (:x, :y) &&
            edge.fixed_sides == (:high, :low)
        ),
        shell.edges,
    )
    isnothing(face_index) &&
        throw(ArgumentError("old White-Lindsey seed face slice was not found"))
    isnothing(edge_index) &&
        throw(ArgumentError("old White-Lindsey seed edge slice was not found"))

    return (;
        face_global_range = seed.inventory.retained_ranges.faces[face_index],
        edge_global_range = seed.inventory.retained_ranges.edges[edge_index],
    )
end

function _global_position_oracle_factors(doside_source_1d)
    pgdg = doside_source_1d.pgdg_intermediate
    return (;
        overlap_1d = (;
            x = pgdg.overlap,
            y = pgdg.overlap,
            z = pgdg.overlap,
        ),
        position_1d = (;
            x = pgdg.position,
            y = pgdg.position,
            z = pgdg.position,
        ),
    )
end

function _global_position_oracle_result_with_ranges(
    result,
    left_column_range::UnitRange{Int},
    right_column_range::UnitRange{Int},
    pair_index::Int,
)
    return CPBMGlobalPositionOracle.PairBlockMaterializationResult(
        result.term,
        result.pair_key,
        result.block,
        result.materialized,
        result.source_operator_blocks_materialized,
        result.final_pair_blocks_materialized,
        result.operator_blocks_materialized,
        result.hamiltonian_data_materialized,
        result.artifacts_materialized,
        merge(
            result.metadata,
            (;
                pair_index,
                selector_family = :white_lindsey_boundary_stratum,
                left_column_range,
                right_column_range,
            ),
        ),
    )
end

function _global_position_oracle_entry(result; block_set_term = :position_x)
    return CPBMGlobalPositionOracle._one_body_local_block_collection_entry(
        result;
        block_set_term,
    )
end

function _global_position_oracle_collection(entries::Tuple)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        terms = (:position_x,),
        requested_terms = (:position_x,),
        materialized_terms = (:position_x,),
        deferred_terms = (),
        entries,
        materialized_entries = entries,
        skipped_entries = (),
    )
end

function _global_position_oracle_max_abs_error(left, right)
    return maximum(abs.(left .- right))
end

@testset "CartesianPairBlockMaterialization global position_x oracle equivalence" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    seed_report =
        GaussletBases._white_lindsey_low_order_materialized_seed_report(;
            expansion,
        )
    @test seed_report.retained_dimension == 223

    seed = seed_report.fixture
    ranges = _global_position_oracle_ranges(seed_report)
    @test ranges.face_global_range == 162:170
    @test ranges.edge_global_range == 210:212

    doside_source_1d =
        _global_position_oracle_seed_source_1d(seed; expansion)
    factors = _global_position_oracle_factors(doside_source_1d)
    real_units = _lw_adapter_real_units(
        doside_source_1d;
        prefix = "global_position_x_oracle",
    )

    face_face_coefficients =
        CPBMGlobalPositionOracle.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_facet_unit,
                real_units.real_facet_unit,
                1,
            ),
        )
    face_edge_coefficients =
        CPBMGlobalPositionOracle.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_facet_unit,
                real_units.real_edge_unit,
                2,
            ),
        )

    face_face_position_x =
        CPBMGlobalPositionOracle.white_lindsey_boundary_stratum_position_block(
            face_face_coefficients;
            axis = :x,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            position_1d = factors.position_1d,
        )
    face_edge_position_x =
        CPBMGlobalPositionOracle.white_lindsey_boundary_stratum_position_block(
            face_edge_coefficients;
            axis = :x,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            position_1d = factors.position_1d,
        )

    @test face_face_position_x.term === :position_x
    @test face_edge_position_x.term === :position_x
    @test face_face_position_x.metadata.position_axis === :x
    @test face_edge_position_x.metadata.position_axis === :x

    face_face_result = _global_position_oracle_result_with_ranges(
        face_face_position_x,
        ranges.face_global_range,
        ranges.face_global_range,
        1,
    )
    face_edge_result = _global_position_oracle_result_with_ranges(
        face_edge_position_x,
        ranges.face_global_range,
        ranges.edge_global_range,
        2,
    )
    entries = (
        _global_position_oracle_entry(face_face_result),
        _global_position_oracle_entry(face_edge_result),
    )
    collection = _global_position_oracle_collection(entries)
    placement_plan =
        CPBMGlobalPositionOracle.one_body_position_x_placement_plan(
            collection;
            global_dimension = seed_report.retained_dimension,
        )
    global_result =
        CPBMGlobalPositionOracle.one_body_global_position_x_matrix(
            placement_plan,
        )

    @test placement_plan.record_count == 2
    @test placement_plan.placeable_count == 2
    @test placement_plan.blocked_count == 0
    @test all(
        record ->
            record.selector_family === :white_lindsey_boundary_stratum &&
            record.block_space === :final_local_space,
        placement_plan.placeable_records,
    )

    @test global_result.status === :materialized_global_position_x_matrix
    @test global_result.global_dimension == 223
    @test size(global_result.matrix) == (223, 223)
    @test global_result.placeable_record_count == 2
    @test global_result.blocked_record_count == 0
    @test global_result.placed_block_count == 3

    oracle_position_x = seed.fixed_block.position_x
    face_face_error = _global_position_oracle_max_abs_error(
        global_result.matrix[
            ranges.face_global_range,
            ranges.face_global_range,
        ],
        oracle_position_x[
            ranges.face_global_range,
            ranges.face_global_range,
        ],
    )
    face_edge_error = _global_position_oracle_max_abs_error(
        global_result.matrix[
            ranges.face_global_range,
            ranges.edge_global_range,
        ],
        oracle_position_x[ranges.face_global_range, ranges.edge_global_range],
    )
    edge_face_error = _global_position_oracle_max_abs_error(
        global_result.matrix[
            ranges.edge_global_range,
            ranges.face_global_range,
        ],
        oracle_position_x[ranges.edge_global_range, ranges.face_global_range],
    )

    @test face_face_error <= 1.0e-10
    @test face_edge_error <= 1.0e-10
    @test edge_face_error <= 1.0e-10

    @test !any(record -> record.selector_family === :pqs_source_pair,
        placement_plan.placement_records)
    @test global_result.operator_matrix_materialized
    @test global_result.global_position_x_matrix_materialized
    @test global_result.global_operator_assembled
    @test !global_result.operator_blocks_materialized
    @test !global_result.hamiltonian_data_materialized
    @test !global_result.coulomb_materialized
    @test !global_result.ida_mwg_data_materialized
    @test !global_result.pqs_lowdin_materialized
    @test !global_result.pqs_shell_projection_materialized
    @test !global_result.exports_materialized
    @test !global_result.artifacts_materialized

    @info "global position_x oracle max abs errors" face_face_error face_edge_error edge_face_error
end
