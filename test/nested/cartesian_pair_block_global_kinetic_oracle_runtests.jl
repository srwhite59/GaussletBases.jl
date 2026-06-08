# Runtime role: focused oracle / global kinetic equivalence.
#
# This compares selected White-Lindsey fixed-block kinetic slices against the
# new local-block placement and global kinetic assembly pilot. It is narrower
# than the broad White-Lindsey oracle comparison gate and covers kinetic only.

using Test

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const CPBMGlobalKineticOracle = CPBMForLWAdapter

function _global_kinetic_oracle_seed_source_1d(seed; expansion)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        seed.basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
end

function _global_kinetic_oracle_ranges(seed_report)
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

function _global_kinetic_oracle_factors(doside_source_1d)
    pgdg = doside_source_1d.pgdg_intermediate
    return (;
        overlap_1d = (;
            x = pgdg.overlap,
            y = pgdg.overlap,
            z = pgdg.overlap,
        ),
        kinetic_1d = (;
            x = pgdg.kinetic,
            y = pgdg.kinetic,
            z = pgdg.kinetic,
        ),
    )
end

function _global_kinetic_oracle_result_with_ranges(
    result,
    left_column_range::UnitRange{Int},
    right_column_range::UnitRange{Int},
    pair_index::Int,
)
    return CPBMGlobalKineticOracle.PairBlockMaterializationResult(
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

function _global_kinetic_oracle_entry(result; block_set_term = :kinetic)
    return CPBMGlobalKineticOracle._one_body_local_block_collection_entry(
        result;
        block_set_term,
    )
end

function _global_kinetic_oracle_collection(entries::Tuple)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        terms = (:kinetic,),
        requested_terms = (:kinetic,),
        materialized_terms = (:kinetic,),
        deferred_terms = (),
        entries,
        materialized_entries = entries,
        skipped_entries = (),
    )
end

function _global_kinetic_oracle_max_abs_error(left, right)
    return maximum(abs.(left .- right))
end

@testset "CartesianPairBlockMaterialization global kinetic oracle equivalence" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    seed_report =
        GaussletBases._white_lindsey_low_order_materialized_seed_report(;
            expansion,
        )
    @test seed_report.retained_dimension == 223

    seed = seed_report.fixture
    ranges = _global_kinetic_oracle_ranges(seed_report)
    @test ranges.face_global_range == 162:170
    @test ranges.edge_global_range == 210:212

    doside_source_1d =
        _global_kinetic_oracle_seed_source_1d(seed; expansion)
    factors = _global_kinetic_oracle_factors(doside_source_1d)
    real_units = _lw_adapter_real_units(
        doside_source_1d;
        prefix = "global_kinetic_oracle",
    )

    face_face_coefficients =
        CPBMGlobalKineticOracle.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_facet_unit,
                real_units.real_facet_unit,
                1,
            ),
        )
    face_edge_coefficients =
        CPBMGlobalKineticOracle.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_facet_unit,
                real_units.real_edge_unit,
                2,
            ),
        )

    face_face_kinetic =
        CPBMGlobalKineticOracle.white_lindsey_boundary_stratum_kinetic_block(
            face_face_coefficients;
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            kinetic_1d = factors.kinetic_1d,
        )
    face_edge_kinetic =
        CPBMGlobalKineticOracle.white_lindsey_boundary_stratum_kinetic_block(
            face_edge_coefficients;
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            kinetic_1d = factors.kinetic_1d,
        )

    @test face_face_kinetic.metadata.kinetic_factor_form ===
          :factorized_cartesian_sum_kss_sks_ssk
    @test face_face_kinetic.metadata.kinetic_component_axes == (:x, :y, :z)
    @test face_edge_kinetic.metadata.kinetic_factor_form ===
          :factorized_cartesian_sum_kss_sks_ssk
    @test face_edge_kinetic.metadata.kinetic_component_axes == (:x, :y, :z)

    face_face_result = _global_kinetic_oracle_result_with_ranges(
        face_face_kinetic,
        ranges.face_global_range,
        ranges.face_global_range,
        1,
    )
    face_edge_result = _global_kinetic_oracle_result_with_ranges(
        face_edge_kinetic,
        ranges.face_global_range,
        ranges.edge_global_range,
        2,
    )
    entries = (
        _global_kinetic_oracle_entry(face_face_result),
        _global_kinetic_oracle_entry(face_edge_result),
    )
    collection = _global_kinetic_oracle_collection(entries)
    placement_plan = CPBMGlobalKineticOracle.one_body_kinetic_placement_plan(
        collection;
        global_dimension = seed_report.retained_dimension,
    )
    global_result =
        CPBMGlobalKineticOracle.one_body_global_kinetic_matrix(placement_plan)

    @test placement_plan.record_count == 2
    @test placement_plan.placeable_count == 2
    @test placement_plan.blocked_count == 0
    @test all(
        record ->
            record.selector_family === :white_lindsey_boundary_stratum &&
            record.block_space === :final_local_space,
        placement_plan.placeable_records,
    )

    @test global_result.status === :materialized_global_kinetic_matrix
    @test global_result.global_dimension == 223
    @test size(global_result.matrix) == (223, 223)
    @test global_result.placeable_record_count == 2
    @test global_result.blocked_record_count == 0
    @test global_result.placed_block_count == 3

    oracle_kinetic = seed.fixed_block.kinetic
    face_face_error = _global_kinetic_oracle_max_abs_error(
        global_result.matrix[
            ranges.face_global_range,
            ranges.face_global_range,
        ],
        oracle_kinetic[ranges.face_global_range, ranges.face_global_range],
    )
    face_edge_error = _global_kinetic_oracle_max_abs_error(
        global_result.matrix[
            ranges.face_global_range,
            ranges.edge_global_range,
        ],
        oracle_kinetic[ranges.face_global_range, ranges.edge_global_range],
    )
    edge_face_error = _global_kinetic_oracle_max_abs_error(
        global_result.matrix[
            ranges.edge_global_range,
            ranges.face_global_range,
        ],
        oracle_kinetic[ranges.edge_global_range, ranges.face_global_range],
    )

    @test face_face_error <= 1.0e-10
    @test face_edge_error <= 1.0e-10
    @test edge_face_error <= 1.0e-10

    @test !any(record -> record.selector_family === :pqs_source_pair,
        placement_plan.placement_records)
    @test global_result.operator_matrix_materialized
    @test global_result.global_kinetic_matrix_materialized
    @test global_result.global_operator_assembled
    @test !global_result.operator_blocks_materialized
    @test !global_result.hamiltonian_data_materialized
    @test !global_result.coulomb_materialized
    @test !global_result.ida_mwg_data_materialized
    @test !global_result.pqs_lowdin_materialized
    @test !global_result.pqs_shell_projection_materialized
    @test !global_result.exports_materialized
    @test !global_result.artifacts_materialized

    @info "global kinetic oracle max abs errors" face_face_error face_edge_error edge_face_error
end
