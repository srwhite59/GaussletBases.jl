# Runtime role: synthetic provider-level CPB overlap placement pilot.
#
# This test validates one local dense CPB overlap block and a tiny local
# overlap collection placement into dense provider-level matrices. It does not
# route-wire, adopt real reports, or build Hamiltonian/operator physics.

using Test
using GaussletBases

const CPBPilot = GaussletBases.CartesianCPBBlockProviders
const _PILOT_BLOCK_KEY = (:tiny_left, :tiny_right)
const _PILOT_ORDERING = :parent_compatible_x_slowest_z_fastest

function _pilot_dense_block(matrix)
    dense = Matrix{Float64}(matrix)
    return CPBPilot.CPBOverlapDenseBlock(
        nothing,
        dense,
        (;
            object_kind = :test_cpb_overlap_dense_block_summary,
            status = :materialized_cpb_overlap_dense_block,
            blocker = nothing,
            dense_block_available = true,
            dense_block_shape = size(dense),
            dense_block_eltype = eltype(dense),
            left_shape = (x = size(dense, 1), y = 1, z = 1),
            right_shape = (x = size(dense, 2), y = 1, z = 1),
            left_support_count = size(dense, 1),
            right_support_count = size(dense, 2),
            local_ordering = _PILOT_ORDERING,
            factor_space = :parent_axis_bundle_pgdg_intermediate,
            factor_convention = :axis_bundle_one_body_overlap,
            normalization_convention = :test_overlap_normalization,
            index_domain = :parent_axis_indices,
            index_domain_source = :axis_bundle_contract,
            index_domain_status =
                :assumed_parent_axis_indexed_by_current_axis_bundle_contract,
            dense_local_block_materialized = true,
            route_driver_wiring = false,
            global_matrix_materialized = false,
            hamiltonian_data_materialized = false,
            coulomb_data_materialized = false,
            ida_mwg_semantics = false,
            exports_or_artifacts = false,
        ),
    )
end

function _pilot_record(block; block_key = _PILOT_BLOCK_KEY)
    dense_summary = CPBPilot.summary(block)
    return CPBPilot.CPBLocalOverlapBlockRecord(
        block_key,
        block,
        (;
            object_kind = :test_cpb_local_overlap_block_record_summary,
            term = :overlap,
            block_key,
            source_kind = :cpb_overlap_dense_block,
            status = :available_cpb_local_overlap_block_record,
            blocker = nothing,
            left_cpb_summary = (;
                support_count = size(block.dense_block, 1),
                shape = dense_summary.left_shape,
            ),
            right_cpb_summary = (;
                support_count = size(block.dense_block, 2),
                shape = dense_summary.right_shape,
            ),
            interval_pair_summary = :test_interval_pair_summary,
            axis_block_summary = :test_axis_block_summary,
            dense_block_summary = dense_summary,
            dense_block_available = true,
            dense_block_shape = size(block.dense_block),
            local_ordering = _PILOT_ORDERING,
            factor_space = dense_summary.factor_space,
            factor_convention = dense_summary.factor_convention,
            normalization_convention = dense_summary.normalization_convention,
            index_domain = dense_summary.index_domain,
            index_domain_source = dense_summary.index_domain_source,
            index_domain_status = dense_summary.index_domain_status,
            placement_status = :unassigned,
            retained_transform_status = :unassigned,
            global_matrix_materialized = false,
            route_driver_wiring = false,
            hamiltonian_data_materialized = false,
            coulomb_data_materialized = false,
            ida_mwg_semantics = false,
            exports_or_artifacts = false,
        ),
    )
end

function _pilot_non_dense_record(block_key)
    return CPBPilot.CPBLocalOverlapBlockRecord(
        block_key,
        nothing,
        (;
            object_kind = :test_cpb_local_overlap_block_record_summary,
            term = :overlap,
            block_key,
            source_kind = :test_non_dense_source,
            status = :available_cpb_local_overlap_block_record,
            blocker = nothing,
            left_cpb_summary = (; support_count = 1, shape = (x = 1, y = 1, z = 1)),
            right_cpb_summary = (; support_count = 1, shape = (x = 1, y = 1, z = 1)),
            interval_pair_summary = :test_interval_pair_summary,
            axis_block_summary = :test_axis_block_summary,
            dense_block_summary = :unavailable,
            dense_block_available = false,
            dense_block_shape = :not_materialized,
            local_ordering = _PILOT_ORDERING,
            factor_space = :parent_axis_bundle_pgdg_intermediate,
            factor_convention = :axis_bundle_one_body_overlap,
            normalization_convention = :test_overlap_normalization,
            index_domain = :parent_axis_indices,
            index_domain_source = :axis_bundle_contract,
            index_domain_status =
                :assumed_parent_axis_indexed_by_current_axis_bundle_contract,
            placement_status = :unassigned,
            retained_transform_status = :unassigned,
            global_matrix_materialized = false,
            route_driver_wiring = false,
            hamiltonian_data_materialized = false,
            coulomb_data_materialized = false,
            ida_mwg_semantics = false,
            exports_or_artifacts = false,
        ),
    )
end

function _pilot_transform(
    side::Symbol,
    matrix;
    block_key = _PILOT_BLOCK_KEY,
    source_support_count = size(matrix, 1),
    target_range = side === :left ? (1:size(matrix, 2)) : (2:1 + size(matrix, 2)),
)
    return CPBPilot.cpb_retained_transform_carry(
        side,
        block_key,
        (; object_kind = :test_cpb_summary, support_count = source_support_count),
        source_support_count,
        _PILOT_ORDERING,
        target_range;
        transform_object = Matrix{Float64}(matrix),
        transform_convention = :test_local_to_retained_columns,
        transform_provenance = :test_fixture,
    )
end

function _pilot_range(
    left_transform,
    right_transform;
    block_key = _PILOT_BLOCK_KEY,
    left_column_range = CPBPilot.summary(left_transform).target_retained_column_range,
    right_column_range = CPBPilot.summary(right_transform).target_retained_column_range,
    global_dimension = max(last(left_column_range), last(right_column_range)),
    global_dimension_source = :test_retained_layout,
)
    return CPBPilot.cpb_source_pair_placement_range(
        block_key;
        left_column_range,
        right_column_range,
        global_dimension,
        global_dimension_source,
        range_source = :test_source_pair_ranges,
        range_provenance = :test_fixture,
        left_transform_carry = left_transform,
        right_transform_carry = right_transform,
    )
end

function _pilot_plan(;
    accepted_block_keys = (_PILOT_BLOCK_KEY,),
    accumulation_rule = :add_explicit_blocks_into_ranges,
)
    return CPBPilot.cpb_reviewed_overlap_placement_plan(;
        placement_plan_kind = :test_reviewed_overlap_placement_plan,
        accumulation_rule,
        accepted_block_keys,
        required_global_dimension_source = :test_retained_layout,
    )
end

function _pilot_facts(block, left_transform, right_transform, placement_range; plan = _pilot_plan())
    collection = CPBPilot.cpb_local_overlap_block_collection((_pilot_record(block),))
    return CPBPilot.cpb_overlap_placement_facts(
        collection;
        transform_carries = (left_transform, right_transform),
        placement_ranges = (placement_range,),
        placement_plan = plan,
    )
end

function _pilot_collection_facts(
    collection,
    transform_carries,
    placement_ranges;
    plan,
    accumulation_rule = nothing,
)
    return CPBPilot.cpb_overlap_placement_facts(
        collection;
        transform_carries,
        placement_ranges,
        placement_plan = plan,
        accumulation_rule,
    )
end

function _test_blocked_pilot_nonclaims(placed)
    placed_summary = CPBPilot.summary(placed)
    @test placed.global_overlap_matrix === nothing
    @test placed_summary.provider_level_matrix_materialized === false
    @test placed_summary.provider_level_overlap_matrix_materialized === false
    @test placed_summary.global_matrix_materialized === false
    @test placed_summary.global_overlap_matrix_materialized === false
    @test placed_summary.route_global_matrix_materialized === false
    @test placed_summary.route_global_overlap_matrix_materialized === false
    @test placed_summary.route_global_overlap_available === false
    @test placed_summary.route_driver_wiring === false
    return nothing
end

function _test_provider_collection_nonclaims(placed)
    placed_summary = CPBPilot.summary(placed)
    @test placed_summary.provider_level_matrix_materialized === true
    @test placed_summary.provider_level_overlap_matrix_materialized === true
    @test placed_summary.global_matrix_materialized === false
    @test placed_summary.global_overlap_matrix_materialized === false
    @test placed_summary.route_global_matrix_materialized === false
    @test placed_summary.route_global_overlap_matrix_materialized === false
    @test placed_summary.route_global_overlap_available === false
    @test placed_summary.route_driver_wiring === false
    return nothing
end

@testset "CPB overlap placement pilot" begin
    dense_block = _pilot_dense_block([2.0 0.5; 0.5 3.0])
    left_transform = _pilot_transform(:left, reshape([1.0, 0.0], 2, 1); target_range = 1:1)
    right_transform = _pilot_transform(:right, reshape([0.0, 1.0], 2, 1); target_range = 2:2)
    placement_range = _pilot_range(
        left_transform,
        right_transform;
        left_column_range = 1:1,
        right_column_range = 2:2,
        global_dimension = 2,
    )
    facts = _pilot_facts(dense_block, left_transform, right_transform, placement_range)

    placed = CPBPilot.cpb_place_overlap_block(
        dense_block,
        left_transform,
        right_transform,
        placement_range,
        facts,
    )
    placed_summary = CPBPilot.summary(placed)
    @test placed_summary.status === :materialized_cpb_overlap_placement_pilot
    @test placed_summary.blocker === nothing
    @test placed_summary.retained_block_shape == (1, 1)
    @test placed_summary.provider_level_global_overlap_matrix_shape == (2, 2)
    @test placed_summary.global_overlap_matrix_shape ===
          :route_global_matrix_not_materialized
    @test placed.global_overlap_matrix !== nothing
    @test placed.global_overlap_matrix == [0.0 0.5; 0.0 0.0]
    @test placed_summary.provider_level_matrix_materialized === true
    @test placed_summary.provider_level_overlap_matrix_materialized === true
    @test placed_summary.provider_level_pilot === true
    @test placed_summary.synthetic_fixture_only === true
    @test placed_summary.global_matrix_materialized === false
    @test placed_summary.global_overlap_matrix_materialized === false
    @test placed_summary.route_global_matrix_materialized === false
    @test placed_summary.route_global_overlap_matrix_materialized === false
    @test placed_summary.route_driver_wiring === false
    @test placed_summary.route_global_overlap_stage_source === false
    @test placed_summary.route_global_overlap_available === false
    @test placed_summary.hamiltonian_data_materialized === false
    @test placed_summary.coulomb_data_materialized === false
    @test placed_summary.ida_mwg_semantics === false
    @test placed_summary.exports_or_artifacts === false

    identity_block = _pilot_dense_block([2.0 0.5; 0.5 3.0])
    identity_left = _pilot_transform(
        :left,
        [1.0 0.0; 0.0 1.0];
        target_range = 1:2,
    )
    identity_right = _pilot_transform(
        :right,
        [1.0 0.0; 0.0 1.0];
        target_range = 3:4,
    )
    identity_range = _pilot_range(
        identity_left,
        identity_right;
        left_column_range = 1:2,
        right_column_range = 3:4,
        global_dimension = 4,
    )
    identity_facts = _pilot_facts(
        identity_block,
        identity_left,
        identity_right,
        identity_range,
    )
    identity_placed = CPBPilot.cpb_place_overlap_block(
        identity_block,
        identity_left,
        identity_right,
        identity_range,
        identity_facts,
    )
    @test CPBPilot.summary(identity_placed).status ===
          :materialized_cpb_overlap_placement_pilot
    @test identity_placed.global_overlap_matrix[1:2, 3:4] ==
          identity_block.dense_block
    @test identity_placed.global_overlap_matrix[3:4, 1:2] == zeros(2, 2)

    mismatch_block = _pilot_dense_block([1.0 0.0; 0.0 1.0; 0.0 0.0])
    mismatch_left = _pilot_transform(
        :left,
        [1.0 0.0; 0.0 1.0];
        source_support_count = 2,
        target_range = 1:2,
    )
    mismatch_right = _pilot_transform(
        :right,
        [1.0 0.0; 0.0 1.0];
        source_support_count = 2,
        target_range = 3:4,
    )
    mismatch_range = _pilot_range(
        mismatch_left,
        mismatch_right;
        left_column_range = 1:2,
        right_column_range = 3:4,
        global_dimension = 4,
    )
    mismatch_facts = _pilot_facts(
        mismatch_block,
        mismatch_left,
        mismatch_right,
        mismatch_range,
    )
    mismatch_placed = CPBPilot.cpb_place_overlap_block(
        mismatch_block,
        mismatch_left,
        mismatch_right,
        mismatch_range,
        mismatch_facts,
    )
    @test CPBPilot.summary(mismatch_placed).status ===
          :blocked_cpb_overlap_placement_pilot
    @test CPBPilot.summary(mismatch_placed).blocker ===
          :retained_transform_shape_mismatch
    _test_blocked_pilot_nonclaims(mismatch_placed)

    placeholder_facts = CPBPilot.cpb_overlap_placement_facts(
        CPBPilot.cpb_local_overlap_block_collection((_pilot_record(dense_block),));
        transform_carries = (left_transform, right_transform),
        placement_ranges = (placement_range,),
        placement_plan = (; kind = :placeholder_overlap_placement_plan),
        accumulation_rule = :add_explicit_blocks_into_ranges,
    )
    placeholder_placed = CPBPilot.cpb_place_overlap_block(
        dense_block,
        left_transform,
        right_transform,
        placement_range,
        placeholder_facts,
    )
    @test CPBPilot.summary(placeholder_placed).blocker ===
          :placement_facts_not_reviewed
    _test_blocked_pilot_nonclaims(placeholder_placed)

    missing_facts = CPBPilot.cpb_overlap_placement_facts(
        CPBPilot.cpb_local_overlap_block_collection((_pilot_record(dense_block),));
        placement_plan = _pilot_plan(),
    )
    missing_placed = CPBPilot.cpb_place_overlap_block(
        dense_block,
        left_transform,
        right_transform,
        placement_range,
        missing_facts,
    )
    @test CPBPilot.summary(missing_placed).blocker ===
          :placement_requirements_missing
    _test_blocked_pilot_nonclaims(missing_placed)

    source_blocked = CPBPilot.cpb_place_overlap_block(
        nothing,
        left_transform,
        right_transform,
        placement_range,
        facts,
    )
    @test CPBPilot.summary(source_blocked).blocker ===
          :local_overlap_source_not_dense
    _test_blocked_pilot_nonclaims(source_blocked)

    key_a = (:a_left, :a_right)
    key_b = (:b_left, :b_right)
    block_a = _pilot_dense_block(fill(2.0, 1, 1))
    block_b = _pilot_dense_block(fill(3.0, 1, 1))
    left_a = _pilot_transform(
        :left,
        fill(1.0, 1, 1);
        block_key = key_a,
        target_range = 1:1,
    )
    right_a = _pilot_transform(
        :right,
        fill(1.0, 1, 1);
        block_key = key_a,
        target_range = 2:2,
    )
    left_b = _pilot_transform(
        :left,
        fill(1.0, 1, 1);
        block_key = key_b,
        target_range = 3:3,
    )
    right_b = _pilot_transform(
        :right,
        fill(1.0, 1, 1);
        block_key = key_b,
        target_range = 4:4,
    )
    range_a = _pilot_range(
        left_a,
        right_a;
        block_key = key_a,
        left_column_range = 1:1,
        right_column_range = 2:2,
        global_dimension = 4,
    )
    range_b = _pilot_range(
        left_b,
        right_b;
        block_key = key_b,
        left_column_range = 3:3,
        right_column_range = 4:4,
        global_dimension = 4,
    )
    collection = CPBPilot.cpb_local_overlap_block_collection((
        _pilot_record(block_a; block_key = key_a),
        _pilot_record(block_b; block_key = key_b),
    ))
    collection_plan = _pilot_plan(accepted_block_keys = (key_a, key_b))
    collection_facts = _pilot_collection_facts(
        collection,
        (left_a, right_a, left_b, right_b),
        (range_a, range_b);
        plan = collection_plan,
    )
    placed_collection = CPBPilot.cpb_place_overlap_collection(
        collection,
        (left_a, right_a, left_b, right_b),
        (range_a, range_b),
        collection_facts,
    )
    placed_collection_summary = CPBPilot.summary(placed_collection)
    expected_collection_matrix = zeros(4, 4)
    expected_collection_matrix[1, 2] = 2.0
    expected_collection_matrix[3, 4] = 3.0
    @test placed_collection_summary.status ===
          :materialized_cpb_overlap_collection_placement_pilot
    @test placed_collection_summary.blocker === nothing
    @test placed_collection_summary.record_count == 2
    @test placed_collection_summary.placed_record_count == 2
    @test placed_collection_summary.blocked_record_count == 0
    @test placed_collection_summary.block_keys == (key_a, key_b)
    @test placed_collection.global_overlap_matrix == expected_collection_matrix
    @test placed_collection_summary.provider_level_global_overlap_matrix_shape == (4, 4)
    _test_provider_collection_nonclaims(placed_collection)

    key_c = (:c_left, :c_right)
    key_d = (:d_left, :d_right)
    block_c = _pilot_dense_block(fill(0.25, 1, 1))
    block_d = _pilot_dense_block(fill(0.75, 1, 1))
    left_c = _pilot_transform(
        :left,
        fill(1.0, 1, 1);
        block_key = key_c,
        target_range = 1:1,
    )
    right_c = _pilot_transform(
        :right,
        fill(1.0, 1, 1);
        block_key = key_c,
        target_range = 2:2,
    )
    left_d = _pilot_transform(
        :left,
        fill(1.0, 1, 1);
        block_key = key_d,
        target_range = 1:1,
    )
    right_d = _pilot_transform(
        :right,
        fill(1.0, 1, 1);
        block_key = key_d,
        target_range = 2:2,
    )
    range_c = _pilot_range(
        left_c,
        right_c;
        block_key = key_c,
        left_column_range = 1:1,
        right_column_range = 2:2,
        global_dimension = 2,
    )
    range_d = _pilot_range(
        left_d,
        right_d;
        block_key = key_d,
        left_column_range = 1:1,
        right_column_range = 2:2,
        global_dimension = 2,
    )
    additive_collection = CPBPilot.cpb_local_overlap_block_collection((
        _pilot_record(block_c; block_key = key_c),
        _pilot_record(block_d; block_key = key_d),
    ))
    additive_plan = _pilot_plan(accepted_block_keys = (key_c, key_d))
    additive_facts = _pilot_collection_facts(
        additive_collection,
        (left_c, right_c, left_d, right_d),
        (range_c, range_d);
        plan = additive_plan,
    )
    additive_placed = CPBPilot.cpb_place_overlap_collection(
        additive_collection,
        (left_c, right_c, left_d, right_d),
        (range_c, range_d),
        additive_facts,
    )
    expected_additive_matrix = zeros(2, 2)
    expected_additive_matrix[1, 2] = 1.0
    @test CPBPilot.summary(additive_placed).status ===
          :materialized_cpb_overlap_collection_placement_pilot
    @test additive_placed.global_overlap_matrix == expected_additive_matrix
    _test_provider_collection_nonclaims(additive_placed)

    duplicate_collection = CPBPilot.cpb_local_overlap_block_collection((
        _pilot_record(block_c; block_key = key_c),
        _pilot_record(block_d; block_key = key_c),
    ))
    duplicate_facts = _pilot_collection_facts(
        duplicate_collection,
        (left_c, right_c),
        (range_c,);
        plan = _pilot_plan(accepted_block_keys = (key_c,)),
    )
    duplicate_placed = CPBPilot.cpb_place_overlap_collection(
        duplicate_collection,
        (left_c, right_c),
        (range_c,),
        duplicate_facts,
    )
    @test CPBPilot.summary(duplicate_placed).blocker ===
          :duplicate_overlap_placement_record
    _test_blocked_pilot_nonclaims(duplicate_placed)

    missing_transform_facts = _pilot_collection_facts(
        collection,
        (left_a, right_a, left_b),
        (range_a, range_b);
        plan = collection_plan,
    )
    missing_transform_collection = CPBPilot.cpb_place_overlap_collection(
        collection,
        (left_a, right_a, left_b),
        (range_a, range_b),
        missing_transform_facts,
    )
    @test CPBPilot.summary(missing_transform_collection).status ===
          :blocked_cpb_overlap_collection_placement_pilot
    @test CPBPilot.summary(missing_transform_collection).blocker ===
          :placement_requirements_missing
    _test_blocked_pilot_nonclaims(missing_transform_collection)

    non_dense_key = (:axis_only_left, :axis_only_right)
    non_dense_left = _pilot_transform(
        :left,
        fill(1.0, 1, 1);
        block_key = non_dense_key,
        target_range = 1:1,
    )
    non_dense_right = _pilot_transform(
        :right,
        fill(1.0, 1, 1);
        block_key = non_dense_key,
        target_range = 2:2,
    )
    non_dense_range = _pilot_range(
        non_dense_left,
        non_dense_right;
        block_key = non_dense_key,
        left_column_range = 1:1,
        right_column_range = 2:2,
        global_dimension = 2,
    )
    non_dense_collection = CPBPilot.cpb_local_overlap_block_collection((
        _pilot_non_dense_record(non_dense_key),
    ))
    non_dense_facts = _pilot_collection_facts(
        non_dense_collection,
        (non_dense_left, non_dense_right),
        (non_dense_range,);
        plan = _pilot_plan(accepted_block_keys = (non_dense_key,)),
    )
    non_dense_placed = CPBPilot.cpb_place_overlap_collection(
        non_dense_collection,
        (non_dense_left, non_dense_right),
        (non_dense_range,),
        non_dense_facts,
    )
    @test CPBPilot.summary(non_dense_placed).blocker ===
          :local_overlap_source_not_dense
    _test_blocked_pilot_nonclaims(non_dense_placed)

    placeholder_collection_facts = _pilot_collection_facts(
        collection,
        (left_a, right_a, left_b, right_b),
        (range_a, range_b);
        plan = (; kind = :placeholder_overlap_placement_plan),
        accumulation_rule = :add_explicit_blocks_into_ranges,
    )
    placeholder_collection = CPBPilot.cpb_place_overlap_collection(
        collection,
        (left_a, right_a, left_b, right_b),
        (range_a, range_b),
        placeholder_collection_facts,
    )
    @test CPBPilot.summary(placeholder_collection).blocker ===
          :placement_facts_not_reviewed
    _test_blocked_pilot_nonclaims(placeholder_collection)
end
