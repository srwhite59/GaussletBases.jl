# Runtime role: metadata-only CPB overlap placement facts contract.
#
# This test validates placement-fact coherence summaries only. It does not
# apply transforms, place local CPB overlap blocks, or assemble route/global
# overlap.

using Test
using GaussletBases

const CPBPlacementFacts = GaussletBases.CartesianCPBBlockProviders

const _PLACEMENT_FACTS_BLOCK_KEY = (:tiny_left, :tiny_right)
const _PLACEMENT_FACTS_SOURCE_SHAPE = (x = 1, y = 2, z = 1)

function _placement_facts_record(;
    status = :available_cpb_local_overlap_block_record,
    blocker = nothing,
    local_ordering = :parent_compatible_x_slowest_z_fastest,
)
    metadata = (;
        object_kind = :test_cpb_local_overlap_block_record_summary,
        term = :overlap,
        block_key = _PLACEMENT_FACTS_BLOCK_KEY,
        source_kind = :cpb_overlap_dense_block,
        status,
        blocker,
        dense_block_available = status === :available_cpb_local_overlap_block_record,
        dense_block_shape =
            status === :available_cpb_local_overlap_block_record ? (2, 2) :
            :not_materialized,
        left_cpb_summary = (;
            object_kind = :test_left_cpb_summary,
            role = :left_fixture,
            shape = _PLACEMENT_FACTS_SOURCE_SHAPE,
        ),
        right_cpb_summary = (;
            object_kind = :test_right_cpb_summary,
            role = :right_fixture,
            shape = _PLACEMENT_FACTS_SOURCE_SHAPE,
        ),
        local_ordering,
        placement_status = :unassigned,
        retained_transform_status = :unassigned,
        global_matrix_materialized = false,
        route_driver_wiring = false,
    )
    return CPBPlacementFacts.CPBLocalOverlapBlockRecord(
        _PLACEMENT_FACTS_BLOCK_KEY,
        nothing,
        metadata,
    )
end

function _placement_facts_collection(; record = _placement_facts_record())
    return CPBPlacementFacts.cpb_local_overlap_block_collection((record,))
end

function _placement_facts_collection(records::Tuple)
    return CPBPlacementFacts.cpb_local_overlap_block_collection(records)
end

function _placement_facts_transform_carry(;
    side = :left,
    source_shape = _PLACEMENT_FACTS_SOURCE_SHAPE,
    target_retained_column_range = side === :left ? (1:2) : (3:4),
    transform_object = [1.0 0.0; 0.0 1.0],
)
    return CPBPlacementFacts.cpb_retained_transform_carry(
        side,
        _PLACEMENT_FACTS_BLOCK_KEY,
        (; object_kind = :test_cpb_source_summary, shape = source_shape),
        source_shape,
        :parent_compatible_x_slowest_z_fastest,
        target_retained_column_range;
        transform_object,
        transform_convention = :test_local_to_retained_columns,
        transform_provenance = :test_fixture,
    )
end

function _placement_facts_range(;
    left_transform_carry = _placement_facts_transform_carry(; side = :left),
    right_transform_carry = _placement_facts_transform_carry(; side = :right),
    left_column_range = 1:2,
    right_column_range = 3:4,
    global_dimension = 4,
    global_dimension_source = :test_retained_layout,
)
    return CPBPlacementFacts.cpb_source_pair_placement_range(
        _PLACEMENT_FACTS_BLOCK_KEY;
        left_column_range,
        right_column_range,
        global_dimension,
        global_dimension_source,
        range_source = :test_source_pair_ranges,
        range_provenance = :test_fixture,
        left_transform_carry,
        right_transform_carry,
    )
end

function _placement_facts_reviewed_plan(;
    accepted_block_keys = (_PLACEMENT_FACTS_BLOCK_KEY,),
    accumulation_rule = :add_explicit_blocks_into_ranges,
)
    return CPBPlacementFacts.cpb_reviewed_overlap_placement_plan(;
        placement_plan_kind = :test_reviewed_overlap_placement_plan,
        accumulation_rule,
        accepted_block_keys,
        required_global_dimension_source = :test_retained_layout,
    )
end

function _complete_placement_facts(; kwargs...)
    left_transform = _placement_facts_transform_carry(; side = :left)
    right_transform = _placement_facts_transform_carry(; side = :right)
    placement_range = _placement_facts_range(;
        left_transform_carry = left_transform,
        right_transform_carry = right_transform,
    )
    return CPBPlacementFacts.cpb_overlap_placement_facts(
        _placement_facts_collection();
        transform_carries = (left_transform, right_transform),
        placement_ranges = (placement_range,),
        kwargs...,
    )
end

@testset "CPB overlap placement facts metadata" begin
    collection = _placement_facts_collection()

    missing_carries = CPBPlacementFacts.cpb_overlap_placement_facts(collection)
    missing_summary = CPBPlacementFacts.summary(missing_carries)
    @test missing_summary.object_kind ===
          :cartesian_cpb_overlap_placement_facts_summary
    @test missing_summary.status === :blocked_cpb_overlap_placement_facts
    @test missing_summary.blocker === :missing_placement_or_retained_transform
    @test missing_summary.collection_available === true
    @test missing_summary.record_count == 1
    @test missing_summary.block_keys === (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test missing_summary.placement_plan_status === :missing_placement_plan
    @test missing_summary.accumulation_rule_status === :missing_accumulation_rule
    @test missing_summary.available_requirements ===
          (:local_cpb_overlap_collection,)
    @test :missing_retained_transform in missing_summary.missing_requirements
    @test :missing_left_column_range in missing_summary.missing_requirements
    @test :missing_right_column_range in missing_summary.missing_requirements
    @test :missing_global_dimension in missing_summary.missing_requirements
    @test :missing_placement_plan in missing_summary.missing_requirements
    @test :missing_accumulation_rule in missing_summary.missing_requirements
    @test missing_summary.global_overlap_status === :blocked
    @test missing_summary.global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test missing_summary.placement_engine_implemented === false
    @test missing_summary.transform_application_implemented === false
    @test missing_summary.global_matrix_materialized === false
    @test missing_summary.route_driver_wiring === false

    missing_record_summary = only(missing_summary.record_fact_summaries)
    @test missing_record_summary.left_transform_status ===
          :missing_retained_transform
    @test missing_record_summary.right_transform_status ===
          :missing_retained_transform
    @test missing_record_summary.placement_range_status ===
          :missing_source_pair_placement_range
    @test missing_record_summary.left_column_range === nothing
    @test missing_record_summary.right_column_range === nothing
    @test missing_record_summary.global_dimension === nothing
    @test missing_record_summary.global_dimension_source === :unavailable
    @test missing_record_summary.left_cpb_summary.object_kind ===
          :test_left_cpb_summary
    @test missing_record_summary.right_cpb_summary.object_kind ===
          :test_right_cpb_summary
    @test missing_record_summary.dense_block_available === true
    @test missing_record_summary.dense_block_shape == (2, 2)
    @test missing_record_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest

    missing_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            missing_carries,
        )
    @test missing_skeleton.object_kind ===
          :cartesian_route_driver_private_global_overlap_placement_plan_skeleton
    @test missing_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test missing_skeleton.blocker ===
          :missing_placement_or_retained_transform
    @test missing_skeleton.local_cpb_overlap_collection_available === true
    @test missing_skeleton.record_count == 1
    @test missing_skeleton.block_keys === (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test length(missing_skeleton.record_placement_summaries) == 1
    @test missing_skeleton.available_requirements ===
          missing_summary.available_requirements
    @test missing_skeleton.missing_requirements ===
          missing_summary.missing_requirements
    @test missing_skeleton.placement_plan_status === :missing_placement_plan
    @test missing_skeleton.placement_plan_kind === :unavailable
    @test missing_skeleton.accumulation_rule_status ===
          :missing_accumulation_rule
    @test missing_skeleton.global_overlap_status === :blocked
    @test missing_skeleton.global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test missing_skeleton.global_matrix_materialized === false
    @test missing_skeleton.global_overlap_matrix_materialized === false
    @test missing_skeleton.route_driver_wiring === false
    @test missing_skeleton.private_global_overlap_input_facts_available === false
    @test missing_skeleton.route_global_overlap_stage_source === false

    missing_skeleton_record =
        only(missing_skeleton.record_placement_summaries)
    @test missing_skeleton_record.block_key === _PLACEMENT_FACTS_BLOCK_KEY
    @test missing_skeleton_record.source_kind ===
          :cpb_overlap_placement_facts_record
    @test missing_skeleton_record.retained_transform_status ===
          :missing_retained_transform
    @test missing_skeleton_record.left_column_range_status ===
          :missing_left_column_range
    @test missing_skeleton_record.right_column_range_status ===
          :missing_right_column_range
    @test missing_skeleton_record.left_cpb_summary.object_kind ===
          :test_left_cpb_summary
    @test missing_skeleton_record.right_cpb_summary.object_kind ===
          :test_right_cpb_summary
    @test missing_skeleton_record.left_transform_status ===
          :missing_retained_transform
    @test missing_skeleton_record.placement_range_status ===
          :missing_source_pair_placement_range

    complete = _complete_placement_facts(;
        placement_plan = (; kind = :test_overlap_placement_plan),
        accumulation_rule = :test_accumulation_rule,
    )
    complete_summary = CPBPlacementFacts.summary(complete)
    @test complete_summary.status === :blocked_cpb_overlap_placement_facts
    @test complete_summary.blocker === :placement_not_implemented
    @test complete_summary.missing_requirements === ()
    @test complete_summary.available_requirements === (
        :local_cpb_overlap_collection,
        :retained_transform,
        :left_column_range,
        :right_column_range,
        :global_dimension,
        :placement_plan,
        :accumulation_rule,
    )
    @test complete_summary.placement_plan_status === :available_placement_plan
    @test complete_summary.placement_plan_kind === :test_overlap_placement_plan
    @test complete_summary.placement_record_inventory_status ===
          :not_checked_cpb_overlap_placement_record_inventory
    @test complete_summary.placement_record_inventory_blocker === nothing
    @test complete_summary.accepted_block_keys === ()
    @test complete_summary.provided_block_keys === (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test complete_summary.rejected_block_keys === ()
    @test complete_summary.duplicate_block_keys === ()
    @test complete_summary.local_ordering_contract_status ===
          :not_checked_cpb_overlap_local_ordering_contract
    @test complete_summary.local_ordering_contract_blocker === nothing
    @test complete_summary.local_ordering_contract === :unavailable
    @test complete_summary.provided_local_orderings ===
          (:parent_compatible_x_slowest_z_fastest,)
    @test complete_summary.mismatched_local_ordering_block_keys === ()
    @test complete_summary.global_dimension_source_contract_status ===
          :not_checked_cpb_overlap_global_dimension_source_contract
    @test complete_summary.global_dimension_source_contract_blocker === nothing
    @test complete_summary.required_global_dimension_source === :unavailable
    @test complete_summary.provided_global_dimension_sources ===
          (:test_retained_layout,)
    @test complete_summary.mismatched_global_dimension_source_block_keys === ()
    @test complete_summary.accumulation_rule_status ===
          :available_accumulation_rule
    @test complete_summary.accumulation_rule === :test_accumulation_rule
    @test complete_summary.global_overlap_status === :blocked
    @test complete_summary.global_overlap_blocker === :placement_not_implemented
    @test complete_summary.placement_engine_implemented === false
    @test complete_summary.transform_application_implemented === false
    @test complete_summary.global_matrix_materialized === false
    @test complete_summary.route_driver_wiring === false

    complete_record_summary = only(complete_summary.record_fact_summaries)
    @test complete_record_summary.left_transform_status ===
          :available_cpb_retained_transform_carry
    @test complete_record_summary.right_transform_status ===
          :available_cpb_retained_transform_carry
    @test complete_record_summary.placement_range_status ===
          :available_cpb_source_pair_placement_range
    @test complete_record_summary.left_column_range == 1:2
    @test complete_record_summary.right_column_range == 3:4
    @test complete_record_summary.global_dimension == 4
    @test complete_record_summary.global_dimension_source ===
          :test_retained_layout

    complete_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            complete,
        )
    @test complete_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test complete_skeleton.blocker === :placement_not_implemented
    @test complete_skeleton.missing_requirements === ()
    @test complete_skeleton.available_requirements ===
          complete_summary.available_requirements
    @test complete_skeleton.placement_plan_status === :available_placement_plan
    @test complete_skeleton.placement_plan_kind === :test_overlap_placement_plan
    @test complete_skeleton.accumulation_rule_status ===
          :available_accumulation_rule
    @test complete_skeleton.global_dimension == 4
    @test complete_skeleton.global_dimension_source === :test_retained_layout
    @test complete_skeleton.global_overlap_status === :blocked
    @test complete_skeleton.global_overlap_blocker === :placement_not_implemented
    @test complete_skeleton.global_matrix_materialized === false
    @test complete_skeleton.global_overlap_matrix_materialized === false
    @test complete_skeleton.route_driver_wiring === false
    @test complete_skeleton.private_global_overlap_input_facts_available === false
    @test complete_skeleton.route_global_overlap_stage_source === false

    complete_skeleton_record =
        only(complete_skeleton.record_placement_summaries)
    @test complete_skeleton_record.retained_transform_status ===
          :available_retained_transform
    @test complete_skeleton_record.left_column_range_status ===
          :available_left_column_range
    @test complete_skeleton_record.right_column_range_status ===
          :available_right_column_range
    @test complete_skeleton_record.left_column_range == 1:2
    @test complete_skeleton_record.right_column_range == 3:4
    @test complete_skeleton_record.global_dimension == 4
    @test complete_skeleton_record.left_cpb_summary.object_kind ===
          :test_left_cpb_summary
    @test complete_skeleton_record.right_cpb_summary.object_kind ===
          :test_right_cpb_summary

    missing_plan = _complete_placement_facts(;
        accumulation_rule = :test_accumulation_rule,
    )
    missing_plan_summary = CPBPlacementFacts.summary(missing_plan)
    @test missing_plan_summary.status === :blocked_cpb_overlap_placement_facts
    @test missing_plan_summary.blocker ===
          :missing_placement_or_retained_transform
    @test missing_plan_summary.placement_plan_status === :missing_placement_plan
    @test :missing_placement_plan in missing_plan_summary.missing_requirements
    @test !(:missing_accumulation_rule in missing_plan_summary.missing_requirements)
    @test missing_plan_summary.global_matrix_materialized === false
    @test missing_plan_summary.route_driver_wiring === false

    reviewed_plan = _placement_facts_reviewed_plan()
    reviewed_plan_facts = _complete_placement_facts(;
        placement_plan = reviewed_plan,
    )
    reviewed_plan_summary = CPBPlacementFacts.summary(reviewed_plan_facts)
    @test reviewed_plan_summary.status === :blocked_cpb_overlap_placement_facts
    @test reviewed_plan_summary.blocker === :placement_not_implemented
    @test reviewed_plan_summary.placement_plan_status === :available_placement_plan
    @test reviewed_plan_summary.placement_plan_kind ===
          :test_reviewed_overlap_placement_plan
    @test reviewed_plan_summary.accumulation_rule_status ===
          :available_accumulation_rule
    @test reviewed_plan_summary.accumulation_rule ===
          :add_explicit_blocks_into_ranges
    @test reviewed_plan_summary.placement_record_inventory_status ===
          :available_cpb_overlap_placement_record_inventory
    @test reviewed_plan_summary.placement_record_inventory_blocker === nothing
    @test reviewed_plan_summary.accepted_block_keys ===
          (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test reviewed_plan_summary.provided_block_keys ===
          (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test reviewed_plan_summary.rejected_block_keys === ()
    @test reviewed_plan_summary.duplicate_block_keys === ()
    @test reviewed_plan_summary.duplicate_record_policy ===
          :reject_duplicate_block_keys
    @test reviewed_plan_summary.local_ordering_contract_status ===
          :available_cpb_overlap_local_ordering_contract
    @test reviewed_plan_summary.local_ordering_contract_blocker === nothing
    @test reviewed_plan_summary.local_ordering_contract ===
          :parent_compatible_x_slowest_z_fastest
    @test reviewed_plan_summary.provided_local_orderings ===
          (:parent_compatible_x_slowest_z_fastest,)
    @test reviewed_plan_summary.mismatched_local_ordering_block_keys === ()
    @test reviewed_plan_summary.global_dimension_source_contract_status ===
          :available_cpb_overlap_global_dimension_source_contract
    @test reviewed_plan_summary.global_dimension_source_contract_blocker ===
          nothing
    @test reviewed_plan_summary.required_global_dimension_source ===
          :test_retained_layout
    @test reviewed_plan_summary.provided_global_dimension_sources ===
          (:test_retained_layout,)
    @test reviewed_plan_summary.mismatched_global_dimension_source_block_keys ===
          ()

    ordering_mismatch_collection = _placement_facts_collection(;
        record = _placement_facts_record(;
            local_ordering = :test_other_local_ordering,
        ),
    )
    ordering_mismatch_left_transform =
        _placement_facts_transform_carry(; side = :left)
    ordering_mismatch_right_transform =
        _placement_facts_transform_carry(; side = :right)
    ordering_mismatch_facts = CPBPlacementFacts.cpb_overlap_placement_facts(
        ordering_mismatch_collection;
        transform_carries = (
            ordering_mismatch_left_transform,
            ordering_mismatch_right_transform,
        ),
        placement_ranges = (
            _placement_facts_range(;
                left_transform_carry = ordering_mismatch_left_transform,
                right_transform_carry = ordering_mismatch_right_transform,
            ),
        ),
        placement_plan = reviewed_plan,
    )
    ordering_mismatch_summary =
        CPBPlacementFacts.summary(ordering_mismatch_facts)
    @test ordering_mismatch_summary.status ===
          :blocked_cpb_overlap_placement_facts
    @test ordering_mismatch_summary.blocker ===
          :overlap_local_ordering_contract_mismatch
    @test ordering_mismatch_summary.local_ordering_contract_status ===
          :blocked_cpb_overlap_local_ordering_contract
    @test ordering_mismatch_summary.local_ordering_contract_blocker ===
          :overlap_local_ordering_contract_mismatch
    @test ordering_mismatch_summary.provided_local_orderings ===
          (:test_other_local_ordering,)
    @test ordering_mismatch_summary.mismatched_local_ordering_block_keys ===
          (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test ordering_mismatch_summary.global_matrix_materialized === false
    @test ordering_mismatch_summary.route_driver_wiring === false

    source_mismatch_left_transform =
        _placement_facts_transform_carry(; side = :left)
    source_mismatch_right_transform =
        _placement_facts_transform_carry(; side = :right)
    source_mismatch_facts = CPBPlacementFacts.cpb_overlap_placement_facts(
        collection;
        transform_carries = (
            source_mismatch_left_transform,
            source_mismatch_right_transform,
        ),
        placement_ranges = (
            _placement_facts_range(;
                left_transform_carry = source_mismatch_left_transform,
                right_transform_carry = source_mismatch_right_transform,
                global_dimension_source = :test_other_retained_layout,
            ),
        ),
        placement_plan = reviewed_plan,
    )
    source_mismatch_summary = CPBPlacementFacts.summary(source_mismatch_facts)
    @test source_mismatch_summary.status ===
          :blocked_cpb_overlap_placement_facts
    @test source_mismatch_summary.blocker ===
          :global_dimension_source_mismatch
    @test source_mismatch_summary.global_dimension_source_contract_status ===
          :blocked_cpb_overlap_global_dimension_source_contract
    @test source_mismatch_summary.global_dimension_source_contract_blocker ===
          :global_dimension_source_mismatch
    @test source_mismatch_summary.required_global_dimension_source ===
          :test_retained_layout
    @test source_mismatch_summary.provided_global_dimension_sources ===
          (:test_other_retained_layout,)
    @test source_mismatch_summary.mismatched_global_dimension_source_block_keys ===
          (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test source_mismatch_summary.global_matrix_materialized === false
    @test source_mismatch_summary.route_driver_wiring === false

    reviewed_missing_ranges = CPBPlacementFacts.cpb_overlap_placement_facts(
        collection;
        placement_plan = reviewed_plan,
    )
    reviewed_missing_ranges_summary =
        CPBPlacementFacts.summary(reviewed_missing_ranges)
    @test reviewed_missing_ranges_summary.status ===
          :blocked_cpb_overlap_placement_facts
    @test reviewed_missing_ranges_summary.blocker ===
          :missing_placement_or_retained_transform
    @test :missing_left_column_range in
          reviewed_missing_ranges_summary.missing_requirements
    @test :missing_right_column_range in
          reviewed_missing_ranges_summary.missing_requirements
    @test :missing_global_dimension in
          reviewed_missing_ranges_summary.missing_requirements
    @test reviewed_missing_ranges_summary.global_dimension_source_contract_status ===
          :not_checked_cpb_overlap_global_dimension_source_contract
    @test reviewed_missing_ranges_summary.global_dimension_source_contract_blocker ===
          nothing
    @test reviewed_missing_ranges_summary.provided_global_dimension_sources ===
          ()
    @test reviewed_missing_ranges_summary.mismatched_global_dimension_source_block_keys ===
          ()

    unaccepted_plan = _placement_facts_reviewed_plan(;
        accepted_block_keys = ((:other_left, :other_right),),
    )
    unaccepted_facts = _complete_placement_facts(;
        placement_plan = unaccepted_plan,
    )
    unaccepted_summary = CPBPlacementFacts.summary(unaccepted_facts)
    @test unaccepted_summary.status === :blocked_cpb_overlap_placement_facts
    @test unaccepted_summary.blocker === :unaccepted_overlap_placement_record
    @test unaccepted_summary.placement_record_inventory_status ===
          :blocked_cpb_overlap_placement_record_inventory
    @test unaccepted_summary.placement_record_inventory_blocker ===
          :unaccepted_overlap_placement_record
    @test unaccepted_summary.rejected_block_keys ===
          (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test unaccepted_summary.duplicate_block_keys === ()
    @test unaccepted_summary.global_matrix_materialized === false
    @test unaccepted_summary.route_driver_wiring === false

    duplicate_collection = _placement_facts_collection((
        _placement_facts_record(),
        _placement_facts_record(),
    ))
    duplicate_facts = CPBPlacementFacts.cpb_overlap_placement_facts(
        duplicate_collection;
        placement_plan = reviewed_plan,
    )
    duplicate_summary = CPBPlacementFacts.summary(duplicate_facts)
    @test duplicate_summary.status === :blocked_cpb_overlap_placement_facts
    @test duplicate_summary.blocker === :duplicate_overlap_placement_record
    @test duplicate_summary.placement_record_inventory_status ===
          :blocked_cpb_overlap_placement_record_inventory
    @test duplicate_summary.placement_record_inventory_blocker ===
          :duplicate_overlap_placement_record
    @test duplicate_summary.duplicate_block_keys === (_PLACEMENT_FACTS_BLOCK_KEY,)
    @test duplicate_summary.rejected_block_keys === ()
    @test duplicate_summary.global_matrix_materialized === false
    @test duplicate_summary.route_driver_wiring === false

    blocked_plan = _placement_facts_reviewed_plan(;
        accumulation_rule = nothing,
    )
    blocked_plan_facts = _complete_placement_facts(;
        placement_plan = blocked_plan,
        accumulation_rule = :test_accumulation_rule,
    )
    blocked_plan_summary = CPBPlacementFacts.summary(blocked_plan_facts)
    @test blocked_plan_summary.status === :blocked_cpb_overlap_placement_facts
    @test blocked_plan_summary.blocker === :missing_accumulation_rule
    @test blocked_plan_summary.placement_plan_status === :blocked_placement_plan
    @test blocked_plan_summary.placement_record_inventory_status ===
          :blocked_cpb_overlap_placement_record_inventory
    @test blocked_plan_summary.placement_record_inventory_blocker ===
          :missing_accumulation_rule
    @test blocked_plan_summary.rejected_block_keys === ()
    @test blocked_plan_summary.duplicate_block_keys === ()

    left_transform = _placement_facts_transform_carry(; side = :left)
    right_transform = _placement_facts_transform_carry(; side = :right)
    blocked_range = _placement_facts_range(;
        left_transform_carry = left_transform,
        right_transform_carry = right_transform,
        left_column_range = 1:3,
        global_dimension = 5,
    )
    range_blocked_facts = CPBPlacementFacts.cpb_overlap_placement_facts(
        collection;
        transform_carries = (left_transform, right_transform),
        placement_ranges = (blocked_range,),
        placement_plan = (; kind = :test_overlap_placement_plan),
        accumulation_rule = :test_accumulation_rule,
    )
    range_blocked_summary = CPBPlacementFacts.summary(range_blocked_facts)
    @test range_blocked_summary.status ===
          :blocked_cpb_overlap_placement_facts
    @test range_blocked_summary.blocker ===
          :left_column_range_dimension_mismatch
    @test only(range_blocked_summary.record_fact_summaries).placement_range_status ===
          :blocked_cpb_source_pair_placement_range
    @test only(range_blocked_summary.record_fact_summaries).placement_range_blocker ===
          :left_column_range_dimension_mismatch
    @test range_blocked_summary.global_matrix_materialized === false
    @test range_blocked_summary.route_driver_wiring === false

    range_blocked_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            range_blocked_facts,
        )
    @test range_blocked_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test range_blocked_skeleton.blocker ===
          :left_column_range_dimension_mismatch
    @test only(range_blocked_skeleton.record_placement_summaries).left_column_range_status ===
          :blocked_left_column_range
    @test range_blocked_skeleton.global_overlap_status === :blocked
    @test range_blocked_skeleton.global_matrix_materialized === false
    @test range_blocked_skeleton.route_driver_wiring === false

    blocked_left_transform = _placement_facts_transform_carry(;
        side = :left,
        source_shape = (x = 1, y = 3, z = 1),
    )
    transform_blocked_facts = CPBPlacementFacts.cpb_overlap_placement_facts(
        collection;
        transform_carries = (
            blocked_left_transform,
            _placement_facts_transform_carry(; side = :right),
        ),
        placement_ranges = (_placement_facts_range(),),
        placement_plan = (; kind = :test_overlap_placement_plan),
        accumulation_rule = :test_accumulation_rule,
    )
    transform_blocked_summary = CPBPlacementFacts.summary(transform_blocked_facts)
    @test transform_blocked_summary.status ===
          :blocked_cpb_overlap_placement_facts
    @test transform_blocked_summary.blocker ===
          :retained_transform_source_shape_mismatch
    @test only(transform_blocked_summary.record_fact_summaries).left_transform_status ===
          :blocked_cpb_retained_transform_carry
    @test only(transform_blocked_summary.record_fact_summaries).left_transform_blocker ===
          :retained_transform_source_shape_mismatch
    @test transform_blocked_summary.global_matrix_materialized === false
    @test transform_blocked_summary.route_driver_wiring === false

    transform_blocked_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            transform_blocked_facts,
        )
    @test transform_blocked_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test transform_blocked_skeleton.blocker ===
          :retained_transform_source_shape_mismatch
    @test only(transform_blocked_skeleton.record_placement_summaries).retained_transform_status ===
          :blocked_retained_transform
    @test transform_blocked_skeleton.global_overlap_status === :blocked
    @test transform_blocked_skeleton.global_matrix_materialized === false
    @test transform_blocked_skeleton.route_driver_wiring === false

    blocked_collection = _placement_facts_collection(;
        record = _placement_facts_record(;
            status = :blocked_cpb_local_overlap_block_record,
            blocker = :test_local_overlap_blocker,
        ),
    )
    blocked_collection_facts =
        CPBPlacementFacts.cpb_overlap_placement_facts(blocked_collection)
    blocked_collection_summary =
        CPBPlacementFacts.summary(blocked_collection_facts)
    @test blocked_collection_summary.status ===
          :blocked_cpb_overlap_placement_facts
    @test blocked_collection_summary.blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_collection_summary.collection_available === false
    @test :missing_local_overlap_collection in
          blocked_collection_summary.missing_requirements
    @test blocked_collection_summary.global_matrix_materialized === false
    @test blocked_collection_summary.route_driver_wiring === false
end
