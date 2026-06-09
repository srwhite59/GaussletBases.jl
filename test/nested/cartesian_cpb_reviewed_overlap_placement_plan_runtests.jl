# Runtime role: metadata-only reviewed CPB overlap placement plan contract.
#
# This test validates reviewed placement-plan metadata only. It does not apply
# transforms, place local CPB overlap blocks, or assemble route/global overlap.

using Test
using GaussletBases

const CPBReviewedPlan = GaussletBases.CartesianCPBBlockProviders

const _REVIEWED_PLAN_BLOCK_KEY = (:tiny_left, :tiny_right)
const _REVIEWED_PLAN_SOURCE_SHAPE = (x = 1, y = 2, z = 1)

function _reviewed_plan(;
    placement_plan_kind = :reviewed_overlap_placement_plan,
    accumulation_rule = :add_explicit_blocks_into_ranges,
    symmetry_policy = :explicit_blocks_only,
    duplicate_record_policy = :reject_duplicate_block_keys,
    accepted_block_keys = (_REVIEWED_PLAN_BLOCK_KEY,),
    required_global_dimension_source = :test_retained_layout,
    local_ordering_contract = :parent_compatible_x_slowest_z_fastest,
)
    return CPBReviewedPlan.cpb_reviewed_overlap_placement_plan(;
        placement_plan_kind,
        accumulation_rule,
        symmetry_policy,
        duplicate_record_policy,
        accepted_block_keys,
        required_global_dimension_source,
        local_ordering_contract,
    )
end

function _reviewed_plan_record()
    metadata = (;
        object_kind = :test_cpb_local_overlap_block_record_summary,
        term = :overlap,
        block_key = _REVIEWED_PLAN_BLOCK_KEY,
        source_kind = :cpb_overlap_dense_block,
        status = :available_cpb_local_overlap_block_record,
        blocker = nothing,
        dense_block_available = true,
        dense_block_shape = (2, 2),
        left_cpb_summary = (;
            object_kind = :test_left_cpb_summary,
            shape = _REVIEWED_PLAN_SOURCE_SHAPE,
        ),
        right_cpb_summary = (;
            object_kind = :test_right_cpb_summary,
            shape = _REVIEWED_PLAN_SOURCE_SHAPE,
        ),
        local_ordering = :parent_compatible_x_slowest_z_fastest,
        placement_status = :unassigned,
        retained_transform_status = :unassigned,
        global_matrix_materialized = false,
        route_driver_wiring = false,
    )
    return CPBReviewedPlan.CPBLocalOverlapBlockRecord(
        _REVIEWED_PLAN_BLOCK_KEY,
        nothing,
        metadata,
    )
end

function _reviewed_plan_collection()
    return CPBReviewedPlan.cpb_local_overlap_block_collection((
        _reviewed_plan_record(),
    ))
end

function _reviewed_plan_transform_carry(side)
    return CPBReviewedPlan.cpb_retained_transform_carry(
        side,
        _REVIEWED_PLAN_BLOCK_KEY,
        (; object_kind = :test_cpb_source_summary, shape = _REVIEWED_PLAN_SOURCE_SHAPE),
        _REVIEWED_PLAN_SOURCE_SHAPE,
        :parent_compatible_x_slowest_z_fastest,
        side === :left ? (1:2) : (3:4);
        transform_object = [1.0 0.0; 0.0 1.0],
        transform_convention = :test_local_to_retained_columns,
        transform_provenance = :test_fixture,
    )
end

function _reviewed_plan_placement_range(left_transform, right_transform)
    return CPBReviewedPlan.cpb_source_pair_placement_range(
        _REVIEWED_PLAN_BLOCK_KEY;
        left_column_range = 1:2,
        right_column_range = 3:4,
        global_dimension = 4,
        global_dimension_source = :test_retained_layout,
        range_source = :test_source_pair_ranges,
        range_provenance = :test_fixture,
        left_transform_carry = left_transform,
        right_transform_carry = right_transform,
    )
end

@testset "CPB reviewed overlap placement plan metadata" begin
    available = _reviewed_plan()
    available_summary = CPBReviewedPlan.summary(available)
    @test available_summary.object_kind ===
          :cartesian_cpb_reviewed_overlap_placement_plan_summary
    @test available_summary.status ===
          :available_cpb_reviewed_overlap_placement_plan
    @test available_summary.blocker === nothing
    @test available_summary.placement_plan_kind ===
          :reviewed_overlap_placement_plan
    @test available_summary.accumulation_rule ===
          :add_explicit_blocks_into_ranges
    @test available_summary.symmetry_policy === :explicit_blocks_only
    @test available_summary.duplicate_record_policy ===
          :reject_duplicate_block_keys
    @test available_summary.accepted_block_keys === (_REVIEWED_PLAN_BLOCK_KEY,)
    @test available_summary.accepted_record_count == 1
    @test available_summary.required_global_dimension_source ===
          :test_retained_layout
    @test available_summary.local_ordering_contract ===
          :parent_compatible_x_slowest_z_fastest
    @test available_summary.transform_application_implemented === false
    @test available_summary.numerical_placement_implemented === false
    @test available_summary.global_matrix_materialized === false
    @test available_summary.global_overlap_matrix_materialized === false
    @test available_summary.route_driver_wiring === false
    @test !hasproperty(available_summary, :dense_block)

    missing_accumulation =
        _reviewed_plan(; accumulation_rule = nothing)
    @test CPBReviewedPlan.summary(missing_accumulation).status ===
          :blocked_cpb_reviewed_overlap_placement_plan
    @test CPBReviewedPlan.summary(missing_accumulation).blocker ===
          :missing_accumulation_rule

    empty_inventory = _reviewed_plan(; accepted_block_keys = ())
    @test CPBReviewedPlan.summary(empty_inventory).status ===
          :blocked_cpb_reviewed_overlap_placement_plan
    @test CPBReviewedPlan.summary(empty_inventory).blocker ===
          :missing_accepted_record_inventory

    missing_inventory = _reviewed_plan(; accepted_block_keys = nothing)
    @test CPBReviewedPlan.summary(missing_inventory).status ===
          :blocked_cpb_reviewed_overlap_placement_plan
    @test CPBReviewedPlan.summary(missing_inventory).blocker ===
          :missing_accepted_record_inventory

    unsupported_symmetry =
        _reviewed_plan(; symmetry_policy = :implicit_transpose_fill)
    @test CPBReviewedPlan.summary(unsupported_symmetry).status ===
          :blocked_cpb_reviewed_overlap_placement_plan
    @test CPBReviewedPlan.summary(unsupported_symmetry).blocker ===
          :unsupported_overlap_symmetry_policy

    unsupported_duplicate =
        _reviewed_plan(; duplicate_record_policy = :sum_duplicate_block_keys)
    @test CPBReviewedPlan.summary(unsupported_duplicate).status ===
          :blocked_cpb_reviewed_overlap_placement_plan
    @test CPBReviewedPlan.summary(unsupported_duplicate).blocker ===
          :unsupported_overlap_duplicate_record_policy

    missing_dimension_source =
        _reviewed_plan(; required_global_dimension_source = :unavailable)
    @test CPBReviewedPlan.summary(missing_dimension_source).status ===
          :blocked_cpb_reviewed_overlap_placement_plan
    @test CPBReviewedPlan.summary(missing_dimension_source).blocker ===
          :missing_required_global_dimension_source

    left_transform = _reviewed_plan_transform_carry(:left)
    right_transform = _reviewed_plan_transform_carry(:right)
    placement_range =
        _reviewed_plan_placement_range(left_transform, right_transform)
    facts = CPBReviewedPlan.cpb_overlap_placement_facts(
        _reviewed_plan_collection();
        transform_carries = (left_transform, right_transform),
        placement_ranges = (placement_range,),
        placement_plan = available,
    )
    facts_summary = CPBReviewedPlan.summary(facts)
    @test facts_summary.status === :blocked_cpb_overlap_placement_facts
    @test facts_summary.blocker === :placement_not_implemented
    @test facts_summary.placement_plan_status === :available_placement_plan
    @test facts_summary.placement_plan_kind ===
          :reviewed_overlap_placement_plan
    @test facts_summary.accumulation_rule_status ===
          :available_accumulation_rule
    @test facts_summary.accumulation_rule ===
          :add_explicit_blocks_into_ranges
    @test facts_summary.missing_requirements === ()
    @test facts_summary.global_overlap_status === :blocked
    @test facts_summary.global_overlap_blocker === :placement_not_implemented
    @test facts_summary.placement_engine_implemented === false
    @test facts_summary.transform_application_implemented === false
    @test facts_summary.global_matrix_materialized === false
    @test facts_summary.route_driver_wiring === false
end
