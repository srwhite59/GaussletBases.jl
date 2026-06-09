# Runtime role: tiny private-overlap audit fingerprint for CPB local overlap.
#
# This test verifies that the private overlap side can recognize the new
# CPB-provider local overlap source as local structured data. It does not make
# the route-global overlap input-facts path consume that local source.

using Test
using GaussletBases

const CPBLocalOverlapFingerprint = GaussletBases.CartesianCPB
const CPGBLocalOverlapFingerprint = GaussletBases.CartesianParentGaussletBases
const CBPLocalOverlapFingerprint = GaussletBases.CartesianCPBBlockProviders

function _local_overlap_fingerprint_parent()
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 2,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBLocalOverlapFingerprint.CartesianParentGaussletBasis3D(axis)
end

function _local_overlap_fingerprint_overlap_1d()
    return (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    )
end

function _local_overlap_fingerprint_axis_bundle()
    overlap_1d = _local_overlap_fingerprint_overlap_1d()
    return (;
        x = (; pgdg_intermediate = (; overlap = overlap_1d.x)),
        y = (; pgdg_intermediate = (; overlap = overlap_1d.y)),
        z = (; pgdg_intermediate = (; overlap = overlap_1d.z)),
    )
end

function _local_overlap_fingerprint_source()
    parent = _local_overlap_fingerprint_parent()
    packet = CPGBLocalOverlapFingerprint.parent_overlap_axis_factor_packet(
        parent,
        _local_overlap_fingerprint_axis_bundle(),
    )
    left = CPBLocalOverlapFingerprint.cpb(
        1:1,
        1:2,
        1:1;
        role = :private_overlap_local_left,
    )
    right = CPBLocalOverlapFingerprint.cpb(
        1:1,
        1:2,
        1:1;
        role = :private_overlap_local_right,
    )
    interval_pair =
        CBPLocalOverlapFingerprint.cpb_interval_pair(parent, left, right)
    axis_blocks =
        CBPLocalOverlapFingerprint.cpb_overlap_axis_blocks(packet, interval_pair)
    dense_block =
        CBPLocalOverlapFingerprint.cpb_overlap_dense_block(axis_blocks)
    record = CBPLocalOverlapFingerprint.cpb_local_overlap_block_record(
        dense_block;
        block_key = (:tiny_left, :tiny_right),
    )
    collection =
        CBPLocalOverlapFingerprint.cpb_local_overlap_block_collection((record,))
    blocked_interval_pair = CBPLocalOverlapFingerprint.cpb_interval_pair(
        parent,
        CPBLocalOverlapFingerprint.cpb(
            0:1,
            1:2,
            1:1;
            role = :private_overlap_local_outside_left,
        ),
        right,
    )
    blocked_axis_blocks =
        CBPLocalOverlapFingerprint.cpb_overlap_axis_blocks(
            packet,
            blocked_interval_pair,
        )
    blocked_record =
        CBPLocalOverlapFingerprint.cpb_local_overlap_block_record(
            blocked_axis_blocks;
            block_key = (:outside_left, :tiny_right),
        )
    blocked_collection =
        CBPLocalOverlapFingerprint.cpb_local_overlap_block_collection((
            blocked_record,
        ))
    return (;
        parent,
        packet,
        interval_pair,
        axis_blocks,
        dense_block,
        collection,
        blocked_collection,
    )
end

@testset "Private overlap CPB local source fingerprint" begin
    source = _local_overlap_fingerprint_source()
    axis_block_summary =
        CBPLocalOverlapFingerprint.summary(source.axis_blocks)
    dense_summary =
        CBPLocalOverlapFingerprint.summary(source.dense_block)
    fingerprint =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_local_source_fingerprint(
            source.dense_block,
        )

    @test fingerprint.status ===
          :available_private_global_overlap_local_source_fingerprint
    @test fingerprint.blocker === nothing
    @test fingerprint.source_kind === :cpb_provider_local_overlap
    @test fingerprint.axis_blocks_available
    @test fingerprint.dense_local_block_materialized
    @test fingerprint.global_matrix_materialized === false
    @test fingerprint.global_overlap_matrix_materialized === false
    @test fingerprint.route_driver_wiring === false
    @test fingerprint.private_global_overlap_input_facts_available === false
    @test fingerprint.route_global_overlap_stage_source === false
    @test fingerprint.factor_space === :parent_axis_bundle_pgdg_intermediate
    @test fingerprint.factor_convention === :axis_bundle_one_body_overlap
    @test fingerprint.normalization_convention ===
          :not_separate_from_axis_bundle_one_body_overlap
    @test fingerprint.index_domain === :parent_axis_indices
    @test fingerprint.index_domain_source === :axis_bundle_contract
    @test fingerprint.index_domain_status ===
          :assumed_parent_axis_indexed_by_current_axis_bundle_contract
    @test fingerprint.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test fingerprint.left_shape == (x = 1, y = 2, z = 1)
    @test fingerprint.right_shape == (x = 1, y = 2, z = 1)
    @test fingerprint.dense_block_shape == (2, 2)
    @test fingerprint.dense_block_shape == dense_summary.dense_block_shape

    axis_fingerprint =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_local_source_fingerprint(
            source.axis_blocks,
        )
    @test axis_fingerprint.status ===
          :available_private_global_overlap_local_source_fingerprint
    @test axis_fingerprint.axis_blocks_available
    @test axis_fingerprint.dense_local_block_materialized === false
    @test axis_fingerprint.dense_block_shape === :not_materialized
    @test axis_fingerprint.left_shape ==
          axis_block_summary.interval_pair_summary.left_shape

    collection_adapter =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_local_collection_adapter(
            source.collection,
        )
    @test collection_adapter.object_kind ===
          :cartesian_route_driver_private_global_overlap_local_collection_adapter
    @test collection_adapter.status ===
          :blocked_private_global_overlap_local_collection_adapter
    @test collection_adapter.blocker ===
          :missing_placement_or_retained_transform
    @test collection_adapter.local_cpb_overlap_collection_status ===
          :available_cpb_local_overlap_block_collection
    @test collection_adapter.local_cpb_overlap_collection_source ===
          :cpb_provider_local_overlap_collection
    @test collection_adapter.local_cpb_overlap_collection_available
    @test collection_adapter.record_count == 1
    @test collection_adapter.block_keys === ((:tiny_left, :tiny_right),)
    @test collection_adapter.dense_block_count == 1
    @test collection_adapter.placement_status === :unassigned
    @test collection_adapter.retained_transform_status === :unassigned
    @test collection_adapter.global_overlap_status === :blocked
    @test collection_adapter.global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test collection_adapter.global_matrix_materialized === false
    @test collection_adapter.global_overlap_matrix_materialized === false
    @test collection_adapter.route_driver_wiring === false
    @test collection_adapter.private_global_overlap_input_facts_available === false
    @test collection_adapter.route_global_overlap_stage_source === false

    placement_requirements =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_requirements_fingerprint(
            collection_adapter,
        )
    @test placement_requirements.object_kind ===
          :cartesian_route_driver_private_global_overlap_placement_requirements_fingerprint
    @test placement_requirements.status ===
          :blocked_private_global_overlap_placement_requirements
    @test placement_requirements.blocker ===
          :missing_placement_or_retained_transform
    @test placement_requirements.local_cpb_overlap_collection_available
    @test placement_requirements.placement_requirements_status ===
          :blocked_missing_placement_requirements
    @test placement_requirements.missing_requirements === (
        :missing_retained_transform,
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_global_dimension,
        :missing_placement_plan,
        :missing_accumulation_rule,
    )
    @test placement_requirements.retained_transform_status ===
          :missing_retained_transform
    @test placement_requirements.left_column_range_status ===
          :missing_left_column_range
    @test placement_requirements.right_column_range_status ===
          :missing_right_column_range
    @test placement_requirements.global_dimension_status ===
          :missing_global_dimension
    @test placement_requirements.placement_plan_status ===
          :missing_placement_plan
    @test placement_requirements.accumulation_rule_status ===
          :missing_accumulation_rule
    @test placement_requirements.global_overlap_status === :blocked
    @test placement_requirements.global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test placement_requirements.global_matrix_materialized === false
    @test placement_requirements.route_driver_wiring === false
    @test placement_requirements.private_global_overlap_input_facts_available === false
    @test placement_requirements.route_global_overlap_stage_source === false

    collection_placement_requirements =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_requirements_fingerprint(
            source.collection,
        )
    @test collection_placement_requirements == placement_requirements

    placement_candidate =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            placement_requirements,
        )
    @test placement_candidate.object_kind ===
          :cartesian_route_driver_private_global_overlap_placement_candidate
    @test placement_candidate.status ===
          :blocked_private_global_overlap_placement_candidate
    @test placement_candidate.blocker ===
          :missing_placement_or_retained_transform
    @test placement_candidate.local_cpb_overlap_collection_available
    @test placement_candidate.placement_candidate_status ===
          :blocked_missing_placement_requirements
    @test placement_candidate.available_requirements ===
          (:local_cpb_overlap_collection,)
    @test placement_candidate.missing_requirements === (
        :missing_retained_transform,
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_global_dimension,
        :missing_placement_plan,
        :missing_accumulation_rule,
    )
    @test placement_candidate.retained_transform_status ===
          :missing_retained_transform
    @test placement_candidate.left_column_range_status ===
          :missing_left_column_range
    @test placement_candidate.right_column_range_status ===
          :missing_right_column_range
    @test placement_candidate.global_dimension_status ===
          :missing_global_dimension
    @test placement_candidate.placement_plan_status ===
          :missing_placement_plan
    @test placement_candidate.accumulation_rule_status ===
          :missing_accumulation_rule
    @test placement_candidate.global_overlap_status === :blocked
    @test placement_candidate.global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test placement_candidate.global_matrix_materialized === false
    @test placement_candidate.global_overlap_matrix_materialized === false
    @test placement_candidate.route_driver_wiring === false
    @test placement_candidate.private_global_overlap_input_facts_available === false
    @test placement_candidate.route_global_overlap_stage_source === false

    collection_placement_candidate =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            source.collection,
        )
    @test collection_placement_candidate == placement_candidate

    partial_placement_candidate =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            source.collection;
            retained_transform = (; kind = :test_retained_transform),
            global_dimension = 2,
        )
    @test partial_placement_candidate.status ===
          :blocked_private_global_overlap_placement_candidate
    @test partial_placement_candidate.blocker ===
          :missing_placement_or_retained_transform
    @test partial_placement_candidate.local_cpb_overlap_collection_available
    @test partial_placement_candidate.placement_candidate_status ===
          :blocked_missing_placement_requirements
    @test partial_placement_candidate.available_requirements === (
        :local_cpb_overlap_collection,
        :retained_transform,
        :global_dimension,
    )
    @test partial_placement_candidate.missing_requirements === (
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_placement_plan,
        :missing_accumulation_rule,
    )
    @test partial_placement_candidate.retained_transform_status ===
          :available_retained_transform
    @test partial_placement_candidate.left_column_range_status ===
          :missing_left_column_range
    @test partial_placement_candidate.right_column_range_status ===
          :missing_right_column_range
    @test partial_placement_candidate.global_dimension_status ===
          :available_global_dimension
    @test partial_placement_candidate.placement_plan_status ===
          :missing_placement_plan
    @test partial_placement_candidate.accumulation_rule_status ===
          :missing_accumulation_rule
    @test partial_placement_candidate.global_overlap_status === :blocked
    @test partial_placement_candidate.global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test partial_placement_candidate.global_matrix_materialized === false
    @test partial_placement_candidate.global_overlap_matrix_materialized === false
    @test partial_placement_candidate.route_driver_wiring === false
    @test partial_placement_candidate.private_global_overlap_input_facts_available === false
    @test partial_placement_candidate.route_global_overlap_stage_source === false

    all_facts_placement_candidate =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            source.collection;
            retained_transform = (; kind = :test_retained_transform),
            left_column_ranges = (; left = 1:2),
            right_column_ranges = (; right = 1:2),
            global_dimension = 2,
            placement_plan = (; kind = :test_placement_plan),
            accumulation_rule = :test_accumulation_rule,
        )
    @test all_facts_placement_candidate.status ===
          :blocked_private_global_overlap_placement_candidate
    @test all_facts_placement_candidate.blocker === :placement_not_implemented
    @test all_facts_placement_candidate.local_cpb_overlap_collection_available
    @test all_facts_placement_candidate.placement_candidate_status ===
          :blocked_placement_not_implemented
    @test all_facts_placement_candidate.available_requirements === (
        :local_cpb_overlap_collection,
        :retained_transform,
        :left_column_range,
        :right_column_range,
        :global_dimension,
        :placement_plan,
        :accumulation_rule,
    )
    @test all_facts_placement_candidate.missing_requirements === ()
    @test all_facts_placement_candidate.retained_transform_status ===
          :available_retained_transform
    @test all_facts_placement_candidate.left_column_range_status ===
          :available_left_column_range
    @test all_facts_placement_candidate.right_column_range_status ===
          :available_right_column_range
    @test all_facts_placement_candidate.global_dimension_status ===
          :available_global_dimension
    @test all_facts_placement_candidate.placement_plan_status ===
          :available_placement_plan
    @test all_facts_placement_candidate.accumulation_rule_status ===
          :available_accumulation_rule
    @test all_facts_placement_candidate.global_overlap_status === :blocked
    @test all_facts_placement_candidate.global_overlap_blocker ===
          :placement_not_implemented
    @test all_facts_placement_candidate.global_matrix_materialized === false
    @test all_facts_placement_candidate.global_overlap_matrix_materialized === false
    @test all_facts_placement_candidate.route_driver_wiring === false
    @test all_facts_placement_candidate.private_global_overlap_input_facts_available === false
    @test all_facts_placement_candidate.route_global_overlap_stage_source === false

    placement_plan_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            source.collection,
        )
    @test placement_plan_skeleton.object_kind ===
          :cartesian_route_driver_private_global_overlap_placement_plan_skeleton
    @test placement_plan_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test placement_plan_skeleton.blocker ===
          :missing_placement_or_retained_transform
    @test placement_plan_skeleton.local_cpb_overlap_collection_available
    @test placement_plan_skeleton.record_count == 1
    @test placement_plan_skeleton.block_keys === ((:tiny_left, :tiny_right),)
    @test placement_plan_skeleton.placement_plan_status ===
          :missing_placement_plan
    @test placement_plan_skeleton.placement_plan_kind === :unavailable
    @test placement_plan_skeleton.retained_transform_status ===
          :missing_retained_transform
    @test placement_plan_skeleton.left_column_range_status ===
          :missing_left_column_range
    @test placement_plan_skeleton.right_column_range_status ===
          :missing_right_column_range
    @test placement_plan_skeleton.global_dimension_status ===
          :missing_global_dimension
    @test placement_plan_skeleton.global_dimension === nothing
    @test placement_plan_skeleton.global_dimension_source === :unavailable
    @test placement_plan_skeleton.accumulation_rule_status ===
          :missing_accumulation_rule
    @test placement_plan_skeleton.accumulation_rule === :unavailable
    @test placement_plan_skeleton.available_requirements ===
          (:local_cpb_overlap_collection,)
    @test placement_plan_skeleton.missing_requirements ===
          placement_candidate.missing_requirements
    @test length(placement_plan_skeleton.record_placement_summaries) == 1
    record_placement_summary =
        only(placement_plan_skeleton.record_placement_summaries)
    @test record_placement_summary.block_key === (:tiny_left, :tiny_right)
    @test record_placement_summary.source_kind === :cpb_overlap_dense_block
    @test record_placement_summary.dense_block_available
    @test record_placement_summary.dense_block_shape == (2, 2)
    @test record_placement_summary.left_cpb_summary.shape ==
          (x = 1, y = 2, z = 1)
    @test record_placement_summary.right_cpb_summary.shape ==
          (x = 1, y = 2, z = 1)
    @test record_placement_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test record_placement_summary.placement_status === :unassigned
    @test record_placement_summary.retained_transform_status ===
          :missing_retained_transform
    @test record_placement_summary.left_column_range_status ===
          :missing_left_column_range
    @test record_placement_summary.right_column_range_status ===
          :missing_right_column_range
    @test placement_plan_skeleton.global_overlap_status === :blocked
    @test placement_plan_skeleton.global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test placement_plan_skeleton.global_matrix_materialized === false
    @test placement_plan_skeleton.global_overlap_matrix_materialized === false
    @test placement_plan_skeleton.route_driver_wiring === false
    @test placement_plan_skeleton.private_global_overlap_input_facts_available === false
    @test placement_plan_skeleton.route_global_overlap_stage_source === false

    partial_placement_plan_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            source.collection;
            retained_transform = (; kind = :test_retained_transform),
            global_dimension = 2,
            global_dimension_source = :test_global_dimension,
        )
    @test partial_placement_plan_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test partial_placement_plan_skeleton.blocker ===
          :missing_placement_or_retained_transform
    @test partial_placement_plan_skeleton.available_requirements ===
          partial_placement_candidate.available_requirements
    @test partial_placement_plan_skeleton.missing_requirements ===
          partial_placement_candidate.missing_requirements
    @test partial_placement_plan_skeleton.retained_transform_status ===
          :available_retained_transform
    @test partial_placement_plan_skeleton.global_dimension_status ===
          :available_global_dimension
    @test partial_placement_plan_skeleton.global_dimension == 2
    @test partial_placement_plan_skeleton.global_dimension_source ===
          :test_global_dimension
    @test only(
        partial_placement_plan_skeleton.record_placement_summaries,
    ).retained_transform_status === :available_retained_transform
    @test partial_placement_plan_skeleton.global_overlap_status === :blocked
    @test partial_placement_plan_skeleton.route_driver_wiring === false
    @test partial_placement_plan_skeleton.global_matrix_materialized === false

    all_facts_placement_plan_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            source.collection;
            retained_transform = (; kind = :test_retained_transform),
            left_column_ranges = (; left = 1:2),
            right_column_ranges = (; right = 1:2),
            global_dimension = 2,
            global_dimension_source = :test_global_dimension,
            placement_plan = (; kind = :test_placement_plan),
            accumulation_rule = :test_accumulation_rule,
        )
    @test all_facts_placement_plan_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test all_facts_placement_plan_skeleton.blocker ===
          :placement_not_implemented
    @test all_facts_placement_plan_skeleton.placement_plan_status ===
          :available_placement_plan
    @test all_facts_placement_plan_skeleton.placement_plan_kind ===
          :test_placement_plan
    @test all_facts_placement_plan_skeleton.accumulation_rule_status ===
          :available_accumulation_rule
    @test all_facts_placement_plan_skeleton.accumulation_rule ===
          :test_accumulation_rule
    @test all_facts_placement_plan_skeleton.available_requirements ===
          all_facts_placement_candidate.available_requirements
    @test all_facts_placement_plan_skeleton.missing_requirements === ()
    @test all_facts_placement_plan_skeleton.global_overlap_status === :blocked
    @test all_facts_placement_plan_skeleton.global_overlap_blocker ===
          :placement_not_implemented
    @test all_facts_placement_plan_skeleton.global_matrix_materialized === false
    @test all_facts_placement_plan_skeleton.route_driver_wiring === false
    @test all_facts_placement_plan_skeleton.route_global_overlap_stage_source === false

    reviewed_left_transform =
        CBPLocalOverlapFingerprint.cpb_retained_transform_carry(
            :left,
            (:tiny_left, :tiny_right),
            (; object_kind = :test_left_cpb_summary, shape = (x = 1, y = 2, z = 1)),
            (x = 1, y = 2, z = 1),
            :parent_compatible_x_slowest_z_fastest,
            1:2;
            transform_object = [1.0 0.0; 0.0 1.0],
            transform_convention = :test_local_to_retained_columns,
            transform_provenance = :test_fixture,
        )
    reviewed_right_transform =
        CBPLocalOverlapFingerprint.cpb_retained_transform_carry(
            :right,
            (:tiny_left, :tiny_right),
            (; object_kind = :test_right_cpb_summary, shape = (x = 1, y = 2, z = 1)),
            (x = 1, y = 2, z = 1),
            :parent_compatible_x_slowest_z_fastest,
            1:2;
            transform_object = [1.0 0.0; 0.0 1.0],
            transform_convention = :test_local_to_retained_columns,
            transform_provenance = :test_fixture,
        )
    reviewed_range =
        CBPLocalOverlapFingerprint.cpb_source_pair_placement_range(
            (:tiny_left, :tiny_right);
            left_column_range = 1:2,
            right_column_range = 1:2,
            global_dimension = 2,
            global_dimension_source = :test_retained_layout,
            range_source = :test_source_pair_ranges,
            range_provenance = :test_fixture,
            left_transform_carry = reviewed_left_transform,
            right_transform_carry = reviewed_right_transform,
        )
    reviewed_plan =
        CBPLocalOverlapFingerprint.cpb_reviewed_overlap_placement_plan(;
            placement_plan_kind = :test_reviewed_overlap_placement_plan,
            accumulation_rule = :add_explicit_blocks_into_ranges,
            accepted_block_keys = ((:tiny_left, :tiny_right),),
            required_global_dimension_source = :test_retained_layout,
        )
    reviewed_plan_facts =
        CBPLocalOverlapFingerprint.cpb_overlap_placement_facts(
            source.collection;
            transform_carries = (reviewed_left_transform, reviewed_right_transform),
            placement_ranges = (reviewed_range,),
            placement_plan = reviewed_plan,
        )
    reviewed_plan_facts_summary =
        CBPLocalOverlapFingerprint.summary(reviewed_plan_facts)
    @test reviewed_plan_facts_summary.status ===
          :blocked_cpb_overlap_placement_facts
    @test reviewed_plan_facts_summary.blocker === :placement_not_implemented
    @test reviewed_plan_facts_summary.placement_plan_status ===
          :available_placement_plan
    @test reviewed_plan_facts_summary.placement_plan_kind ===
          :test_reviewed_overlap_placement_plan
    @test reviewed_plan_facts_summary.accumulation_rule_status ===
          :available_accumulation_rule
    @test reviewed_plan_facts_summary.accumulation_rule ===
          :add_explicit_blocks_into_ranges
    @test reviewed_plan_facts_summary.missing_requirements === ()

    reviewed_plan_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            reviewed_plan_facts,
        )
    @test reviewed_plan_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test reviewed_plan_skeleton.blocker === :placement_not_implemented
    @test reviewed_plan_skeleton.placement_plan_status ===
          :available_placement_plan
    @test reviewed_plan_skeleton.placement_plan_kind ===
          :test_reviewed_overlap_placement_plan
    @test reviewed_plan_skeleton.accumulation_rule_status ===
          :available_accumulation_rule
    @test reviewed_plan_skeleton.accumulation_rule ===
          :add_explicit_blocks_into_ranges
    @test reviewed_plan_skeleton.global_overlap_status === :blocked
    @test reviewed_plan_skeleton.global_overlap_blocker ===
          :placement_not_implemented
    @test reviewed_plan_skeleton.route_driver_wiring === false
    @test reviewed_plan_skeleton.transform_application_implemented === false
    @test reviewed_plan_skeleton.global_matrix_materialized === false
    @test reviewed_plan_skeleton.route_global_overlap_stage_source === false

    blocked_collection_adapter =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_local_collection_adapter(
            source.blocked_collection,
        )
    @test blocked_collection_adapter.status ===
          :blocked_private_global_overlap_local_collection_adapter
    @test blocked_collection_adapter.blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_collection_adapter.local_cpb_overlap_collection_status ===
          :blocked_cpb_local_overlap_block_collection
    @test blocked_collection_adapter.local_cpb_overlap_collection_available === false
    @test blocked_collection_adapter.record_count == 1
    @test blocked_collection_adapter.block_keys === ((:outside_left, :tiny_right),)
    @test blocked_collection_adapter.dense_block_count == 0
    @test blocked_collection_adapter.global_overlap_status === :blocked
    @test blocked_collection_adapter.global_overlap_blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_collection_adapter.global_matrix_materialized === false
    @test blocked_collection_adapter.route_driver_wiring === false

    blocked_placement_requirements =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_requirements_fingerprint(
            blocked_collection_adapter,
        )
    @test blocked_placement_requirements.status ===
          :blocked_private_global_overlap_placement_requirements
    @test blocked_placement_requirements.blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_placement_requirements.local_cpb_overlap_collection_available === false
    @test blocked_placement_requirements.placement_requirements_status ===
          :blocked_missing_local_overlap_collection
    @test blocked_placement_requirements.missing_requirements ===
          (:missing_local_overlap_collection,)
    @test blocked_placement_requirements.retained_transform_status === :unavailable
    @test blocked_placement_requirements.global_overlap_status === :blocked
    @test blocked_placement_requirements.global_overlap_blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_placement_requirements.global_matrix_materialized === false
    @test blocked_placement_requirements.route_driver_wiring === false

    blocked_placement_candidate =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            blocked_placement_requirements,
        )
    @test blocked_placement_candidate.status ===
          :blocked_private_global_overlap_placement_candidate
    @test blocked_placement_candidate.blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_placement_candidate.local_cpb_overlap_collection_available === false
    @test blocked_placement_candidate.placement_candidate_status ===
          :blocked_missing_local_overlap_collection
    @test blocked_placement_candidate.available_requirements === ()
    @test blocked_placement_candidate.missing_requirements ===
          (:missing_local_overlap_collection,)
    @test blocked_placement_candidate.retained_transform_status === :unavailable
    @test blocked_placement_candidate.global_overlap_status === :blocked
    @test blocked_placement_candidate.global_overlap_blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_placement_candidate.global_matrix_materialized === false
    @test blocked_placement_candidate.route_driver_wiring === false

    blocked_placement_plan_skeleton =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            source.blocked_collection,
        )
    @test blocked_placement_plan_skeleton.status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test blocked_placement_plan_skeleton.blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_placement_plan_skeleton.local_cpb_overlap_collection_available === false
    @test blocked_placement_plan_skeleton.record_count == 1
    @test blocked_placement_plan_skeleton.block_keys ===
          ((:outside_left, :tiny_right),)
    @test blocked_placement_plan_skeleton.available_requirements === ()
    @test blocked_placement_plan_skeleton.missing_requirements ===
          (:missing_local_overlap_collection,)
    @test length(blocked_placement_plan_skeleton.record_placement_summaries) == 1
    blocked_record_placement_summary =
        only(blocked_placement_plan_skeleton.record_placement_summaries)
    @test blocked_record_placement_summary.block_key ===
          (:outside_left, :tiny_right)
    @test blocked_record_placement_summary.source_kind ===
          :cpb_overlap_axis_block_set
    @test blocked_record_placement_summary.dense_block_available === false
    @test blocked_record_placement_summary.dense_block_shape === :not_materialized
    @test blocked_record_placement_summary.placement_status === :unassigned
    @test blocked_record_placement_summary.retained_transform_status === :unavailable
    @test blocked_record_placement_summary.left_column_range_status === :unavailable
    @test blocked_record_placement_summary.right_column_range_status === :unavailable
    @test blocked_placement_plan_skeleton.retained_transform_status === :unavailable
    @test blocked_placement_plan_skeleton.global_overlap_status === :blocked
    @test blocked_placement_plan_skeleton.global_overlap_blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_placement_plan_skeleton.global_matrix_materialized === false
    @test blocked_placement_plan_skeleton.route_driver_wiring === false

    facts_report = (;
        retained_dimension = 2,
        parent_axis_counts = (2, 2, 2),
        cpb_provider_local_overlap = source.dense_block,
    )
    facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            facts_report,
        )
    @test facts.status === :blocked_private_global_overlap_input_facts
    @test facts.blocker === :missing_parent_axis_bundle_overlap_factors
    @test facts.input_source === :unavailable
    @test facts.overlap_1d === nothing
    @test facts.global_dimension == 2
    @test facts.factor_space === :unavailable
    @test facts.factor_convention === :unavailable
end
