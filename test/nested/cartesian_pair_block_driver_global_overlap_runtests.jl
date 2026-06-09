# Runtime role: tiny smoke / private driver-facing global overlap hook test.
#
# This verifies the overlap-only driver-facing hook over structured
# pair-block materialization state. It does not cover Hamiltonians, Coulomb,
# IDA/MWG, PQS Lowdin/projection, exports, route-driver wiring, or
# White-Lindsey oracle fixtures.

using Test
using GaussletBases

const CPBMDriverOverlap = GaussletBases.CartesianPairBlockMaterialization
const CPBDriverOverlap = GaussletBases.CartesianCPB
const CPBProviderDriverOverlap = GaussletBases.CartesianCPBBlockProviders
const CPGBDriverOverlap = GaussletBases.CartesianParentGaussletBases
const CTLDriverOverlap = GaussletBases.CartesianTerminalLowering
const CRUDriverOverlap = GaussletBases.CartesianRetainedUnits
const CRTCDriverOverlap =
    GaussletBases.CartesianRetainedUnitTransformContracts
const CUPDriverOverlap = GaussletBases.CartesianUnitPairs
const CPOPDriverOverlap = GaussletBases.CartesianPairOperatorPlans

function _driver_overlap_pair_operator_plan()
    lowering_plan = CTLDriverOverlap.TerminalLoweringPlan(
        CTLDriverOverlap.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
    retained_plan = CRUDriverOverlap.RetainedUnitPlan(
        CRUDriverOverlap.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
    unit_pair_plan = CUPDriverOverlap.UnitPairPlan(
        CUPDriverOverlap.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
    transform_plan = CRTCDriverOverlap.RetainedUnitTransformContractPlan(
        CRTCDriverOverlap.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
    return CPOPDriverOverlap.PairOperatorPlan(
        CPOPDriverOverlap.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
end

function _driver_overlap_plan()
    direct_source = CPBDriverOverlap.cpb(
        1:1,
        1:2,
        1:1;
        role = :driver_global_overlap_direct_source,
    )
    record = CPBMDriverOverlap.PairBlockMaterializationRecord(
        (:direct_diag, :direct_diag),
        1,
        :direct_direct,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
        (:overlap,),
        :direct_direct_pair_block_materialization_pilot,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        (;
            left_source_cpbs = (direct_source,),
            right_source_cpbs = (direct_source,),
            pair_index = 1,
            selector_family = :direct_direct,
            left_column_range = 1:2,
            right_column_range = 1:2,
        ),
    )
    return CPBMDriverOverlap.PairBlockMaterializationPlan(
        CPBMDriverOverlap.MetadataOnlyPairBlockMaterialization(),
        _driver_overlap_pair_operator_plan(),
        (record,),
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :driver_global_overlap_hook),
    )
end

function _driver_overlap_inputs()
    return (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = (;
            x = [1.0 0.2; 0.2 1.1],
            y = [1.2 0.3; 0.3 1.3],
            z = [1.4 0.4; 0.4 1.5],
        ),
    )
end

function _driver_overlap_expected_matrix()
    overlap = _driver_overlap_inputs().overlap_1d
    scale = overlap.x[1, 1] * overlap.z[1, 1]
    return scale .* overlap.y
end

function _driver_overlap_stage_report(plan = _driver_overlap_plan())
    return (;
        low_order_route_summary = (;
            terminal_route_state = (; pair_block_materialization_plan = plan),
        ),
    )
end

function _driver_overlap_axis_bundle_object()
    overlap = _driver_overlap_inputs().overlap_1d
    return (;
        x = (; pgdg_intermediate = (; overlap = overlap.x)),
        y = (; pgdg_intermediate = (; overlap = overlap.y)),
        z = (; pgdg_intermediate = (; overlap = overlap.z)),
    )
end

function _driver_overlap_parent_object()
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 2,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBDriverOverlap.CartesianParentGaussletBasis3D(axis)
end

function _driver_overlap_facts_report(plan = _driver_overlap_plan())
    return merge(
        _driver_overlap_stage_report(plan),
        (;
            retained_dimension = 2,
            parent_axis_counts = (2, 2, 2),
            parent_axis_bundle_object = _driver_overlap_axis_bundle_object(),
        ),
    )
end

function _driver_overlap_example_overrides()
    example_module = Module(:PrivateGlobalOverlapOptionExample)
    example_path =
        joinpath(@__DIR__, "..", "..", "examples", "private_global_overlap_option.jl")
    Base.include(example_module, example_path)
    return (;
        private_global_overlap_requested =
            Core.eval(example_module, :private_global_overlap_requested),
        private_global_overlap_global_dimension =
            Core.eval(example_module, :private_global_overlap_global_dimension),
        private_global_overlap_inputs =
            Core.eval(example_module, :private_global_overlap_inputs),
    )
end

function _driver_overlap_real_report_fingerprint()
    return GaussletBases._pqs_source_box_route_driver_dry_run(
        route_family = :pqs_source_box,
        route_kind = :overlap_facts_fingerprint_probe,
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0,
        parent_axis_counts = (x = 5, y = 5, z = 5),
        map_backend = :pgdg_localized_experimental,
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = nothing,
        probe_parent_axis_construction = false,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
        probe_raw_product_box_plans = false,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (:overlap,),
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities =
            (:support_row_oracle, :dense_parent_projection),
    )
end

function _driver_overlap_probe_enabled_real_report_fingerprint()
    return GaussletBases._pqs_source_box_route_driver_dry_run(
        route_family = :pqs_source_box,
        route_kind = :overlap_facts_probe_enabled_fingerprint,
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0,
        parent_axis_counts = (x = 5, y = 5, z = 5),
        map_backend = :pgdg_localized_experimental,
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = nothing,
        probe_parent_axis_construction = :auto,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
        probe_raw_product_box_plans = false,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (:overlap,),
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities =
            (:support_row_oracle, :dense_parent_projection),
    )
end

function _driver_overlap_source_box_cpb(report, source_key::Symbol)
    hasproperty(report, :source_boxes) || return nothing
    hasproperty(report.source_boxes, source_key) || return nothing
    source_box = getproperty(report.source_boxes, source_key)
    return CPBDriverOverlap.cpb(
        source_box.x,
        source_box.y,
        source_box.z;
        role = Symbol("real_report_", source_key, "_source_cpb"),
    )
end

function _driver_overlap_smallest_source_pair_key(report)
    hasproperty(report, :pair_entries) || return nothing
    hasproperty(report, :source_boxes) || return nothing
    candidates = filter(report.pair_entries) do entry
        left_key, right_key = entry.pair_key
        return hasproperty(report.source_boxes, left_key) &&
               hasproperty(report.source_boxes, right_key)
    end
    isempty(candidates) && return nothing
    return argmin(candidates) do entry
        left_key, right_key = entry.pair_key
        left_cpb = _driver_overlap_source_box_cpb(report, left_key)
        right_cpb = _driver_overlap_source_box_cpb(report, right_key)
        return CPBDriverOverlap.support_count(left_cpb) *
               CPBDriverOverlap.support_count(right_cpb)
    end
end

function _driver_overlap_nested_property(object, path::Tuple)
    current = object
    for key in path
        current =
            !isnothing(current) && hasproperty(current, key) ?
            getproperty(current, key) :
            nothing
        isnothing(current) && return nothing
    end
    return current
end

function _driver_overlap_first_report_source(report, candidates)
    for candidate in candidates
        value = _driver_overlap_nested_property(report, candidate.path)
        isnothing(value) || return (; source = candidate.source, value)
    end
    return (; source = :unavailable, value = nothing)
end

function _driver_overlap_recorded_range(value)
    isnothing(value) && return nothing
    if hasproperty(value, :column_range)
        return value.column_range
    elseif hasproperty(value, :retained_range)
        return value.retained_range
    end
    return nothing
end

function _driver_overlap_retained_units_dimension(units)
    isnothing(units) && return nothing
    ranges = UnitRange{Int}[]
    for unit in units
        range = _driver_overlap_recorded_range(unit)
        isnothing(range) && return nothing
        push!(ranges, range)
    end
    isempty(ranges) && return nothing
    return maximum(last, ranges)
end

function _driver_overlap_global_dimension_source_audit(report)
    retained_unit_sources = (
        (;
            source = :low_order_terminal_route_state_retained_units,
            path = (:low_order_route_summary, :terminal_route_state, :retained_units),
        ),
        (;
            source = :terminal_route_state_retained_units,
            path = (:terminal_route_state, :retained_units),
        ),
        (; source = :report_retained_units, path = (:retained_units,)),
    )
    for candidate in retained_unit_sources
        units = _driver_overlap_nested_property(report, candidate.path)
        dimension = _driver_overlap_retained_units_dimension(units)
        isnothing(dimension) || return (;
            status = :available_retained_layout_global_dimension,
            source = candidate.source,
            global_dimension = dimension,
        )
    end
    if hasproperty(report, :retained_dimension) && !isnothing(report.retained_dimension)
        return (;
            status = :available_compatibility_global_dimension,
            source = :report_retained_dimension,
            global_dimension = report.retained_dimension,
        )
    end
    return (;
        status = :missing_global_dimension_source,
        source = :unavailable,
        global_dimension = nothing,
    )
end

function _driver_overlap_pair_entry_column_range_source_status(pair_entry, side::Symbol)
    isnothing(pair_entry) && return (;
        status = :missing_structured_source_box_pair,
        source = :unavailable,
    )
    keys =
        side === :left ?
        (:left_column_range, :left_final_column_range, :left_retained_column_range) :
        (:right_column_range, :right_final_column_range, :right_retained_column_range)
    for key in keys
        if hasproperty(pair_entry, key) && !isnothing(getproperty(pair_entry, key))
            return (;
                status =
                    :ambiguous_pair_entry_column_range_without_reviewed_placement_plan,
                source = Symbol(:report_pair_entry_, key),
            )
        end
    end
    return (;
        status = :missing_source_pair_retained_column_range,
        source = :unavailable,
    )
end

function _driver_overlap_real_report_overlap_placement_source_audit(
    report,
    local_block_collection,
    pair_entry,
)
    retained_transform_source = _driver_overlap_first_report_source(
        report,
        (
            (;
                source = :route_materializer_payload_retained_transform,
                path = (:route_materializer_payload, :retained_transform),
            ),
            (;
                source = :route_materializer_payload_retained_transforms,
                path = (:route_materializer_payload, :retained_transforms),
            ),
            (;
                source = :terminal_route_state_retained_transform,
                path = (:low_order_route_summary, :terminal_route_state, :retained_transform),
            ),
            (;
                source = :terminal_route_state_retained_transforms,
                path = (:low_order_route_summary, :terminal_route_state, :retained_transforms),
            ),
        ),
    )
    left_column_range_source =
        _driver_overlap_pair_entry_column_range_source_status(pair_entry, :left)
    right_column_range_source =
        _driver_overlap_pair_entry_column_range_source_status(pair_entry, :right)
    global_dimension_source = _driver_overlap_global_dimension_source_audit(report)
    placement_plan_source = _driver_overlap_first_report_source(
        report,
        (
            (;
                source = :terminal_route_state_overlap_placement_plan,
                path = (:low_order_route_summary, :terminal_route_state, :overlap_placement_plan),
            ),
            (;
                source = :terminal_route_state_pair_block_materialization_plan,
                path = (:low_order_route_summary, :terminal_route_state, :pair_block_materialization_plan),
            ),
        ),
    )
    placement_plan_source_status =
        isnothing(placement_plan_source.value) ?
        :missing_reviewed_overlap_placement_plan :
        placement_plan_source.source === :terminal_route_state_pair_block_materialization_plan ?
        :ambiguous_pair_block_materialization_plan_not_overlap_placement_plan :
        :available_reviewed_overlap_placement_plan
    accumulation_rule_source = _driver_overlap_first_report_source(
        report,
        (
            (;
                source = :terminal_route_state_overlap_accumulation_rule,
                path = (:low_order_route_summary, :terminal_route_state, :overlap_accumulation_rule),
            ),
            (;
                source = :route_materializer_payload_overlap_accumulation_rule,
                path = (:route_materializer_payload, :overlap_accumulation_rule),
            ),
        ),
    )
    accumulation_rule_source_status =
        isnothing(accumulation_rule_source.value) ?
        :missing_overlap_accumulation_rule :
        :available_overlap_accumulation_rule
    retained_transform =
        isnothing(retained_transform_source.value) ?
        nothing :
        retained_transform_source.value
    global_dimension =
        global_dimension_source.status in (
            :available_retained_layout_global_dimension,
            :available_compatibility_global_dimension,
        ) ?
        global_dimension_source.global_dimension :
        nothing
    candidate =
        isnothing(local_block_collection) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            local_block_collection;
            retained_transform,
            global_dimension,
        )
    skeleton =
        isnothing(local_block_collection) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            local_block_collection;
            retained_transform,
            global_dimension,
            global_dimension_source = global_dimension_source.source,
        )
    placement_facts =
        isnothing(local_block_collection) ?
        nothing :
        CPBProviderDriverOverlap.cpb_overlap_placement_facts(
            local_block_collection;
            transform_carries = (),
            placement_ranges = (),
            placement_plan = nothing,
            accumulation_rule = nothing,
        )
    placement_facts_summary =
        isnothing(placement_facts) ?
        nothing :
        CPBProviderDriverOverlap.summary(placement_facts)
    placement_facts_record_summary =
        isnothing(placement_facts_summary) ||
        isempty(placement_facts_summary.record_fact_summaries) ?
        nothing :
        only(placement_facts_summary.record_fact_summaries)
    placement_facts_skeleton =
        isnothing(placement_facts) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            placement_facts,
        )
    placement_facts_skeleton_record_summary =
        isnothing(placement_facts_skeleton) ||
        isempty(placement_facts_skeleton.record_placement_summaries) ?
        nothing :
        only(placement_facts_skeleton.record_placement_summaries)
    reviewed_placement_plan =
        isnothing(local_block_collection) || isnothing(pair_entry) ?
        nothing :
        CPBProviderDriverOverlap.cpb_reviewed_overlap_placement_plan(;
            placement_plan_kind =
                :real_report_reviewed_overlap_placement_plan_fingerprint,
            accumulation_rule = :add_explicit_blocks_into_ranges,
            accepted_block_keys = (pair_entry.pair_key,),
            required_global_dimension_source = global_dimension_source.source,
        )
    reviewed_placement_facts =
        isnothing(reviewed_placement_plan) || isnothing(local_block_collection) ?
        nothing :
        CPBProviderDriverOverlap.cpb_overlap_placement_facts(
            local_block_collection;
            transform_carries = (),
            placement_ranges = (),
            placement_plan = reviewed_placement_plan,
        )
    reviewed_placement_facts_summary =
        isnothing(reviewed_placement_facts) ?
        nothing :
        CPBProviderDriverOverlap.summary(reviewed_placement_facts)
    reviewed_placement_facts_skeleton =
        isnothing(reviewed_placement_facts) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            reviewed_placement_facts,
        )
    return (;
        retained_transform_source_status =
            isnothing(retained_transform_source.value) ?
            :missing_actual_retained_transform_source :
            :available_actual_retained_transform_source,
        retained_transform_source = retained_transform_source.source,
        left_column_range_source_status = left_column_range_source.status,
        left_column_range_source = left_column_range_source.source,
        right_column_range_source_status = right_column_range_source.status,
        right_column_range_source = right_column_range_source.source,
        global_dimension_source_status = global_dimension_source.status,
        global_dimension_source = global_dimension_source.source,
        global_dimension = global_dimension_source.global_dimension,
        placement_plan_source_status,
        placement_plan_source = placement_plan_source.source,
        accumulation_rule_source_status,
        accumulation_rule_source = accumulation_rule_source.source,
        candidate_status_with_real_sources =
            isnothing(candidate) ?
            :not_attempted_missing_local_overlap_collection :
            candidate.status,
        candidate_blocker_with_real_sources =
            isnothing(candidate) ?
            :not_attempted_missing_local_overlap_collection :
            candidate.blocker,
        skeleton_status_with_real_sources =
            isnothing(skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            skeleton.status,
        skeleton_blocker_with_real_sources =
            isnothing(skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            skeleton.blocker,
        missing_requirements_with_real_sources =
            isnothing(skeleton) ? () : skeleton.missing_requirements,
        available_requirements_with_real_sources =
            isnothing(skeleton) ? () : skeleton.available_requirements,
        placement_facts_status =
            isnothing(placement_facts_summary) ?
            :not_attempted_missing_local_overlap_collection :
            placement_facts_summary.status,
        placement_facts_blocker =
            isnothing(placement_facts_summary) ?
            :not_attempted_missing_local_overlap_collection :
            placement_facts_summary.blocker,
        placement_facts_available_requirements =
            isnothing(placement_facts_summary) ?
            () :
            placement_facts_summary.available_requirements,
        placement_facts_missing_requirements =
            isnothing(placement_facts_summary) ?
            () :
            placement_facts_summary.missing_requirements,
        placement_facts_record_left_cpb_summary =
            isnothing(placement_facts_record_summary) ?
            :unavailable :
            placement_facts_record_summary.left_cpb_summary,
        placement_facts_record_right_cpb_summary =
            isnothing(placement_facts_record_summary) ?
            :unavailable :
            placement_facts_record_summary.right_cpb_summary,
        placement_facts_skeleton_status =
            isnothing(placement_facts_skeleton) ?
            :not_attempted_missing_placement_facts :
            placement_facts_skeleton.status,
        placement_facts_skeleton_blocker =
            isnothing(placement_facts_skeleton) ?
            :not_attempted_missing_placement_facts :
            placement_facts_skeleton.blocker,
        placement_facts_skeleton_missing_requirements =
            isnothing(placement_facts_skeleton) ?
            () :
            placement_facts_skeleton.missing_requirements,
        placement_facts_skeleton_global_overlap_status =
            isnothing(placement_facts_skeleton) ?
            :not_attempted_missing_placement_facts :
            placement_facts_skeleton.global_overlap_status,
        placement_facts_skeleton_route_driver_wiring =
            !isnothing(placement_facts_skeleton) &&
            placement_facts_skeleton.route_driver_wiring,
        placement_facts_skeleton_global_matrix_materialized =
            !isnothing(placement_facts_skeleton) &&
            placement_facts_skeleton.global_matrix_materialized,
        placement_facts_skeleton_private_input_facts_available =
            !isnothing(placement_facts_skeleton) &&
            placement_facts_skeleton.private_global_overlap_input_facts_available,
        placement_facts_skeleton_route_global_stage_source =
            !isnothing(placement_facts_skeleton) &&
            placement_facts_skeleton.route_global_overlap_stage_source,
        placement_facts_skeleton_record_left_cpb_summary =
            isnothing(placement_facts_skeleton_record_summary) ?
            :unavailable :
            placement_facts_skeleton_record_summary.left_cpb_summary,
        placement_facts_skeleton_record_right_cpb_summary =
            isnothing(placement_facts_skeleton_record_summary) ?
            :unavailable :
            placement_facts_skeleton_record_summary.right_cpb_summary,
        reviewed_placement_facts_status =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.status,
        reviewed_placement_facts_blocker =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.blocker,
        reviewed_placement_facts_available_requirements =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.available_requirements,
        reviewed_placement_facts_missing_requirements =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.missing_requirements,
        reviewed_placement_facts_placement_plan_status =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.placement_plan_status,
        reviewed_placement_facts_placement_plan_kind =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.placement_plan_kind,
        reviewed_placement_facts_accumulation_rule_status =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.accumulation_rule_status,
        reviewed_placement_facts_accumulation_rule =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.accumulation_rule,
        reviewed_placement_facts_inventory_status =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.placement_record_inventory_status,
        reviewed_placement_facts_inventory_blocker =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.placement_record_inventory_blocker,
        reviewed_placement_facts_accepted_block_keys =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.accepted_block_keys,
        reviewed_placement_facts_provided_block_keys =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.provided_block_keys,
        reviewed_placement_facts_rejected_block_keys =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.rejected_block_keys,
        reviewed_placement_facts_duplicate_block_keys =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.duplicate_block_keys,
        reviewed_placement_facts_duplicate_record_policy =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.duplicate_record_policy,
        reviewed_placement_facts_local_ordering_contract_status =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.local_ordering_contract_status,
        reviewed_placement_facts_local_ordering_contract_blocker =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.local_ordering_contract_blocker,
        reviewed_placement_facts_local_ordering_contract =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.local_ordering_contract,
        reviewed_placement_facts_provided_local_orderings =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.provided_local_orderings,
        reviewed_placement_facts_mismatched_local_ordering_block_keys =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.mismatched_local_ordering_block_keys,
        reviewed_placement_facts_global_dimension_source_contract_status =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.global_dimension_source_contract_status,
        reviewed_placement_facts_global_dimension_source_contract_blocker =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.global_dimension_source_contract_blocker,
        reviewed_placement_facts_required_global_dimension_source =
            isnothing(reviewed_placement_facts_summary) ?
            :not_attempted_missing_reviewed_placement_plan :
            reviewed_placement_facts_summary.required_global_dimension_source,
        reviewed_placement_facts_provided_global_dimension_sources =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.provided_global_dimension_sources,
        reviewed_placement_facts_mismatched_global_dimension_source_block_keys =
            isnothing(reviewed_placement_facts_summary) ?
            () :
            reviewed_placement_facts_summary.mismatched_global_dimension_source_block_keys,
        reviewed_placement_facts_skeleton_status =
            isnothing(reviewed_placement_facts_skeleton) ?
            :not_attempted_missing_reviewed_placement_facts :
            reviewed_placement_facts_skeleton.status,
        reviewed_placement_facts_skeleton_blocker =
            isnothing(reviewed_placement_facts_skeleton) ?
            :not_attempted_missing_reviewed_placement_facts :
            reviewed_placement_facts_skeleton.blocker,
        reviewed_placement_facts_skeleton_global_overlap_status =
            isnothing(reviewed_placement_facts_skeleton) ?
            :not_attempted_missing_reviewed_placement_facts :
            reviewed_placement_facts_skeleton.global_overlap_status,
        reviewed_placement_facts_skeleton_route_driver_wiring =
            !isnothing(reviewed_placement_facts_skeleton) &&
            reviewed_placement_facts_skeleton.route_driver_wiring,
        reviewed_placement_facts_skeleton_global_matrix_materialized =
            !isnothing(reviewed_placement_facts_skeleton) &&
            reviewed_placement_facts_skeleton.global_matrix_materialized,
        reviewed_placement_facts_skeleton_private_input_facts_available =
            !isnothing(reviewed_placement_facts_skeleton) &&
            reviewed_placement_facts_skeleton.private_global_overlap_input_facts_available,
        reviewed_placement_facts_skeleton_route_global_stage_source =
            !isnothing(reviewed_placement_facts_skeleton) &&
            reviewed_placement_facts_skeleton.route_global_overlap_stage_source,
        global_overlap_status =
            isnothing(skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            skeleton.global_overlap_status,
        global_overlap_blocker =
            isnothing(skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            skeleton.global_overlap_blocker,
        global_matrix_materialized = false,
        route_driver_wiring = false,
        route_global_overlap_stage_source = false,
    )
end

function _driver_overlap_real_report_local_cpb_provider_fingerprint(report)
    payload =
        hasproperty(report, :route_materializer_payload) ?
        report.route_materializer_payload :
        nothing
    parent_qw_basis =
        !isnothing(payload) && hasproperty(payload, :parent_qw_basis_object) ?
        payload.parent_qw_basis_object :
        nothing
    parent =
        isnothing(parent_qw_basis) ?
        nothing :
        CPGBDriverOverlap.cartesian_parent_gausslet_basis(parent_qw_basis)
    axis_bundle =
        !isnothing(payload) && hasproperty(payload, :parent_axis_bundle_object) ?
        payload.parent_axis_bundle_object :
        nothing
    packet =
        isnothing(parent) || isnothing(axis_bundle) ?
        nothing :
        CPGBDriverOverlap.parent_overlap_axis_factor_packet(parent, axis_bundle)
    packet_summary = isnothing(packet) ? nothing : CPGBDriverOverlap.summary(packet)
    pair_entry =
        isnothing(parent) ?
        nothing :
        _driver_overlap_smallest_source_pair_key(report)
    left_cpb =
        isnothing(pair_entry) ?
        nothing :
        _driver_overlap_source_box_cpb(report, pair_entry.pair_key[1])
    right_cpb =
        isnothing(pair_entry) ?
        nothing :
        _driver_overlap_source_box_cpb(report, pair_entry.pair_key[2])
    interval_pair =
        isnothing(parent) || isnothing(left_cpb) || isnothing(right_cpb) ?
        nothing :
        CPBProviderDriverOverlap.cpb_interval_pair(parent, left_cpb, right_cpb)
    interval_summary =
        isnothing(interval_pair) ?
        nothing :
        CPBProviderDriverOverlap.summary(interval_pair)
    axis_blocks =
        isnothing(packet) || isnothing(interval_pair) ?
        nothing :
        CPBProviderDriverOverlap.cpb_overlap_axis_blocks(packet, interval_pair)
    axis_block_summary =
        isnothing(axis_blocks) ?
        nothing :
        CPBProviderDriverOverlap.summary(axis_blocks)
    dense_block =
        !isnothing(axis_block_summary) &&
        axis_block_summary.status === :available_cpb_overlap_axis_blocks ?
        CPBProviderDriverOverlap.cpb_overlap_dense_block(axis_blocks) :
        nothing
    dense_summary =
        isnothing(dense_block) ?
        nothing :
        CPBProviderDriverOverlap.summary(dense_block)
    local_source_fingerprint =
        isnothing(dense_block) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_local_source_fingerprint(
            dense_block,
        )
    local_block_record =
        isnothing(dense_block) || isnothing(pair_entry) ?
        nothing :
        CPBProviderDriverOverlap.cpb_local_overlap_block_record(
            dense_block;
            block_key = pair_entry.pair_key,
        )
    local_block_record_summary =
        isnothing(local_block_record) ?
        nothing :
        CPBProviderDriverOverlap.summary(local_block_record)
    local_block_collection =
        isnothing(local_block_record) ?
        nothing :
        CPBProviderDriverOverlap.cpb_local_overlap_block_collection((
            local_block_record,
        ))
    local_block_collection_summary =
        isnothing(local_block_collection) ?
        nothing :
        CPBProviderDriverOverlap.summary(local_block_collection)
    local_block_collection_adapter =
        isnothing(local_block_collection) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_local_collection_adapter(
            local_block_collection,
        )
    placement_requirements =
        isnothing(local_block_collection_adapter) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_requirements_fingerprint(
            local_block_collection_adapter,
        )
    placement_candidate =
        isnothing(placement_requirements) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            placement_requirements,
        )
    placement_candidate_global_dimension =
        hasproperty(report, :retained_dimension) ? report.retained_dimension : 2
    partial_placement_candidate =
        isnothing(local_block_collection) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            local_block_collection;
            retained_transform = (; kind = :test_retained_transform),
            global_dimension = placement_candidate_global_dimension,
        )
    all_facts_placement_candidate =
        isnothing(local_block_collection) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_candidate(
            local_block_collection;
            retained_transform = (; kind = :test_retained_transform),
            left_column_ranges = (; product = 1:25),
            right_column_ranges = (; product = 1:25),
            global_dimension = placement_candidate_global_dimension,
            placement_plan = (; kind = :test_placement_plan),
            accumulation_rule = :test_accumulation_rule,
        )
    all_facts_placement_plan_skeleton =
        isnothing(local_block_collection) ?
        nothing :
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_placement_plan_skeleton(
            local_block_collection;
            retained_transform = (; kind = :test_retained_transform),
            left_column_ranges = (; product = 1:25),
            right_column_ranges = (; product = 1:25),
            global_dimension = placement_candidate_global_dimension,
            global_dimension_source = :report_retained_dimension,
            placement_plan = (; kind = :test_placement_plan),
            accumulation_rule = :test_accumulation_rule,
        )
    placement_source_audit =
        _driver_overlap_real_report_overlap_placement_source_audit(
            report,
            local_block_collection,
            pair_entry,
        )

    return (;
        parent_source_status =
            isnothing(parent) ?
            :missing_route_materializer_payload_parent_qw_basis_object :
            :available_route_materializer_payload_parent_qw_basis_object,
        parent_axis_counts =
            isnothing(parent) ?
            nothing :
            CPGBDriverOverlap.parent_axis_counts(parent),
        parent_overlap_packet_source_status =
            isnothing(packet_summary) ?
            :missing_parent_overlap_packet_source :
            packet_summary.status,
        parent_overlap_packet_blocker =
            isnothing(packet_summary) ? :missing_parent_overlap_packet_source :
            packet_summary.blocker,
        cpb_source_pair_status =
            isnothing(pair_entry) ?
            :missing_structured_source_box_pair :
            :available_report_pair_entry_source_boxes,
        cpb_source_pair_key =
            isnothing(pair_entry) ? nothing : pair_entry.pair_key,
        interval_pair_status =
            isnothing(interval_summary) ?
            :not_attempted_missing_cpb_source_pair :
            interval_summary.status,
        interval_pair_blocker =
            isnothing(interval_summary) ? :not_attempted_missing_cpb_source_pair :
            interval_summary.blocker,
        overlap_axis_blocks_status =
            isnothing(axis_block_summary) ?
            :not_attempted_missing_parent_packet_or_interval_pair :
            axis_block_summary.status,
        overlap_axis_blocks_blocker =
            isnothing(axis_block_summary) ?
            :not_attempted_missing_parent_packet_or_interval_pair :
            axis_block_summary.blocker,
        dense_local_overlap_status =
            isnothing(dense_summary) ?
            :not_attempted_missing_available_axis_blocks :
            dense_summary.status,
        dense_local_overlap_blocker =
            isnothing(dense_summary) ?
            :not_attempted_missing_available_axis_blocks :
            dense_summary.blocker,
        dense_block_shape =
            isnothing(dense_summary) ? :not_materialized : dense_summary.dense_block_shape,
        local_source_fingerprint_status =
            isnothing(local_source_fingerprint) ?
            :not_attempted_missing_dense_local_overlap :
            local_source_fingerprint.status,
        local_source_fingerprint_blocker =
            isnothing(local_source_fingerprint) ?
            :not_attempted_missing_dense_local_overlap :
            local_source_fingerprint.blocker,
        local_overlap_record_status =
            isnothing(local_block_record_summary) ?
            :not_attempted_missing_dense_local_overlap :
            local_block_record_summary.status,
        local_overlap_record_blocker =
            isnothing(local_block_record_summary) ?
            :not_attempted_missing_dense_local_overlap :
            local_block_record_summary.blocker,
        local_overlap_record_source_kind =
            isnothing(local_block_record_summary) ?
            :not_attempted_missing_dense_local_overlap :
            local_block_record_summary.source_kind,
        local_overlap_record_term =
            isnothing(local_block_record_summary) ?
            :not_attempted_missing_dense_local_overlap :
            local_block_record_summary.term,
        local_overlap_record_dense_block_available =
            !isnothing(local_block_record_summary) &&
            local_block_record_summary.dense_block_available,
        local_overlap_record_dense_block_shape =
            isnothing(local_block_record_summary) ?
            :not_materialized :
            local_block_record_summary.dense_block_shape,
        local_overlap_record_placement_status =
            isnothing(local_block_record_summary) ?
            :not_attempted_missing_dense_local_overlap :
            local_block_record_summary.placement_status,
        local_overlap_record_retained_transform_status =
            isnothing(local_block_record_summary) ?
            :not_attempted_missing_dense_local_overlap :
            local_block_record_summary.retained_transform_status,
        local_overlap_collection_status =
            isnothing(local_block_collection_summary) ?
            :not_attempted_missing_local_overlap_record :
            local_block_collection_summary.status,
        local_overlap_collection_blocker =
            isnothing(local_block_collection_summary) ?
            :not_attempted_missing_local_overlap_record :
            local_block_collection_summary.blocker,
        local_overlap_collection_record_count =
            isnothing(local_block_collection_summary) ?
            0 :
            local_block_collection_summary.record_count,
        local_overlap_collection_terms =
            isnothing(local_block_collection_summary) ?
            () :
            local_block_collection_summary.terms,
        local_overlap_collection_block_keys =
            isnothing(local_block_collection_summary) ?
            () :
            local_block_collection_summary.block_keys,
        local_overlap_collection_dense_block_count =
            isnothing(local_block_collection_summary) ?
            0 :
            local_block_collection_summary.dense_block_count,
        local_overlap_collection_placement_status =
            isnothing(local_block_collection_summary) ?
            :not_attempted_missing_local_overlap_record :
            local_block_collection_summary.placement_status,
        local_overlap_collection_retained_transform_status =
            isnothing(local_block_collection_summary) ?
            :not_attempted_missing_local_overlap_record :
            local_block_collection_summary.retained_transform_status,
        local_overlap_collection_global_matrix_materialized =
            !isnothing(local_block_collection_summary) &&
            local_block_collection_summary.global_matrix_materialized,
        local_overlap_collection_route_driver_wiring =
            !isnothing(local_block_collection_summary) &&
            local_block_collection_summary.route_driver_wiring,
        local_overlap_collection_adapter_status =
            isnothing(local_block_collection_adapter) ?
            :not_attempted_missing_local_overlap_collection :
            local_block_collection_adapter.status,
        local_overlap_collection_adapter_blocker =
            isnothing(local_block_collection_adapter) ?
            :not_attempted_missing_local_overlap_collection :
            local_block_collection_adapter.blocker,
        local_overlap_collection_adapter_global_overlap_status =
            isnothing(local_block_collection_adapter) ?
            :not_attempted_missing_local_overlap_collection :
            local_block_collection_adapter.global_overlap_status,
        local_overlap_collection_adapter_global_overlap_blocker =
            isnothing(local_block_collection_adapter) ?
            :not_attempted_missing_local_overlap_collection :
            local_block_collection_adapter.global_overlap_blocker,
        local_overlap_collection_adapter_private_input_facts_available =
            !isnothing(local_block_collection_adapter) &&
            local_block_collection_adapter.private_global_overlap_input_facts_available,
        local_overlap_collection_adapter_route_global_stage_source =
            !isnothing(local_block_collection_adapter) &&
            local_block_collection_adapter.route_global_overlap_stage_source,
        placement_requirements_status =
            isnothing(placement_requirements) ?
            :not_attempted_missing_local_overlap_collection :
            placement_requirements.status,
        placement_requirements_blocker =
            isnothing(placement_requirements) ?
            :not_attempted_missing_local_overlap_collection :
            placement_requirements.blocker,
        placement_missing_requirements =
            isnothing(placement_requirements) ?
            () :
            placement_requirements.missing_requirements,
        placement_global_overlap_status =
            isnothing(placement_requirements) ?
            :not_attempted_missing_local_overlap_collection :
            placement_requirements.global_overlap_status,
        placement_global_overlap_blocker =
            isnothing(placement_requirements) ?
            :not_attempted_missing_local_overlap_collection :
            placement_requirements.global_overlap_blocker,
        placement_candidate_status =
            isnothing(placement_candidate) ?
            :not_attempted_missing_placement_requirements :
            placement_candidate.status,
        placement_candidate_blocker =
            isnothing(placement_candidate) ?
            :not_attempted_missing_placement_requirements :
            placement_candidate.blocker,
        placement_candidate_available_requirements =
            isnothing(placement_candidate) ?
            () :
            placement_candidate.available_requirements,
        placement_candidate_missing_requirements =
            isnothing(placement_candidate) ?
            () :
            placement_candidate.missing_requirements,
        placement_candidate_global_overlap_status =
            isnothing(placement_candidate) ?
            :not_attempted_missing_placement_requirements :
            placement_candidate.global_overlap_status,
        placement_candidate_global_overlap_blocker =
            isnothing(placement_candidate) ?
            :not_attempted_missing_placement_requirements :
            placement_candidate.global_overlap_blocker,
        placement_candidate_route_driver_wiring =
            !isnothing(placement_candidate) &&
            placement_candidate.route_driver_wiring,
        placement_candidate_global_matrix_materialized =
            !isnothing(placement_candidate) &&
            placement_candidate.global_matrix_materialized,
        partial_placement_candidate_global_dimension_source =
            hasproperty(report, :retained_dimension) ?
            :report_retained_dimension :
            :fallback_test_dimension,
        partial_placement_candidate_global_dimension =
            placement_candidate_global_dimension,
        partial_placement_candidate_status =
            isnothing(partial_placement_candidate) ?
            :not_attempted_missing_local_overlap_collection :
            partial_placement_candidate.status,
        partial_placement_candidate_blocker =
            isnothing(partial_placement_candidate) ?
            :not_attempted_missing_local_overlap_collection :
            partial_placement_candidate.blocker,
        partial_placement_candidate_available_requirements =
            isnothing(partial_placement_candidate) ?
            () :
            partial_placement_candidate.available_requirements,
        partial_placement_candidate_missing_requirements =
            isnothing(partial_placement_candidate) ?
            () :
            partial_placement_candidate.missing_requirements,
        partial_placement_candidate_global_overlap_status =
            isnothing(partial_placement_candidate) ?
            :not_attempted_missing_local_overlap_collection :
            partial_placement_candidate.global_overlap_status,
        partial_placement_candidate_route_driver_wiring =
            !isnothing(partial_placement_candidate) &&
            partial_placement_candidate.route_driver_wiring,
        partial_placement_candidate_global_matrix_materialized =
            !isnothing(partial_placement_candidate) &&
            partial_placement_candidate.global_matrix_materialized,
        all_facts_placement_candidate_status =
            isnothing(all_facts_placement_candidate) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_candidate.status,
        all_facts_placement_candidate_status_detail =
            isnothing(all_facts_placement_candidate) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_candidate.placement_candidate_status,
        all_facts_placement_candidate_blocker =
            isnothing(all_facts_placement_candidate) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_candidate.blocker,
        all_facts_placement_candidate_missing_requirements =
            isnothing(all_facts_placement_candidate) ?
            () :
            all_facts_placement_candidate.missing_requirements,
        all_facts_placement_candidate_global_overlap_status =
            isnothing(all_facts_placement_candidate) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_candidate.global_overlap_status,
        all_facts_placement_candidate_global_overlap_blocker =
            isnothing(all_facts_placement_candidate) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_candidate.global_overlap_blocker,
        all_facts_placement_candidate_route_driver_wiring =
            !isnothing(all_facts_placement_candidate) &&
            all_facts_placement_candidate.route_driver_wiring,
        all_facts_placement_candidate_global_matrix_materialized =
            !isnothing(all_facts_placement_candidate) &&
            all_facts_placement_candidate.global_matrix_materialized,
        all_facts_placement_candidate_route_global_stage_source =
            !isnothing(all_facts_placement_candidate) &&
            all_facts_placement_candidate.route_global_overlap_stage_source,
        all_facts_placement_plan_skeleton_status =
            isnothing(all_facts_placement_plan_skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_plan_skeleton.status,
        all_facts_placement_plan_skeleton_blocker =
            isnothing(all_facts_placement_plan_skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_plan_skeleton.blocker,
        all_facts_placement_plan_skeleton_record_count =
            isnothing(all_facts_placement_plan_skeleton) ?
            0 :
            all_facts_placement_plan_skeleton.record_count,
        all_facts_placement_plan_skeleton_block_keys =
            isnothing(all_facts_placement_plan_skeleton) ?
            () :
            all_facts_placement_plan_skeleton.block_keys,
        all_facts_placement_plan_skeleton_global_dimension =
            isnothing(all_facts_placement_plan_skeleton) ?
            nothing :
            all_facts_placement_plan_skeleton.global_dimension,
        all_facts_placement_plan_skeleton_global_dimension_source =
            isnothing(all_facts_placement_plan_skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_plan_skeleton.global_dimension_source,
        all_facts_placement_plan_skeleton_missing_requirements =
            isnothing(all_facts_placement_plan_skeleton) ?
            () :
            all_facts_placement_plan_skeleton.missing_requirements,
        all_facts_placement_plan_skeleton_placement_plan_status =
            isnothing(all_facts_placement_plan_skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_plan_skeleton.placement_plan_status,
        all_facts_placement_plan_skeleton_accumulation_rule_status =
            isnothing(all_facts_placement_plan_skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_plan_skeleton.accumulation_rule_status,
        all_facts_placement_plan_skeleton_global_overlap_status =
            isnothing(all_facts_placement_plan_skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_plan_skeleton.global_overlap_status,
        all_facts_placement_plan_skeleton_global_overlap_blocker =
            isnothing(all_facts_placement_plan_skeleton) ?
            :not_attempted_missing_local_overlap_collection :
            all_facts_placement_plan_skeleton.global_overlap_blocker,
        all_facts_placement_plan_skeleton_route_driver_wiring =
            !isnothing(all_facts_placement_plan_skeleton) &&
            all_facts_placement_plan_skeleton.route_driver_wiring,
        all_facts_placement_plan_skeleton_global_matrix_materialized =
            !isnothing(all_facts_placement_plan_skeleton) &&
            all_facts_placement_plan_skeleton.global_matrix_materialized,
        all_facts_placement_plan_skeleton_route_global_stage_source =
            !isnothing(all_facts_placement_plan_skeleton) &&
            all_facts_placement_plan_skeleton.route_global_overlap_stage_source,
        placement_source_audit,
        route_driver_wiring = false,
        global_matrix_materialized = false,
        route_global_overlap_stage_source = false,
    )
end

function _test_driver_overlap_nonclaim_flags(result)
    @test !result.route_driver_wiring
    @test !result.hamiltonian_data_materialized
    @test !result.global_hamiltonian_data_materialized
    @test !result.coulomb_materialized
    @test !result.ida_mwg_data_materialized
    @test !result.pqs_lowdin_materialized
    @test !result.pqs_shell_projection_materialized
    @test !result.artifacts_materialized
    @test !result.exports_materialized
    @test !result.full_white_lindsey_route_assembled
end

@testset "CartesianPairBlockMaterialization driver global overlap hook" begin
    plan = _driver_overlap_plan()
    expected = CPBMDriverOverlap.route_state_global_overlap_matrix(
        plan;
        global_dimension = 2,
        inputs = _driver_overlap_inputs(),
    )
    direct = CPBMDriverOverlap.driver_global_overlap_result(
        plan;
        global_dimension = 2,
        inputs = _driver_overlap_inputs(),
    )
    wrapped = CPBMDriverOverlap.driver_global_overlap_result(
        (; pair_block_materialization_plan = plan);
        global_dimension = 2,
        factors = _driver_overlap_inputs(),
    )
    terminal_route_state_source =
        (; terminal_route_state = (; pair_block_materialization_plan = plan))
    terminal_route_state = CPBMDriverOverlap.driver_global_overlap_result(
        terminal_route_state_source;
        global_dimension = 2,
        inputs = _driver_overlap_inputs(),
    )
    expected_matrix = _driver_overlap_expected_matrix()

    @test direct.status === :materialized_route_global_overlap_matrix
    @test direct.global_one_body_term_matrix_materialized
    @test direct.global_overlap_matrix_materialized
    @test size(terminal_route_state.global_matrix_result.matrix) == (2, 2)
    @test terminal_route_state.status === :materialized_route_global_overlap_matrix
    @test terminal_route_state.global_overlap_matrix_materialized
    @test terminal_route_state.global_one_body_term_matrix_materialized
    @test terminal_route_state.global_matrix_result.matrix ≈ expected_matrix
    @test terminal_route_state.global_matrix_result.matrix[1, 1] ≈ 1.68
    @test terminal_route_state.global_matrix_result.matrix[1, 2] ≈ 0.42
    @test terminal_route_state.global_matrix_result.matrix[2, 1] ≈ 0.42
    @test terminal_route_state.global_matrix_result.matrix[2, 2] ≈ 1.82
    @test direct.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test direct.global_matrix_result.matrix ≈ expected_matrix
    @test wrapped.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test terminal_route_state.global_matrix_result.matrix ≈
          expected.global_matrix_result.matrix
    @test direct.pair_block_materialization_plan === plan

    missing_plan = CPBMDriverOverlap.driver_global_overlap_result(
        (; route_stage = :missing_pair_block_plan);
        global_dimension = 2,
        inputs = _driver_overlap_inputs(),
    )
    @test missing_plan.status === :blocked_route_global_overlap_matrix
    @test missing_plan.blocker === :missing_pair_block_materialization_plan
    @test !missing_plan.global_one_body_term_matrix_materialized
    @test !missing_plan.global_overlap_matrix_materialized

    missing_dimension = CPBMDriverOverlap.driver_global_overlap_result(
        (; pair_block_materialization_plan = plan);
        inputs = _driver_overlap_inputs(),
    )
    @test missing_dimension.status === :blocked_route_global_overlap_matrix
    @test missing_dimension.blocker === :missing_global_dimension
    @test !missing_dimension.global_one_body_term_matrix_materialized

    _test_driver_overlap_nonclaim_flags(direct)
    _test_driver_overlap_nonclaim_flags(missing_plan)
end

@testset "PQS route driver private global overlap option" begin
    plan = _driver_overlap_plan()
    report = _driver_overlap_stage_report(plan)

    off_stage =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            report,
        )
    @test Tuple(keys(off_stage)) ==
          (:private_global_overlap_result, :private_global_overlap_summary)
    @test isnothing(off_stage.private_global_overlap_result)
    @test isnothing(off_stage.private_global_overlap_summary)

    missing_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            (;);
            private_global_overlap_requested = true,
        )
    @test Tuple(keys(missing_facts)) ==
          (:private_global_overlap_result, :private_global_overlap_summary)
    @test isnothing(missing_facts.private_global_overlap_result)
    @test missing_facts.private_global_overlap_summary.blocker ===
          :missing_final_retained_column_layout

    missing_plan =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            (;);
            private_global_overlap_requested = true,
            private_global_overlap_global_dimension = 2,
            private_global_overlap_inputs = _driver_overlap_inputs(),
        )
    @test missing_plan.private_global_overlap_result.blocker ===
          :missing_pair_block_materialization_plan
    @test missing_plan.private_global_overlap_summary.blocker ===
          :missing_pair_block_materialization_plan

    missing_dimension =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            report;
            private_global_overlap_requested = true,
            private_global_overlap_inputs = _driver_overlap_inputs(),
        )
    @test missing_dimension.private_global_overlap_result.blocker ===
          :missing_global_dimension
    @test missing_dimension.private_global_overlap_summary.blocker ===
          :missing_global_dimension

    missing_structured_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            report;
            private_global_overlap_requested = true,
            private_global_overlap_global_dimension = 2,
        )
    @test isnothing(missing_structured_facts.private_global_overlap_result)
    @test missing_structured_facts.private_global_overlap_summary.status ===
          :blocked_private_global_overlap
    @test missing_structured_facts.private_global_overlap_summary.blocker ===
          :missing_final_retained_column_layout

    structured_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            _driver_overlap_facts_report(plan),
        )
    @test structured_facts.status ===
          :available_private_global_overlap_input_facts
    @test structured_facts.global_dimension == 2
    @test structured_facts.parent_axis_counts == (2, 2, 2)
    @test structured_facts.final_layout_source === :retained_dimension_compatibility
    @test structured_facts.parent_axis_counts_source === :report_parent_axis_counts
    @test structured_facts.axis_bundle_source === :report_parent_axis_bundle_object
    @test structured_facts.factor_space === :parent_axis_bundle_pgdg_intermediate
    @test structured_facts.factor_convention === :axis_bundle_one_body_overlap
    @test structured_facts.overlap_1d_source === :report_parent_axis_bundle_object
    @test structured_facts.overlap_1d.x ≈ _driver_overlap_inputs().overlap_1d.x

    carried_parent_facts_report = merge(
        _driver_overlap_stage_report(plan),
        (;
            retained_dimension = 2,
            parent = _driver_overlap_parent_object(),
            parent_axis_bundle_object = _driver_overlap_axis_bundle_object(),
        ),
    )
    carried_parent_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            carried_parent_facts_report,
        )
    @test carried_parent_facts.status ===
          :available_private_global_overlap_input_facts
    @test carried_parent_facts.global_dimension == 2
    @test carried_parent_facts.parent_axis_counts == (2, 2, 2)
    @test carried_parent_facts.parent_axis_counts_source ===
          :parent_object_parent_axis_counts
    @test carried_parent_facts.axis_bundle_source ===
          :report_parent_axis_bundle_object
    @test carried_parent_facts.factor_space ===
          :parent_axis_bundle_pgdg_intermediate
    @test carried_parent_facts.factor_convention ===
          :axis_bundle_one_body_overlap

    materialized_from_carried_parent_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            carried_parent_facts_report;
            private_global_overlap_requested = true,
        )
    @test materialized_from_carried_parent_facts.private_global_overlap_result.status ===
          :materialized_route_global_overlap_matrix
    @test materialized_from_carried_parent_facts.private_global_overlap_result.global_matrix_result.matrix ≈
          _driver_overlap_expected_matrix()

    retained_units_facts_report = (;
        low_order_route_summary = (;
            terminal_route_state = (;
                pair_block_materialization_plan = plan,
                retained_units = ((; column_range = 1:2),),
            ),
        ),
        parent = _driver_overlap_parent_object(),
        parent_axis_bundle_object = _driver_overlap_axis_bundle_object(),
    )
    retained_units_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            retained_units_facts_report,
        )
    @test retained_units_facts.status ===
          :available_private_global_overlap_input_facts
    @test retained_units_facts.global_dimension == 2
    @test retained_units_facts.final_layout_source ===
          :terminal_route_state_retained_units
    @test retained_units_facts.parent_axis_counts == (2, 2, 2)
    @test retained_units_facts.parent_axis_counts_source ===
          :parent_object_parent_axis_counts
    @test retained_units_facts.axis_bundle_source ===
          :report_parent_axis_bundle_object
    @test retained_units_facts.factor_space ===
          :parent_axis_bundle_pgdg_intermediate
    @test retained_units_facts.factor_convention ===
          :axis_bundle_one_body_overlap

    materialized_from_retained_units_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            retained_units_facts_report;
            private_global_overlap_requested = true,
        )
    @test materialized_from_retained_units_facts.private_global_overlap_result.status ===
          :materialized_route_global_overlap_matrix
    @test materialized_from_retained_units_facts.private_global_overlap_result.global_matrix_result.matrix ≈
          _driver_overlap_expected_matrix()

    missing_parent_axis_counts_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            (;
                low_order_route_summary = (;
                    terminal_route_state =
                        (; pair_block_materialization_plan = plan),
                ),
                retained_dimension = 2,
                parent_axis_bundle_object = _driver_overlap_axis_bundle_object(),
            ),
        )
    @test missing_parent_axis_counts_facts.status ===
          :blocked_private_global_overlap_input_facts
    @test missing_parent_axis_counts_facts.blocker ===
          :missing_parent_axis_counts
    @test missing_parent_axis_counts_facts.global_dimension == 2
    @test missing_parent_axis_counts_facts.parent_axis_counts === nothing
    @test missing_parent_axis_counts_facts.final_layout_source ===
          :retained_dimension_compatibility
    @test missing_parent_axis_counts_facts.parent_axis_counts_source ===
          :unavailable
    @test missing_parent_axis_counts_facts.factor_space === :unavailable
    @test missing_parent_axis_counts_facts.factor_convention === :unavailable
    @test missing_parent_axis_counts_facts.overlap_1d === nothing

    missing_axis_bundle_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            (;
                terminal_route_state = (; pair_block_materialization_plan = plan),
                retained_dimension = 2,
                parent_axis_counts = (2, 2, 2),
            ),
        )
    @test missing_axis_bundle_facts.status ===
          :blocked_private_global_overlap_input_facts
    @test missing_axis_bundle_facts.blocker ===
          :missing_parent_axis_bundle_overlap_factors
    @test missing_axis_bundle_facts.global_dimension == 2
    @test missing_axis_bundle_facts.parent_axis_counts == (2, 2, 2)
    @test missing_axis_bundle_facts.final_layout_source ===
          :retained_dimension_compatibility
    @test missing_axis_bundle_facts.parent_axis_counts_source ===
          :report_parent_axis_counts
    @test missing_axis_bundle_facts.axis_bundle_source === :unavailable
    @test missing_axis_bundle_facts.factor_space === :unavailable
    @test missing_axis_bundle_facts.factor_convention === :unavailable
    @test missing_axis_bundle_facts.overlap_1d === nothing

    materialized_from_facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            _driver_overlap_facts_report(plan);
            private_global_overlap_requested = true,
        )
    @test materialized_from_facts.private_global_overlap_result.status ===
          :materialized_route_global_overlap_matrix
    @test materialized_from_facts.private_global_overlap_result.global_matrix_result.matrix ≈
          _driver_overlap_expected_matrix()

    materialized =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            report;
            private_global_overlap_requested = true,
            private_global_overlap_global_dimension = 2,
            private_global_overlap_inputs = _driver_overlap_inputs(),
        )
    @test materialized.private_global_overlap_result.status ===
          :materialized_route_global_overlap_matrix
    @test materialized.private_global_overlap_summary.status ===
          :materialized_route_global_overlap_matrix
    @test materialized.private_global_overlap_summary.result_available
    @test materialized.private_global_overlap_summary.global_overlap_matrix_materialized
    @test materialized.private_global_overlap_result.global_matrix_result.matrix ≈
          _driver_overlap_expected_matrix()
    _test_driver_overlap_nonclaim_flags(materialized.private_global_overlap_result)
    @test !materialized.private_global_overlap_summary.route_driver_wiring
    @test !materialized.private_global_overlap_summary.hamiltonian_data_materialized
    @test !materialized.private_global_overlap_summary.global_hamiltonian_data_materialized
    @test !materialized.private_global_overlap_summary.coulomb_materialized
    @test !materialized.private_global_overlap_summary.ida_mwg_data_materialized
    @test !materialized.private_global_overlap_summary.pqs_lowdin_materialized
    @test !materialized.private_global_overlap_summary.pqs_shell_projection_materialized
    @test !materialized.private_global_overlap_summary.artifacts_materialized
    @test !materialized.private_global_overlap_summary.exports_materialized
    @test !materialized.private_global_overlap_summary.full_white_lindsey_route_assembled
end

@testset "PQS route driver private global overlap real report fingerprint" begin
    report = _driver_overlap_real_report_fingerprint()
    facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            report,
        )

    @test hasproperty(report, :low_order_route_summary)
    @test facts.status === :blocked_private_global_overlap_input_facts
    @test facts.blocker === :missing_parent_axis_bundle_overlap_factors
    @test facts.global_dimension == 221
    @test facts.final_layout_source === :retained_dimension_compatibility
    @test facts.parent_axis_counts == (5, 5, 5)
    @test facts.parent_axis_counts_source === :report_route_axis_counts
    @test facts.axis_bundle_source === :unavailable
    @test facts.factor_space === :unavailable
    @test facts.factor_convention === :unavailable
end

@testset "PQS route driver private global overlap probe-enabled report fingerprint" begin
    report = _driver_overlap_probe_enabled_real_report_fingerprint()
    facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            report,
        )

    @test report.route_materializer_payload.parent_axis_bundle_object_available
    @test facts.status === :available_private_global_overlap_input_facts
    @test facts.global_dimension == report.retained_dimension
    @test facts.final_layout_source === :retained_dimension_compatibility
    @test facts.parent_axis_counts == (
        report.route_axis_counts.parent_axis_counts.x,
        report.route_axis_counts.parent_axis_counts.y,
        report.route_axis_counts.parent_axis_counts.z,
    )
    @test facts.parent_axis_counts_source === :report_route_axis_counts
    @test facts.axis_bundle_source ===
          :route_materializer_payload_parent_axis_bundle_object
    @test facts.factor_space === :parent_axis_bundle_pgdg_intermediate
    @test facts.factor_convention === :axis_bundle_one_body_overlap
    @test !isnothing(facts.overlap_1d)

    local_cpb_overlap_fingerprint =
        _driver_overlap_real_report_local_cpb_provider_fingerprint(report)
    @test local_cpb_overlap_fingerprint.parent_source_status ===
          :available_route_materializer_payload_parent_qw_basis_object
    @test local_cpb_overlap_fingerprint.parent_axis_counts == (31, 17, 17)
    @test local_cpb_overlap_fingerprint.parent_overlap_packet_source_status ===
          :available_parent_overlap_axis_factors
    @test local_cpb_overlap_fingerprint.parent_overlap_packet_blocker === nothing
    @test local_cpb_overlap_fingerprint.cpb_source_pair_status ===
          :available_report_pair_entry_source_boxes
    @test local_cpb_overlap_fingerprint.cpb_source_pair_key ===
          (:product, :product)
    @test local_cpb_overlap_fingerprint.interval_pair_status ===
          :available_cpb_interval_pair
    @test local_cpb_overlap_fingerprint.interval_pair_blocker === nothing
    @test local_cpb_overlap_fingerprint.overlap_axis_blocks_status ===
          :available_cpb_overlap_axis_blocks
    @test local_cpb_overlap_fingerprint.overlap_axis_blocks_blocker === nothing
    @test local_cpb_overlap_fingerprint.dense_local_overlap_status ===
          :materialized_cpb_overlap_dense_block
    @test local_cpb_overlap_fingerprint.dense_local_overlap_blocker === nothing
    @test local_cpb_overlap_fingerprint.dense_block_shape == (25, 25)
    @test local_cpb_overlap_fingerprint.local_source_fingerprint_status ===
          :available_private_global_overlap_local_source_fingerprint
    @test local_cpb_overlap_fingerprint.local_source_fingerprint_blocker === nothing
    @test local_cpb_overlap_fingerprint.local_overlap_record_status ===
          :available_cpb_local_overlap_block_record
    @test local_cpb_overlap_fingerprint.local_overlap_record_blocker === nothing
    @test local_cpb_overlap_fingerprint.local_overlap_record_source_kind ===
          :cpb_overlap_dense_block
    @test local_cpb_overlap_fingerprint.local_overlap_record_term === :overlap
    @test local_cpb_overlap_fingerprint.local_overlap_record_dense_block_available
    @test local_cpb_overlap_fingerprint.local_overlap_record_dense_block_shape ==
          (25, 25)
    @test local_cpb_overlap_fingerprint.local_overlap_record_placement_status ===
          :unassigned
    @test local_cpb_overlap_fingerprint.local_overlap_record_retained_transform_status ===
          :unassigned
    @test local_cpb_overlap_fingerprint.local_overlap_collection_status ===
          :available_cpb_local_overlap_block_collection
    @test local_cpb_overlap_fingerprint.local_overlap_collection_blocker === nothing
    @test local_cpb_overlap_fingerprint.local_overlap_collection_record_count == 1
    @test local_cpb_overlap_fingerprint.local_overlap_collection_terms ===
          (:overlap,)
    @test local_cpb_overlap_fingerprint.local_overlap_collection_block_keys ===
          ((:product, :product),)
    @test local_cpb_overlap_fingerprint.local_overlap_collection_dense_block_count == 1
    @test local_cpb_overlap_fingerprint.local_overlap_collection_placement_status ===
          :unassigned
    @test local_cpb_overlap_fingerprint.local_overlap_collection_retained_transform_status ===
          :unassigned
    @test !local_cpb_overlap_fingerprint.local_overlap_collection_global_matrix_materialized
    @test !local_cpb_overlap_fingerprint.local_overlap_collection_route_driver_wiring
    @test local_cpb_overlap_fingerprint.local_overlap_collection_adapter_status ===
          :blocked_private_global_overlap_local_collection_adapter
    @test local_cpb_overlap_fingerprint.local_overlap_collection_adapter_blocker ===
          :missing_placement_or_retained_transform
    @test local_cpb_overlap_fingerprint.local_overlap_collection_adapter_global_overlap_status ===
          :blocked
    @test local_cpb_overlap_fingerprint.local_overlap_collection_adapter_global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test !local_cpb_overlap_fingerprint.local_overlap_collection_adapter_private_input_facts_available
    @test !local_cpb_overlap_fingerprint.local_overlap_collection_adapter_route_global_stage_source
    @test local_cpb_overlap_fingerprint.placement_requirements_status ===
          :blocked_private_global_overlap_placement_requirements
    @test local_cpb_overlap_fingerprint.placement_requirements_blocker ===
          :missing_placement_or_retained_transform
    @test local_cpb_overlap_fingerprint.placement_missing_requirements === (
        :missing_retained_transform,
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_global_dimension,
        :missing_placement_plan,
        :missing_accumulation_rule,
    )
    @test local_cpb_overlap_fingerprint.placement_global_overlap_status ===
          :blocked
    @test local_cpb_overlap_fingerprint.placement_global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test local_cpb_overlap_fingerprint.placement_candidate_status ===
          :blocked_private_global_overlap_placement_candidate
    @test local_cpb_overlap_fingerprint.placement_candidate_blocker ===
          :missing_placement_or_retained_transform
    @test local_cpb_overlap_fingerprint.placement_candidate_available_requirements ===
          (:local_cpb_overlap_collection,)
    @test local_cpb_overlap_fingerprint.placement_candidate_missing_requirements ===
          local_cpb_overlap_fingerprint.placement_missing_requirements
    @test local_cpb_overlap_fingerprint.placement_candidate_global_overlap_status ===
          :blocked
    @test local_cpb_overlap_fingerprint.placement_candidate_global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test !local_cpb_overlap_fingerprint.placement_candidate_route_driver_wiring
    @test !local_cpb_overlap_fingerprint.placement_candidate_global_matrix_materialized
    @test local_cpb_overlap_fingerprint.partial_placement_candidate_global_dimension_source ===
          :report_retained_dimension
    @test local_cpb_overlap_fingerprint.partial_placement_candidate_global_dimension ==
          report.retained_dimension
    @test local_cpb_overlap_fingerprint.partial_placement_candidate_status ===
          :blocked_private_global_overlap_placement_candidate
    @test local_cpb_overlap_fingerprint.partial_placement_candidate_blocker ===
          :missing_placement_or_retained_transform
    @test local_cpb_overlap_fingerprint.partial_placement_candidate_available_requirements === (
        :local_cpb_overlap_collection,
        :retained_transform,
        :global_dimension,
    )
    @test local_cpb_overlap_fingerprint.partial_placement_candidate_missing_requirements === (
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_placement_plan,
        :missing_accumulation_rule,
    )
    @test local_cpb_overlap_fingerprint.partial_placement_candidate_global_overlap_status ===
          :blocked
    @test !local_cpb_overlap_fingerprint.partial_placement_candidate_route_driver_wiring
    @test !local_cpb_overlap_fingerprint.partial_placement_candidate_global_matrix_materialized
    @test local_cpb_overlap_fingerprint.all_facts_placement_candidate_status ===
          :blocked_private_global_overlap_placement_candidate
    @test local_cpb_overlap_fingerprint.all_facts_placement_candidate_status_detail ===
          :blocked_placement_not_implemented
    @test local_cpb_overlap_fingerprint.all_facts_placement_candidate_blocker ===
          :placement_not_implemented
    @test local_cpb_overlap_fingerprint.all_facts_placement_candidate_missing_requirements ===
          ()
    @test local_cpb_overlap_fingerprint.all_facts_placement_candidate_global_overlap_status ===
          :blocked
    @test local_cpb_overlap_fingerprint.all_facts_placement_candidate_global_overlap_blocker ===
          :placement_not_implemented
    @test !local_cpb_overlap_fingerprint.all_facts_placement_candidate_route_driver_wiring
    @test !local_cpb_overlap_fingerprint.all_facts_placement_candidate_global_matrix_materialized
    @test !local_cpb_overlap_fingerprint.all_facts_placement_candidate_route_global_stage_source
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_blocker ===
          :placement_not_implemented
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_record_count ==
          1
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_block_keys ===
          ((:product, :product),)
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_global_dimension ==
          report.retained_dimension
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_global_dimension_source ===
          :report_retained_dimension
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_missing_requirements ===
          ()
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_placement_plan_status ===
          :available_placement_plan
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_accumulation_rule_status ===
          :available_accumulation_rule
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_global_overlap_status ===
          :blocked
    @test local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_global_overlap_blocker ===
          :placement_not_implemented
    @test !local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_route_driver_wiring
    @test !local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_global_matrix_materialized
    @test !local_cpb_overlap_fingerprint.all_facts_placement_plan_skeleton_route_global_stage_source
    placement_source_audit =
        local_cpb_overlap_fingerprint.placement_source_audit
    @test placement_source_audit.retained_transform_source_status ===
          :missing_actual_retained_transform_source
    @test placement_source_audit.retained_transform_source === :unavailable
    @test placement_source_audit.left_column_range_source_status ===
          :missing_source_pair_retained_column_range
    @test placement_source_audit.left_column_range_source === :unavailable
    @test placement_source_audit.right_column_range_source_status ===
          :missing_source_pair_retained_column_range
    @test placement_source_audit.right_column_range_source === :unavailable
    @test placement_source_audit.global_dimension_source_status ===
          :available_retained_layout_global_dimension
    @test placement_source_audit.global_dimension_source === :report_retained_units
    @test placement_source_audit.global_dimension == report.retained_dimension
    @test placement_source_audit.placement_plan_source_status ===
          :missing_reviewed_overlap_placement_plan
    @test placement_source_audit.placement_plan_source === :unavailable
    @test placement_source_audit.accumulation_rule_source_status ===
          :missing_overlap_accumulation_rule
    @test placement_source_audit.accumulation_rule_source === :unavailable
    @test placement_source_audit.candidate_status_with_real_sources ===
          :blocked_private_global_overlap_placement_candidate
    @test placement_source_audit.candidate_blocker_with_real_sources ===
          :missing_placement_or_retained_transform
    @test placement_source_audit.skeleton_status_with_real_sources ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test placement_source_audit.skeleton_blocker_with_real_sources ===
          :missing_placement_or_retained_transform
    @test placement_source_audit.available_requirements_with_real_sources === (
        :local_cpb_overlap_collection,
        :global_dimension,
    )
    @test placement_source_audit.missing_requirements_with_real_sources === (
        :missing_retained_transform,
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_placement_plan,
        :missing_accumulation_rule,
    )
    @test placement_source_audit.placement_facts_status ===
          :blocked_cpb_overlap_placement_facts
    @test placement_source_audit.placement_facts_blocker ===
          :missing_placement_or_retained_transform
    @test placement_source_audit.placement_facts_available_requirements ===
          (:local_cpb_overlap_collection,)
    @test placement_source_audit.placement_facts_missing_requirements === (
        :missing_retained_transform,
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_global_dimension,
        :missing_placement_plan,
        :missing_accumulation_rule,
    )
    @test placement_source_audit.placement_facts_record_left_cpb_summary !==
          :unavailable
    @test placement_source_audit.placement_facts_record_right_cpb_summary !==
          :unavailable
    @test placement_source_audit.placement_facts_skeleton_status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test placement_source_audit.placement_facts_skeleton_blocker ===
          :missing_placement_or_retained_transform
    @test placement_source_audit.placement_facts_skeleton_missing_requirements ===
          placement_source_audit.placement_facts_missing_requirements
    @test placement_source_audit.placement_facts_skeleton_global_overlap_status ===
          :blocked
    @test !placement_source_audit.placement_facts_skeleton_route_driver_wiring
    @test !placement_source_audit.placement_facts_skeleton_global_matrix_materialized
    @test !placement_source_audit.placement_facts_skeleton_private_input_facts_available
    @test !placement_source_audit.placement_facts_skeleton_route_global_stage_source
    @test placement_source_audit.placement_facts_skeleton_record_left_cpb_summary ==
          placement_source_audit.placement_facts_record_left_cpb_summary
    @test placement_source_audit.placement_facts_skeleton_record_right_cpb_summary ==
          placement_source_audit.placement_facts_record_right_cpb_summary
    @test placement_source_audit.reviewed_placement_facts_status ===
          :blocked_cpb_overlap_placement_facts
    @test placement_source_audit.reviewed_placement_facts_blocker ===
          :missing_placement_or_retained_transform
    @test placement_source_audit.reviewed_placement_facts_available_requirements === (
        :local_cpb_overlap_collection,
        :placement_plan,
        :accumulation_rule,
    )
    @test placement_source_audit.reviewed_placement_facts_missing_requirements === (
        :missing_retained_transform,
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_global_dimension,
    )
    @test placement_source_audit.reviewed_placement_facts_placement_plan_status ===
          :available_placement_plan
    @test placement_source_audit.reviewed_placement_facts_placement_plan_kind ===
          :real_report_reviewed_overlap_placement_plan_fingerprint
    @test placement_source_audit.reviewed_placement_facts_accumulation_rule_status ===
          :available_accumulation_rule
    @test placement_source_audit.reviewed_placement_facts_accumulation_rule ===
          :add_explicit_blocks_into_ranges
    @test placement_source_audit.reviewed_placement_facts_inventory_status ===
          :available_cpb_overlap_placement_record_inventory
    @test placement_source_audit.reviewed_placement_facts_inventory_blocker ===
          nothing
    @test placement_source_audit.reviewed_placement_facts_accepted_block_keys ===
          ((:product, :product),)
    @test placement_source_audit.reviewed_placement_facts_provided_block_keys ===
          ((:product, :product),)
    @test placement_source_audit.reviewed_placement_facts_rejected_block_keys ===
          ()
    @test placement_source_audit.reviewed_placement_facts_duplicate_block_keys ===
          ()
    @test placement_source_audit.reviewed_placement_facts_duplicate_record_policy ===
          :reject_duplicate_block_keys
    @test placement_source_audit.reviewed_placement_facts_local_ordering_contract_status ===
          :available_cpb_overlap_local_ordering_contract
    @test placement_source_audit.reviewed_placement_facts_local_ordering_contract_blocker ===
          nothing
    @test placement_source_audit.reviewed_placement_facts_local_ordering_contract ===
          :parent_compatible_x_slowest_z_fastest
    @test placement_source_audit.reviewed_placement_facts_provided_local_orderings ===
          (:parent_compatible_x_slowest_z_fastest,)
    @test placement_source_audit.reviewed_placement_facts_mismatched_local_ordering_block_keys ===
          ()
    @test placement_source_audit.reviewed_placement_facts_global_dimension_source_contract_status ===
          :not_checked_cpb_overlap_global_dimension_source_contract
    @test placement_source_audit.reviewed_placement_facts_global_dimension_source_contract_blocker ===
          nothing
    @test placement_source_audit.reviewed_placement_facts_required_global_dimension_source ===
          :report_retained_units
    @test placement_source_audit.reviewed_placement_facts_provided_global_dimension_sources ===
          ()
    @test placement_source_audit.reviewed_placement_facts_mismatched_global_dimension_source_block_keys ===
          ()
    @test placement_source_audit.reviewed_placement_facts_skeleton_status ===
          :blocked_private_global_overlap_placement_plan_skeleton
    @test placement_source_audit.reviewed_placement_facts_skeleton_blocker ===
          :missing_placement_or_retained_transform
    @test placement_source_audit.reviewed_placement_facts_skeleton_global_overlap_status ===
          :blocked
    @test !placement_source_audit.reviewed_placement_facts_skeleton_route_driver_wiring
    @test !placement_source_audit.reviewed_placement_facts_skeleton_global_matrix_materialized
    @test !placement_source_audit.reviewed_placement_facts_skeleton_private_input_facts_available
    @test !placement_source_audit.reviewed_placement_facts_skeleton_route_global_stage_source
    @test placement_source_audit.global_overlap_status === :blocked
    @test placement_source_audit.global_overlap_blocker ===
          :missing_placement_or_retained_transform
    @test !placement_source_audit.global_matrix_materialized
    @test !placement_source_audit.route_driver_wiring
    @test !placement_source_audit.route_global_overlap_stage_source
    @test local_cpb_overlap_fingerprint.route_driver_wiring === false
    @test local_cpb_overlap_fingerprint.global_matrix_materialized === false
    @test local_cpb_overlap_fingerprint.route_global_overlap_stage_source === false

    stage =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            report;
            private_global_overlap_requested = true,
        )
    @test !isnothing(stage.private_global_overlap_result)
    @test stage.private_global_overlap_result.status ===
          :blocked_route_global_overlap_matrix
    @test stage.private_global_overlap_result.blocker ===
          :missing_pair_block_materialization_plan
    @test stage.private_global_overlap_result.global_dimension === nothing
    @test !stage.private_global_overlap_result.global_overlap_matrix_materialized
    @test !stage.private_global_overlap_result.global_one_body_term_matrix_materialized
    @test stage.private_global_overlap_summary.status ===
          :blocked_route_global_overlap_matrix
    @test stage.private_global_overlap_summary.blocker ===
          :missing_pair_block_materialization_plan
    @test stage.private_global_overlap_summary.result_available
    @test stage.private_global_overlap_summary.global_dimension === nothing
    @test !stage.private_global_overlap_summary.global_overlap_matrix_materialized
    @test !stage.private_global_overlap_summary.global_one_body_term_matrix_materialized
    _test_driver_overlap_nonclaim_flags(stage.private_global_overlap_result)
end

@testset "PQS route driver private global overlap option config" begin
    overrides = _driver_overlap_example_overrides()
    inputs = overrides.private_global_overlap_inputs

    @test overrides.private_global_overlap_requested === true
    @test overrides.private_global_overlap_global_dimension == 2
    @test haskey(inputs, :parent_axis_counts)
    @test haskey(inputs, :overlap_1d)
    @test inputs.parent_axis_counts == (2, 2, 2)
    @test all(key -> haskey(inputs.overlap_1d, key), (:x, :y, :z))
    @test size(inputs.overlap_1d.x) == (2, 2)
    @test size(inputs.overlap_1d.y) == (2, 2)
    @test size(inputs.overlap_1d.z) == (2, 2)

    off_stage =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            _driver_overlap_stage_report(),
        )
    @test isnothing(off_stage.private_global_overlap_result)
    @test isnothing(off_stage.private_global_overlap_summary)

    missing_inputs =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            _driver_overlap_stage_report();
            private_global_overlap_requested =
                overrides.private_global_overlap_requested,
            private_global_overlap_global_dimension =
                overrides.private_global_overlap_global_dimension,
        )
    @test missing_inputs.private_global_overlap_summary.blocker ===
          :missing_final_retained_column_layout

    materialized =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_stage(
            _driver_overlap_stage_report();
            private_global_overlap_requested =
                overrides.private_global_overlap_requested,
            private_global_overlap_global_dimension =
                overrides.private_global_overlap_global_dimension,
            private_global_overlap_inputs = inputs,
        )
    @test materialized.private_global_overlap_result.status ===
          :materialized_route_global_overlap_matrix
    @test materialized.private_global_overlap_result.global_matrix_result.matrix ≈
          _driver_overlap_expected_matrix()
    @test materialized.private_global_overlap_summary.global_overlap_matrix_materialized
    _test_driver_overlap_nonclaim_flags(materialized.private_global_overlap_result)
end
