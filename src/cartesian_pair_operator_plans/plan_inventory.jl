# Retained-unit pair to metadata-only pair-operator plan rules.

const _PAIR_OPERATOR_PLAN_DEFAULT_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

const _DIRECT_RETAINED_UNIT_KINDS = (:direct_cpb_retained_unit,)

"""
    pair_operator_plan(unit_pair_plan; policy = MetadataOnlyPairOperatorPlans())

Build metadata-only pair-operator plans from retained-unit pairs.

The result classifies source-space, realization, and final-block paths. It does
not materialize source operator blocks, final pair blocks, Hamiltonian data, or
artifacts.
"""
function pair_operator_plan(
    unit_pair_plan::CUP.UnitPairPlan;
    policy::PairOperatorPlanPolicy = MetadataOnlyPairOperatorPlans(),
    supported_terms = _PAIR_OPERATOR_PLAN_DEFAULT_TERMS,
    metadata = (;),
)
    terms = _supported_terms_tuple(supported_terms)
    route_core_plan = _route_core_pair_operator_plan_or_nothing(
        unit_pair_plan;
        supported_terms = terms,
    )
    route_core_plans =
        isnothing(route_core_plan.inventory) ?
        nothing :
        CRC.pair_operator_plans(route_core_plan.inventory)

    records = Tuple(
        _pair_operator_plan_record(
            pair,
            terms,
            isnothing(route_core_plans) ? nothing : route_core_plans[pair.pair_index],
            route_core_plan,
        ) for pair in CUP.unit_pairs(unit_pair_plan)
    )
    plan_summary =
        _pair_operator_plan_summary(policy, unit_pair_plan, records, route_core_plan)
    return PairOperatorPlan(
        policy,
        unit_pair_plan,
        records,
        route_core_plan.inventory,
        plan_summary,
        NamedTuple(metadata),
    )
end

function _supported_terms_tuple(supported_terms)
    terms = Tuple(supported_terms)
    all(term -> term isa Symbol, terms) ||
        throw(ArgumentError("supported_terms must contain Symbols"))
    return terms
end

function _route_core_pair_operator_plan_or_nothing(
    unit_pair_plan::CUP.UnitPairPlan;
    supported_terms,
)
    route_core_pairs = CUP.route_core_pair_inventory(unit_pair_plan)
    if isnothing(route_core_pairs)
        pair_summary = CUP.summary(unit_pair_plan)
        return (;
            inventory = nothing,
            status = :blocked_missing_route_core_pair_inventory,
            blocker = pair_summary.route_core_pair_inventory_blocker,
        )
    end

    try
        inventory = CRC.pair_operator_plan_inventory(
            route_core_pairs;
            supported_terms,
            metadata = (; source = :cartesian_pair_operator_plans),
        )
        return (;
            inventory,
            status = :available_route_core_pair_operator_plan_inventory,
            blocker = nothing,
        )
    catch err
        return (;
            inventory = nothing,
            status = :blocked_route_core_pair_operator_plan_error,
            blocker = sprint(showerror, err),
        )
    end
end

function _pair_operator_plan_record(
    pair::CUP.UnitPairRecord,
    supported_terms,
    route_core_sidecar,
    route_core_plan,
)
    source_path, source_blocker = _source_operator_path(pair)
    realization = (;
        left = _realization_path(pair.left_unit_kind),
        right = _realization_path(pair.right_unit_kind),
    )
    blocker =
        isnothing(route_core_plan.inventory) ?
        _first_blocker(route_core_plan.blocker, source_blocker) :
        source_blocker
    final_path = _final_block_path(pair, blocker)
    status =
        isnothing(blocker) ?
        :metadata_only_not_materialized :
        :blocked_metadata_only_not_materialized

    return PairOperatorPlanRecord(
        pair.pair_key,
        pair.pair_index,
        pair.pair_family,
        pair.left_unit_key,
        pair.right_unit_key,
        pair.left_unit_kind,
        pair.right_unit_kind,
        source_path,
        realization,
        final_path,
        supported_terms,
        status,
        blocker,
        route_core_sidecar,
        false,
        (;
            route_core_pair_operator_sidecar_status =
                isnothing(route_core_sidecar) ?
                route_core_plan.status :
                :available_route_core_pair_operator_plan,
            route_core_pair_operator_sidecar_blocker =
                isnothing(route_core_sidecar) ? route_core_plan.blocker : nothing,
        ),
    )
end

function _source_operator_path(pair::CUP.UnitPairRecord)
    left = pair.left_unit_kind
    right = pair.right_unit_kind
    if _unit_kind_is_pqs(left) || _unit_kind_is_pqs(right)
        return :pqs_source_cpb_1d_factor_path, nothing
    end
    if _unit_kind_is_distorted(left) || _unit_kind_is_distorted(right)
        return :distorted_product_box_operator_path, nothing
    end
    if _unit_kind_is_white_lindsey(left) || _unit_kind_is_white_lindsey(right)
        return :white_lindsey_boundary_stratum_adapter_path, nothing
    end
    if _unit_kind_is_direct(left) && _unit_kind_is_direct(right)
        return :direct_identity_cpb_path, nothing
    end
    if _unit_kind_is_direct(left) || _unit_kind_is_direct(right)
        return :direct_source_cpb_1d_factor_path, nothing
    end
    return :pending_pair_operator_path, :unclassified_pair_operator_path
end

_unit_kind_is_direct(unit_kind::Symbol) = unit_kind in _DIRECT_RETAINED_UNIT_KINDS
_unit_kind_is_pqs(unit_kind::Symbol) = unit_kind === :pqs_shell_retained_unit
_unit_kind_is_white_lindsey(unit_kind::Symbol) =
    unit_kind === :white_lindsey_boundary_stratum_retained_unit
_unit_kind_is_distorted(unit_kind::Symbol) =
    unit_kind === :distorted_product_box_retained_unit

function _realization_path(unit_kind::Symbol)
    if _unit_kind_is_pqs(unit_kind)
        return :shell_projection_lowdin_planned
    end
    if _unit_kind_is_distorted(unit_kind)
        return :distorted_product_realization_planned
    end
    if _unit_kind_is_direct(unit_kind) || _unit_kind_is_white_lindsey(unit_kind)
        return :identity_or_trivial_embedding
    end
    return :planned_shell_realization_path
end

function _final_block_path(pair::CUP.UnitPairRecord, blocker)
    isnothing(blocker) || return :blocked_final_pair_block_path
    if _unit_kind_is_pqs(pair.left_unit_kind) || _unit_kind_is_pqs(pair.right_unit_kind)
        return :source_block_realization_then_final_block
    end
    return :source_block_direct_to_final_block
end

function _first_blocker(blockers...)
    for blocker in blockers
        isnothing(blocker) || return blocker
    end
    return nothing
end
