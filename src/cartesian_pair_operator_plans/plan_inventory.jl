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
    pair_operator_plan(unit_pair_plan, transform_contract_plan; policy)

Build metadata-only pair-operator plans from retained-unit pairs.

The result classifies source-space, realization, and final-block paths. It does
not materialize source operator blocks, final pair blocks, Hamiltonian data, or
artifacts.
"""
function pair_operator_plan(
    unit_pair_plan::CartesianUnitPairs.UnitPairPlan;
    policy::PairOperatorPlanPolicy = MetadataOnlyPairOperatorPlans(),
    supported_terms = _PAIR_OPERATOR_PLAN_DEFAULT_TERMS,
    route_core_sidecars::Bool = true,
    require_route_core_crosscheck::Bool = false,
    metadata = (;),
)
    transform_contract_plan =
        CartesianRetainedUnitTransformContracts.retained_unit_transform_contract_plan(
            unit_pair_plan.retained_unit_plan,
        )
    return pair_operator_plan(
        unit_pair_plan,
        transform_contract_plan;
        policy,
        supported_terms,
        route_core_sidecars,
        require_route_core_crosscheck,
        metadata,
    )
end

function pair_operator_plan(
    unit_pair_plan::CartesianUnitPairs.UnitPairPlan,
    transform_contract_plan::CartesianRetainedUnitTransformContracts.RetainedUnitTransformContractPlan;
    policy::PairOperatorPlanPolicy = MetadataOnlyPairOperatorPlans(),
    supported_terms = _PAIR_OPERATOR_PLAN_DEFAULT_TERMS,
    route_core_sidecars::Bool = true,
    require_route_core_crosscheck::Bool = false,
    metadata = (;),
)
    terms = _supported_terms_tuple(supported_terms)
    transform_lookup = _retained_unit_transform_contract_lookup(transform_contract_plan)
    route_core_plan =
        route_core_sidecars ?
        _route_core_pair_operator_plan_or_nothing(
        unit_pair_plan;
        supported_terms = terms,
        ) :
        _route_core_pair_operator_plan_not_requested()
    route_core_plans =
        isnothing(route_core_plan.inventory) ?
        nothing :
        CartesianRouteCore.pair_operator_plans(route_core_plan.inventory)

    records = Tuple(
        _pair_operator_plan_record(
            pair,
            terms,
            isnothing(route_core_plans) ? nothing : route_core_plans[pair.pair_index],
            route_core_plan,
            require_route_core_crosscheck,
            transform_lookup,
        ) for pair in CartesianUnitPairs.unit_pairs(unit_pair_plan)
    )
    plan_summary =
        _pair_operator_plan_summary(policy, unit_pair_plan, records, route_core_plan)
    return PairOperatorPlan(
        policy,
        unit_pair_plan,
        transform_contract_plan,
        records,
        route_core_plan.inventory,
        plan_summary,
        NamedTuple(metadata),
    )
end

function _route_core_pair_operator_plan_not_requested()
    return (;
        inventory = nothing,
        status = :not_requested_route_core_pair_operator_plan_inventory,
        blocker = nothing,
    )
end

function _retained_unit_transform_contract_lookup(
    transform_contract_plan::CartesianRetainedUnitTransformContracts.RetainedUnitTransformContractPlan,
)
    contracts_by_key = Dict{
        Symbol,
        CartesianRetainedUnitTransformContracts.RetainedUnitTransformContract,
    }()
    duplicate_keys = Symbol[]
    for contract in
            CartesianRetainedUnitTransformContracts.transform_contracts(transform_contract_plan)
        if haskey(contracts_by_key, contract.unit_key)
            push!(duplicate_keys, contract.unit_key)
        else
            contracts_by_key[contract.unit_key] = contract
        end
    end
    blocker =
        isempty(duplicate_keys) ?
        nothing :
        :duplicate_retained_unit_transform_contract
    return (;
        contracts_by_key,
        duplicate_unit_keys = Tuple(duplicate_keys),
        blocker,
    )
end

function _supported_terms_tuple(supported_terms)
    terms = Tuple(supported_terms)
    all(term -> term isa Symbol, terms) ||
        throw(ArgumentError("supported_terms must contain Symbols"))
    return terms
end

function _route_core_pair_operator_plan_or_nothing(
    unit_pair_plan::CartesianUnitPairs.UnitPairPlan;
    supported_terms,
)
    route_core_pairs = CartesianUnitPairs.route_core_pair_inventory(unit_pair_plan)
    if isnothing(route_core_pairs)
        pair_summary = CartesianUnitPairs.summary(unit_pair_plan)
        return (;
            inventory = nothing,
            status = :blocked_missing_route_core_pair_inventory,
            blocker = pair_summary.route_core_pair_inventory_blocker,
        )
    end

    try
        inventory = CartesianRouteCore.pair_operator_plan_inventory(
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
    pair::CartesianUnitPairs.UnitPairRecord,
    supported_terms,
    route_core_sidecar,
    route_core_plan,
    require_route_core_crosscheck::Bool,
    transform_lookup,
)
    source_path, source_blocker = _source_operator_path(pair)
    left_transform = _transform_contract_for_pair_unit(pair, :left, transform_lookup)
    right_transform = _transform_contract_for_pair_unit(pair, :right, transform_lookup)
    transform = (;
        left = left_transform.transform_path,
        right = right_transform.transform_path,
    )
    realization = (;
        left = left_transform.realization_path,
        right = right_transform.realization_path,
    )
    blocker =
        require_route_core_crosscheck && isnothing(route_core_plan.inventory) ?
        _first_blocker(
            transform_lookup.blocker,
            left_transform.blocker,
            right_transform.blocker,
            route_core_plan.blocker,
            source_blocker,
        ) :
        _first_blocker(
            transform_lookup.blocker,
            left_transform.blocker,
            right_transform.blocker,
            source_blocker,
        )
    final_path = _final_block_path(realization, blocker)
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
        transform,
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
            duplicate_transform_contract_unit_keys =
                transform_lookup.duplicate_unit_keys,
        ),
    )
end

function _transform_contract_for_pair_unit(
    pair::CartesianUnitPairs.UnitPairRecord,
    side::Symbol,
    transform_lookup,
)
    side in (:left, :right) ||
        throw(ArgumentError("pair transform contract side must be :left or :right"))
    unit_key = side === :left ? pair.left_unit_key : pair.right_unit_key
    unit_index = side === :left ? pair.left_index : pair.right_index
    unit_kind = side === :left ? pair.left_unit_kind : pair.right_unit_kind
    missing_blocker =
        side === :left ?
        :missing_left_transform_contract :
        :missing_right_transform_contract
    contract = get(transform_lookup.contracts_by_key, unit_key, nothing)
    if isnothing(contract)
        return (;
            transform_path = :missing_retained_unit_transform_contract,
            realization_path = :missing_retained_unit_realization_contract,
            blocker = missing_blocker,
        )
    end
    if contract.unit_index != unit_index ||
       contract.unit_key != unit_key ||
       contract.unit_kind != unit_kind
        return (;
            transform_path = contract.transform_path,
            realization_path = contract.realization_path,
            blocker = :transform_contract_unit_mismatch,
        )
    end
    return (;
        transform_path = contract.transform_path,
        realization_path = contract.realization_path,
        blocker = contract.blocker,
    )
end

function _source_operator_path(pair::CartesianUnitPairs.UnitPairRecord)
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

function _final_block_path(realization::NamedTuple, blocker)
    isnothing(blocker) || return :blocked_final_pair_block_path
    if realization.left === :shell_projection_lowdin_planned ||
       realization.right === :shell_projection_lowdin_planned
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
