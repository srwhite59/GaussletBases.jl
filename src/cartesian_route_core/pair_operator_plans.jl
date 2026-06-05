# Metadata-only pair-operator planning records.
#
# These records describe the mathematical path for a pair of final retained
# units before any numerical source block, realization matrix, final block, or
# Hamiltonian data is materialized.

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

const _PAIR_OPERATOR_DIRECT_SOURCE_READY_PATHS = (
    :direct_identity_cpb_path,
    :direct_source_cpb_1d_factor_path,
    :pqs_source_cpb_1d_factor_path,
)

function _pair_operator_supported_terms(terms)
    term_tuple = Tuple(terms)
    all(term -> term isa Symbol, term_tuple) ||
        throw(ArgumentError("supported_terms must contain Symbols"))
    return term_tuple
end

function _pair_operator_status(blocker)
    return isnothing(blocker) ?
        :metadata_only_not_materialized :
        :blocked_metadata_only_not_materialized
end

function _pair_operator_first_blocker(blockers...)
    for candidate in blockers
        isnothing(candidate) || return candidate
    end
    return nothing
end

function _pair_operator_unit_is_pqs(unit::FinalRetainedUnit)
    return lowering_recipe(unit) === :pqs_filled_source_cpb ||
        unit.intermediate.retained_rule === :pqs_boundary_comx_product_modes ||
        unit.shell_realization.realization_kind === :shell_projection_lowdin
end

function _pair_operator_unit_is_direct_identity(unit::FinalRetainedUnit)
    return lowering_recipe(unit) === :direct_identity_cpb &&
        unit.intermediate.retained_rule === :identity_source_modes
end

function _pair_operator_unit_is_direct_boundary_slab_set(unit::FinalRetainedUnit)
    return lowering_recipe(unit) === :direct_boundary_slab_set &&
        unit.intermediate.retained_rule === :direct_slab_set_identity_modes
end

function _pair_operator_unit_is_white_lindsey_boundary_strata(unit::FinalRetainedUnit)
    return lowering_recipe(unit) === :white_lindsey_boundary_strata &&
        unit.intermediate.retained_rule === :white_lindsey_boundary_stratum_product
end

function _pair_operator_unit_is_aggregate_subtree(unit::FinalRetainedUnit)
    return lowering_recipe(unit) === :white_lindsey_atom_local_child_shellification ||
        unit.intermediate.retained_rule === :atom_local_child_shellification_sequence
end

function _infer_source_operator_path(left::FinalRetainedUnit, right::FinalRetainedUnit)
    if _pair_operator_unit_is_aggregate_subtree(left) ||
       _pair_operator_unit_is_aggregate_subtree(right)
        return :aggregate_subtree_adapter_required,
            :aggregate_subtree_operator_plan_required
    end
    if _pair_operator_unit_is_pqs(left) || _pair_operator_unit_is_pqs(right)
        return :pqs_source_cpb_1d_factor_path, nothing
    end
    if _pair_operator_unit_is_direct_identity(left) &&
       _pair_operator_unit_is_direct_identity(right)
        return :direct_identity_cpb_path, nothing
    end
    if _pair_operator_unit_is_direct_boundary_slab_set(left) ||
       _pair_operator_unit_is_direct_boundary_slab_set(right)
        return :direct_boundary_slab_set_adapter_path, nothing
    end
    if _pair_operator_unit_is_white_lindsey_boundary_strata(left) ||
       _pair_operator_unit_is_white_lindsey_boundary_strata(right)
        return :white_lindsey_boundary_stratum_adapter_path, nothing
    end
    return :pending_pair_operator_path, :unclassified_pair_operator_path
end

function _infer_realization_path(unit::FinalRetainedUnit)
    if unit.shell_realization.realization_kind === :shell_projection_lowdin
        return :shell_projection_lowdin_planned
    end
    if unit.shell_realization.realization_kind === :direct_or_trivial_embedding
        return :identity_or_trivial_embedding
    end
    return :planned_shell_realization_path
end

function _infer_final_block_path(
    source_plan,
    left_realization_plan,
    right_realization_plan,
)
    source_plan.blocker === nothing || return :blocked_final_pair_block_path
    left_realization_plan.blocker === nothing ||
        return :blocked_final_pair_block_path
    right_realization_plan.blocker === nothing ||
        return :blocked_final_pair_block_path
    if realization_path(left_realization_plan) === :shell_projection_lowdin_planned ||
       realization_path(right_realization_plan) === :shell_projection_lowdin_planned
        return :source_block_realization_then_final_block
    end
    return :source_block_direct_to_final_block
end

struct SourceOperatorPlan
    source_operator_path::Symbol
    left::FinalRetainedUnit
    right::FinalRetainedUnit
    left_lowering_recipe::Symbol
    right_lowering_recipe::Symbol
    left_source_cpbs::Tuple{Vararg{CoordinateProductBox}}
    right_source_cpbs::Tuple{Vararg{CoordinateProductBox}}
    supported_terms::Tuple{Vararg{Symbol}}
    materialization_status::Symbol
    materialized::Bool
    blocker::Union{Symbol,Nothing}
    metadata::NamedTuple
end

struct RealizationApplicationPlan
    side::Symbol
    realization_path::Symbol
    final_unit::FinalRetainedUnit
    intermediate::IntermediateRetainedSpace
    shell_realization::ShellRealization
    supported_terms::Tuple{Vararg{Symbol}}
    materialization_status::Symbol
    materialized::Bool
    blocker::Union{Symbol,Nothing}
    metadata::NamedTuple
end

struct FinalPairBlockPlan
    final_block_path::Symbol
    left::FinalRetainedUnit
    right::FinalRetainedUnit
    source_operator_plan::SourceOperatorPlan
    left_realization_plan::RealizationApplicationPlan
    right_realization_plan::RealizationApplicationPlan
    supported_terms::Tuple{Vararg{Symbol}}
    materialization_status::Symbol
    materialized::Bool
    blocker::Union{Symbol,Nothing}
    metadata::NamedTuple
end

struct PairOperatorPlan
    pair::UnitPair
    left::FinalRetainedUnit
    right::FinalRetainedUnit
    source_operator_plan::SourceOperatorPlan
    left_realization_plan::RealizationApplicationPlan
    right_realization_plan::RealizationApplicationPlan
    final_pair_block_plan::FinalPairBlockPlan
    supported_terms::Tuple{Vararg{Symbol}}
    materialization_status::Symbol
    materialized::Bool
    blocker::Union{Symbol,Nothing}
    metadata::NamedTuple
end

struct PairOperatorPlanInventory
    pair_inventory::UnitPairInventory
    plans::Tuple{Vararg{PairOperatorPlan}}
    supported_terms::Tuple{Vararg{Symbol}}
    materialization_status::Symbol
    materialized::Bool
    blocker::Union{Symbol,Nothing}
    metadata::NamedTuple
end

function source_operator_plan(
    left::FinalRetainedUnit,
    right::FinalRetainedUnit;
    source_operator_path = nothing,
    supported_terms = _PAIR_OPERATOR_PLAN_DEFAULT_TERMS,
    materialization_status = nothing,
    materialized::Bool = false,
    blocker = nothing,
    metadata = (;),
)
    inferred_path, inferred_blocker = _infer_source_operator_path(left, right)
    path = isnothing(source_operator_path) ? inferred_path : source_operator_path
    plan_blocker = isnothing(blocker) ? inferred_blocker : blocker
    if (_pair_operator_unit_is_aggregate_subtree(left) ||
        _pair_operator_unit_is_aggregate_subtree(right)) &&
       path in _PAIR_OPERATOR_DIRECT_SOURCE_READY_PATHS
        throw(ArgumentError("aggregate/subtree pair plans cannot claim direct source readiness"))
    end
    status =
        isnothing(materialization_status) ?
        _pair_operator_status(plan_blocker) :
        materialization_status
    return SourceOperatorPlan(
        path,
        left,
        right,
        lowering_recipe(left),
        lowering_recipe(right),
        source_cpbs(left),
        source_cpbs(right),
        _pair_operator_supported_terms(supported_terms),
        status,
        materialized,
        plan_blocker,
        NamedTuple(metadata),
    )
end

function realization_application_plan(
    side::Symbol,
    unit::FinalRetainedUnit;
    realization_path = nothing,
    supported_terms = _PAIR_OPERATOR_PLAN_DEFAULT_TERMS,
    materialization_status = nothing,
    materialized::Bool = false,
    blocker = nothing,
    metadata = (;),
)
    side in (:left, :right) ||
        throw(ArgumentError("realization side must be :left or :right"))
    path = isnothing(realization_path) ? _infer_realization_path(unit) : realization_path
    status =
        isnothing(materialization_status) ?
        _pair_operator_status(blocker) :
        materialization_status
    return RealizationApplicationPlan(
        side,
        path,
        unit,
        unit.intermediate,
        unit.shell_realization,
        _pair_operator_supported_terms(supported_terms),
        status,
        materialized,
        blocker,
        NamedTuple(metadata),
    )
end

function final_pair_block_plan(
    left::FinalRetainedUnit,
    right::FinalRetainedUnit,
    source_plan::SourceOperatorPlan,
    left_realization_plan::RealizationApplicationPlan,
    right_realization_plan::RealizationApplicationPlan;
    final_block_path = nothing,
    supported_terms = _PAIR_OPERATOR_PLAN_DEFAULT_TERMS,
    materialization_status = nothing,
    materialized::Bool = false,
    blocker = nothing,
    metadata = (;),
)
    source_plan.left === left && source_plan.right === right ||
        throw(ArgumentError("source plan units do not match final pair block units"))
    left_realization_plan.final_unit === left ||
        throw(ArgumentError("left realization plan does not match left unit"))
    right_realization_plan.final_unit === right ||
        throw(ArgumentError("right realization plan does not match right unit"))

    inferred_blocker = _pair_operator_first_blocker(
        source_plan.blocker,
        left_realization_plan.blocker,
        right_realization_plan.blocker,
    )
    plan_blocker = isnothing(blocker) ? inferred_blocker : blocker
    path =
        isnothing(final_block_path) ?
        _infer_final_block_path(
            source_plan,
            left_realization_plan,
            right_realization_plan,
        ) :
        final_block_path
    status =
        isnothing(materialization_status) ?
        _pair_operator_status(plan_blocker) :
        materialization_status

    return FinalPairBlockPlan(
        path,
        left,
        right,
        source_plan,
        left_realization_plan,
        right_realization_plan,
        _pair_operator_supported_terms(supported_terms),
        status,
        materialized,
        plan_blocker,
        NamedTuple(metadata),
    )
end

function pair_operator_plan(
    pair::UnitPair;
    supported_terms = _PAIR_OPERATOR_PLAN_DEFAULT_TERMS,
    materialization_status = nothing,
    materialized::Bool = false,
    metadata = (;),
)
    terms = _pair_operator_supported_terms(supported_terms)
    source_plan = source_operator_plan(
        pair.left,
        pair.right;
        supported_terms = terms,
        materialized,
        metadata = (; plan_layer = :source_operator, pair_key = pair.pair_key),
    )
    left_realization = realization_application_plan(
        :left,
        pair.left;
        supported_terms = terms,
        materialized,
        metadata = (; plan_layer = :left_realization, pair_key = pair.pair_key),
    )
    right_realization = realization_application_plan(
        :right,
        pair.right;
        supported_terms = terms,
        materialized,
        metadata = (; plan_layer = :right_realization, pair_key = pair.pair_key),
    )
    final_block = final_pair_block_plan(
        pair.left,
        pair.right,
        source_plan,
        left_realization,
        right_realization;
        supported_terms = terms,
        materialized,
        metadata = (; plan_layer = :final_pair_block, pair_key = pair.pair_key),
    )
    plan_blocker = _pair_operator_first_blocker(
        blocker(source_plan),
        blocker(left_realization),
        blocker(right_realization),
        blocker(final_block),
    )
    status =
        isnothing(materialization_status) ?
        _pair_operator_status(plan_blocker) :
        materialization_status

    return PairOperatorPlan(
        pair,
        pair.left,
        pair.right,
        source_plan,
        left_realization,
        right_realization,
        final_block,
        terms,
        status,
        materialized,
        plan_blocker,
        NamedTuple(metadata),
    )
end

function pair_operator_plan_inventory(
    pair_inventory::UnitPairInventory;
    supported_terms = _PAIR_OPERATOR_PLAN_DEFAULT_TERMS,
    materialization_status = nothing,
    materialized::Bool = false,
    metadata = (;),
)
    terms = _pair_operator_supported_terms(supported_terms)
    plans = Tuple(
        pair_operator_plan(
            pair;
            supported_terms = terms,
            materialized,
            metadata = (; inventory_source = :unit_pair_inventory),
        ) for pair in pair_entries(pair_inventory)
    )
    plan_blocker = _pair_operator_first_blocker((blocker(plan) for plan in plans)...)
    status =
        isnothing(materialization_status) ?
        _pair_operator_status(plan_blocker) :
        materialization_status
    return PairOperatorPlanInventory(
        pair_inventory,
        plans,
        terms,
        status,
        materialized,
        plan_blocker,
        NamedTuple(metadata),
    )
end

source_operator_path(plan::SourceOperatorPlan) = plan.source_operator_path
source_operator_path(plan::PairOperatorPlan) =
    source_operator_path(plan.source_operator_plan)

realization_path(plan::RealizationApplicationPlan) = plan.realization_path
realization_path(plan::PairOperatorPlan) = (;
    left = realization_path(plan.left_realization_plan),
    right = realization_path(plan.right_realization_plan),
)
left_realization_path(plan::PairOperatorPlan) =
    realization_path(plan.left_realization_plan)
right_realization_path(plan::PairOperatorPlan) =
    realization_path(plan.right_realization_plan)

final_block_path(plan::FinalPairBlockPlan) = plan.final_block_path
final_block_path(plan::PairOperatorPlan) =
    final_block_path(plan.final_pair_block_plan)

supported_terms(plan::SourceOperatorPlan) = plan.supported_terms
supported_terms(plan::RealizationApplicationPlan) = plan.supported_terms
supported_terms(plan::FinalPairBlockPlan) = plan.supported_terms
supported_terms(plan::PairOperatorPlan) = plan.supported_terms
supported_terms(inventory::PairOperatorPlanInventory) = inventory.supported_terms

materialization_status(plan::SourceOperatorPlan) = plan.materialization_status
materialization_status(plan::RealizationApplicationPlan) =
    plan.materialization_status
materialization_status(plan::FinalPairBlockPlan) = plan.materialization_status
materialization_status(plan::PairOperatorPlan) = plan.materialization_status
materialization_status(inventory::PairOperatorPlanInventory) =
    inventory.materialization_status

blocker(plan::SourceOperatorPlan) = plan.blocker
blocker(plan::RealizationApplicationPlan) = plan.blocker
blocker(plan::FinalPairBlockPlan) = plan.blocker
blocker(plan::PairOperatorPlan) = plan.blocker
blocker(inventory::PairOperatorPlanInventory) = inventory.blocker

pair_operator_plans(inventory::PairOperatorPlanInventory) = inventory.plans
pair_operator_plan_count(inventory::PairOperatorPlanInventory) =
    length(inventory.plans)

function pair_operator_materialization_readiness_requirements()
    return (
        :nonempty_pair_operator_plan_inventory,
        :no_blocked_typed_pair_operator_plans,
        :no_pending_pair_operator_source_paths,
        :no_already_materialized_typed_pair_operator_plans,
    )
end

function pair_operator_materialization_readiness(
    inventory::PairOperatorPlanInventory,
)
    plans = pair_operator_plans(inventory)
    plan_count = pair_operator_plan_count(inventory)
    blocked_count = count(plan -> !isnothing(blocker(plan)), plans)
    materialized_count = count(plan -> plan.materialized, plans)
    pending_count =
        count(plan -> source_operator_path(plan) === :pending_pair_operator_path, plans)
    requirements = pair_operator_materialization_readiness_requirements()

    if plan_count == 0
        return (;
            ready = false,
            status = :blocked_empty_pair_operator_plan_inventory,
            blocker = :empty_pair_operator_plan_inventory,
            requirements,
            plan_count,
            blocked_count,
            materialized_count,
        )
    end

    if blocked_count > 0
        first_blocked_plan =
            first(plan for plan in plans if !isnothing(blocker(plan)))
        return (;
            ready = false,
            status = :blocked_pair_operator_materialization,
            blocker = blocker(first_blocked_plan),
            requirements,
            plan_count,
            blocked_count,
            materialized_count,
        )
    end

    if pending_count > 0
        return (;
            ready = false,
            status = :blocked_pending_pair_operator_source_path,
            blocker = :pending_pair_operator_path,
            requirements,
            plan_count,
            blocked_count,
            materialized_count,
        )
    end

    if materialized_count > 0
        return (;
            ready = false,
            status = :blocked_already_materialized_pair_operator_plans,
            blocker = :typed_pair_operator_plans_already_materialized,
            requirements,
            plan_count,
            blocked_count,
            materialized_count,
        )
    end

    return (;
        ready = true,
        status = :ready_pair_operator_materialization,
        blocker = nothing,
        requirements,
        plan_count,
        blocked_count,
        materialized_count,
    )
end

function _pair_operator_count_entries(values, value_field::Symbol)
    counts = Dict{Any,Int}()
    order = Any[]
    for value in values
        if !haskey(counts, value)
            counts[value] = 0
            push!(order, value)
        end
        counts[value] += 1
    end
    return Tuple(
        NamedTuple{(value_field, :pair_count)}((value, counts[value]))
        for value in order
    )
end

function pair_operator_source_path_counts(inventory::PairOperatorPlanInventory)
    return _pair_operator_count_entries(
        (source_operator_path(plan) for plan in pair_operator_plans(inventory)),
        :source_operator_path,
    )
end

function pair_operator_final_block_path_counts(inventory::PairOperatorPlanInventory)
    return _pair_operator_count_entries(
        (final_block_path(plan) for plan in pair_operator_plans(inventory)),
        :final_block_path,
    )
end

function pair_operator_materialization_status_counts(
    inventory::PairOperatorPlanInventory,
)
    return _pair_operator_count_entries(
        (materialization_status(plan) for plan in pair_operator_plans(inventory)),
        :materialization_status,
    )
end

function pair_operator_blocker_counts(inventory::PairOperatorPlanInventory)
    return _pair_operator_count_entries(
        (blocker(plan) for plan in pair_operator_plans(inventory)),
        :blocker,
    )
end

function pair_operator_plan_family_counts(inventory::PairOperatorPlanInventory)
    counts = Dict{Any,Int}()
    order = Any[]
    for plan in pair_operator_plans(inventory)
        key = (
            source_operator_path(plan),
            final_block_path(plan),
            materialization_status(plan),
            blocker(plan),
            plan.materialized,
        )
        if !haskey(counts, key)
            counts[key] = 0
            push!(order, key)
        end
        counts[key] += 1
    end
    return Tuple(
        (;
            source_operator_path = key[1],
            final_block_path = key[2],
            materialization_status = key[3],
            blocker = key[4],
            materialized = key[5],
            pair_count = counts[key],
        ) for key in order
    )
end
