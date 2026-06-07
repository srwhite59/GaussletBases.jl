# Compact pair-block materialization summaries.

function _count_by_value(values, field::Symbol)
    counts = Dict{Any,Int}()
    order = Any[]
    for value in values
        if !haskey(counts, value)
            counts[value] = 0
            push!(order, value)
        end
        counts[value] += 1
    end
    return Tuple((; field => value, count = counts[value]) for value in order)
end

function _pair_block_materialization_plan_summary(
    policy::PairBlockMaterializationPolicy,
    pair_operator_plan::CPOP.PairOperatorPlan,
    records,
)
    blockers = Tuple(record.blocker for record in records)
    ready_count =
        count(record -> record.readiness_status === :ready_metadata_only_not_materialized, records)
    blocked_count = length(records) - ready_count

    return (;
        object_kind = :cartesian_pair_block_materialization_plan_summary,
        status =
            blocked_count == 0 ?
            :ready_pair_block_materialization_plan :
            :blocked_pair_block_materialization_plan,
        blocker =
            blocked_count == 0 ?
            nothing :
            first(blocker for blocker in blockers if !isnothing(blocker)),
        policy_kind = policy_kind(policy),
        pair_operator_plan_count = length(CPOP.pair_operator_records(pair_operator_plan)),
        pair_block_record_count = length(records),
        ready_record_count = ready_count,
        blocked_record_count = blocked_count,
        materialization_path_counts =
            _count_by_value(
                (record.materialization_path for record in records),
                :materialization_path,
            ),
        readiness_status_counts =
            _count_by_value(
                (record.readiness_status for record in records),
                :readiness_status,
            ),
        blocker_counts = _count_by_value(blockers, :blocker),
        materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function unavailable_summary(status::Symbol, blocker = nothing)
    return (;
        object_kind = :cartesian_pair_block_materialization_plan_summary,
        status,
        blocker,
        pair_operator_plan_count = 0,
        pair_block_record_count = 0,
        ready_record_count = 0,
        blocked_record_count = 0,
        materialization_path_counts = (),
        readiness_status_counts = (),
        blocker_counts = (),
        materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end
