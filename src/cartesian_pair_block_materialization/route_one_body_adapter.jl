# Private route-shaped local one-body block adapter.
#
# This bridge accepts structured pair-block state, delegates numerical work to
# the existing mixed one-body block-set consumer, and returns the existing local
# block collection vocabulary. It does not assemble global operators,
# Hamiltonians, Coulomb data, IDA/MWG data, exports, artifacts, or PQS
# shell/Lowdin realizations.

function route_local_one_body_block_collection(
    plan::PairBlockMaterializationPlan;
    terms = _ONE_BODY_TERMS,
    inputs = (;),
    provider = nothing,
    materialize_terms = (),
)
    block_set_consumption = _one_body_pair_block_set_consumption(
        plan;
        terms,
        inputs,
        provider,
        materialize_terms,
    )
    block_set_consumption_summary =
        _one_body_pair_block_set_consumption_summary(block_set_consumption)
    local_block_collection =
        _one_body_local_block_collection(block_set_consumption)
    local_block_collection_summary =
        _one_body_local_block_collection_summary(local_block_collection)

    return (;
        object_kind = :cartesian_pair_block_route_local_one_body_block_adapter,
        status = local_block_collection_summary.status,
        blocker = local_block_collection_summary.blocker,
        pair_block_materialization_plan = plan,
        block_set_consumption,
        block_set_consumption_summary,
        local_block_collection,
        local_block_collection_summary,
        terms = local_block_collection_summary.terms,
        requested_terms = local_block_collection_summary.requested_terms,
        requested_materialize_terms =
            local_block_collection_summary.requested_materialize_terms,
        materialized_terms = local_block_collection_summary.materialized_terms,
        deferred_terms = local_block_collection_summary.deferred_terms,
        entry_count = local_block_collection_summary.entry_count,
        materialized_entry_count =
            local_block_collection_summary.materialized_entry_count,
        skipped_entry_count =
            local_block_collection_summary.skipped_entry_count,
        total_materialized_count =
            block_set_consumption_summary.total_materialized_count,
        total_skipped_count = block_set_consumption_summary.total_skipped_count,
        source_space_entry_count =
            local_block_collection_summary.source_space_entry_count,
        final_local_entry_count =
            local_block_collection_summary.final_local_entry_count,
        term_separated_entries =
            local_block_collection_summary.term_separated_entries,
        pair_separated_entries =
            local_block_collection_summary.pair_separated_entries,
        block_set_results_summed =
            local_block_collection_summary.block_set_results_summed,
        result_terms_remain_separated =
            block_set_consumption_summary.result_terms_remain_separated,
        block_matrices_copied_into_collection =
            local_block_collection.block_matrices_copied_into_collection,
        factors_constructed = block_set_consumption_summary.factors_constructed,
        source_operator_blocks_materialized =
            local_block_collection_summary.source_operator_blocks_materialized,
        final_pair_blocks_materialized =
            local_block_collection_summary.final_pair_blocks_materialized,
        _route_local_one_body_adapter_nonclaim_flags(
            local_block_collection_summary,
            block_set_consumption_summary,
        )...,
    )
end

function route_local_one_body_block_collection(source; kwargs...)
    return route_local_one_body_block_collection(
        _route_one_body_pair_block_materialization_plan(source);
        kwargs...,
    )
end

function _route_one_body_pair_block_materialization_plan(
    plan::PairBlockMaterializationPlan,
)
    return plan
end

function _route_one_body_pair_block_materialization_plan(
    pair_operator_plan::CPOP.PairOperatorPlan,
)
    return pair_block_materialization_plan(pair_operator_plan)
end

function _route_one_body_pair_block_materialization_plan(source)
    if hasproperty(source, :pair_block_materialization_plan)
        plan = getproperty(source, :pair_block_materialization_plan)
        plan isa PairBlockMaterializationPlan && return plan
        throw(
            ArgumentError(
                "route local one-body adapter requires pair_block_materialization_plan to be a PairBlockMaterializationPlan",
            ),
        )
    end

    if hasproperty(source, :terminal_route_state)
        return _route_one_body_terminal_route_state_plan(
            getproperty(source, :terminal_route_state),
        )
    end

    throw(
        ArgumentError(
            "route local one-body adapter requires a PairBlockMaterializationPlan, PairOperatorPlan, or object carrying pair_block_materialization_plan",
        ),
    )
end

function _route_one_body_terminal_route_state_plan(route_state)
    if hasproperty(route_state, :pair_block_materialization_plan)
        plan = getproperty(route_state, :pair_block_materialization_plan)
        plan isa PairBlockMaterializationPlan && return plan
        throw(
            ArgumentError(
                "terminal_route_state.pair_block_materialization_plan must be a PairBlockMaterializationPlan",
            ),
        )
    end

    if hasproperty(route_state, :pair_operator_plan)
        pair_operator_plan = getproperty(route_state, :pair_operator_plan)
        pair_operator_plan isa CPOP.PairOperatorPlan &&
            return pair_block_materialization_plan(pair_operator_plan)
        throw(
            ArgumentError(
                "terminal_route_state.pair_operator_plan must be a PairOperatorPlan",
            ),
        )
    end

    throw(
        ArgumentError(
            "terminal_route_state must carry pair_block_materialization_plan or pair_operator_plan",
        ),
    )
end

function _route_local_one_body_adapter_nonclaim_flags(
    local_summary,
    block_set_summary,
)
    return (;
        local_operator_assembled = local_summary.local_operator_assembled,
        global_operator_assembled = local_summary.global_operator_assembled,
        route_driver_wiring = local_summary.route_driver_wiring,
        operator_blocks_materialized =
            local_summary.operator_blocks_materialized,
        hamiltonian_data_materialized =
            local_summary.hamiltonian_data_materialized,
        artifacts_materialized = local_summary.artifacts_materialized,
        global_operator_blocks_materialized =
            local_summary.global_operator_blocks_materialized,
        global_hamiltonian_data_materialized =
            local_summary.global_hamiltonian_data_materialized,
        global_artifacts_materialized =
            local_summary.global_artifacts_materialized,
        exports_materialized = false,
        coulomb_materialized = local_summary.coulomb_materialized,
        density_density_materialized =
            local_summary.density_density_materialized,
        ida_mwg_data_materialized = local_summary.ida_mwg_data_materialized,
        pqs_lowdin_materialized = local_summary.pqs_lowdin_materialized,
        pqs_shell_projection_materialized =
            local_summary.pqs_shell_projection_materialized,
        full_white_lindsey_route_assembled =
            local_summary.full_white_lindsey_route_assembled,
        mixed_dispatcher_materialized =
            block_set_summary.mixed_dispatcher_materialized,
        numerical_blocks_materialized =
            block_set_summary.numerical_blocks_materialized,
    )
end
