# PQS/PQS raw source-space one-body term selector.

"""
    pqs_source_pair_one_body_block(record, term; overlap_1d, ...)

Materialize one supported PQS/PQS raw source-space one-body block by delegating
to the term-specific source helper. The accepted terms are `:overlap`,
`:position_x/y/z`, `:x2_x/y/z`, and `:kinetic`.
"""
function pqs_source_pair_one_body_block(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    overlap_1d,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    descriptor = _supported_pqs_source_safe_term_descriptor(term)
    if descriptor.family === :overlap
        return pqs_source_pair_overlap_block(record; overlap_1d)
    end

    if descriptor.family === :position
        isnothing(position_1d) &&
            throw(ArgumentError(_pqs_source_required_factor_message(descriptor)))
        return pqs_source_pair_position_block(
            record;
            axis = descriptor.axis,
            overlap_1d,
            position_1d,
        )
    end

    if descriptor.family === :x2
        isnothing(x2_1d) &&
            throw(ArgumentError(_pqs_source_required_factor_message(descriptor)))
        return pqs_source_pair_x2_block(
            record;
            axis = descriptor.axis,
            overlap_1d,
            x2_1d,
        )
    end

    if descriptor.family === :kinetic
        isnothing(kinetic_1d) &&
            throw(ArgumentError(_pqs_source_required_factor_message(descriptor)))
        return pqs_source_pair_kinetic_block(
            record;
            overlap_1d,
            kinetic_1d,
        )
    end

    throw(ArgumentError("unsupported PQS source one-body term: $(term)"))
end

"""
    pqs_source_pair_one_body_blocks(plan, term; overlap_1d, ...)

Materialize supported PQS/PQS raw source-space one-body blocks for a plan by
delegating to the term-specific source plan helper.
"""
function pqs_source_pair_one_body_blocks(
    plan::PairBlockMaterializationPlan,
    term::Symbol;
    overlap_1d,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    descriptor = _supported_pqs_source_safe_term_descriptor(term)
    if descriptor.family === :overlap
        return pqs_source_pair_overlap_blocks(plan; overlap_1d)
    end

    if descriptor.family === :position
        isnothing(position_1d) &&
            throw(ArgumentError(_pqs_source_required_factor_message(descriptor)))
        return pqs_source_pair_position_blocks(
            plan;
            axis = descriptor.axis,
            overlap_1d,
            position_1d,
        )
    end

    if descriptor.family === :x2
        isnothing(x2_1d) &&
            throw(ArgumentError(_pqs_source_required_factor_message(descriptor)))
        return pqs_source_pair_x2_blocks(
            plan;
            axis = descriptor.axis,
            overlap_1d,
            x2_1d,
        )
    end

    if descriptor.family === :kinetic
        isnothing(kinetic_1d) &&
            throw(ArgumentError(_pqs_source_required_factor_message(descriptor)))
        return pqs_source_pair_kinetic_blocks(
            plan;
            overlap_1d,
            kinetic_1d,
        )
    end

    throw(ArgumentError("unsupported PQS source one-body term: $(term)"))
end

"""
    pqs_source_pair_retained_one_body_block(record, term; overlap_1d, ...)

Materialize one retained PQS source-mode one-body block by delegating through
the raw source-space one-body selector and then applying the retained
source-mode boundary selector. Supported terms are currently `:overlap` and
`:kinetic`.
"""
function pqs_source_pair_retained_one_body_block(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    overlap_1d,
    kinetic_1d = nothing,
)
    descriptor = _supported_pqs_source_retained_safe_term_descriptor(term)
    if descriptor.family === :overlap
        return pqs_source_pair_retained_overlap_block(record; overlap_1d)
    end

    if descriptor.family === :kinetic
        isnothing(kinetic_1d) &&
            throw(ArgumentError(_pqs_source_required_factor_message(descriptor)))
        return pqs_source_pair_retained_kinetic_block(
            record;
            overlap_1d,
            kinetic_1d,
        )
    end

    throw(ArgumentError("unsupported retained PQS source one-body term: $(term)"))
end

"""
    pqs_source_pair_retained_one_body_blocks(plan, term; overlap_1d, ...)

Materialize retained PQS source-mode one-body blocks for ready PQS/PQS
source-pair records in a plan. Supported terms are currently `:overlap` and
`:kinetic`.
"""
function pqs_source_pair_retained_one_body_blocks(
    plan::PairBlockMaterializationPlan,
    term::Symbol;
    overlap_1d,
    kinetic_1d = nothing,
)
    descriptor = _supported_pqs_source_retained_safe_term_descriptor(term)
    return _pqs_source_pair_batch_results(
        record -> pqs_source_pair_retained_one_body_block(
            record,
            descriptor.requested_term;
            overlap_1d,
            kinetic_1d,
        ),
        plan,
        _retained_pqs_source_term(descriptor.source_term),
        _pqs_source_retained_batch_materialization_path(descriptor),
        _pqs_source_retained_unsupported_record_blocker(descriptor),
        (; retained_transform_kind = :source_mode_column_selector),
    )
end

"""
    pqs_retained_source_one_body_matrix(batch_result)
    pqs_retained_source_one_body_matrix(plan, term; overlap_1d, ...)

Build the first one-unit retained PQS source-mode dense matrix from a retained
overlap/kinetic batch. This helper only accepts one materialized retained
self-pair block. It does not perform multi-unit placement, shell realization,
Lowdin cleanup, IDA, Hamiltonian assembly, driver adoption, exports, or
artifacts.
"""
function pqs_retained_source_one_body_matrix(
    batch_result::PairBlockMaterializationBatchResult,
)
    term = batch_result.term
    _pqs_retained_source_matrix_supported_term(term) ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :unsupported_pqs_retained_source_matrix_term,
            nothing,
            0,
            nothing,
        )

    batch_result.materialized_count == 1 ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :requires_exactly_one_materialized_retained_self_pair,
            nothing,
            0,
            nothing,
        )

    length(batch_result.materialized_results) == 1 ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :materialized_result_count_mismatch,
            nothing,
            0,
            nothing,
        )

    _pqs_retained_source_skipped_records_allowed(batch_result.skipped_records) ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :skipped_ready_retained_self_pair,
            nothing,
            0,
            nothing,
        )

    result = only(batch_result.materialized_results)
    result.pair_key[1] == result.pair_key[2] ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :requires_retained_self_pair_block,
            nothing,
            0,
            result.pair_key,
        )
    _pqs_retained_source_result_is_ready(result) ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :missing_retained_source_mode_block,
            nothing,
            0,
            result.pair_key,
        )

    matrix = result.block
    size(matrix, 1) == size(matrix, 2) ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :retained_source_matrix_not_square,
            nothing,
            0,
            result.pair_key,
        )
    all(isfinite, matrix) ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :retained_source_matrix_not_finite,
            nothing,
            0,
            result.pair_key,
        )
    isapprox(matrix, transpose(matrix); rtol = 1e-12, atol = 1e-12) ||
        return _pqs_retained_source_one_body_matrix_result(
            batch_result,
            :blocked_pqs_retained_source_one_body_matrix,
            :retained_source_matrix_not_symmetric,
            nothing,
            0,
            result.pair_key,
        )

    return _pqs_retained_source_one_body_matrix_result(
        batch_result,
        :materialized_pqs_retained_source_one_body_matrix,
        nothing,
        matrix,
        size(matrix, 1),
        result.pair_key,
    )
end

function pqs_retained_source_one_body_matrix(
    plan::PairBlockMaterializationPlan,
    term::Symbol;
    overlap_1d,
    kinetic_1d = nothing,
)
    batch_result = pqs_source_pair_retained_one_body_blocks(
        plan,
        term;
        overlap_1d,
        kinetic_1d,
    )
    return pqs_retained_source_one_body_matrix(batch_result)
end

function _pqs_source_required_factor_message(descriptor)
    term_text =
        descriptor.requested_term === :kinetic ?
        ":kinetic" :
        string(descriptor.requested_term)
    return "$(descriptor.required_factor_name) is required for $(term_text)"
end

function _supported_pqs_source_retained_safe_term_descriptor(term::Symbol)
    descriptor = _supported_pqs_source_safe_term_descriptor(term)
    descriptor.family in (:overlap, :kinetic) ||
        throw(ArgumentError("unsupported retained PQS source one-body term: $(term)"))
    descriptor.family === :kinetic && isnothing(descriptor.required_factor_name) &&
        throw(ArgumentError("retained PQS kinetic term is missing factor metadata"))
    return descriptor
end

function _pqs_source_retained_batch_materialization_path(descriptor)
    return Symbol("ready_pqs_source_retained_", String(descriptor.family), "_blocks_only")
end

function _pqs_source_retained_unsupported_record_blocker(descriptor)
    return Symbol("unsupported_pqs_source_retained_", String(descriptor.family), "_materialization_record")
end

function _pqs_retained_source_matrix_supported_term(term::Symbol)
    return term in (:retained_source_overlap, :retained_source_kinetic)
end

function _pqs_retained_source_skipped_records_allowed(skipped_records)
    for skipped in skipped_records
        materialization_path =
            haskey(skipped, :materialization_path) ? skipped.materialization_path : nothing
        readiness_status =
            haskey(skipped, :readiness_status) ? skipped.readiness_status : nothing
        pair_key = haskey(skipped, :pair_key) ? skipped.pair_key : nothing
        if materialization_path === :pqs_source_pair_preflight &&
           readiness_status === :ready_metadata_only_not_materialized &&
           pair_key isa Tuple &&
           length(pair_key) == 2 &&
           pair_key[1] == pair_key[2]
            return false
        end
    end
    return true
end

function _pqs_retained_source_result_is_ready(
    result::PairBlockMaterializationResult,
)
    return result.materialized &&
           result.source_operator_blocks_materialized &&
           haskey(result.metadata, :block_space) &&
           result.metadata.block_space === :retained_pqs_source_modes
end

function _pqs_retained_source_one_body_matrix_result(
    batch_result::PairBlockMaterializationBatchResult,
    status::Symbol,
    blocker,
    matrix,
    retained_dimension::Int,
    pair_key,
)
    materialized = isnothing(blocker)
    return (;
        object_kind = :pqs_retained_source_one_body_matrix,
        status,
        blocker,
        term = batch_result.term,
        pair_key,
        matrix,
        matrix_space = :retained_pqs_source_modes,
        retained_dimension,
        materialized_result_count = batch_result.materialized_count,
        skipped_record_count = batch_result.skipped_count,
        matrix_materialized = materialized,
        retained_source_matrix_materialized = materialized,
        source_space_input_used = materialized,
        shell_realization_materialized = false,
        lowdin_cleanup_used = false,
        nonclaim_capabilities = (
            :electron_nuclear,
            :density_density,
            :ida_data,
            :hamiltonian_data,
            :driver_route,
            :exports,
            :artifacts,
        ),
    )
end
