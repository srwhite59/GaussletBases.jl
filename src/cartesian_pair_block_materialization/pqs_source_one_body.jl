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
    source_result = pqs_source_pair_one_body_block(
        record,
        descriptor.requested_term;
        overlap_1d,
        kinetic_1d,
    )
    return pqs_source_pair_retained_one_body_block(source_result)
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
