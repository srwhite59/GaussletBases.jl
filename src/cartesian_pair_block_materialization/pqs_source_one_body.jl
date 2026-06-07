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

function _pqs_source_required_factor_message(descriptor)
    term_text =
        descriptor.requested_term === :kinetic ?
        ":kinetic" :
        string(descriptor.requested_term)
    return "$(descriptor.required_factor_name) is required for $(term_text)"
end
