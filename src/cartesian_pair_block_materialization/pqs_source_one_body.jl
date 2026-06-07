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
    term === :overlap && return pqs_source_pair_overlap_block(
        record;
        overlap_1d,
    )

    position_axis = _position_axis_for_term(term)
    if !isnothing(position_axis)
        isnothing(position_1d) &&
            throw(ArgumentError("position_1d is required for $(term)"))
        return pqs_source_pair_position_block(
            record;
            axis = position_axis,
            overlap_1d,
            position_1d,
        )
    end

    x2_axis = _x2_axis_for_term(term)
    if !isnothing(x2_axis)
        isnothing(x2_1d) &&
            throw(ArgumentError("x2_1d is required for $(term)"))
        return pqs_source_pair_x2_block(
            record;
            axis = x2_axis,
            overlap_1d,
            x2_1d,
        )
    end

    if term === :kinetic
        isnothing(kinetic_1d) &&
            throw(ArgumentError("kinetic_1d is required for :kinetic"))
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
    term === :overlap && return pqs_source_pair_overlap_blocks(
        plan;
        overlap_1d,
    )

    position_axis = _position_axis_for_term(term)
    if !isnothing(position_axis)
        isnothing(position_1d) &&
            throw(ArgumentError("position_1d is required for $(term)"))
        return pqs_source_pair_position_blocks(
            plan;
            axis = position_axis,
            overlap_1d,
            position_1d,
        )
    end

    x2_axis = _x2_axis_for_term(term)
    if !isnothing(x2_axis)
        isnothing(x2_1d) &&
            throw(ArgumentError("x2_1d is required for $(term)"))
        return pqs_source_pair_x2_blocks(
            plan;
            axis = x2_axis,
            overlap_1d,
            x2_1d,
        )
    end

    if term === :kinetic
        isnothing(kinetic_1d) &&
            throw(ArgumentError("kinetic_1d is required for :kinetic"))
        return pqs_source_pair_kinetic_blocks(
            plan;
            overlap_1d,
            kinetic_1d,
        )
    end

    throw(ArgumentError("unsupported PQS source one-body term: $(term)"))
end
