# Direct/direct one-body term selector.

"""
    direct_direct_one_body_block(record, term; parent_axis_counts, overlap_1d, ...)

Materialize one supported direct/direct one-body pair block by delegating to the
term-specific helper.
"""
function direct_direct_one_body_block(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    parent_axis_counts,
    overlap_1d,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    term === :overlap && return direct_direct_overlap_block(
        record;
        parent_axis_counts,
        overlap_1d,
    )

    position_axis = _position_axis_for_term(term)
    if !isnothing(position_axis)
        isnothing(position_1d) &&
            throw(ArgumentError("position_1d is required for $(term)"))
        return direct_direct_position_block(
            record;
            axis = position_axis,
            parent_axis_counts,
            overlap_1d,
            position_1d,
        )
    end

    x2_axis = _x2_axis_for_term(term)
    if !isnothing(x2_axis)
        isnothing(x2_1d) &&
            throw(ArgumentError("x2_1d is required for $(term)"))
        return direct_direct_x2_block(
            record;
            axis = x2_axis,
            parent_axis_counts,
            overlap_1d,
            x2_1d,
        )
    end

    if term === :kinetic
        isnothing(kinetic_1d) &&
            throw(ArgumentError("kinetic_1d is required for :kinetic"))
        return direct_direct_kinetic_block(
            record;
            parent_axis_counts,
            overlap_1d,
            kinetic_1d,
        )
    end

    throw(ArgumentError("unsupported direct/direct one-body term: $(term)"))
end

"""
    direct_direct_one_body_blocks(plan, term; parent_axis_counts, overlap_1d, ...)

Materialize supported direct/direct one-body blocks for a plan by delegating to
the term-specific plan helper.
"""
function direct_direct_one_body_blocks(
    plan::PairBlockMaterializationPlan,
    term::Symbol;
    parent_axis_counts,
    overlap_1d,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    term === :overlap && return direct_direct_overlap_blocks(
        plan;
        parent_axis_counts,
        overlap_1d,
    )

    position_axis = _position_axis_for_term(term)
    if !isnothing(position_axis)
        isnothing(position_1d) &&
            throw(ArgumentError("position_1d is required for $(term)"))
        return direct_direct_position_blocks(
            plan;
            axis = position_axis,
            parent_axis_counts,
            overlap_1d,
            position_1d,
        )
    end

    x2_axis = _x2_axis_for_term(term)
    if !isnothing(x2_axis)
        isnothing(x2_1d) &&
            throw(ArgumentError("x2_1d is required for $(term)"))
        return direct_direct_x2_blocks(
            plan;
            axis = x2_axis,
            parent_axis_counts,
            overlap_1d,
            x2_1d,
        )
    end

    if term === :kinetic
        isnothing(kinetic_1d) &&
            throw(ArgumentError("kinetic_1d is required for :kinetic"))
        return direct_direct_kinetic_blocks(
            plan;
            parent_axis_counts,
            overlap_1d,
            kinetic_1d,
        )
    end

    throw(ArgumentError("unsupported direct/direct one-body term: $(term)"))
end

function _position_axis_for_term(term::Symbol)
    term === :position_x && return :x
    term === :position_y && return :y
    term === :position_z && return :z
    return nothing
end

function _x2_axis_for_term(term::Symbol)
    term === :x2_x && return :x
    term === :x2_y && return :y
    term === :x2_z && return :z
    return nothing
end
