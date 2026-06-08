# Selector over White--Lindsey one-body pair-block pilots.

const _WHITE_LINDSEY_ONE_BODY_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

"""
    white_lindsey_boundary_stratum_one_body_block(pair_unit_coefficients, term; ...)
    white_lindsey_boundary_stratum_one_body_block(unit_pair, term; ...)

Dispatch one supported White--Lindsey boundary-stratum one-body safe-term
pilot. Supported terms are overlap, position_x/y/z, x2_x/y/z, and kinetic.
This selector does not build plan-level operators, Hamiltonian data, exports,
artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_one_body_block(
    pair_unit_coefficients,
    term;
    parent_axis_counts,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    term = _white_lindsey_one_body_term(term)
    if term === :overlap
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        return white_lindsey_boundary_stratum_overlap_block(
            pair_unit_coefficients;
            parent_axis_counts,
            overlap_1d,
        )
    elseif term in (:position_x, :position_y, :position_z)
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :position_1d, position_1d)
        return white_lindsey_boundary_stratum_position_block(
            pair_unit_coefficients;
            axis = _white_lindsey_term_axis(term, "position_"),
            parent_axis_counts,
            overlap_1d,
            position_1d,
        )
    elseif term in (:x2_x, :x2_y, :x2_z)
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :x2_1d, x2_1d)
        return white_lindsey_boundary_stratum_x2_block(
            pair_unit_coefficients;
            axis = _white_lindsey_term_axis(term, "x2_"),
            parent_axis_counts,
            overlap_1d,
            x2_1d,
        )
    elseif term === :kinetic
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :kinetic_1d, kinetic_1d)
        return white_lindsey_boundary_stratum_kinetic_block(
            pair_unit_coefficients;
            parent_axis_counts,
            overlap_1d,
            kinetic_1d,
        )
    end
    throw(ArgumentError("unsupported White--Lindsey one-body term $(term)"))
end

function white_lindsey_boundary_stratum_one_body_block(
    unit_pair::CUP.UnitPairRecord,
    term;
    parent_axis_counts,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    return white_lindsey_boundary_stratum_one_body_block(
        white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pair),
        term;
        parent_axis_counts,
        overlap_1d,
        position_1d,
        x2_1d,
        kinetic_1d,
    )
end

function _white_lindsey_one_body_term(term)
    term isa Symbol || throw(
        ArgumentError("White--Lindsey one-body term must be a Symbol"),
    )
    term in _WHITE_LINDSEY_ONE_BODY_TERMS || throw(
        ArgumentError("unsupported White--Lindsey one-body term $(term)"),
    )
    return term
end

function _white_lindsey_require_one_body_factor(term::Symbol, name::Symbol, value)
    isnothing(value) && throw(
        ArgumentError("White--Lindsey $(term) requires $(name)"),
    )
    return nothing
end

function _white_lindsey_term_axis(term::Symbol, prefix::AbstractString)
    text = String(term)
    startswith(text, prefix) || throw(
        ArgumentError("White--Lindsey term $(term) does not have prefix $(prefix)"),
    )
    axis = Symbol(text[(lastindex(prefix) + 1):end])
    axis in (:x, :y, :z) || throw(
        ArgumentError("White--Lindsey term $(term) does not name axis x, y, or z"),
    )
    return axis
end
