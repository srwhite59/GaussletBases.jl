# Support-space density inputs for complete core/shell multi-layer PQS probes.

function _pqs_multilayer_support_states(plan)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer support density helpers require a pqs_multilayer_shell_source_plan"))
    plan.status === :available_pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer support density helpers require an available source plan"))
    return vcat(plan.core_support_states, plan.shell_support_states)
end

function _pqs_multilayer_common_or_axis_tuple(axis_data, name::AbstractString)
    axis_data isa AbstractVector{<:Real} && return (axis_data, axis_data, axis_data)
    return _pqs_multilayer_axis_tuple(axis_data, name)
end

"""
    pqs_multilayer_support_weights(plan; axis_weights)

Build product support weights in the complete core/shell support ordering:
`core_support_states` followed by `shell_support_states`.
"""
function pqs_multilayer_support_weights(plan; axis_weights)
    states = _pqs_multilayer_support_states(plan)
    weights_x, weights_y, weights_z =
        _pqs_multilayer_common_or_axis_tuple(axis_weights, "axis_weights")
    weights = Vector{Float64}(undef, length(states))
    @inbounds for (index, (ix, iy, iz)) in pairs(states)
        weights[index] =
            Float64(weights_x[ix]) * Float64(weights_y[iy]) * Float64(weights_z[iz])
    end
    return weights
end

function _pqs_multilayer_raw_pair_factor_terms(raw_pair_factor_terms)
    if raw_pair_factor_terms isa AbstractArray{<:Real,3}
        return (
            raw_pair_factor_terms,
            raw_pair_factor_terms,
            raw_pair_factor_terms,
        )
    end
    return _pqs_multilayer_axis_tuple(
        raw_pair_factor_terms,
        "raw_pair_factor_terms",
    )
end

function _pqs_multilayer_validate_raw_pair_factor_terms(axis_terms, term_count::Int)
    for axis in 1:3
        ndims(axis_terms[axis]) == 3 ||
            throw(ArgumentError("PQS multi-layer raw pair factor terms must be term-first 3D arrays"))
        size(axis_terms[axis], 1) == term_count ||
            throw(ArgumentError("PQS multi-layer raw pair factor term count mismatch"))
    end
    return axis_terms
end

"""
    pqs_multilayer_support_pair_raw_numerator_matrix(plan; raw_pair_factor_terms, coulomb_expansion)

Build the complete core/shell support raw pair numerator matrix from raw
PGDG pair-factor terms. This helper intentionally does not divide by weights
and does not consume density-normalized `pair_factor_terms` as authority.
"""
function pqs_multilayer_support_pair_raw_numerator_matrix(
    plan;
    raw_pair_factor_terms,
    coulomb_expansion,
)
    states = _pqs_multilayer_support_states(plan)
    coefficients = Float64.(coulomb_expansion.coefficients)
    all(>(0.0), coefficients) ||
        throw(ArgumentError("PQS multi-layer raw pair numerator requires positive Coulomb expansion coefficients"))
    axis_terms =
        _pqs_multilayer_validate_raw_pair_factor_terms(
            _pqs_multilayer_raw_pair_factor_terms(raw_pair_factor_terms),
            length(coefficients),
        )
    result = zeros(Float64, length(states), length(states))
    for term_index in eachindex(coefficients)
        result .+=
            coefficients[term_index] *
            _pqs_multilayer_support_product_matrix(
                states,
                states,
                @view(axis_terms[1][term_index, :, :]),
                @view(axis_terms[2][term_index, :, :]),
                @view(axis_terms[3][term_index, :, :]),
            )
    end
    return result
end
