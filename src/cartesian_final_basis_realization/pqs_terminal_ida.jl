function _check_terminal_ida_weights(weights, n, name)
    length(weights) >= n ||
        throw(DimensionMismatch("terminal IDA $name weights are too short"))
    all(isfinite, weights) ||
        throw(ArgumentError("terminal IDA $name weights are not finite"))
    return weights
end

function _terminal_ida_support_weights(block, weights)
    wx, wy, wz = weights
    return Float64[
        Float64(wx[state[1]]) * Float64(wy[state[2]]) * Float64(wz[state[3]])
        for state in block.support_states
    ]
end

function _terminal_ida_final_weights(block, weights, weight_atol)
    support_weights = _terminal_ida_support_weights(block, weights)
    final_weights = isnothing(block.coefficients) ?
        support_weights :
        vec(transpose(block.coefficients) * support_weights)
    all(weight -> isfinite(weight) && weight > weight_atol, final_weights) ||
        throw(ArgumentError("terminal IDA final weights must be finite and positive"))
    return final_weights
end

function _normalize_terminal_ida_block!(block, left_weights, right_weights)
    size(block) == (length(left_weights), length(right_weights)) ||
        throw(DimensionMismatch("terminal IDA block and weight dimensions differ"))
    @inbounds for col in axes(block, 2), row in axes(block, 1)
        block[row, col] /= left_weights[row] * right_weights[col]
    end
    return block
end

function assemble_terminal_ida_interaction!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    raw_pair_terms_x,
    raw_pair_terms_y,
    raw_pair_terms_z,
    weights_x,
    weights_y,
    weights_z;
    weight_atol = 1.0e-12,
    symmetry_atol = 1.0e-10,
)
    size(destination) == (basis.final_dimension, basis.final_dimension) ||
        throw(DimensionMismatch("terminal IDA destination size is wrong"))
    coeffs = Float64.(coefficients)
    !isempty(coeffs) && all(isfinite, coeffs) ||
        throw(ArgumentError("terminal IDA coefficients must be finite"))
    max_state = (0, 0, 0)
    for block in basis.blocks, state in block.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    nterms = length(coeffs)
    raw_terms = (
        _check_terminal_factor_terms(_terminal_factor_terms(raw_pair_terms_x), nterms, max_state[1], "raw pair x"),
        _check_terminal_factor_terms(_terminal_factor_terms(raw_pair_terms_y), nterms, max_state[2], "raw pair y"),
        _check_terminal_factor_terms(_terminal_factor_terms(raw_pair_terms_z), nterms, max_state[3], "raw pair z"),
    )
    weights = (
        _check_terminal_ida_weights(weights_x, max_state[1], "x"),
        _check_terminal_ida_weights(weights_y, max_state[2], "y"),
        _check_terminal_ida_weights(weights_z, max_state[3], "z"),
    )
    block_weights = [
        _terminal_ida_final_weights(block, weights, Float64(weight_atol))
        for block in basis.blocks
    ]
    for right in eachindex(basis.blocks), left in 1:right
        left_block = basis.blocks[left]
        right_block = basis.blocks[right]
        rows = left_block.column_range
        cols = right_block.column_range
        block = _terminal_gaussian_sum_action(left_block, right_block, coeffs, raw_terms)
        _normalize_terminal_ida_block!(block, block_weights[left], block_weights[right])
        if left == right
            norm(block - transpose(block), Inf) <= symmetry_atol ||
                throw(ArgumentError("terminal IDA diagonal block is not symmetric"))
        end
        destination[rows, cols] .+= block
        left != right && (destination[cols, rows] .+= transpose(block))
    end
    norm(destination - transpose(destination), Inf) <= symmetry_atol ||
        throw(ArgumentError("terminal IDA destination is not symmetric"))
    return destination
end
