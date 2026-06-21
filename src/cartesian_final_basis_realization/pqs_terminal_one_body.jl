function _check_terminal_axis_factor(axis, n, name)
    size(axis, 1) >= n && size(axis, 2) >= n ||
        throw(DimensionMismatch("terminal one-body $name factor is too small"))
    all(isfinite, axis) || throw(ArgumentError("terminal one-body $name factor is not finite"))
    norm(axis - transpose(axis), Inf) <= 1.0e-10 ||
        throw(ArgumentError("terminal one-body $name factor is not symmetric"))
    return axis
end

function _terminal_product_action(left, right, axes)
    ncols = isnothing(right.coefficients) ?
        length(right.support_states) :
        size(right.coefficients, 2)
    action = zeros(Float64, length(left.support_states), ncols)
    maxcols = max(1, _TERMINAL_WORKSPACE_BYTES ÷ (8 * max(length(left.support_states), 1)))
    for firstcol in 1:maxcols:length(right.support_states)
        cols = firstcol:min(firstcol + maxcols - 1, length(right.support_states))
        tile = _support_cross(left.support_states, @view(right.support_states[cols]), axes)
        if isnothing(right.coefficients)
            action[:, cols] .+= tile
        else
            action .+= tile * @view(right.coefficients[cols, :])
        end
    end
    return isnothing(left.coefficients) ? action : transpose(left.coefficients) * action
end

function assemble_terminal_product_operator!(
    destination,
    basis::CartesianTerminalBasisRealization,
    axis_x,
    axis_y,
    axis_z;
    scale = 1.0,
)
    size(destination) == (basis.final_dimension, basis.final_dimension) ||
        throw(DimensionMismatch("terminal product destination size is wrong"))
    max_state = (0, 0, 0)
    for block in basis.blocks, state in block.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    axes = (
        _check_terminal_axis_factor(axis_x, max_state[1], "x"),
        _check_terminal_axis_factor(axis_y, max_state[2], "y"),
        _check_terminal_axis_factor(axis_z, max_state[3], "z"),
    )
    factor = Float64(scale)
    for right in eachindex(basis.blocks), left in 1:right
        left_block = basis.blocks[left]
        right_block = basis.blocks[right]
        rows = left_block.column_range
        cols = right_block.column_range
        block = factor .* _terminal_product_action(left_block, right_block, axes)
        destination[rows, cols] .+= block
        left != right && (destination[cols, rows] .+= transpose(block))
    end
    return destination
end
