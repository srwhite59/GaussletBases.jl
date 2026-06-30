_terminal_face_axis_index(axis::Symbol) =
    axis === :x ? 1 : axis === :y ? 2 : axis === :z ? 3 :
    throw(ArgumentError("terminal face normal axis must be :x, :y, or :z"))

_terminal_face_active_axes(axis::Symbol) =
    axis === :x ? (:y, :z) : axis === :y ? (:x, :z) : axis === :z ? (:x, :y) :
    throw(ArgumentError("terminal face normal axis must be :x, :y, or :z"))

_terminal_face_kind(axis::Symbol) =
    axis === :x ? :yz : axis === :y ? :xz : axis === :z ? :xy :
    throw(ArgumentError("terminal face normal axis must be :x, :y, or :z"))

_terminal_face_fixed_side(side::Symbol) = side === :midpoint ? :low :
    side in (:low, :high) ? side : throw(ArgumentError("terminal face fixed side must be :low, :high, or :midpoint"))
function _terminal_face_product_block(
    source_cpb, bundles; normal_axis::Symbol, fixed_indices,
    retained_count::Integer, fixed_side::Symbol = :low,
)
    retained_count > 0 ||
        throw(ArgumentError("terminal face retained count must be positive"))
    fixed_side = _terminal_face_fixed_side(fixed_side)
    intervals = CartesianCPB.intervals(source_cpb)
    dims = _nested_axis_lengths(bundles)
    fixed_axis = _terminal_face_axis_index(normal_axis)
    active_axes = _terminal_face_active_axes(normal_axis)
    side(axis) = _nested_doside_1d(
        _nested_axis_pgdg(bundles, axis), intervals[_terminal_face_axis_index(axis)],
        retained_count; enforce_symmetric_odd = false,
    )
    first_side = side(active_axes[1])
    second_side = side(active_axes[2])
    indices = _nested_box_support_indices(intervals..., dims)
    fixed = Int.(collect(fixed_indices))
    all(index -> index in intervals[fixed_axis], fixed) ||
        throw(ArgumentError("terminal face fixed indices must be inside source interval"))
    coefficients = zeros(length(indices), length(fixed) * retained_count^2)
    col = 1
    for fixed_index in fixed
        product = _nested_face_product(_terminal_face_kind(normal_axis), fixed_side,
            first_side, second_side, fixed_index, dims)
        ncol = size(product.coefficient_matrix, 2)
        coefficients[:, col:(col + ncol - 1)] .=
            Matrix{Float64}(product.coefficient_matrix[indices, :])
        col += ncol
    end
    states = NTuple{3,Int}[_cartesian_unflat_index(index, dims) for index in indices]
    return indices, states, coefficients
end

function _terminal_compact_thin_slab_block(source_cpbs, metadata, bundles)
    source_cpb = only(source_cpbs)
    axis = metadata.slab_normal_axis
    q = metadata.thin_slab_retained_count_1d
    intervals = CartesianCPB.intervals(source_cpb)
    fixed_axis = _terminal_face_axis_index(axis)
    return _terminal_face_product_block(source_cpb, bundles; normal_axis = axis,
        fixed_indices = intervals[fixed_axis], retained_count = q,
        fixed_side = get(metadata, :slab_side, :low))
end
