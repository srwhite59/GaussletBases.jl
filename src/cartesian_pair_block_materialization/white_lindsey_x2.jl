# X2-only White--Lindsey pair-block pilot.

"""
    white_lindsey_boundary_stratum_x2_block(pair_unit_coefficients; axis, parent_axis_counts, overlap_1d, x2_1d)
    white_lindsey_boundary_stratum_x2_block(unit_pair; axis, parent_axis_counts, overlap_1d, x2_1d)

Materialize one White--Lindsey boundary-stratum x2 block from prepared
left/right unit coefficient maps. This uses the same support-restricted
contraction as the overlap and position pilots and does not build kinetic,
Hamiltonian data, exports, artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_x2_block(
    pair_unit_coefficients;
    axis,
    parent_axis_counts,
    overlap_1d,
    x2_1d,
)
    axis = _white_lindsey_x2_axis(axis)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    overlap_axes = _overlap_1d_tuple(overlap_1d)
    x2_axes = _operator_1d_tuple(x2_1d, "x2_1d")
    _assert_overlap_axis_sizes(overlap_axes, axis_counts)
    _assert_operator_axis_sizes(x2_axes, axis_counts, "x2_1d")

    operator_axes =
        axis === :x ? (x2_axes[1], overlap_axes[2], overlap_axes[3]) :
        axis === :y ? (overlap_axes[1], x2_axes[2], overlap_axes[3]) :
        (overlap_axes[1], overlap_axes[2], x2_axes[3])
    term = Symbol("x2_", String(axis))
    return _white_lindsey_boundary_stratum_product_block(
        pair_unit_coefficients,
        term,
        axis_counts,
        operator_axes,
        :white_lindsey_boundary_stratum_x2_adapter,
        (;
            x2_axis = axis,
            operator_factor_form = :x2_on_selected_axis_overlap_on_inactive_axes,
            supported_terms = (:x2_x, :x2_y, :x2_z),
        ),
    )
end

function white_lindsey_boundary_stratum_x2_block(
    unit_pair::CUP.UnitPairRecord;
    axis,
    parent_axis_counts,
    overlap_1d,
    x2_1d,
)
    return white_lindsey_boundary_stratum_x2_block(
        white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pair);
        axis,
        parent_axis_counts,
        overlap_1d,
        x2_1d,
    )
end

function _white_lindsey_x2_axis(axis)
    axis isa Symbol || throw(
        ArgumentError("White--Lindsey x2 axis must be :x, :y, or :z"),
    )
    axis in (:x, :y, :z) || throw(
        ArgumentError("White--Lindsey x2 axis must be :x, :y, or :z"),
    )
    return axis
end
