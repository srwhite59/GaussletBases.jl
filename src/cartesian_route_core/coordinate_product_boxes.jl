# Coordinate Product Box primitives.
#
# A CoordinateProductBox (CPB) is an axis-aligned product of three integer
# coordinate intervals. Singleton intervals are allowed, so codimension-1,
# codimension-2, and codimension-3 boundary strata can be represented by the
# same object as filled boxes. Shells such as B_outer \ B_inner must be
# represented as owned support, not as a CPB.

const _CPB_AXES = (:x, :y, :z)
const _CPB_SIDES = (:low, :high)

"""
    CoordinateProductBox

Axis-aligned product of three integer coordinate intervals. This is the typed
source/support object used by lowering; shells such as `B_outer \\ B_inner` are
represented by `OwnedSupport`, not by this type.
"""
struct CoordinateProductBox
    intervals::NTuple{3,UnitRange{Int}}
    role::Symbol
    codimension::Int
    metadata::NamedTuple
end

function _as_unit_range(interval)
    interval isa AbstractUnitRange{<:Integer} ||
        throw(ArgumentError("CPB intervals must be integer unit ranges"))
    isempty(interval) && throw(ArgumentError("CPB intervals must be nonempty"))
    return Int(first(interval)):Int(last(interval))
end

function _normalize_intervals(intervals)
    length(intervals) == 3 ||
        throw(ArgumentError("a CoordinateProductBox requires exactly three intervals"))
    return (_as_unit_range(intervals[1]), _as_unit_range(intervals[2]), _as_unit_range(intervals[3]))
end

function _cpb_codimension(intervals::NTuple{3,UnitRange{Int}})
    return count(interval -> length(interval) == 1, intervals)
end

"""
    cpb(ix, iy, iz; role = :coordinate_product_box, metadata = (;))
    cpb((ix, iy, iz); role = :coordinate_product_box, metadata = (;))

Construct a Coordinate Product Box (CPB), an axis-aligned product of three
integer coordinate intervals. Singleton intervals are allowed, so the same type
can represent filled boxes, slabs, facets, edges, and corners.

This constructor does not represent shells such as `B_outer \\ B_inner`; use
`complete_shell_support` for shell-owned support.
"""
function cpb(intervals; role::Symbol = :coordinate_product_box, metadata = (;))
    normalized = _normalize_intervals(intervals)
    return CoordinateProductBox(
        normalized,
        role,
        _cpb_codimension(normalized),
        NamedTuple(metadata),
    )
end

function cpb(ix, iy, iz; role::Symbol = :coordinate_product_box, metadata = (;))
    return cpb((ix, iy, iz); role, metadata)
end

"""
    filled_cpb(ix, iy, iz; role = :filled_cpb, metadata = (;))

Construct a codimension-0 CPB. All three intervals must be non-singleton.

Use this for PQS filled source boxes and ordinary volume/source boxes.
"""
function filled_cpb(ix, iy, iz; role::Symbol = :filled_cpb, metadata = (;))
    box = cpb(ix, iy, iz; role, metadata)
    box.codimension == 0 ||
        throw(ArgumentError("filled_cpb requires three nonsingleton intervals"))
    return box
end

"""
    slab_cpb(ix, iy, iz; role = :slab_cpb, metadata = (;))

Construct a codimension-1 CPB. Exactly one interval must be singleton.

Use this for direct slabs, midpoint/contact slabs, boundary slabs, or other
slab-like source/support pieces.
"""
function slab_cpb(ix, iy, iz; role::Symbol = :slab_cpb, metadata = (;))
    box = cpb(ix, iy, iz; role, metadata)
    box.codimension == 1 ||
        throw(ArgumentError("slab_cpb requires exactly one singleton interval"))
    return box
end

"""
    intervals(box)

Return the three coordinate intervals of a CPB.
"""
intervals(box::CoordinateProductBox) = box.intervals

"""
    shape(box)

Return the three interval lengths of a CPB.
"""
shape(box::CoordinateProductBox) = Tuple(length(interval) for interval in box.intervals)

"""
    codimension(box)

Return the number of singleton axes in a CPB.
"""
codimension(box::CoordinateProductBox) = box.codimension

"""
    role(object)

Return the symbolic role attached to a CPB, shellification region, or final
retained unit.
"""
role(box::CoordinateProductBox) = box.role

"""
    support_count(object)

Return the number of parent rows/sites represented by the object when that
count is known from geometry.
"""
support_count(box::CoordinateProductBox) = prod(shape(box))

_axis_index(axis::Symbol) =
    axis === :x ? 1 :
    axis === :y ? 2 :
    axis === :z ? 3 :
    throw(ArgumentError("axis must be :x, :y, or :z"))

function _axis_symbol(axis_index::Int)
    1 <= axis_index <= 3 || throw(ArgumentError("axis index must be 1, 2, or 3"))
    return _CPB_AXES[axis_index]
end

function _side_point(interval::UnitRange{Int}, side::Symbol)
    side === :low && return first(interval):first(interval)
    side === :high && return last(interval):last(interval)
    throw(ArgumentError("side must be :low or :high"))
end

function _side_label(side::Symbol)
    side in _CPB_SIDES || throw(ArgumentError("side must be :low or :high"))
    return side === :low ? "low" : "high"
end

function _role_symbol(parts::AbstractVector{String}, suffix::String)
    return Symbol(join(parts, "_") * "_" * suffix)
end

function _assert_complete_shell_boxes(
    outer_box::CoordinateProductBox,
    inner_box::CoordinateProductBox,
)
    outer_box.codimension == 0 ||
        throw(ArgumentError("complete-shell outer box must be a filled CPB"))
    inner_box.codimension == 0 ||
        throw(ArgumentError("complete-shell inner box must be a filled CPB"))
    for axis_index in 1:3
        outer = outer_box.intervals[axis_index]
        inner = inner_box.intervals[axis_index]
        length(outer) >= 3 ||
            throw(ArgumentError("complete-shell outer intervals must have length at least 3"))
        first(inner) == first(outer) + 1 ||
            throw(ArgumentError("inner interval must remove exactly the low boundary point"))
        last(inner) == last(outer) - 1 ||
            throw(ArgumentError("inner interval must remove exactly the high boundary point"))
    end
    return nothing
end

function _facet_cpb(
    outer_box::CoordinateProductBox,
    inner_box::CoordinateProductBox,
    axis_index::Int,
    side::Symbol,
)
    intervals = collect(inner_box.intervals)
    intervals[axis_index] = _side_point(outer_box.intervals[axis_index], side)
    axis = _axis_symbol(axis_index)
    role = _role_symbol([String(axis), _side_label(side)], "facet")
    return cpb(Tuple(intervals); role, metadata = (; stratum_kind = :facet_cpb, axis, side))
end

function _edge_cpb(
    outer_box::CoordinateProductBox,
    inner_box::CoordinateProductBox,
    axis_indices::Tuple{Int,Int},
    sides::Tuple{Symbol,Symbol},
)
    intervals = collect(inner_box.intervals)
    role_parts = String[]
    fixed_axes = Symbol[]
    for (axis_index, side) in zip(axis_indices, sides)
        intervals[axis_index] = _side_point(outer_box.intervals[axis_index], side)
        axis = _axis_symbol(axis_index)
        push!(fixed_axes, axis)
        push!(role_parts, String(axis), _side_label(side))
    end
    role = _role_symbol(role_parts, "edge")
    free_axis = only(setdiff(collect(1:3), collect(axis_indices)))
    return cpb(
        Tuple(intervals);
        role,
        metadata = (;
            stratum_kind = :edge_cpb,
            fixed_axes = Tuple(fixed_axes),
            free_axis = _axis_symbol(free_axis),
            sides,
        ),
    )
end

function _corner_cpb(
    outer_box::CoordinateProductBox,
    sides::NTuple{3,Symbol},
)
    intervals = Vector{UnitRange{Int}}(undef, 3)
    role_parts = String[]
    for axis_index in 1:3
        side = sides[axis_index]
        axis = _axis_symbol(axis_index)
        intervals[axis_index] = _side_point(outer_box.intervals[axis_index], side)
        push!(role_parts, String(axis), _side_label(side))
    end
    role = _role_symbol(role_parts, "corner")
    return cpb(
        Tuple(intervals);
        role,
        metadata = (; stratum_kind = :corner_cpb, sides),
    )
end

"""
    complete_shell_boundary_strata(outer_box, inner_box)

Return the disjoint CPB boundary-stratum decomposition of a complete one-layer
shell `outer_box \\ inner_box`.

The result is a named tuple with `facets`, `edges`, `corners`, and `all_strata`.
The complete-shell invariant is checked: the inner box must remove exactly one
low and one high boundary point on every axis.

This is the White--Lindsey shell-decomposition geometry. PQS shell lowering uses
a filled source CPB instead of this facet/edge/corner breakdown.
"""
function complete_shell_boundary_strata(
    outer_box::CoordinateProductBox,
    inner_box::CoordinateProductBox,
)
    _assert_complete_shell_boxes(outer_box, inner_box)

    facets = CoordinateProductBox[
        _facet_cpb(outer_box, inner_box, axis_index, side)
        for axis_index in 1:3 for side in _CPB_SIDES
    ]

    edge_axis_pairs = ((1, 2), (1, 3), (2, 3))
    edges = CoordinateProductBox[]
    for axis_pair in edge_axis_pairs, side_a in _CPB_SIDES, side_b in _CPB_SIDES
        push!(edges, _edge_cpb(outer_box, inner_box, axis_pair, (side_a, side_b)))
    end

    corners = CoordinateProductBox[
        _corner_cpb(outer_box, (side_x, side_y, side_z))
        for side_x in _CPB_SIDES for side_y in _CPB_SIDES for side_z in _CPB_SIDES
    ]

    all_strata = Tuple(vcat(facets, edges, corners))
    shell_count = support_count(outer_box) - support_count(inner_box)
    stratum_count = sum(support_count, all_strata; init = 0)
    stratum_count == shell_count ||
        throw(ArgumentError("complete-shell CPB strata do not cover the shell support count"))

    return (;
        object_kind = :complete_shell_boundary_strata,
        outer_box,
        inner_box,
        facets = Tuple(facets),
        edges = Tuple(edges),
        corners = Tuple(corners),
        all_strata,
        shell_support_count = shell_count,
        stratum_support_count = stratum_count,
    )
end
