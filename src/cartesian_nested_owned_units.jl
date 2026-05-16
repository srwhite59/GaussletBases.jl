"""
    _CartesianNestedOwnedUnit3D

Internal experimental scaffolding for future nested owned-unit/endcap-panel
construction. An owned unit records one declared support region and the local
contraction map on that support. It is not a public OPCU API and is not wired
into nested source construction.
"""
struct _CartesianNestedOwnedUnit3D{M}
    role::Symbol
    support_indices::Vector{Int}
    coefficient_matrix::_CartesianCoefficientMap
    metadata::M
end

function _CartesianNestedOwnedUnit3D(
    role::Symbol,
    support_indices::AbstractVector{<:Integer},
    coefficient_matrix::AbstractMatrix{<:Real};
    metadata = (;),
)
    support = Int.(support_indices)
    coefficients = _cartesian_coefficient_map_storage(coefficient_matrix)
    size(coefficients, 1) == length(support) || throw(
        DimensionMismatch("nested owned-unit coefficient rows must match support-index count"),
    )
    return _CartesianNestedOwnedUnit3D{typeof(metadata)}(
        role,
        support,
        coefficients,
        metadata,
    )
end

"""
    _CartesianNestedOwnedUnitCoverageAudit

Compact report for owned-unit support coverage. `owned_support_count` counts
unique owned support indices. `duplicate_count` counts extra support-owner
occurrences beyond the first owner of each support index.
"""
struct _CartesianNestedOwnedUnitCoverageAudit
    expected_support_count::Int
    owned_support_count::Int
    duplicate_count::Int
    missing_count::Int
    outside_count::Int
    retained_count::Int
    coverage_ok::Bool
end

function _nested_owned_unit_coefficient_values(unit::_CartesianNestedOwnedUnit3D)
    matrix = unit.coefficient_matrix
    return matrix isa SparseArrays.SparseMatrixCSC ? SparseArrays.nonzeros(matrix) : vec(matrix)
end

function _nested_validate_owned_unit_coefficients(unit::_CartesianNestedOwnedUnit3D)
    retained_count = size(unit.coefficient_matrix, 2)
    retained_count >= 1 || throw(
        ArgumentError("nested owned-unit $(unit.role) must retain at least one column"),
    )
    all(isfinite, _nested_owned_unit_coefficient_values(unit)) || throw(
        ArgumentError("nested owned-unit $(unit.role) coefficient map must contain only finite values"),
    )
    return retained_count
end

"""
    _nested_owned_unit_coverage_audit(units, expected_support_indices)

Audit exact support ownership for future endcap/panel shell units. The helper
checks only declared support coverage and local contraction-map sanity; it does
not assemble packets or mutate source builders.
"""
function _nested_owned_unit_coverage_audit(
    units::AbstractVector{<:_CartesianNestedOwnedUnit3D},
    expected_support_indices::AbstractVector{<:Integer},
)
    expected_counts = Dict{Int,Int}()
    for index in expected_support_indices
        expected_index = Int(index)
        expected_counts[expected_index] = get(expected_counts, expected_index, 0) + 1
    end
    all(count == 1 for count in values(expected_counts)) || throw(
        ArgumentError("nested owned-unit coverage audit requires unique expected support indices"),
    )
    expected = Set(keys(expected_counts))

    owned_counts = Dict{Int,Int}()
    retained_count = 0
    for unit in units
        retained_count += _nested_validate_owned_unit_coefficients(unit)
        for index in unit.support_indices
            owned_counts[index] = get(owned_counts, index, 0) + 1
        end
    end
    owned = Set(keys(owned_counts))
    duplicate_count = sum(max(count - 1, 0) for count in values(owned_counts))
    missing_count = length(setdiff(expected, owned))
    outside_count = length(setdiff(owned, expected))
    coverage_ok = duplicate_count == 0 && missing_count == 0 && outside_count == 0
    return _CartesianNestedOwnedUnitCoverageAudit(
        length(expected),
        length(owned),
        duplicate_count,
        missing_count,
        outside_count,
        retained_count,
        coverage_ok,
    )
end

struct _CartesianNestedEndcapPanelOwnedUnits3D{U<:Tuple}
    units::U
    expected_support_indices::Vector{Int}
    audit::_CartesianNestedOwnedUnitCoverageAudit
    current_box::NTuple{3,UnitRange{Int}}
    inner_box::NTuple{3,UnitRange{Int}}
    bond_axis::Symbol
    support_contract::Symbol
    q::Int
    L::Int
end

function _nested_axis_symbol(axis_index::Int)
    axis_index == 1 && return :x
    axis_index == 2 && return :y
    axis_index == 3 && return :z
    throw(ArgumentError("nested axis symbol lookup requires axis index 1, 2, or 3"))
end

function _nested_owned_unit_side_indices(
    current::UnitRange{Int},
    inner::UnitRange{Int},
    side::Symbol,
)
    if side == :low
        return collect(first(current):(first(inner) - 1))
    elseif side == :high
        return collect((last(inner) + 1):last(current))
    end
    throw(ArgumentError("nested owned-unit side lookup requires side = :low or :high"))
end

function _nested_validate_endcap_panel_boxes(
    dims::NTuple{3,Int},
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
)
    for axis in 1:3
        current = current_box[axis]
        inner = inner_box[axis]
        first(current) >= 1 && last(current) <= dims[axis] || throw(
            ArgumentError("nested endcap/panel current box must lie inside the Cartesian parent dimensions"),
        )
        first(current) < first(inner) <= last(inner) < last(current) || throw(
            ArgumentError("nested endcap/panel inner box must be a strict interior box on every axis"),
        )
        first(inner) == first(current) + 1 && last(inner) == last(current) - 1 || throw(
            ArgumentError("nested endcap/panel producer currently supports only one-cell-thick endcap/perimeter shells"),
        )
    end
    return nothing
end

function _nested_endcap_panel_support_indices(
    dims::NTuple{3,Int},
    x_indices::AbstractVector{<:Integer},
    y_indices::AbstractVector{<:Integer},
    z_indices::AbstractVector{<:Integer},
)
    support = Int[]
    sizehint!(support, length(x_indices) * length(y_indices) * length(z_indices))
    for ix in x_indices, iy in y_indices, iz in z_indices
        push!(support, _cartesian_flat_index(Int(ix), Int(iy), Int(iz), dims))
    end
    sort!(support)
    return support
end

function _nested_endcap_panel_support_indices(
    dims::NTuple{3,Int},
    axis_indices::NTuple{3,AbstractVector{<:Integer}},
)
    return _nested_endcap_panel_support_indices(dims, axis_indices...)
end

function _nested_owned_unit_selector_coefficients(
    nsupport::Int,
    retained_count::Int,
    role::Symbol,
)
    retained_count >= 1 || throw(
        ArgumentError("nested endcap/panel unit $role requires positive retained count"),
    )
    nsupport >= retained_count || throw(
        ArgumentError("nested endcap/panel unit $role requires support count >= retained count"),
    )
    row_indices = Vector{Int}(undef, retained_count)
    @inbounds for column in 1:retained_count
        row_indices[column] = fld((column - 1) * nsupport, retained_count) + 1
    end
    return _nested_sparse_coefficient_map(
        row_indices,
        collect(1:retained_count),
        ones(Float64, retained_count),
        nsupport,
        retained_count,
    )
end

function _nested_endcap_panel_unit(
    role::Symbol,
    support_indices::AbstractVector{Int},
    retained_count::Int;
    q::Int,
    L::Int,
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    enforce_symmetric_odd::Bool,
)
    coefficients = _nested_owned_unit_selector_coefficients(
        length(support_indices),
        retained_count,
        role,
    )
    return _CartesianNestedOwnedUnit3D(
        role,
        support_indices,
        coefficients;
        metadata = (
            q = q,
            L = L,
            current_box = current_box,
            inner_box = inner_box,
            bond_axis = bond_axis,
            support_contract = :thin_endcap_box_perimeter,
            enforce_symmetric_odd = enforce_symmetric_odd,
            support_count = length(support_indices),
            retained_count = retained_count,
        ),
    )
end

"""
    _nested_endcap_panel_owned_units(dims, current_box, inner_box; bond_axis=:z, q, L)

Internal experimental support/count producer for the validated thin
endcap-box/perimeter shared-shell contract. The support is `current_box \\
inner_box`, but only for a one-cell-thick shell: two full transverse endcaps on
the bond-axis ends plus four side panels covering the transverse perimeter over
the inner bond-axis span. Panel corners are assigned asymmetrically so each
support index has exactly one owner. The coefficient maps are sparse
direct-selector scaffolds sized to the intended retained counts; they are not
yet the physical high-order endcap/panel contraction maps.
"""
function _nested_endcap_panel_owned_units(
    dims::NTuple{3,Int},
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    q::Int,
    L::Int,
    enforce_symmetric_odd::Bool = false,
)
    q >= 1 || throw(ArgumentError("nested endcap/panel producer requires q >= 1"))
    L >= 1 || throw(ArgumentError("nested endcap/panel producer requires L >= 1"))
    _nested_validate_endcap_panel_boxes(dims, current_box, inner_box)
    bond_axis_index = _nested_axis_index(bond_axis)
    transverse_axes = Tuple(axis for axis in 1:3 if axis != bond_axis_index)
    first_transverse, second_transverse = transverse_axes

    inner_indices = ntuple(axis -> collect(inner_box[axis]), 3)
    low_sides = ntuple(
        axis -> _nested_owned_unit_side_indices(current_box[axis], inner_box[axis], :low),
        3,
    )
    high_sides = ntuple(
        axis -> _nested_owned_unit_side_indices(current_box[axis], inner_box[axis], :high),
        3,
    )

    function endcap_support(side_indices)
        axis_indices = ntuple(axis -> axis == bond_axis_index ? side_indices : collect(current_box[axis]), 3)
        return _nested_endcap_panel_support_indices(dims, axis_indices)
    end

    function perimeter_panel_support(fixed_axis::Int, fixed_value::Int, trace_axis::Int, trace_range)
        axis_indices = ntuple(axis -> begin
            if axis == bond_axis_index
                inner_indices[axis]
            elseif axis == fixed_axis
                [fixed_value]
            elseif axis == trace_axis
                collect(trace_range)
            else
                throw(ArgumentError("nested endcap/panel perimeter support received inconsistent transverse axes"))
            end
        end, 3)
        return _nested_endcap_panel_support_indices(dims, axis_indices)
    end

    first_axis = _nested_axis_symbol(first_transverse)
    second_axis = _nested_axis_symbol(second_transverse)
    first_low = first(current_box[first_transverse])
    first_high = last(current_box[first_transverse])
    second_low = first(current_box[second_transverse])
    second_high = last(current_box[second_transverse])
    units = (
        _nested_endcap_panel_unit(
            :endcap_low,
            endcap_support(low_sides[bond_axis_index]),
            q * q;
            q,
            L,
            current_box,
            inner_box,
            bond_axis,
            enforce_symmetric_odd,
        ),
        _nested_endcap_panel_unit(
            :endcap_high,
            endcap_support(high_sides[bond_axis_index]),
            q * q;
            q,
            L,
            current_box,
            inner_box,
            bond_axis,
            enforce_symmetric_odd,
        ),
        _nested_endcap_panel_unit(
            Symbol(:panel_, second_axis, :_low),
            perimeter_panel_support(
                second_transverse,
                second_low,
                first_transverse,
                first_low:(first_high - 1),
            ),
            q * L;
            q,
            L,
            current_box,
            inner_box,
            bond_axis,
            enforce_symmetric_odd,
        ),
        _nested_endcap_panel_unit(
            Symbol(:panel_, first_axis, :_high),
            perimeter_panel_support(
                first_transverse,
                first_high,
                second_transverse,
                second_low:(second_high - 1),
            ),
            q * L;
            q,
            L,
            current_box,
            inner_box,
            bond_axis,
            enforce_symmetric_odd,
        ),
        _nested_endcap_panel_unit(
            Symbol(:panel_, second_axis, :_high),
            perimeter_panel_support(
                second_transverse,
                second_high,
                first_transverse,
                (first_low + 1):first_high,
            ),
            q * L;
            q,
            L,
            current_box,
            inner_box,
            bond_axis,
            enforce_symmetric_odd,
        ),
        _nested_endcap_panel_unit(
            Symbol(:panel_, first_axis, :_low),
            perimeter_panel_support(
                first_transverse,
                first_low,
                second_transverse,
                (second_low + 1):second_high,
            ),
            q * L;
            q,
            L,
            current_box,
            inner_box,
            bond_axis,
            enforce_symmetric_odd,
        ),
    )
    expected_support_indices = setdiff(
        _nested_box_support_indices(current_box..., dims),
        _nested_box_support_indices(inner_box..., dims),
    )
    sort!(expected_support_indices)
    audit = _nested_owned_unit_coverage_audit(collect(units), expected_support_indices)
    return _CartesianNestedEndcapPanelOwnedUnits3D(
        units,
        expected_support_indices,
        audit,
        current_box,
        inner_box,
        bond_axis,
        :thin_endcap_box_perimeter,
        q,
        L,
    )
end

function _nested_endcap_panel_owned_units(
    bundles::_CartesianNestedAxisBundles3D,
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}};
    kwargs...,
)
    return _nested_endcap_panel_owned_units(
        _nested_axis_lengths(bundles),
        current_box,
        inner_box;
        kwargs...,
    )
end
