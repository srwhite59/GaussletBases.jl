"""
    _CartesianNestedOwnedUnit3D

Internal experimental scaffolding for owned-unit/endcap-panel construction.
An owned unit records one declared support region and the local contraction map
on that support. It is not a public OPCU API; the bounded endcap/panel route is
wired into the bond-aligned diatomic nested source only through the guarded
`:endcap_panel_owned` shared-shell policy.
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
    coefficient_contract::Symbol
    q::Int
    L::Int
end

struct _CartesianNestedEndcapPanelShellLayer3D{O,P} <: _AbstractCartesianNestedShellLayer3D
    owned_units::O
    unit_column_ranges::Vector{UnitRange{Int}}
    coefficient_matrix::_CartesianCoefficientMap
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    packet::_CartesianNestedShellPacket3D
    provenance::P
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

function _nested_owned_unit_direct_selector_coefficients(
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
    coefficient_matrix::AbstractMatrix{<:Real};
    q::Int,
    L::Int,
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    coefficient_contract::Symbol,
    enforce_symmetric_odd::Bool,
    product_metadata = (;),
)
    retained_count = size(coefficient_matrix, 2)
    metadata = merge(
        (
            q = q,
            L = L,
            current_box = current_box,
            inner_box = inner_box,
            bond_axis = bond_axis,
            support_contract = :thin_endcap_box_perimeter,
            coefficient_contract = coefficient_contract,
            enforce_symmetric_odd = enforce_symmetric_odd,
            support_count = length(support_indices),
            retained_count = retained_count,
        ),
        product_metadata,
    )
    return _CartesianNestedOwnedUnit3D(
        role,
        support_indices,
        coefficient_matrix;
        metadata,
    )
end

function _nested_endcap_panel_direct_selector_unit(
    spec;
    q::Int,
    L::Int,
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    enforce_symmetric_odd::Bool,
)
    coefficients = _nested_owned_unit_direct_selector_coefficients(
        length(spec.support_indices),
        spec.retained_count,
        spec.role,
    )
    return _nested_endcap_panel_unit(
        spec.role,
        spec.support_indices,
        coefficients;
        q,
        L,
        current_box,
        inner_box,
        bond_axis,
        coefficient_contract = :direct_selector,
        enforce_symmetric_odd,
    )
end

"""
    _nested_endcap_panel_owned_units(bundles, current_box, inner_box; bond_axis=:z, q, L)

Internal experimental producer for the validated thin endcap-box/perimeter
shared-shell contract. The support is `current_box \\ inner_box`, but only for
a one-cell-thick shell: two full transverse endcaps on the bond-axis ends plus
four side panels covering the transverse perimeter over the inner bond-axis
span. Panel corners are assigned asymmetrically so each support index has
exactly one owner.

The bundle-based route builds product `doside` contraction maps over each owned
support. The dims-only route is retained as explicit `:direct_selector` debug
scaffolding and is not the intended integration producer.
"""
function _nested_endcap_panel_unit_specs(
    dims::NTuple{3,Int},
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    q::Int,
    L::Int,
)
    q >= 1 || throw(ArgumentError("nested endcap/panel producer requires q >= 1"))
    L >= 1 || throw(ArgumentError("nested endcap/panel producer requires L >= 1"))
    _nested_validate_endcap_panel_boxes(dims, current_box, inner_box)
    bond_axis_index = _nested_axis_index(bond_axis)
    transverse_axes = Tuple(axis for axis in 1:3 if axis != bond_axis_index)
    first_transverse, second_transverse = transverse_axes

    low_sides = ntuple(
        axis -> _nested_owned_unit_side_indices(current_box[axis], inner_box[axis], :low),
        3,
    )
    high_sides = ntuple(
        axis -> _nested_owned_unit_side_indices(current_box[axis], inner_box[axis], :high),
        3,
    )

    function support_indices(
        fixed_axis::Int,
        fixed_value::Int,
        first_axis::Int,
        first_interval::UnitRange{Int},
        second_axis::Int,
        second_interval::UnitRange{Int},
    )
        axis_indices = ntuple(axis -> begin
            if axis == fixed_axis
                [fixed_value]
            elseif axis == first_axis
                collect(first_interval)
            elseif axis == second_axis
                collect(second_interval)
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

    function spec(
        role::Symbol,
        fixed_axis::Int,
        fixed_index::Int,
        first_axis::Int,
        first_interval::UnitRange{Int},
        first_retained_count::Int,
        second_axis::Int,
        second_interval::UnitRange{Int},
        second_retained_count::Int,
    )
        return (
            role = role,
            support_indices = support_indices(
                fixed_axis,
                fixed_index,
                first_axis,
                first_interval,
                second_axis,
                second_interval,
            ),
            fixed_axis = fixed_axis,
            fixed_index = fixed_index,
            first_axis = first_axis,
            first_interval = first_interval,
            first_retained_count = first_retained_count,
            second_axis = second_axis,
            second_interval = second_interval,
            second_retained_count = second_retained_count,
            retained_count = first_retained_count * second_retained_count,
        )
    end

    return (
        spec(
            :endcap_low,
            bond_axis_index,
            only(low_sides[bond_axis_index]),
            first_transverse,
            current_box[first_transverse],
            q,
            second_transverse,
            current_box[second_transverse],
            q,
        ),
        spec(
            :endcap_high,
            bond_axis_index,
            only(high_sides[bond_axis_index]),
            first_transverse,
            current_box[first_transverse],
            q,
            second_transverse,
            current_box[second_transverse],
            q,
        ),
        spec(
            Symbol(:panel_, second_axis, :_low),
            second_transverse,
            second_low,
            first_transverse,
            first_low:(first_high - 1),
            q,
            bond_axis_index,
            inner_box[bond_axis_index],
            L,
        ),
        spec(
            Symbol(:panel_, first_axis, :_high),
            first_transverse,
            first_high,
            second_transverse,
            second_low:(second_high - 1),
            q,
            bond_axis_index,
            inner_box[bond_axis_index],
            L,
        ),
        spec(
            Symbol(:panel_, second_axis, :_high),
            second_transverse,
            second_high,
            first_transverse,
            (first_low + 1):first_high,
            q,
            bond_axis_index,
            inner_box[bond_axis_index],
            L,
        ),
        spec(
            Symbol(:panel_, first_axis, :_low),
            first_transverse,
            first_low,
            second_transverse,
            (second_low + 1):second_high,
            q,
            bond_axis_index,
            inner_box[bond_axis_index],
            L,
        ),
    )
end

function _nested_endcap_panel_product_coefficients(
    dims::NTuple{3,Int},
    support_indices::AbstractVector{Int},
    fixed_axis::Int,
    fixed_index::Int,
    first_axis::Int,
    first_side::_CartesianNestedDoSide1D,
    second_axis::Int,
    second_side::_CartesianNestedDoSide1D,
)
    length(unique((fixed_axis, first_axis, second_axis))) == 3 || throw(
        ArgumentError("nested endcap/panel product contraction requires three distinct axes"),
    )
    row_lookup = Dict{Int,Int}(index => row for (row, index) in enumerate(support_indices))
    nfirst = size(first_side.coefficient_matrix, 2)
    nsecond = size(second_side.coefficient_matrix, 2)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    column = 0
    state = Vector{Int}(undef, 3)
    state[fixed_axis] = fixed_index
    for ifirst in 1:nfirst, isecond in 1:nsecond
        column += 1
        for (local_first, index_first) in enumerate(first_side.interval)
            value_first = Float64(first_side.local_coefficients[local_first, ifirst])
            iszero(value_first) && continue
            state[first_axis] = index_first
            for (local_second, index_second) in enumerate(second_side.interval)
                value_second = Float64(second_side.local_coefficients[local_second, isecond])
                iszero(value_second) && continue
                state[second_axis] = index_second
                flat_index = _cartesian_flat_index(state[1], state[2], state[3], dims)
                row = get(row_lookup, flat_index, 0)
                row > 0 || throw(
                    ArgumentError("nested endcap/panel product contraction generated support outside its owned unit"),
                )
                push!(row_indices, row)
                push!(col_indices, column)
                push!(values, value_first * value_second)
            end
        end
    end
    return _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        length(support_indices),
        nfirst * nsecond,
    )
end

function _nested_endcap_panel_product_unit(
    spec,
    bundles::_CartesianNestedAxisBundles3D,
    dims::NTuple{3,Int};
    q::Int,
    L::Int,
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    enforce_symmetric_odd::Bool,
)
    first_side = _nested_doside_1d(
        _nested_axis_bundle(bundles, _nested_axis_symbol(spec.first_axis)),
        spec.first_interval,
        spec.first_retained_count;
        enforce_symmetric_odd,
    )
    second_side = _nested_doside_1d(
        _nested_axis_bundle(bundles, _nested_axis_symbol(spec.second_axis)),
        spec.second_interval,
        spec.second_retained_count;
        enforce_symmetric_odd,
    )
    coefficients = _nested_endcap_panel_product_coefficients(
        dims,
        spec.support_indices,
        spec.fixed_axis,
        spec.fixed_index,
        spec.first_axis,
        first_side,
        spec.second_axis,
        second_side,
    )
    return _nested_endcap_panel_unit(
        spec.role,
        spec.support_indices,
        coefficients;
        q,
        L,
        current_box,
        inner_box,
        bond_axis,
        coefficient_contract = :product_doside,
        enforce_symmetric_odd,
        product_metadata = (
            fixed_axis = spec.fixed_axis,
            fixed_index = spec.fixed_index,
            first_axis = spec.first_axis,
            first_interval = spec.first_interval,
            first_retained_count = spec.first_retained_count,
            first_coefficients = Matrix{Float64}(first_side.local_coefficients),
            second_axis = spec.second_axis,
            second_interval = spec.second_interval,
            second_retained_count = spec.second_retained_count,
            second_coefficients = Matrix{Float64}(second_side.local_coefficients),
        ),
    )
end

function _nested_endcap_panel_owned_units_result(
    units,
    dims::NTuple{3,Int},
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    coefficient_contract::Symbol,
    q::Int,
    L::Int,
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
        coefficient_contract,
        q,
        L,
    )
end

function _nested_endcap_panel_owned_units(
    dims::NTuple{3,Int},
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    q::Int,
    L::Int,
    coefficient_contract::Symbol = :direct_selector,
    enforce_symmetric_odd::Bool = false,
)
    coefficient_contract == :direct_selector || throw(
        ArgumentError("nested endcap/panel dims-only producer supports only coefficient_contract = :direct_selector; pass axis bundles for :product_doside"),
    )
    specs = _nested_endcap_panel_unit_specs(dims, current_box, inner_box; bond_axis, q, L)
    units = Tuple(
        _nested_endcap_panel_direct_selector_unit(
            spec;
            q,
            L,
            current_box,
            inner_box,
            bond_axis,
            enforce_symmetric_odd,
        ) for spec in specs
    )
    return _nested_endcap_panel_owned_units_result(
        units,
        dims,
        current_box,
        inner_box,
        bond_axis,
        coefficient_contract,
        q,
        L,
    )
end

function _nested_endcap_panel_owned_units(
    bundles::_CartesianNestedAxisBundles3D,
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    q::Int,
    L::Int,
    coefficient_contract::Symbol = :product_doside,
    enforce_symmetric_odd::Bool = false,
)
    dims = _nested_axis_lengths(bundles)
    coefficient_contract in (:product_doside, :direct_selector) || throw(
        ArgumentError("nested endcap/panel producer supports coefficient_contract = :product_doside or :direct_selector"),
    )
    coefficient_contract == :direct_selector && return _nested_endcap_panel_owned_units(
        dims,
        current_box,
        inner_box;
        bond_axis,
        q,
        L,
        coefficient_contract,
        enforce_symmetric_odd,
    )

    specs = _nested_endcap_panel_unit_specs(dims, current_box, inner_box; bond_axis, q, L)
    units = Tuple(
        _nested_endcap_panel_product_unit(
            spec,
            bundles,
            dims;
            q,
            L,
            current_box,
            inner_box,
            bond_axis,
            enforce_symmetric_odd,
        ) for spec in specs
    )
    return _nested_endcap_panel_owned_units_result(
        units,
        dims,
        current_box,
        inner_box,
        bond_axis,
        coefficient_contract,
        q,
        L,
    )
end

function _nested_endcap_panel_parent_coefficient_block(
    unit::_CartesianNestedOwnedUnit3D,
    parent_row_count::Int,
)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    matrix = unit.coefficient_matrix
    if SparseArrays.issparse(matrix)
        local_rows, cols, nzvals = SparseArrays.findnz(matrix)
        sizehint!(row_indices, length(nzvals))
        sizehint!(col_indices, length(nzvals))
        sizehint!(values, length(nzvals))
        for index in eachindex(nzvals)
            push!(row_indices, unit.support_indices[Int(local_rows[index])])
            push!(col_indices, Int(cols[index]))
            push!(values, Float64(nzvals[index]))
        end
    else
        for column in axes(matrix, 2), local_row in axes(matrix, 1)
            value = Float64(matrix[local_row, column])
            iszero(value) && continue
            push!(row_indices, unit.support_indices[local_row])
            push!(col_indices, column)
            push!(values, value)
        end
    end
    return _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        parent_row_count,
        size(matrix, 2),
    )
end

function _nested_endcap_panel_parent_coefficients(
    owned_units::_CartesianNestedEndcapPanelOwnedUnits3D,
    dims::NTuple{3,Int},
)
    parent_row_count = prod(dims)
    blocks = AbstractMatrix{Float64}[]
    unit_column_ranges = UnitRange{Int}[]
    column_start = 1
    for unit in owned_units.units
        block = _nested_endcap_panel_parent_coefficient_block(unit, parent_row_count)
        push!(blocks, block)
        column_count = size(block, 2)
        push!(unit_column_ranges, column_start:(column_start + column_count - 1))
        column_start += column_count
    end
    return _nested_hcat_coefficient_maps(blocks), unit_column_ranges
end

function _nested_endcap_panel_shell_layer(
    owned_units::_CartesianNestedEndcapPanelOwnedUnits3D,
    bundles::_CartesianNestedAxisBundles3D;
    packet_kernel::Symbol = :support_reference,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    verify_factorized_reconstruction::Bool = true,
)
    owned_units.coefficient_contract == :product_doside || throw(
        ArgumentError("nested endcap/panel shell-layer assembly requires product_doside owned-unit coefficients"),
    )
    owned_units.audit.coverage_ok || throw(
        ArgumentError("nested endcap/panel shell-layer assembly requires exact owned-unit support coverage"),
    )
    dims = _nested_axis_lengths(bundles)
    coefficient_matrix, unit_column_ranges =
        _nested_endcap_panel_parent_coefficients(owned_units, dims)
    support_indices = copy(owned_units.expected_support_indices)
    packet_data = _nested_shell_packet(
        bundles,
        coefficient_matrix,
        support_indices;
        packet_kernel,
        term_coefficients,
        verify_factorized_reconstruction,
    )
    provenance = (
        support_contract = owned_units.support_contract,
        coefficient_contract = owned_units.coefficient_contract,
        current_box = owned_units.current_box,
        inner_box = owned_units.inner_box,
        bond_axis = owned_units.bond_axis,
        q = owned_units.q,
        L = owned_units.L,
        packet_kernel = packet_kernel,
    )
    return _CartesianNestedEndcapPanelShellLayer3D(
        owned_units,
        unit_column_ranges,
        coefficient_matrix,
        support_indices,
        packet_data.support_states,
        packet_data.packet,
        provenance,
    )
end

function _nested_endcap_panel_shell_layer(
    bundles::_CartesianNestedAxisBundles3D,
    current_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    q::Int,
    L::Int,
    packet_kernel::Symbol = :support_reference,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    verify_factorized_reconstruction::Bool = true,
    enforce_symmetric_odd::Bool = false,
)
    owned_units = _nested_endcap_panel_owned_units(
        bundles,
        current_box,
        inner_box;
        bond_axis,
        q,
        L,
        coefficient_contract = :product_doside,
        enforce_symmetric_odd,
    )
    return _nested_endcap_panel_shell_layer(
        owned_units,
        bundles;
        packet_kernel,
        term_coefficients,
        verify_factorized_reconstruction,
    )
end
