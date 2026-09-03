"""
    BondAlignedDiatomicGeometryPoint3D

One emitted point in the narrow bond-aligned diatomic geometry surface.

The first public payload stays backend-free and only records stable geometry
data together with lightweight grouping metadata.
"""
struct BondAlignedDiatomicGeometryPoint3D
    x::Float64
    y::Float64
    z::Float64
    kind::Symbol
    label::String
    group_kind::Symbol
    group_id::Int
end

"""
    BondAlignedDiatomicGeometryNucleus3D

One nuclear center in the narrow bond-aligned diatomic geometry surface.
"""
struct BondAlignedDiatomicGeometryNucleus3D
    x::Float64
    y::Float64
    z::Float64
    label::String
    group_id::Int
end

"""
    BondAlignedDiatomicGeometryBox3D

Axis-aligned outline metadata for one bond-aligned diatomic box-like region.

The current payload records parent-grid index ranges together with the
corresponding physical center ranges on each axis.
"""
struct BondAlignedDiatomicGeometryBox3D
    label::String
    group_kind::Symbol
    group_id::Int
    ix::UnitRange{Int}
    iy::UnitRange{Int}
    iz::UnitRange{Int}
    x_range::NTuple{2,Float64}
    y_range::NTuple{2,Float64}
    z_range::NTuple{2,Float64}
end

"""
    BondAlignedDiatomicGeometryShellProvenance3D

Explicit provenance for one shared nested shell layer in the bond-aligned
diatomic geometry/export surface.

`source_box` is the outer box being peeled for that shell, while
`next_inner_box` is the remaining inner box after peeling it. `source_point_count`
records the shell annulus point count, not the full outer-box volume.
"""
struct BondAlignedDiatomicGeometryShellProvenance3D
    label::String
    group_kind::Symbol
    group_id::Int
    source_box::NTuple{3,UnitRange{Int}}
    next_inner_box::NTuple{3,UnitRange{Int}}
    source_point_count::Int
    retained_fixed_count::Int
end

function Base.:(==)(
    left::BondAlignedDiatomicGeometryShellProvenance3D,
    right::BondAlignedDiatomicGeometryShellProvenance3D,
)
    return (
        left.label == right.label &&
        left.group_kind == right.group_kind &&
        left.group_id == right.group_id &&
        left.source_box == right.source_box &&
        left.next_inner_box == right.next_inner_box &&
        left.source_point_count == right.source_point_count &&
        left.retained_fixed_count == right.retained_fixed_count
    )
end

"""
    BondAlignedDiatomicGeometryPayload3D

Backend-free geometry payload for the landed bond-aligned diatomic line.
"""
struct BondAlignedDiatomicGeometryPayload3D
    points::Vector{BondAlignedDiatomicGeometryPoint3D}
    nuclei::Vector{BondAlignedDiatomicGeometryNucleus3D}
    bond_axis::Symbol
    box_outlines::Vector{BondAlignedDiatomicGeometryBox3D}
    shell_provenance::Vector{BondAlignedDiatomicGeometryShellProvenance3D}
end

"""
    BondAlignedDiatomicGeometryPlaneSlice3D

Strict/near-strict representative-plane selection result for a
[`BondAlignedDiatomicGeometryPayload3D`](@ref).
"""
struct BondAlignedDiatomicGeometryPlaneSlice3D
    points::Vector{BondAlignedDiatomicGeometryPoint3D}
    nuclei::Vector{BondAlignedDiatomicGeometryNucleus3D}
    bond_axis::Symbol
    box_outlines::Vector{BondAlignedDiatomicGeometryBox3D}
    shell_provenance::Vector{BondAlignedDiatomicGeometryShellProvenance3D}
    plane_axis::Symbol
    plane_value::Float64
    plane_tol::Float64
    selected_count::Int
    total_count::Int
end

function Base.show(io::IO, payload::BondAlignedDiatomicGeometryPayload3D)
    print(
        io,
        "BondAlignedDiatomicGeometryPayload3D(npoints=",
        length(payload.points),
        ", nnuclei=",
        length(payload.nuclei),
        ", bond_axis=:",
        payload.bond_axis,
        ", nboxes=",
        length(payload.box_outlines),
        ", nshells=",
        length(payload.shell_provenance),
        ")",
    )
end

function Base.show(io::IO, slice::BondAlignedDiatomicGeometryPlaneSlice3D)
    print(
        io,
        "BondAlignedDiatomicGeometryPlaneSlice3D(axis=:",
        slice.plane_axis,
        ", value=",
        slice.plane_value,
        ", tol=",
        slice.plane_tol,
        ", selected=",
        slice.selected_count,
        "/",
        slice.total_count,
        ", nshells=",
        length(slice.shell_provenance),
        ")",
    )
end

function _bond_aligned_geometry_nuclei(
    nuclei::AbstractVector{<:NTuple{3,<:Real}},
)
    return [
        BondAlignedDiatomicGeometryNucleus3D(
            Float64(nucleus[1]),
            Float64(nucleus[2]),
            Float64(nucleus[3]),
            index == 1 ? "A" : index == 2 ? "B" : "N$index",
            index,
        ) for (index, nucleus) in pairs(nuclei)
    ]
end

function _bond_aligned_geometry_box(
    basis::BondAlignedDiatomicQWBasis3D,
    box::NTuple{3,UnitRange{Int}},
    label::AbstractString,
    group_kind::Symbol,
    group_id::Integer,
)
    x_centers = centers(basis.basis_x)
    y_centers = centers(basis.basis_y)
    z_centers = centers(basis.basis_z)
    return BondAlignedDiatomicGeometryBox3D(
        String(label),
        group_kind,
        Int(group_id),
        box[1],
        box[2],
        box[3],
        (Float64(x_centers[first(box[1])]), Float64(x_centers[last(box[1])])),
        (Float64(y_centers[first(box[2])]), Float64(y_centers[last(box[2])])),
        (Float64(z_centers[first(box[3])]), Float64(z_centers[last(box[3])])),
    )
end

function _bond_aligned_full_box(
    basis::BondAlignedDiatomicQWBasis3D,
)
    return (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
end

function _bond_aligned_default_plane_axis(
    bond_axis::Symbol,
)
    bond_axis == :z && return :y
    bond_axis == :x && return :z
    bond_axis == :y && return :z
    throw(ArgumentError("bond-aligned diatomic plane selection requires bond_axis = :x, :y, or :z"))
end

function _bond_aligned_same_basis_geometry(
    left::BondAlignedDiatomicQWBasis3D,
    right::BondAlignedDiatomicQWBasis3D,
)
    left.bond_axis == right.bond_axis || return false
    _qwrg_same_nuclei(left.nuclei, right.nuclei) || return false
    return (
        length(left.basis_x) == length(right.basis_x) &&
        length(left.basis_y) == length(right.basis_y) &&
        length(left.basis_z) == length(right.basis_z)
    )
end

function _bond_aligned_point_coordinate(
    point::BondAlignedDiatomicGeometryPoint3D,
    axis::Symbol,
)
    axis == :x && return point.x
    axis == :y && return point.y
    axis == :z && return point.z
    throw(ArgumentError("plane selection requires plane_axis = :x, :y, or :z"))
end

function _bond_aligned_nucleus_coordinate(
    point::BondAlignedDiatomicGeometryNucleus3D,
    axis::Symbol,
)
    axis == :x && return point.x
    axis == :y && return point.y
    axis == :z && return point.z
    throw(ArgumentError("plane selection requires plane_axis = :x, :y, or :z"))
end

function _bond_aligned_ordinary_points(
    basis::BondAlignedDiatomicQWBasis3D,
)
    orbitals = _mapped_cartesian_orbitals(
        centers(basis.basis_x),
        centers(basis.basis_y),
        centers(basis.basis_z),
    )
    return [
        BondAlignedDiatomicGeometryPoint3D(
            orbital.x,
            orbital.y,
            orbital.z,
            :gausslet,
            "g($(orbital.ix),$(orbital.iy),$(orbital.iz))",
            :gausslet_product,
            1,
        ) for orbital in orbitals
    ]
end

function _bond_aligned_qw_points(
    ops::OrdinaryCartesianOperators3D,
)
    points = BondAlignedDiatomicGeometryPoint3D[]
    residual_group = 0
    for orbital in ops.orbital_data
        if orbital.kind == :residual_gaussian
            residual_group += 1
            group_kind = :residual_gaussian
            group_id = residual_group
        else
            group_kind = :gausslet_product
            group_id = 1
        end
        push!(
            points,
            BondAlignedDiatomicGeometryPoint3D(
                orbital.x,
                orbital.y,
                orbital.z,
                orbital.kind,
                orbital.label,
                group_kind,
                group_id,
            ),
        )
    end
    return points
end

function _bond_aligned_nested_group_specs(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    specs = NamedTuple{(:range, :group_kind, :group_id, :label_prefix)}[]
    if source.geometry.did_split
        push!(
            specs,
            (
                range = source.child_column_ranges[1],
                group_kind = :left_child,
                group_id = 1,
                label_prefix = "left_child",
            ),
        )
        if !isnothing(source.midpoint_slab_column_range)
            push!(
                specs,
                (
                    range = source.midpoint_slab_column_range,
                    group_kind = :shared_midpoint_slab,
                    group_id = 1,
                    label_prefix = "shared_midpoint_slab",
                ),
            )
        end
        push!(
            specs,
            (
                range = source.child_column_ranges[2],
                group_kind = :right_child,
                group_id = 2,
                label_prefix = "right_child",
            ),
        )
    else
        push!(
            specs,
            (
                range = source.sequence.core_column_range,
                group_kind = :shared_child,
                group_id = 1,
                label_prefix = "shared_child",
            ),
        )
    end

    for (layer_index, range) in pairs(source.sequence.layer_column_ranges)
        push!(
            specs,
            (
                range = range,
                group_kind = :shared_shell_layer,
                group_id = layer_index,
                label_prefix = "shared_shell_$(layer_index)",
            ),
        )
    end
    return specs
end

function _bond_aligned_nested_points(
    fixed_centers::AbstractMatrix{<:Real},
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    size(fixed_centers, 2) == 3 || throw(
        ArgumentError("bond-aligned nested geometry requires an n×3 fixed center matrix"),
    )
    nfixed = size(fixed_centers, 1)
    points = BondAlignedDiatomicGeometryPoint3D[]
    covered = falses(nfixed)
    for spec in _bond_aligned_nested_group_specs(source)
        last(spec.range) <= nfixed || throw(
            ArgumentError("nested geometry grouping exceeds the fixed center count"),
        )
        for column in spec.range
            covered[column] = true
            local_index = column - first(spec.range) + 1
            push!(
                points,
                BondAlignedDiatomicGeometryPoint3D(
                    Float64(fixed_centers[column, 1]),
                    Float64(fixed_centers[column, 2]),
                    Float64(fixed_centers[column, 3]),
                    :nested_fixed,
                    string(spec.label_prefix, "_", local_index),
                    spec.group_kind,
                    spec.group_id,
                ),
            )
        end
    end
    all(covered) || throw(
        ArgumentError("nested geometry grouping did not cover all fixed centers"),
    )
    return points
end

function _bond_aligned_nested_box_outlines(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    outlines = BondAlignedDiatomicGeometryBox3D[]
    push!(
        outlines,
        _bond_aligned_geometry_box(
            source.basis,
            source.geometry.parent_box,
            "parent_box",
            :parent_box,
            1,
        ),
    )
    push!(
        outlines,
        _bond_aligned_geometry_box(
            source.basis,
            source.geometry.working_box,
            "working_box",
            :working_box,
            1,
        ),
    )
    for (child_index, child_box) in pairs(source.geometry.child_boxes)
        label = source.geometry.did_split ?
            (child_index == 1 ? "left_child_box" : "right_child_box") :
            "shared_child_box"
        push!(
            outlines,
            _bond_aligned_geometry_box(
                source.basis,
                child_box,
                label,
                :child_box,
                child_index,
            ),
        )
    end
    if !isnothing(source.geometry.shared_midpoint_box)
        push!(
            outlines,
            _bond_aligned_geometry_box(
                source.basis,
                source.geometry.shared_midpoint_box,
                "shared_midpoint_slab_box",
                :shared_midpoint_slab_box,
                1,
            ),
        )
    end
    return outlines
end

function _bond_aligned_nested_shell_provenance(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return [
        BondAlignedDiatomicGeometryShellProvenance3D(
            "shared_shell_$(layer_index)",
            :shared_shell_layer,
            Int(layer_index),
            layer.provenance.source_box,
            layer.provenance.next_inner_box,
            layer.provenance.source_point_count,
            layer.provenance.retained_fixed_count,
        ) for (layer_index, layer) in pairs(source.shared_shell_layers)
    ]
end

function _bond_aligned_support_states_for_export(
    support_indices::AbstractVector{Int},
    dims::NTuple{3,Int},
)
    return [_cartesian_unflat_index(index, dims) for index in support_indices]
end

function _bond_aligned_support_states_for_export(
    support_states::AbstractVector{<:NTuple{3,Int}},
    _support_indices::AbstractVector{Int},
    _dims::NTuple{3,Int},
)
    return support_states
end

function _bond_aligned_support_states_for_export(
    support_states::Nothing,
    support_indices::AbstractVector{Int},
    dims::NTuple{3,Int},
)
    return _bond_aligned_support_states_for_export(support_indices, dims)
end

function _bond_aligned_source_region_points(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    basis = source.basis
    dims = (length(basis.basis_x), length(basis.basis_y), length(basis.basis_z))
    x_centers = centers(basis.basis_x)
    y_centers = centers(basis.basis_y)
    z_centers = centers(basis.basis_z)
    points = BondAlignedDiatomicGeometryPoint3D[]

    for (layer_index, layer) in pairs(source.shared_shell_layers)
        layer_states = _bond_aligned_support_states_for_export(
            layer.support_states,
            layer.support_indices,
            dims,
        )
        for (local_index, state) in pairs(layer_states)
            push!(
                points,
                BondAlignedDiatomicGeometryPoint3D(
                    Float64(x_centers[state[1]]),
                    Float64(y_centers[state[2]]),
                    Float64(z_centers[state[3]]),
                    :source_region,
                    "shared_shell_region_$(layer_index)_$(local_index)",
                    :shared_shell_region,
                    layer_index,
                ),
            )
        end
    end

    for (child_index, child) in pairs(source.child_sequences)
        group_kind = if source.geometry.did_split
            child_index == 1 ? :left_child_region : :right_child_region
        else
            :shared_child_region
        end
        label_prefix = if source.geometry.did_split
            child_index == 1 ? "left_child_region" : "right_child_region"
        else
            "shared_child_region"
        end
        child_states = _bond_aligned_support_states_for_export(
            child.support_states,
            child.support_indices,
            dims,
        )
        for (local_index, state) in pairs(child_states)
            push!(
                points,
                BondAlignedDiatomicGeometryPoint3D(
                    Float64(x_centers[state[1]]),
                    Float64(y_centers[state[2]]),
                    Float64(z_centers[state[3]]),
                    :source_region,
                    "$(label_prefix)_$(local_index)",
                    group_kind,
                    child_index,
                ),
            )
        end
    end

    if !isnothing(source.geometry.shared_midpoint_box)
        midpoint_states = [
            _cartesian_unflat_index(index, dims) for index in
            _nested_box_support_indices(source.geometry.shared_midpoint_box..., dims)
        ]
        for (local_index, state) in pairs(midpoint_states)
            push!(
                points,
                BondAlignedDiatomicGeometryPoint3D(
                    Float64(x_centers[state[1]]),
                    Float64(y_centers[state[2]]),
                    Float64(z_centers[state[3]]),
                    :source_region,
                    "shared_midpoint_slab_region_$(local_index)",
                    :shared_midpoint_slab_region,
                    1,
                ),
            )
        end
    end

    return points
end

function bond_aligned_diatomic_geometry_payload(
    basis::BondAlignedDiatomicQWBasis3D,
)
    return BondAlignedDiatomicGeometryPayload3D(
        _bond_aligned_ordinary_points(basis),
        _bond_aligned_geometry_nuclei(basis.nuclei),
        basis.bond_axis,
        BondAlignedDiatomicGeometryBox3D[
            _bond_aligned_geometry_box(
                basis,
                _bond_aligned_full_box(basis),
                "basis_box",
                :basis_box,
                1,
            ),
        ],
        BondAlignedDiatomicGeometryShellProvenance3D[],
    )
end

function bond_aligned_diatomic_geometry_payload(
    ops::OrdinaryCartesianOperators3D,
)
    basis = ops.basis
    basis isa BondAlignedDiatomicQWBasis3D || throw(
        ArgumentError("bond_aligned_diatomic_geometry_payload(ops) currently supports only ordinary bond-aligned diatomic QW operators"),
    )
    return BondAlignedDiatomicGeometryPayload3D(
        _bond_aligned_qw_points(ops),
        _bond_aligned_geometry_nuclei(basis.nuclei),
        basis.bond_axis,
        BondAlignedDiatomicGeometryBox3D[
            _bond_aligned_geometry_box(
                basis,
                _bond_aligned_full_box(basis),
                "basis_box",
                :basis_box,
                1,
            ),
        ],
        BondAlignedDiatomicGeometryShellProvenance3D[],
    )
end

function bond_aligned_diatomic_geometry_payload(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    fixed_centers = hcat(
        diag(source.sequence.packet.position_x),
        diag(source.sequence.packet.position_y),
        diag(source.sequence.packet.position_z),
    )
    return BondAlignedDiatomicGeometryPayload3D(
        _bond_aligned_nested_points(fixed_centers, source),
        _bond_aligned_geometry_nuclei(source.basis.nuclei),
        source.basis.bond_axis,
        _bond_aligned_nested_box_outlines(source),
        _bond_aligned_nested_shell_provenance(source),
    )
end

"""
    bond_aligned_diatomic_source_geometry_payload(source)

Emit the underlying raw region/source geometry for the landed bond-aligned
diatomic nested line.

Unlike [`bond_aligned_diatomic_geometry_payload(source)`](@ref), which shows
compressed fixed centers, this payload records the original parent-grid rows
assigned to the shared shell region and child regions used to define the split.
"""
function bond_aligned_diatomic_source_geometry_payload(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return BondAlignedDiatomicGeometryPayload3D(
        _bond_aligned_source_region_points(source),
        _bond_aligned_geometry_nuclei(source.basis.nuclei),
        source.basis.bond_axis,
        _bond_aligned_nested_box_outlines(source),
        _bond_aligned_nested_shell_provenance(source),
    )
end

function bond_aligned_diatomic_geometry_payload(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
)
    basis = fixed_block.parent_basis
    return BondAlignedDiatomicGeometryPayload3D(
        [
            BondAlignedDiatomicGeometryPoint3D(
                Float64(fixed_block.fixed_centers[index, 1]),
                Float64(fixed_block.fixed_centers[index, 2]),
                Float64(fixed_block.fixed_centers[index, 3]),
                :nested_fixed,
                "nf$index",
                :nested_fixed_block,
                1,
            ) for index in axes(fixed_block.fixed_centers, 1)
        ],
        _bond_aligned_geometry_nuclei(basis.nuclei),
        basis.bond_axis,
        BondAlignedDiatomicGeometryBox3D[
            _bond_aligned_geometry_box(
                basis,
                _bond_aligned_full_box(basis),
                "basis_box",
                :basis_box,
                1,
            ),
        ],
        BondAlignedDiatomicGeometryShellProvenance3D[],
    )
end

function bond_aligned_diatomic_geometry_payload(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    _bond_aligned_same_basis_geometry(fixed_block.parent_basis, source.basis) || throw(
        ArgumentError("bond-aligned nested geometry payload requires a fixed block and source built from the same bond-aligned parent basis geometry"),
    )
    return BondAlignedDiatomicGeometryPayload3D(
        _bond_aligned_nested_points(fixed_block.fixed_centers, source),
        _bond_aligned_geometry_nuclei(source.basis.nuclei),
        source.basis.bond_axis,
        _bond_aligned_nested_box_outlines(source),
        _bond_aligned_nested_shell_provenance(source),
    )
end

function bond_aligned_diatomic_geometry_payload(
    ops::OrdinaryCartesianOperators3D,
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    fixed_block = ops.basis
    fixed_block isa _NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D} || throw(
        ArgumentError("bond_aligned_diatomic_geometry_payload(ops, source) currently supports only nested bond-aligned diatomic QW operators"),
    )
    fixed_points = _bond_aligned_nested_points(fixed_block.fixed_centers, source)
    residual_points = BondAlignedDiatomicGeometryPoint3D[]
    residual_group = 0
    for orbital in Iterators.drop(ops.orbital_data, ops.gausslet_count)
        orbital.kind == :residual_gaussian || continue
        residual_group += 1
        push!(
            residual_points,
            BondAlignedDiatomicGeometryPoint3D(
                orbital.x,
                orbital.y,
                orbital.z,
                orbital.kind,
                orbital.label,
                :residual_gaussian,
                residual_group,
            ),
        )
    end
    return BondAlignedDiatomicGeometryPayload3D(
        vcat(fixed_points, residual_points),
        _bond_aligned_geometry_nuclei(source.basis.nuclei),
        source.basis.bond_axis,
        _bond_aligned_nested_box_outlines(source),
        _bond_aligned_nested_shell_provenance(source),
    )
end

"""
    bond_aligned_diatomic_plane_slice(
        payload::BondAlignedDiatomicGeometryPayload3D;
        plane_axis = nothing,
        plane_value = 0.0,
        plane_tol = 1.0e-8,
    )

Select one strict/near-strict representative plane from a bond-aligned
diatomic geometry payload.

The helper does not apply any faded off-plane logic; it keeps only points and
nuclei within `plane_tol` of the requested plane.
"""
function bond_aligned_diatomic_plane_slice(
    payload::BondAlignedDiatomicGeometryPayload3D;
    plane_axis::Union{Nothing,Symbol} = nothing,
    plane_value::Real = 0.0,
    plane_tol::Real = 1.0e-8,
)
    axis = plane_axis === nothing ? _bond_aligned_default_plane_axis(payload.bond_axis) : plane_axis
    tol_value = Float64(plane_tol)
    tol_value >= 0.0 || throw(ArgumentError("bond-aligned diatomic plane selection requires plane_tol >= 0"))
    plane_value_float = Float64(plane_value)

    selected_points = [
        point for point in payload.points if
        abs(_bond_aligned_point_coordinate(point, axis) - plane_value_float) <= tol_value
    ]
    selected_nuclei = [
        nucleus for nucleus in payload.nuclei if
        abs(_bond_aligned_nucleus_coordinate(nucleus, axis) - plane_value_float) <= tol_value
    ]
    return BondAlignedDiatomicGeometryPlaneSlice3D(
        selected_points,
        selected_nuclei,
        payload.bond_axis,
        payload.box_outlines,
        payload.shell_provenance,
        axis,
        plane_value_float,
        tol_value,
        length(selected_points),
        length(payload.points),
    )
end
