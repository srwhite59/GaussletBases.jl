function _bond_aligned_projection_axes(
    plane_axis::Symbol,
)
    plane_axis == :x && return (:y, :z)
    plane_axis == :y && return (:x, :z)
    plane_axis == :z && return (:x, :y)
    throw(ArgumentError("bond-aligned projection requires plane_axis = :x, :y, or :z"))
end

function _bond_aligned_project_point(
    point::BondAlignedDiatomicGeometryPoint3D,
    plane_axis::Symbol,
)
    axis1, axis2 = _bond_aligned_projection_axes(plane_axis)
    return (
        _bond_aligned_point_coordinate(point, axis1),
        _bond_aligned_point_coordinate(point, axis2),
    )
end

function _bond_aligned_project_nucleus(
    nucleus::BondAlignedDiatomicGeometryNucleus3D,
    plane_axis::Symbol,
)
    axis1, axis2 = _bond_aligned_projection_axes(plane_axis)
    return (
        _bond_aligned_nucleus_coordinate(nucleus, axis1),
        _bond_aligned_nucleus_coordinate(nucleus, axis2),
    )
end

function _bond_aligned_box_plane_range(
    box::BondAlignedDiatomicGeometryBox3D,
    plane_axis::Symbol,
)
    plane_axis == :x && return box.x_range
    plane_axis == :y && return box.y_range
    plane_axis == :z && return box.z_range
    throw(ArgumentError("bond-aligned box projection requires plane_axis = :x, :y, or :z"))
end

function _bond_aligned_project_box_outline(
    box::BondAlignedDiatomicGeometryBox3D,
    plane_axis::Symbol,
)
    axis1, axis2 = _bond_aligned_projection_axes(plane_axis)
    range1 = axis1 == :x ? box.x_range : axis1 == :y ? box.y_range : box.z_range
    range2 = axis2 == :x ? box.x_range : axis2 == :y ? box.y_range : box.z_range
    low1, high1 = range1
    low2, high2 = range2
    return [
        (low1, low2),
        (high1, low2),
        (high1, high2),
        (low1, high2),
        (low1, low2),
    ]
end

function _bond_aligned_group_keys(
    points::AbstractVector{BondAlignedDiatomicGeometryPoint3D},
)
    keys = Tuple{Symbol,Int}[]
    seen = Set{Tuple{Symbol,Int}}()
    for point in points
        key = (point.group_kind, point.group_id)
        if !(key in seen)
            push!(keys, key)
            push!(seen, key)
        end
    end
    return keys
end

function _bond_aligned_box_metadata_line(
    box::BondAlignedDiatomicGeometryBox3D,
)
    return string(
        "# box label=", box.label,
        " group_kind=", box.group_kind,
        " group_id=", box.group_id,
        " ix=", first(box.ix), ":", last(box.ix),
        " iy=", first(box.iy), ":", last(box.iy),
        " iz=", first(box.iz), ":", last(box.iz),
        " x_range=", box.x_range[1], ":", box.x_range[2],
        " y_range=", box.y_range[1], ":", box.y_range[2],
        " z_range=", box.z_range[1], ":", box.z_range[2],
    )
end

function _bond_aligned_shell_provenance_metadata_line(
    shell::BondAlignedDiatomicGeometryShellProvenance3D,
)
    return string(
        "# shell label=", shell.label,
        " group_kind=", shell.group_kind,
        " group_id=", shell.group_id,
        " source_box=", _nested_box_dimension_string(shell.source_box),
        " source_points=", shell.source_point_count,
        " retained_fixed_count=", shell.retained_fixed_count,
        " next_inner_box=", _nested_box_dimension_string(shell.next_inner_box),
    )
end

function _bond_aligned_points3d_lines(
    payload::BondAlignedDiatomicGeometryPayload3D;
    include_box_metadata::Bool = true,
)
    lines = String[
        "# GaussletBases bond-aligned diatomic 3d geometry",
        "# bond_axis = $(payload.bond_axis)",
        "# point_count = $(length(payload.points))",
        "# nucleus_count = $(length(payload.nuclei))",
        "# columns = x y z role kind group_kind group_id label",
    ]

    for shell in payload.shell_provenance
        push!(lines, _bond_aligned_shell_provenance_metadata_line(shell))
    end

    if include_box_metadata
        for box in payload.box_outlines
            push!(lines, _bond_aligned_box_metadata_line(box))
        end
    end

    for point in payload.points
        push!(
            lines,
            string(
                point.x, "\t",
                point.y, "\t",
                point.z, "\tpoint\t",
                point.kind, "\t",
                point.group_kind, "\t",
                point.group_id, "\t",
                point.label,
            ),
        )
    end

    for nucleus in payload.nuclei
        push!(
            lines,
            string(
                nucleus.x, "\t",
                nucleus.y, "\t",
                nucleus.z, "\tnucleus\tnucleus\tnucleus\t",
                nucleus.group_id, "\t",
                nucleus.label,
            ),
        )
    end

    return lines
end

function _bond_aligned_plane_projection_lines(
    slice::BondAlignedDiatomicGeometryPlaneSlice3D;
    include_box_outlines::Bool = false,
)
    axis1, axis2 = _bond_aligned_projection_axes(slice.plane_axis)
    lines = String[
        "# GaussletBases bond-aligned diatomic plane projection",
        "# bond_axis = $(slice.bond_axis)",
        "# plane_axis = $(slice.plane_axis)",
        "# plane_value = $(slice.plane_value)",
        "# plane_tol = $(slice.plane_tol)",
        "# selected_count = $(slice.selected_count)",
        "# total_count = $(slice.total_count)",
        "# projection_axes = $(axis1) $(axis2)",
        "# columns = $(axis1) $(axis2)",
    ]

    for shell in slice.shell_provenance
        push!(lines, _bond_aligned_shell_provenance_metadata_line(shell))
    end

    dataset_index = 0
    for key in _bond_aligned_group_keys(slice.points)
        group_points = [point for point in slice.points if (point.group_kind, point.group_id) == key]
        isempty(group_points) && continue
        dataset_index += 1
        push!(
            lines,
            "# dataset $(dataset_index) role=point group_kind=$(key[1]) group_id=$(key[2]) kind=$(group_points[1].kind) count=$(length(group_points))",
        )
        for point in group_points
            u, v = _bond_aligned_project_point(point, slice.plane_axis)
            push!(lines, "$(u) $(v)")
        end
        push!(lines, "@")
    end

    for nucleus in slice.nuclei
        dataset_index += 1
        u, v = _bond_aligned_project_nucleus(nucleus, slice.plane_axis)
        push!(
            lines,
            "# dataset $(dataset_index) role=nucleus group_kind=nucleus group_id=$(nucleus.group_id) label=$(nucleus.label) count=1",
        )
        push!(lines, "$(u) $(v)")
        push!(lines, "@")
    end

    if include_box_outlines
        for box in slice.box_outlines
            plane_range = _bond_aligned_box_plane_range(box, slice.plane_axis)
            plane_range[1] - slice.plane_tol <= slice.plane_value <= plane_range[2] + slice.plane_tol || continue
            dataset_index += 1
            push!(
                lines,
                "# dataset $(dataset_index) role=box_outline group_kind=$(box.group_kind) group_id=$(box.group_id) label=$(box.label) count=5",
            )
            for (u, v) in _bond_aligned_project_box_outline(box, slice.plane_axis)
                push!(lines, "$(u) $(v)")
            end
            push!(lines, "@")
        end
    end

    return lines
end

"""
    write_bond_aligned_diatomic_points3d(
        path,
        payload::BondAlignedDiatomicGeometryPayload3D;
        include_box_metadata = true,
    )

Write one plain-text 3D geometry artifact for a bond-aligned diatomic payload.

The file is a single tabular text block with comment metadata header lines and
one row per emitted point or nucleus.
"""
function write_bond_aligned_diatomic_points3d(
    path::AbstractString,
    payload::BondAlignedDiatomicGeometryPayload3D;
    include_box_metadata::Bool = true,
)
    lines = _bond_aligned_points3d_lines(
        payload;
        include_box_metadata = include_box_metadata,
    )
    mkpath(dirname(String(path)))
    open(path, "w") do io
        for line in lines
            write(io, line, "\n")
        end
    end
    return payload
end

"""
    write_bond_aligned_diatomic_plane_projection(
        path,
        payload::BondAlignedDiatomicGeometryPayload3D;
        plane_axis = nothing,
        plane_value = 0.0,
        plane_tol = 1.0e-8,
        include_box_outlines = false,
    )

Write one plain-text 2D projection artifact for a bond-aligned diatomic
geometry payload.

The format is intentionally lightweight:

- header comments with explicit plane metadata
- grouped numeric datasets separated by `@`
- no plotting-backend dependency
"""
function write_bond_aligned_diatomic_plane_projection(
    path::AbstractString,
    payload::BondAlignedDiatomicGeometryPayload3D;
    plane_axis::Union{Nothing,Symbol} = nothing,
    plane_value::Real = 0.0,
    plane_tol::Real = 1.0e-8,
    include_box_outlines::Bool = false,
)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = plane_axis,
        plane_value = plane_value,
        plane_tol = plane_tol,
    )
    return write_bond_aligned_diatomic_plane_projection(
        path,
        slice;
        include_box_outlines = include_box_outlines,
    )
end

function write_bond_aligned_diatomic_plane_projection(
    path::AbstractString,
    slice::BondAlignedDiatomicGeometryPlaneSlice3D;
    include_box_outlines::Bool = false,
)
    lines = _bond_aligned_plane_projection_lines(
        slice;
        include_box_outlines = include_box_outlines,
    )
    mkpath(dirname(String(path)))
    open(path, "w") do io
        for line in lines
            write(io, line, "\n")
        end
    end
    return slice
end
