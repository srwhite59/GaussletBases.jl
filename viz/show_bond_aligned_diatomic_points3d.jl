#!/usr/bin/env julia
struct GeometryPoint3D
    x::Float64
    y::Float64
    z::Float64
    role::Symbol
    kind::String
    group_kind::Symbol
    group_id::Int
    label::String
end

struct GeometryBox3D
    label::String
    group_kind::Symbol
    group_id::Int
    x_range::Tuple{Float64,Float64}
    y_range::Tuple{Float64,Float64}
    z_range::Tuple{Float64,Float64}
end

struct GeometryPath3D
    label::String
    group_kind::Symbol
    group_id::Int
    points::Vector{NTuple{3,Float64}}
end

struct GeometryPayload3D
    bond_axis::Symbol
    point_count::Int
    nucleus_count::Int
    points::Vector{GeometryPoint3D}
    boxes::Vector{GeometryBox3D}
    paths::Vector{GeometryPath3D}
end

const DEFAULT_TITLE_3D = "Bond-aligned diatomic 3D geometry"
const DEFAULT_PATH_TITLE_3D = "Ordering path 3D"

function _parse_header_assignment_3d(line::AbstractString)
    body = strip(line[2:end])
    startswith(body, "box ") && return nothing
    startswith(body, "path ") && return nothing
    startswith(body, "columns =") && return nothing
    occursin("=", body) || return nothing
    key, value = split(body, "="; limit = 2)
    return Symbol(strip(key)), strip(value)
end

function _parse_range_pair(value::AbstractString)
    low, high = split(value, ":"; limit = 2)
    return (parse(Float64, low), parse(Float64, high))
end

function _parse_box_metadata_3d(line::AbstractString)
    startswith(line, "# box ") || throw(ArgumentError("expected box metadata line"))
    tokens = split(strip(line[3:end]))
    metadata = Dict{String,String}()
    for token in tokens[2:end]
        if occursin("=", token)
            key, value = split(token, "="; limit = 2)
            metadata[key] = value
        end
    end
    return GeometryBox3D(
        get(metadata, "label", "box"),
        Symbol(get(metadata, "group_kind", "box_outline")),
        parse(Int, get(metadata, "group_id", "0")),
        _parse_range_pair(get(metadata, "x_range", "0.0:0.0")),
        _parse_range_pair(get(metadata, "y_range", "0.0:0.0")),
        _parse_range_pair(get(metadata, "z_range", "0.0:0.0")),
    )
end

function _parse_path_metadata_3d(line::AbstractString)
    startswith(line, "# path ") || throw(ArgumentError("expected path metadata line"))
    tokens = split(strip(line[3:end]))
    metadata = Dict{String,String}()
    for token in tokens[2:end]
        if occursin("=", token)
            key, value = split(token, "="; limit = 2)
            metadata[key] = value
        end
    end
    return (
        label = get(metadata, "label", "path"),
        group_kind = Symbol(get(metadata, "group_kind", "ordering_path")),
        group_id = parse(Int, get(metadata, "group_id", "0")),
    )
end

function read_bond_aligned_diatomic_points3d(path::AbstractString)
    metadata = Dict{Symbol,Any}()
    points = GeometryPoint3D[]
    boxes = GeometryBox3D[]
    paths = GeometryPath3D[]
    pending_path = nothing
    pending_path_points = NTuple{3,Float64}[]

    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        if startswith(line, "#")
            assignment = _parse_header_assignment_3d(line)
            if !isnothing(assignment)
                key, value = assignment
                if key == :bond_axis
                    metadata[key] = Symbol(value)
                elseif key == :point_count || key == :nucleus_count
                    metadata[key] = parse(Int, value)
                end
                continue
            end
            if startswith(line, "# box ")
                push!(boxes, _parse_box_metadata_3d(line))
            elseif startswith(line, "# path ")
                pending_path = _parse_path_metadata_3d(line)
                empty!(pending_path_points)
            end
            continue
        end

        if line == "@"
            isnothing(pending_path) && continue
            push!(
                paths,
                GeometryPath3D(
                    pending_path.label,
                    pending_path.group_kind,
                    pending_path.group_id,
                    copy(pending_path_points),
                ),
            )
            pending_path = nothing
            empty!(pending_path_points)
            continue
        end

        if !isnothing(pending_path)
            tokens = split(line)
            length(tokens) == 3 || throw(ArgumentError("expected 3-column path row in $path"))
            push!(
                pending_path_points,
                (
                    parse(Float64, tokens[1]),
                    parse(Float64, tokens[2]),
                    parse(Float64, tokens[3]),
                ),
            )
            continue
        end

        tokens = split(line)
        length(tokens) == 8 || throw(ArgumentError("expected 8-column 3D row in $path"))
        push!(
            points,
            GeometryPoint3D(
                parse(Float64, tokens[1]),
                parse(Float64, tokens[2]),
                parse(Float64, tokens[3]),
                Symbol(tokens[4]),
                tokens[5],
                Symbol(tokens[6]),
                parse(Int, tokens[7]),
                tokens[8],
            ),
        )
    end

    return GeometryPayload3D(
        get(metadata, :bond_axis, :z),
        get(metadata, :point_count, count(point -> point.role == :point, points)),
        get(metadata, :nucleus_count, count(point -> point.role == :nucleus, points)),
        points,
        boxes,
        paths,
    )
end

function _group_style_3d(group_kind::Symbol)
    if group_kind == :gausslet_product
        return (color = :black, marker = :circle, markersize = 26, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :bridge_point
        return (color = RGBAf(0.35, 0.35, 0.35, 0.28), marker = :circle, markersize = 18, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :residual_gaussian
        return (color = :magenta4, marker = :star8, markersize = 34, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :left_child
        return (color = :dodgerblue3, marker = :utriangle, markersize = 28, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :shared_midpoint_slab
        return (color = :black, marker = :rect, markersize = 30, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :right_child
        return (color = :firebrick3, marker = :dtriangle, markersize = 28, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :shared_shell_layer
        return (color = :darkgreen, marker = :diamond, markersize = 28, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :left_child_region
        return (color = :dodgerblue3, marker = :utriangle, markersize = 20, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :shared_midpoint_slab_region
        return (color = :black, marker = :rect, markersize = 22, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :right_child_region
        return (color = :firebrick3, marker = :dtriangle, markersize = 20, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :shared_shell_region
        return (color = :darkgreen, marker = :diamond, markersize = 20, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :nucleus
        return (color = :goldenrod2, marker = :star5, markersize = 40, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :parent_box
        return (color = :gray35, marker = :none, markersize = 0, linewidth = 2.5, linestyle = :solid)
    elseif group_kind == :working_box
        return (color = :goldenrod3, marker = :none, markersize = 0, linewidth = 2.5, linestyle = :solid)
    elseif group_kind == :child_box
        return (color = :royalblue3, marker = :none, markersize = 0, linewidth = 2.5, linestyle = :solid)
    elseif group_kind == :shared_midpoint_slab_box
        return (color = :deeppink4, marker = :none, markersize = 0, linewidth = 3.0, linestyle = :solid)
    elseif group_kind == :ordering_path
        return (color = :darkorange2, marker = :none, markersize = 0, linewidth = 4.0, linestyle = :solid)
    else
        return (color = :gray40, marker = :circle, markersize = 24, linewidth = 0.0, linestyle = :solid)
    end
end

function _legend_group_key_3d(point::GeometryPoint3D)
    point.role == :nucleus && return :nucleus
    return point.group_kind
end

function _nucleus_legend_label_3d(payload::GeometryPayload3D)
    nucleus_labels = sort(unique(point.label for point in payload.points if point.role == :nucleus))
    if length(nucleus_labels) == 2 && nucleus_labels == ["A", "B"]
        return "H nuclei"
    elseif length(nucleus_labels) == 1 && nucleus_labels[1] == "A"
        return "H nucleus"
    end
    return "nuclei"
end

function _legend_group_label_3d(payload::GeometryPayload3D, legend_key::Symbol)
    legend_key == :gausslet_product && return "gausslet product"
    legend_key == :bridge_point && return "basis points"
    legend_key == :residual_gaussian && return "residual Gaussian"
    legend_key == :left_child && return "left child"
    legend_key == :shared_midpoint_slab && return "midpoint slab"
    legend_key == :right_child && return "right child"
    legend_key == :shared_shell_layer && return "shared shell"
    legend_key == :left_child_region && return "left child region"
    legend_key == :shared_midpoint_slab_region && return "midpoint slab region"
    legend_key == :right_child_region && return "right child region"
    legend_key == :shared_shell_region && return "shared shell region"
    legend_key == :parent_box && return "parent box"
    legend_key == :working_box && return "working box"
    legend_key == :child_box && return "child box"
    legend_key == :shared_midpoint_slab_box && return "midpoint slab box"
    legend_key == :ordering_path && return "ordering path"
    legend_key == :nucleus && return _nucleus_legend_label_3d(payload)
    return string(legend_key)
end

function _path_color_3d(group_id::Int)
    palette = [
        :darkorange2,
        :dodgerblue3,
        :firebrick3,
        :darkgreen,
        :purple4,
        :goldenrod3,
    ]
    return palette[mod1(group_id, length(palette))]
end

function _primary_ordering_path_payload_3d(payload::GeometryPayload3D)
    selected_path = nothing
    for path in payload.paths
        isempty(path.points) && continue
        selected_path = path
        break
    end
    selected_paths = isnothing(selected_path) ? GeometryPath3D[] : [selected_path]
    nucleus_points = [point for point in payload.points if point.role == :nucleus]
    return GeometryPayload3D(
        payload.bond_axis,
        payload.point_count,
        payload.nucleus_count,
        nucleus_points,
        GeometryBox3D[],
        selected_paths,
    )
end

function _interpolate_rgba_3d(c1, c2, t::Float64, GLM)
    s = clamp(t, 0.0, 1.0)
    return GLM.RGBAf(
        (1 - s) * c1.r + s * c2.r,
        (1 - s) * c1.g + s * c2.g,
        (1 - s) * c1.b + s * c2.b,
        (1 - s) * c1.alpha + s * c2.alpha,
    )
end

function _path_segment_palette_3d(nsegments::Int, GLM)
    nsegments <= 0 && return GLM.RGBAf[]
    anchors = [
        GLM.RGBAf(0.00, 0.38, 1.00, 1.0),
        GLM.RGBAf(0.00, 0.78, 1.00, 1.0),
        GLM.RGBAf(0.00, 0.80, 0.52, 1.0),
        GLM.RGBAf(0.55, 0.82, 0.12, 1.0),
        GLM.RGBAf(0.98, 0.78, 0.05, 1.0),
        GLM.RGBAf(1.00, 0.42, 0.00, 1.0),
        GLM.RGBAf(0.92, 0.12, 0.62, 1.0),
        GLM.RGBAf(0.88, 0.12, 0.18, 1.0),
    ]
    if nsegments == 1
        return [anchors[1]]
    end

    colors = Vector{typeof(anchors[1])}(undef, nsegments)
    nblocks = length(anchors) - 1
    @inbounds for idx in 1:nsegments
        t = (idx - 1) / (nsegments - 1)
        scaled = t * nblocks
        block = min(floor(Int, scaled) + 1, nblocks)
        local_t = scaled - (block - 1)
        colors[idx] = _interpolate_rgba_3d(anchors[block], anchors[block + 1], local_t, GLM)
    end
    return colors
end

function _path_point_palette_3d(npoints::Int, GLM)
    npoints <= 0 && return GLM.RGBAf[]
    npoints == 1 && return [GLM.RGBAf(0.10, 0.35, 0.95, 1.0)]
    segment_colors = _path_segment_palette_3d(npoints - 1, GLM)
    return vcat(segment_colors, [segment_colors[end]])
end

function _legend_elements_3d(payload::GeometryPayload3D, GLM)
    elements = Any[]
    labels = String[]
    seen = Set{Symbol}()

    for point in payload.points
        key = _legend_group_key_3d(point)
        key in seen && continue
        push!(seen, key)
        style = _group_style_3d(key)
        push!(elements, GLM.MarkerElement(color = style.color, marker = style.marker, markersize = style.markersize))
        push!(labels, _legend_group_label_3d(payload, key))
    end

    for box in payload.boxes
        key = box.group_kind
        key in seen && continue
        push!(seen, key)
        style = _group_style_3d(key)
        push!(elements, GLM.LineElement(color = style.color, linewidth = style.linewidth, linestyle = style.linestyle))
        push!(labels, _legend_group_label_3d(payload, key))
    end

    for path in payload.paths
        push!(
            elements,
            GLM.LineElement(
                color = _path_color_3d(path.group_id),
                linewidth = _group_style_3d(path.group_kind).linewidth,
                linestyle = _group_style_3d(path.group_kind).linestyle,
            ),
        )
        push!(labels, path.label)
    end

    return elements, labels
end

function _box_edges_3d(box::GeometryBox3D)
    xlo, xhi = box.x_range
    ylo, yhi = box.y_range
    zlo, zhi = box.z_range
    corners = [
        (xlo, ylo, zlo),
        (xhi, ylo, zlo),
        (xhi, yhi, zlo),
        (xlo, yhi, zlo),
        (xlo, ylo, zhi),
        (xhi, ylo, zhi),
        (xhi, yhi, zhi),
        (xlo, yhi, zhi),
    ]
    edge_pairs = [
        (1, 2), (2, 3), (3, 4), (4, 1),
        (5, 6), (6, 7), (7, 8), (8, 5),
        (1, 5), (2, 6), (3, 7), (4, 8),
    ]
    return [(corners[i], corners[j]) for (i, j) in edge_pairs]
end

function _geometry_bounds_3d(payload::GeometryPayload3D)
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    for point in payload.points
        push!(xs, point.x)
        push!(ys, point.y)
        push!(zs, point.z)
    end
    for box in payload.boxes
        append!(xs, [box.x_range[1], box.x_range[2]])
        append!(ys, [box.y_range[1], box.y_range[2]])
        append!(zs, [box.z_range[1], box.z_range[2]])
    end
    for path in payload.paths
        for (x, y, z) in path.points
            push!(xs, x)
            push!(ys, y)
            push!(zs, z)
        end
    end
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    zmin, zmax = extrema(zs)
    xpad = max((xmax - xmin) * 0.08, 0.5)
    ypad = max((ymax - ymin) * 0.08, 0.5)
    zpad = max((zmax - zmin) * 0.08, 0.5)
    return (
        xmin - xpad,
        xmax + xpad,
        ymin - ypad,
        ymax + ypad,
        zmin - zpad,
        zmax + zpad,
    )
end

function _load_glmakie()
    isdefined(Main, :GLMakie) || Core.eval(Main, :(using GLMakie))
    return getfield(Main, :GLMakie)
end

function render_bond_aligned_diatomic_points3d(
    input_path::AbstractString;
    title::AbstractString = DEFAULT_TITLE_3D,
)
    payload = read_bond_aligned_diatomic_points3d(input_path)
    GLM = _load_glmakie()
    GLM.activate!(; focus_on_show = true)

    limits = _geometry_bounds_3d(payload)
    figure = GLM.Figure(size = (1600, 1100), fontsize = 22)
    GLM.Label(figure[0, 1], title; fontsize = 28, tellwidth = false)
    GLM.Label(
        figure[1, 1],
        "bond=$(payload.bond_axis), final gausslet count=$(payload.point_count), + $(payload.nucleus_count) nuclei";
        fontsize = 20,
        tellwidth = false,
    )

    axis = GLM.Axis3(
        figure[2, 1];
        xlabel = "x",
        ylabel = "y",
        zlabel = "z",
        aspect = :data,
        perspectiveness = 0.9,
        limits = limits,
        protrusions = (40, 40, 80, 20),
    )

    for box in payload.boxes
        style = _group_style_3d(box.group_kind)
        for ((x1, y1, z1), (x2, y2, z2)) in _box_edges_3d(box)
            GLM.lines!(
                axis,
                [x1, x2],
                [y1, y2],
                [z1, z2];
                color = style.color,
                linewidth = style.linewidth,
                linestyle = style.linestyle,
            )
        end
    end

    for point in payload.points
        style = _group_style_3d(_legend_group_key_3d(point))
        GLM.scatter!(
            axis,
            [point.x],
            [point.y],
            [point.z];
            color = style.color,
            marker = style.marker,
            markersize = style.markersize,
        )
    end

    for path in payload.paths
        isempty(path.points) && continue
        style = _group_style_3d(path.group_kind)
        xs = [point[1] for point in path.points]
        ys = [point[2] for point in path.points]
        zs = [point[3] for point in path.points]
        point_colors = _path_point_palette_3d(length(path.points), GLM)
        segment_colors = _path_segment_palette_3d(length(path.points) - 1, GLM)
        @inbounds for idx in 1:(length(path.points) - 1)
            GLM.lines!(
                axis,
                [xs[idx], xs[idx + 1]],
                [ys[idx], ys[idx + 1]],
                [zs[idx], zs[idx + 1]];
                color = segment_colors[idx],
                linewidth = style.linewidth,
                linestyle = style.linestyle,
            )
        end
        GLM.scatter!(
            axis,
            xs,
            ys,
            zs;
            color = point_colors,
            marker = :circle,
            markersize = 12,
        )
        GLM.scatter!(
            axis,
            [xs[1]],
            [ys[1]],
            [zs[1]];
            color = :gold,
            marker = :circle,
            markersize = 26,
        )
        GLM.scatter!(
            axis,
            [xs[end]],
            [ys[end]],
            [zs[end]];
            color = :black,
            marker = :rect,
            markersize = 22,
        )
    end

    elements, labels = _legend_elements_3d(payload, GLM)
    GLM.Legend(figure[2, 2], elements, labels; tellheight = false)

    screen = GLM.display(figure)
    wait(screen)
    return payload
end

function _ordering_path_legend_elements_3d(payload::GeometryPayload3D, GLM)
    elements = Any[]
    labels = String[]
    if any(point -> point.role == :nucleus, payload.points)
        style = _group_style_3d(:nucleus)
        push!(elements, GLM.MarkerElement(color = style.color, marker = style.marker, markersize = style.markersize))
        push!(labels, _nucleus_legend_label_3d(payload))
    end
    for path in payload.paths
        push!(
            elements,
            GLM.LineElement(
                color = isempty(path.points) ? :darkorange2 : _path_segment_palette_3d(max(length(path.points) - 1, 1), GLM)[clamp(length(path.points) ÷ 2, 1, max(length(path.points) - 1, 1))],
                linewidth = _group_style_3d(path.group_kind).linewidth,
                linestyle = _group_style_3d(path.group_kind).linestyle,
            ),
        )
        push!(labels, "ordering path")
    end
    return elements, labels
end

function render_ordering_path3d(
    input_path::AbstractString;
    title::AbstractString = DEFAULT_PATH_TITLE_3D,
)
    payload = _primary_ordering_path_payload_3d(read_bond_aligned_diatomic_points3d(input_path))
    GLM = _load_glmakie()
    GLM.activate!(; focus_on_show = true)

    limits = _geometry_bounds_3d(payload)
    figure = GLM.Figure(size = (1500, 1100), fontsize = 22)
    GLM.Label(figure[0, 1], title; fontsize = 28, tellwidth = false)
    GLM.Label(
        figure[1, 1],
        "showing $(length(payload.paths)) path, nuclei=$(payload.nucleus_count), point cloud hidden";
        fontsize = 20,
        tellwidth = false,
    )

    axis = GLM.Axis3(
        figure[2, 1];
        xlabel = "x",
        ylabel = "y",
        zlabel = "z",
        aspect = :data,
        perspectiveness = 1.0,
        limits = limits,
        protrusions = (40, 40, 80, 20),
    )

    for point in payload.points
        point.role == :nucleus || continue
        style = _group_style_3d(:nucleus)
        GLM.scatter!(
            axis,
            [point.x],
            [point.y],
            [point.z];
            color = style.color,
            marker = style.marker,
            markersize = style.markersize,
        )
    end

    for path in payload.paths
        isempty(path.points) && continue
        style = _group_style_3d(path.group_kind)
        xs = [point[1] for point in path.points]
        ys = [point[2] for point in path.points]
        zs = [point[3] for point in path.points]
        point_colors = _path_point_palette_3d(length(path.points), GLM)
        segment_colors = _path_segment_palette_3d(length(path.points) - 1, GLM)
        path_linewidth = max(style.linewidth, 3.15)
        @inbounds for idx in 1:(length(path.points) - 1)
            GLM.lines!(
                axis,
                [xs[idx], xs[idx + 1]],
                [ys[idx], ys[idx + 1]],
                [zs[idx], zs[idx + 1]];
                color = segment_colors[idx],
                linewidth = path_linewidth,
                linestyle = style.linestyle,
            )
        end
        GLM.scatter!(
            axis,
            xs,
            ys,
            zs;
            color = point_colors,
            marker = :circle,
            markersize = 18,
        )
        GLM.scatter!(
            axis,
            [xs[1]],
            [ys[1]],
            [zs[1]];
            color = :gold,
            marker = :circle,
            markersize = 26,
        )
        GLM.scatter!(
            axis,
            [xs[end]],
            [ys[end]],
            [zs[end]];
            color = :black,
            marker = :rect,
            markersize = 22,
        )
        mid = cld(length(xs), 2)
        GLM.scatter!(
            axis,
            [xs[mid]],
            [ys[mid]],
            [zs[mid]];
            color = :white,
            marker = :diamond,
            markersize = 20,
            strokecolor = :black,
            strokewidth = 2,
        )
    end

    elements, labels = _ordering_path_legend_elements_3d(payload, GLM)
    isempty(elements) || GLM.Legend(figure[2, 2], elements, labels; tellheight = false)

    screen = GLM.display(figure)
    wait(screen)
    return payload
end

function describe_bond_aligned_diatomic_points3d(path::AbstractString)
    payload = read_bond_aligned_diatomic_points3d(path)
    println("input_path=", path)
    println("bond_axis=", payload.bond_axis)
    println("point_count=", payload.point_count)
    println("nucleus_count=", payload.nucleus_count)
    println("box_count=", length(payload.boxes))
    println("path_count=", length(payload.paths))
    for point in payload.points
        point.role == :point || continue
    end
    seen_groups = Set{Tuple{Symbol,Int}}()
    for point in payload.points
        point.role == :point || continue
        key = (point.group_kind, point.group_id)
        key in seen_groups && continue
        push!(seen_groups, key)
        count_in_group = count(p -> p.role == :point && p.group_kind == point.group_kind && p.group_id == point.group_id, payload.points)
        println("group_kind=", point.group_kind, " group_id=", point.group_id, " count=", count_in_group)
    end
    for box in payload.boxes
        println(
            "box label=", box.label,
            " group_kind=", box.group_kind,
            " group_id=", box.group_id,
            " x_range=", box.x_range[1], ":", box.x_range[2],
            " y_range=", box.y_range[1], ":", box.y_range[2],
            " z_range=", box.z_range[1], ":", box.z_range[2],
        )
    end
    for path in payload.paths
        println("path label=", path.label, " group_kind=", path.group_kind, " group_id=", path.group_id, " count=", length(path.points))
    end
    return payload
end

function _usage_3d()
    println("usage: viz/showpoints3d.jl [--describe] input.dat [title]")
end

function main_show_bond_aligned_diatomic_points3d(args = ARGS)
    describe_only = false
    args = copy(args)
    if !isempty(args) && args[1] == "--describe"
        describe_only = true
        popfirst!(args)
    end

    if !(1 <= length(args) <= 2)
        _usage_3d()
        exit(1)
    end

    input_path = args[1]
    title = length(args) == 2 ? args[2] : DEFAULT_TITLE_3D

    if describe_only
        describe_bond_aligned_diatomic_points3d(input_path)
    else
        render_bond_aligned_diatomic_points3d(
            input_path;
            title = title,
        )
    end
end

function _usage_path_3d()
    println("usage: viz/showpath3d.jl [--describe] input.dat [title]")
end

function main_show_ordering_path3d(args = ARGS)
    describe_only = false
    args = copy(args)
    if !isempty(args) && args[1] == "--describe"
        describe_only = true
        popfirst!(args)
    end

    if !(1 <= length(args) <= 2)
        _usage_path_3d()
        exit(1)
    end

    input_path = args[1]
    title = length(args) == 2 ? args[2] : DEFAULT_PATH_TITLE_3D

    if describe_only
        describe_bond_aligned_diatomic_points3d(input_path)
    else
        render_ordering_path3d(
            input_path;
            title = title,
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_show_bond_aligned_diatomic_points3d(ARGS)
end
