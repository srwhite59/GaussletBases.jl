using CairoMakie

struct ProjectionDataset2D
    role::Symbol
    group_kind::Symbol
    group_id::Int
    kind::String
    label::String
    count::Int
    points::Vector{Tuple{Float64,Float64}}
end

struct ProjectionPayload2D
    bond_axis::Symbol
    plane_axis::Symbol
    plane_value::Float64
    plane_tol::Float64
    selected_count::Int
    total_count::Int
    projection_axes::Tuple{Symbol,Symbol}
    datasets::Vector{ProjectionDataset2D}
end

const DEFAULT_TITLE = "Bond-aligned diatomic projection"

function _nucleus_legend_label(
    payload::ProjectionPayload2D,
)
    nucleus_labels = sort(unique(dataset.label for dataset in payload.datasets if dataset.role == :nucleus))
    if length(nucleus_labels) == 2 && nucleus_labels == ["A", "B"]
        return "H"
    elseif length(nucleus_labels) == 1 && nucleus_labels[1] == "A"
        return "H"
    end
    return "nuclei"
end

function _parse_header_assignment(line::AbstractString)
    body = strip(line[2:end])
    startswith(body, "dataset ") && return nothing
    startswith(body, "columns =") && return nothing
    occursin("=", body) || return nothing
    key, value = split(body, "="; limit = 2)
    return Symbol(strip(key)), strip(value)
end

function _parse_dataset_metadata(line::AbstractString)
    startswith(line, "# dataset ") || throw(ArgumentError("expected dataset metadata line"))
    tokens = split(strip(line[3:end]))
    metadata = Dict{String,String}()
    dataset_number = parse(Int, tokens[2])
    metadata["dataset"] = string(dataset_number)
    for token in tokens[3:end]
        if occursin("=", token)
            key, value = split(token, "="; limit = 2)
            metadata[key] = value
        end
    end
    return metadata
end

function read_bond_aligned_diatomic_projection(path::AbstractString)
    metadata = Dict{Symbol,Any}()
    datasets = ProjectionDataset2D[]
    pending_metadata = nothing
    pending_points = Tuple{Float64,Float64}[]

    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        if startswith(line, "#")
            assignment = _parse_header_assignment(line)
            if !isnothing(assignment)
                key, value = assignment
                if key == :bond_axis || key == :plane_axis
                    metadata[key] = Symbol(value)
                elseif key == :plane_value || key == :plane_tol
                    metadata[key] = parse(Float64, value)
                elseif key == :selected_count || key == :total_count
                    metadata[key] = parse(Int, value)
                elseif key == :projection_axes
                    axis_tokens = split(value)
                    metadata[key] = (Symbol(axis_tokens[1]), Symbol(axis_tokens[2]))
                end
                continue
            end
            if startswith(line, "# dataset ")
                pending_metadata = _parse_dataset_metadata(line)
                empty!(pending_points)
            end
            continue
        end

        if line == "@"
            isnothing(pending_metadata) && continue
            push!(
                datasets,
                ProjectionDataset2D(
                    Symbol(get(pending_metadata, "role", "point")),
                    Symbol(get(pending_metadata, "group_kind", "unknown")),
                    parse(Int, get(pending_metadata, "group_id", "0")),
                    get(pending_metadata, "kind", ""),
                    get(pending_metadata, "label", ""),
                    parse(Int, get(pending_metadata, "count", string(length(pending_points)))),
                    copy(pending_points),
                ),
            )
            pending_metadata = nothing
            empty!(pending_points)
            continue
        end

        xy = split(line)
        length(xy) == 2 || throw(ArgumentError("expected 2-column projection row in $path"))
        push!(pending_points, (parse(Float64, xy[1]), parse(Float64, xy[2])))
    end

    return ProjectionPayload2D(
        get(metadata, :bond_axis, :z),
        get(metadata, :plane_axis, :y),
        get(metadata, :plane_value, 0.0),
        get(metadata, :plane_tol, 0.0),
        get(metadata, :selected_count, 0),
        get(metadata, :total_count, 0),
        get(metadata, :projection_axes, (:x, :z)),
        datasets,
    )
end

function _group_style(group_kind::Symbol)
    if group_kind == :gausslet_product
        return (color = :black, marker = :circle, markersize = 20, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :residual_gaussian
        return (color = :magenta4, marker = :star8, markersize = 28, strokewidth = 2.5, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :left_child
        return (color = :dodgerblue3, marker = :utriangle, markersize = 24, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :shared_midpoint_slab
        return (color = :black, marker = :rect, markersize = 24, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :right_child
        return (color = :firebrick3, marker = :dtriangle, markersize = 24, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :shared_shell_layer
        return (color = :darkgreen, marker = :diamond, markersize = 24, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :left_child_region
        return (color = :dodgerblue3, marker = :utriangle, markersize = 16, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :shared_midpoint_slab_region
        return (color = :black, marker = :rect, markersize = 16, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :right_child_region
        return (color = :firebrick3, marker = :dtriangle, markersize = 16, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :shared_shell_region
        return (color = :darkgreen, marker = :diamond, markersize = 16, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :nucleus
        return (color = :black, marker = :star5, markersize = 34, strokewidth = 2.0, linewidth = 0.0, linestyle = :solid)
    elseif group_kind == :parent_box
        return (color = :gray35, marker = :none, markersize = 0, strokewidth = 1.0, linewidth = 3.0, linestyle = :solid)
    elseif group_kind == :working_box
        return (color = :goldenrod3, marker = :none, markersize = 0, strokewidth = 1.0, linewidth = 3.0, linestyle = :solid)
    elseif group_kind == :child_box
        return (color = :royalblue3, marker = :none, markersize = 0, strokewidth = 1.0, linewidth = 3.0, linestyle = :solid)
    elseif group_kind == :shared_midpoint_slab_box
        return (color = :deeppink4, marker = :none, markersize = 0, strokewidth = 1.0, linewidth = 4.0, linestyle = :solid)
    else
        return (color = :gray40, marker = :circle, markersize = 20, strokewidth = 1.0, linewidth = 0.0, linestyle = :solid)
    end
end

function _skip_dataset(dataset::ProjectionDataset2D)
    return dataset.group_kind == :working_box
end

function _reference_grid_payload(
    input_path::AbstractString,
    payload::ProjectionPayload2D,
)
    path = String(input_path)
    if occursin("nested_hybrid", basename(path))
        source_path = joinpath(dirname(path), replace(basename(path), "nested_hybrid" => "nested_source"))
        return isfile(source_path) ? read_bond_aligned_diatomic_projection(source_path) : nothing
    elseif occursin("ordinary_hybrid", basename(path))
        return payload
    end
    return nothing
end

function _reference_grid_coordinates(
    input_path::AbstractString,
    payload::ProjectionPayload2D,
)
    grid_payload = _reference_grid_payload(input_path, payload)
    isnothing(grid_payload) && return nothing
    point_coords = Tuple{Float64,Float64}[]
    for dataset in grid_payload.datasets
        dataset.role == :point || continue
        dataset.group_kind == :residual_gaussian && continue
        append!(point_coords, dataset.points)
    end
    isempty(point_coords) && return nothing
    xs = sort(unique(first.(point_coords)))
    ys = sort(unique(last.(point_coords)))
    return xs, ys
end

function _legend_group_key(dataset::ProjectionDataset2D)
    dataset.role == :nucleus && return :nucleus
    return dataset.group_kind
end

function _legend_group_label(
    payload::ProjectionPayload2D,
    legend_key::Symbol,
)
    legend_key == :gausslet_product && return "gausslet product"
    if legend_key == :residual_gaussian
        rg_count = count(dataset -> dataset.group_kind == :residual_gaussian, payload.datasets)
        return rg_count > 1 ? "residual Gaussian (n=$(rg_count))" : "residual Gaussian"
    end
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
    legend_key == :nucleus && return _nucleus_legend_label(payload)
    return string(legend_key)
end

function _legend_element_for_group(
    legend_key::Symbol,
)
    style = _group_style(legend_key)
    if endswith(String(legend_key), "_box")
        return LineElement(color = style.color, linewidth = style.linewidth, linestyle = style.linestyle)
    end
    return MarkerElement(color = style.color, marker = style.marker, markersize = style.markersize, strokewidth = style.strokewidth)
end

function render_bond_aligned_diatomic_projection(
    input_path::AbstractString,
    output_path::AbstractString;
    title::AbstractString = DEFAULT_TITLE,
    subtitle::Union{Nothing,AbstractString} = nothing,
)
    payload = read_bond_aligned_diatomic_projection(input_path)
    axis1, axis2 = payload.projection_axes
    caption = isnothing(subtitle) ?
        "bond=$(payload.bond_axis), plane: $(payload.plane_axis)=$(payload.plane_value), selected=$(payload.selected_count)/$(payload.total_count)" :
        String(subtitle)

    figure = Figure(size = (1200, 900), fontsize = 24)
    Label(figure[0, 1], title; fontsize = 28, tellwidth = false)
    Label(
        figure[1, 1],
        "★ $(_nucleus_legend_label(payload))";
        fontsize = Int(round(_group_style(:nucleus).markersize)),
        tellwidth = false,
    )
    axis = Axis(
        figure[2, 1];
        xlabel = String(axis1),
        ylabel = String(axis2),
        aspect = DataAspect(),
        title = caption,
        xgridvisible = false,
        ygridvisible = false,
        xminorgridvisible = false,
        yminorgridvisible = false,
    )

    reference_grid = _reference_grid_coordinates(input_path, payload)
    if !isnothing(reference_grid)
        grid_xs, grid_ys = reference_grid
        vlines!(axis, grid_xs; color = RGBAf(0.0, 0.0, 0.0, 0.12), linewidth = 1.0)
        hlines!(axis, grid_ys; color = RGBAf(0.0, 0.0, 0.0, 0.12), linewidth = 1.0)
    end

    for dataset in payload.datasets
        _skip_dataset(dataset) && continue
        style = _group_style(dataset.group_kind)
        if dataset.role == :box_outline
            xs = first.(dataset.points)
            ys = last.(dataset.points)
            lines!(
                axis,
                xs,
                ys;
                color = style.color,
                linewidth = style.linewidth,
                linestyle = style.linestyle,
            )
        else
            xs = first.(dataset.points)
            ys = last.(dataset.points)
            scatter!(
                axis,
                xs,
                ys;
                color = style.color,
                marker = style.marker,
                markersize = style.markersize,
                strokewidth = style.strokewidth,
            )
        end
    end

    mkpath(dirname(String(output_path)))
    save(output_path, figure)
    return payload
end

function _usage()
    println("usage: viz/showpoints2d.jl input.dat output.png [title]")
end

function main_render_bond_aligned_diatomic_projection(args = ARGS)
    args = copy(args)
    if !(2 <= length(args) <= 3)
        _usage()
        exit(1)
    end
    input_path = args[1]
    output_path = args[2]
    title = length(args) == 3 ? args[3] : DEFAULT_TITLE
    payload = render_bond_aligned_diatomic_projection(
        input_path,
        output_path;
        title = title,
    )
    println("input_path=", input_path)
    println("output_path=", output_path)
    println("bond_axis=", payload.bond_axis)
    println("plane_axis=", payload.plane_axis)
    println("plane_value=", payload.plane_value)
    println("plane_tol=", payload.plane_tol)
    println("selected_count=", payload.selected_count)
    println("total_count=", payload.total_count)
    println("projection_axes=", payload.projection_axes[1], ",", payload.projection_axes[2])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_render_bond_aligned_diatomic_projection(ARGS)
end
