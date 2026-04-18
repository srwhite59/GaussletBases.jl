function Base.show(io::IO, timed::TimedNestedFixedBlockBuild)
    print(
        io,
        "TimedNestedFixedBlockBuild(nfixed=",
        size(timed.fixed_block.overlap, 1),
        ", total=",
        _nested_timeg_report_seconds(timed.timings, "fixed_block.total"),
        "s)",
    )
end

function _nested_timeg_report_seconds(
    report::TimeG.TimingReport,
    label::AbstractString,
)
    return sum(_nested_timeg_report_seconds(node, label) for node in report.roots)
end

function _nested_timeg_report_seconds(
    node::TimeG.TimingNode,
    label::AbstractString,
)
    total = node.label == label ? node.elapsed_seconds : 0.0
    total += sum(_nested_timeg_report_seconds(child, label) for child in node.children)
    return total
end

function _nested_timing_requested(timing)
    return timing === true || timing === :report
end

function _nested_capture_timeg_report(
    build,
    timing::Union{Bool,Symbol},
    timing_io::IO,
)
    timing === false && return build()
    _nested_timing_requested(timing) || throw(
        ArgumentError("one-center atomic fixed-block timing must be false, true, or :report"),
    )
    old_config = TimeG._TIMING_CONFIG[]
    state = TimeG._timing_state()
    parent = isempty(state.stack) ? nothing : state.stack[end]
    root_count = length(state.roots)
    child_count = isnothing(parent) ? 0 : length(parent.children)
    try
        set_timing!(true)
        set_timing_live!(false)
        set_timing_thresholds!(expand = 0.0, drop = 0.0)
        fixed_block = @timeg "fixed_block.total" begin
            build()
        end
        node = isnothing(parent) ? state.roots[root_count + 1] : parent.children[child_count + 1]
        timings = TimeG.TimingReport(TimeG.TimingNode[deepcopy(node)])
        timing === :report && nested_fixed_block_timing_report(timing_io, timings)
        return TimedNestedFixedBlockBuild(fixed_block, timings)
    finally
        TimeG._TIMING_CONFIG[] = old_config
    end
end

function nested_fixed_block_timing_report(
    io::IO,
    timings::TimeG.TimingReport,
)
    old_config = TimeG._TIMING_CONFIG[]
    try
        set_timing_thresholds!(expand = 0.0, drop = 0.0)
        return timing_report(io, timings)
    finally
        TimeG._TIMING_CONFIG[] = old_config
    end
end

function nested_fixed_block_timing_report(
    timings::TimeG.TimingReport,
)
    return sprint(io -> nested_fixed_block_timing_report(io, timings))
end

function nested_fixed_block_timing_report(
    io::IO,
    timed::TimedNestedFixedBlockBuild,
)
    return nested_fixed_block_timing_report(io, timed.timings)
end

function nested_fixed_block_timing_report(
    timed::TimedNestedFixedBlockBuild,
)
    return nested_fixed_block_timing_report(timed.timings)
end

"""
    _CartesianNestedDoSideTrace1D

Structured diagnostic record for one local 1D `doside` / COMX contraction
used in the current nested Cartesian source language.

This keeps enough information to make the local center-loss question explicit:

- where the contraction was used
- which parent interval it came from
- the physical centers before contraction
- the localized centers and signed weights after COMX cleanup
- whether the parent interval is symmetric about zero
- whether the localized centers retain a near-zero center
"""
struct _CartesianNestedDoSideTrace1D
    context_label::String
    group_kind::Symbol
    layer_index::Int
    piece_kind::Symbol
    axis::Symbol
    usage_label::String
    interval::UnitRange{Int}
    parent_centers::Vector{Float64}
    retained_count::Int
    localized_centers::Vector{Float64}
    localized_weights::Vector{Float64}
    symmetric_about_zero::Bool
    symmetry_error::Float64
    contains_near_zero_center::Bool
    even_retained_count::Bool
end

function Base.show(io::IO, trace::_CartesianNestedDoSideTrace1D)
    print(
        io,
        "_CartesianNestedDoSideTrace1D(context=",
        trace.context_label,
        ", axis=:",
        trace.axis,
        ", retained=",
        trace.retained_count,
        ", near_zero=",
        trace.contains_near_zero_center,
        ")",
    )
end

function _nested_contains_near_zero(
    values::AbstractVector{<:Real};
    tol::Float64 = 1.0e-8,
)
    return any(abs(Float64(value)) <= tol for value in values)
end

function _nested_doside_trace(
    side::_CartesianNestedDoSide1D;
    context_label::AbstractString,
    group_kind::Symbol,
    layer_index::Integer,
    piece_kind::Symbol,
    axis::Symbol,
    usage_label::AbstractString,
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    symmetric_about_zero, symmetry_error = _nested_is_symmetric_about_zero(
        side.local_centers;
        tol = symmetry_tol,
    )
    contains_near_zero_center = _nested_contains_near_zero(
        side.localized_centers;
        tol = zero_tol,
    )
    return _CartesianNestedDoSideTrace1D(
        String(context_label),
        group_kind,
        Int(layer_index),
        piece_kind,
        axis,
        String(usage_label),
        side.interval,
        copy(side.local_centers),
        side.retained_count,
        copy(side.localized_centers),
        copy(side.localized_weights),
        symmetric_about_zero,
        symmetry_error,
        contains_near_zero_center,
        iseven(side.retained_count),
    )
end

function _nested_first_matching_face(
    shell::_CartesianNestedCompleteShell3D,
    face_kind::Symbol,
    fixed_side::Symbol,
)
    index = findfirst(face -> face.face_kind == face_kind && face.fixed_side == fixed_side, shell.faces)
    isnothing(index) && throw(ArgumentError("nested doside trace requires a $(face_kind) $(fixed_side) face"))
    return shell.faces[index]
end

function _nested_first_matching_edge(
    shell::_CartesianNestedCompleteShell3D,
    free_axis::Symbol,
)
    index = findfirst(edge -> edge.free_axis == free_axis, shell.edges)
    isnothing(index) && throw(ArgumentError("nested doside trace requires a $(free_axis)-edge representative"))
    return shell.edges[index]
end

function _nested_complete_shell_doside_traces(
    shell::_CartesianNestedCompleteShell3D,
    context_prefix::AbstractString,
    group_kind::Symbol,
    layer_index::Integer;
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    traces = _CartesianNestedDoSideTrace1D[]

    face_xy = _nested_first_matching_face(shell, :xy, :low)
    push!(
        traces,
        _nested_doside_trace(
            face_xy.side_first;
            context_label = string(context_prefix, "/face_xy/tangential_x"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :x,
            usage_label = "face_kind=:xy shared_by=low/high tangential_axis=:x",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    push!(
        traces,
        _nested_doside_trace(
            face_xy.side_second;
            context_label = string(context_prefix, "/face_xy/tangential_y"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :y,
            usage_label = "face_kind=:xy shared_by=low/high tangential_axis=:y",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )

    face_xz = _nested_first_matching_face(shell, :xz, :low)
    push!(
        traces,
        _nested_doside_trace(
            face_xz.side_first;
            context_label = string(context_prefix, "/face_xz/tangential_x"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :x,
            usage_label = "face_kind=:xz shared_by=low/high tangential_axis=:x",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    push!(
        traces,
        _nested_doside_trace(
            face_xz.side_second;
            context_label = string(context_prefix, "/face_xz/tangential_z"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :z,
            usage_label = "face_kind=:xz shared_by=low/high tangential_axis=:z",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )

    face_yz = _nested_first_matching_face(shell, :yz, :low)
    push!(
        traces,
        _nested_doside_trace(
            face_yz.side_first;
            context_label = string(context_prefix, "/face_yz/tangential_y"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :y,
            usage_label = "face_kind=:yz shared_by=low/high tangential_axis=:y",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    push!(
        traces,
        _nested_doside_trace(
            face_yz.side_second;
            context_label = string(context_prefix, "/face_yz/tangential_z"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :z,
            usage_label = "face_kind=:yz shared_by=low/high tangential_axis=:z",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )

    edge_x = _nested_first_matching_edge(shell, :x)
    push!(
        traces,
        _nested_doside_trace(
            edge_x.side;
            context_label = string(context_prefix, "/edge_x/free_axis_x"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :edge_free,
            axis = :x,
            usage_label = "free_axis=:x shared_by=all_boundary_sign_pairs",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    edge_y = _nested_first_matching_edge(shell, :y)
    push!(
        traces,
        _nested_doside_trace(
            edge_y.side;
            context_label = string(context_prefix, "/edge_y/free_axis_y"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :edge_free,
            axis = :y,
            usage_label = "free_axis=:y shared_by=all_boundary_sign_pairs",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    edge_z = _nested_first_matching_edge(shell, :z)
    push!(
        traces,
        _nested_doside_trace(
            edge_z.side;
            context_label = string(context_prefix, "/edge_z/free_axis_z"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :edge_free,
            axis = :z,
            usage_label = "free_axis=:z shared_by=all_boundary_sign_pairs",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )

    return traces
end

function _nested_sequence_doside_traces(
    sequence::_CartesianNestedShellSequence3D,
    region_label::AbstractString,
    group_kind::Symbol;
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    traces = _CartesianNestedDoSideTrace1D[]
    for (layer_index, layer) in pairs(sequence.shell_layers)
        layer isa _CartesianNestedCompleteShell3D || continue
        append!(
            traces,
            _nested_complete_shell_doside_traces(
                layer,
                string(region_label, "/layer_", layer_index),
                group_kind,
                layer_index;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
    end
    return traces
end

function _bond_aligned_diatomic_doside_traces(
    source;
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    traces = _CartesianNestedDoSideTrace1D[]
    for (layer_index, layer) in pairs(source.shared_shell_layers)
        append!(
            traces,
            _nested_complete_shell_doside_traces(
                layer,
                string("shared_shell/layer_", layer_index),
                :shared_shell,
                layer_index;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
    end
    if source.geometry.did_split
        append!(
            traces,
            _nested_sequence_doside_traces(
                source.child_sequences[1],
                "left_child",
                :left_child;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
        append!(
            traces,
            _nested_sequence_doside_traces(
                source.child_sequences[2],
                "right_child",
                :right_child;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
    else
        append!(
            traces,
            _nested_sequence_doside_traces(
                source.child_sequences[1],
                "shared_child",
                :shared_child;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
    end
    return traces
end

function _nested_trace_vector_string(values::AbstractVector{<:Real})
    return "[" * join((string(Float64(value)) for value in values), ", ") * "]"
end

function _bond_aligned_diatomic_doside_trace_notes(source)
    notes = String[]
    isempty(source.shared_shell_layers) && push!(
        notes,
        "# note shared_shell has no local side contractions",
    )
    if source.geometry.did_split
        isempty(source.child_sequences[1].shell_layers) && push!(
            notes,
            "# note left_child has no local side contractions; it remains a direct core block",
        )
        isempty(source.child_sequences[2].shell_layers) && push!(
            notes,
            "# note right_child has no local side contractions; it remains a direct core block",
        )
    else
        isempty(source.child_sequences[1].shell_layers) && push!(
            notes,
            "# note shared_child has no local side contractions; it remains a direct core block",
        )
    end
    return notes
end

function _write_bond_aligned_diatomic_doside_trace(
    path::AbstractString,
    source;
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    traces = _bond_aligned_diatomic_doside_traces(
        source;
        symmetry_tol = symmetry_tol,
        zero_tol = zero_tol,
    )
    mkpath(dirname(String(path)))
    open(path, "w") do io
        write(io, "# GaussletBases bond-aligned diatomic doside/COMX trace\n")
        write(io, "# bond_axis = $(source.basis.bond_axis)\n")
        write(io, "# working_box = $(source.geometry.working_box)\n")
        write(io, "# did_split = $(source.geometry.did_split)\n")
        if !isnothing(source.geometry.shared_midpoint_box)
            write(io, "# shared_midpoint_box = $(source.geometry.shared_midpoint_box)\n")
        end
        write(io, "# symmetry_tol = $(symmetry_tol)\n")
        write(io, "# zero_tol = $(zero_tol)\n")
        write(io, "# trace_count = $(length(traces))\n")
        for note in _bond_aligned_diatomic_doside_trace_notes(source)
            write(io, note, "\n")
        end
        for (index, trace) in pairs(traces)
            write(io, "\n[trace $(index)]\n")
            write(io, "context_label = $(trace.context_label)\n")
            write(io, "group_kind = $(trace.group_kind)\n")
            write(io, "layer_index = $(trace.layer_index)\n")
            write(io, "piece_kind = $(trace.piece_kind)\n")
            write(io, "axis = $(trace.axis)\n")
            write(io, "usage = $(trace.usage_label)\n")
            write(io, "interval = $(first(trace.interval)):$(last(trace.interval))\n")
            write(io, "parent_centers = $(_nested_trace_vector_string(trace.parent_centers))\n")
            write(io, "retained_count = $(trace.retained_count)\n")
            write(io, "localized_centers = $(_nested_trace_vector_string(trace.localized_centers))\n")
            write(io, "localized_weights = $(_nested_trace_vector_string(trace.localized_weights))\n")
            write(io, "symmetric_about_zero = $(trace.symmetric_about_zero)\n")
            write(io, "symmetry_error = $(trace.symmetry_error)\n")
            write(io, "contains_near_zero_center = $(trace.contains_near_zero_center)\n")
            write(io, "even_retained_count = $(trace.even_retained_count)\n")
        end
    end
    return traces
end
