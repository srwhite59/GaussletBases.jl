struct TimingConfig
    enabled::Bool
    live_print::Bool
    expand_threshold_seconds::Float64
    drop_threshold_seconds::Float64
end

mutable struct TimingNode
    label::String
    elapsed_seconds::Float64
    self_seconds::Float64
    call_count::Int
    children::Vector{TimingNode}
end

mutable struct _TimingFrame
    label::String
    start_ns::UInt64
    children::Vector{TimingNode}
    child_elapsed_seconds::Float64
end

struct TimingReport
    roots::Vector{TimingNode}
end

function Base.show(io::IO, report::TimingReport)
    print(io, "TimingReport(nroots=", length(report.roots), ")")
end

const _TIMING_CONFIG = Ref(TimingConfig(true, true, 10.0, 0.1))

function _timing_state()
    storage = task_local_storage()
    if !haskey(storage, :gaussletbases_timing_state)
        state = (
            stack = _TimingFrame[],
            roots = TimingNode[],
        )
        task_local_storage(:gaussletbases_timing_state, state)
        return state
    end
    return storage[:gaussletbases_timing_state]
end

timing_enabled() = _TIMING_CONFIG[].enabled
timing_live_enabled() = _TIMING_CONFIG[].live_print

function set_timing!(enabled::Bool)
    config = _TIMING_CONFIG[]
    _TIMING_CONFIG[] = TimingConfig(
        enabled,
        config.live_print,
        config.expand_threshold_seconds,
        config.drop_threshold_seconds,
    )
    return enabled
end

function set_timing_live!(enabled::Bool)
    config = _TIMING_CONFIG[]
    _TIMING_CONFIG[] = TimingConfig(
        config.enabled,
        enabled,
        config.expand_threshold_seconds,
        config.drop_threshold_seconds,
    )
    return enabled
end

function set_timing_thresholds!(; expand::Real = _TIMING_CONFIG[].expand_threshold_seconds, drop::Real = _TIMING_CONFIG[].drop_threshold_seconds)
    expand_value = Float64(expand)
    drop_value = Float64(drop)
    expand_value >= 0.0 || throw(ArgumentError("timing expand threshold must be nonnegative"))
    drop_value >= 0.0 || throw(ArgumentError("timing drop threshold must be nonnegative"))
    _TIMING_CONFIG[] = TimingConfig(
        _TIMING_CONFIG[].enabled,
        _TIMING_CONFIG[].live_print,
        expand_value,
        drop_value,
    )
    return _TIMING_CONFIG[]
end

function reset_timing_report!()
    state = _timing_state()
    empty!(state.stack)
    empty!(state.roots)
    return nothing
end

function current_timing_report()
    state = _timing_state()
    isempty(state.stack) || throw(
        ArgumentError("cannot snapshot timing report while timed scopes are still active"),
    )
    return TimingReport(deepcopy(state.roots))
end

function _timing_attach_node!(node::TimingNode)
    state = _timing_state()
    if isempty(state.stack)
        push!(state.roots, node)
    else
        parent = state.stack[end]
        push!(parent.children, node)
        parent.child_elapsed_seconds += node.elapsed_seconds
    end
    return nothing
end

function _time_scope(f, label)
    !timing_enabled() && return f()

    state = _timing_state()
    frame = _TimingFrame(String(label), time_ns(), TimingNode[], 0.0)
    push!(state.stack, frame)
    try
        return f()
    finally
        finished = pop!(state.stack)
        elapsed_seconds = (time_ns() - finished.start_ns) / 1.0e9
        node = TimingNode(
            finished.label,
            elapsed_seconds,
            max(0.0, elapsed_seconds - finished.child_elapsed_seconds),
            1,
            finished.children,
        )
        _timing_attach_node!(node)
        config = _TIMING_CONFIG[]
        if config.live_print && elapsed_seconds >= config.drop_threshold_seconds
            depth = length(state.stack)
            println(
                "  "^depth,
                finished.label,
                ": ",
                round(elapsed_seconds; sigdigits = 4),
                " seconds",
            )
        end
    end
end

macro timeg(label, expr)
    return esc(:(GaussletBases._time_scope(() -> $(expr), $(label))))
end

function _timing_merge_nodes(nodes::AbstractVector{TimingNode})
    grouped = Dict{String,TimingNode}()
    order = String[]
    for node in nodes
        if !haskey(grouped, node.label)
            grouped[node.label] = TimingNode(node.label, 0.0, 0.0, 0, TimingNode[])
            push!(order, node.label)
        end
        merged = grouped[node.label]
        merged.elapsed_seconds += node.elapsed_seconds
        merged.self_seconds += node.self_seconds
        merged.call_count += node.call_count
        append!(merged.children, node.children)
    end
    return [grouped[label] for label in order]
end

function _timing_report_lines!(
    lines::Vector{String},
    nodes::AbstractVector{TimingNode},
    config::TimingConfig,
    indent::Int,
)
    merged_nodes = _timing_merge_nodes(nodes)
    for node in merged_nodes
        node.elapsed_seconds >= config.drop_threshold_seconds || continue
        push!(
            lines,
            string(
                "  "^indent,
                node.label,
                ": elapsed=",
                node.elapsed_seconds,
                " s, self=",
                node.self_seconds,
                " s, calls=",
                node.call_count,
            ),
        )
        if node.elapsed_seconds >= config.expand_threshold_seconds
            _timing_report_lines!(lines, node.children, config, indent + 1)
        end
    end
    return lines
end

function timing_report(io::IO, report::TimingReport = current_timing_report())
    config = _TIMING_CONFIG[]
    println(io, "GaussletBases timing report")
    _timing_report_lines!(String[], report.roots, config, 0) |> lines -> begin
        for line in lines
            println(io, line)
        end
    end
    return report
end

function timing_report(report::TimingReport = current_timing_report())
    return sprint(io -> timing_report(io, report))
end
