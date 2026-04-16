"""
    TimeG

Internal timing utilities for lightweight scope-level instrumentation inside
`GaussletBases`.

`TimeG` is for:
- hierarchical wall-clock timing of selected code regions
- low-friction live timing during development
- capturing compact timing trees for later rendering

`TimeG` is not for:
- full statistical profiling
- allocation or memory attribution
- persistent benchmarking infrastructure

Main public entry points:
- `@timeg`
- `timing_enabled`, `timing_live_enabled`
- `set_timing!`, `set_timing_live!`, `set_timing_thresholds!`
- `reset_timing_report!`, `current_timing_report`, `timing_report`

Example:
```julia
julia> using GaussletBases

julia> reset_timing_report!(); set_timing!(true); set_timing_live!(false);

julia> @timeg "demo" begin
           sleep(0.01)
       end

julia> println(timing_report())
GaussletBases timing report
demo: elapsed=...
```
"""
module TimeG

export @timeg,
       timing_enabled,
       timing_live_enabled,
       set_timing!,
       set_timing_live!,
       set_timing_thresholds!,
       reset_timing_report!,
       current_timing_report,
       timing_report

"""
    TimingConfig

Internal global timing configuration shared by `@timeg`.

Fields:
- `enabled`: whether timed scopes are recorded at all
- `live_print`: whether completed scopes are printed immediately
- `expand_threshold_seconds`: minimum elapsed time for recursive report expansion
- `drop_threshold_seconds`: minimum elapsed time shown in reports or live output
"""
struct TimingConfig
    enabled::Bool
    live_print::Bool
    expand_threshold_seconds::Float64
    drop_threshold_seconds::Float64
end

"""
    TimingNode

One node in a timing report tree.

This is the public-facing report record returned through `TimingReport.roots`.
Its fields are stable for inspection:
- `label`: scope label passed to `@timeg`
- `elapsed_seconds`: total wall-clock time for the scope
- `self_seconds`: elapsed time excluding recorded children
- `call_count`: merged call count in rendered reports
- `children`: nested timing nodes
"""
mutable struct TimingNode
    label::String
    elapsed_seconds::Float64
    self_seconds::Float64
    call_count::Int
    children::Vector{TimingNode}
end

"""
    _TimingFrame

Internal active-scope frame held on the task-local timing stack while a timed
scope is executing.
"""
mutable struct _TimingFrame
    label::String
    start_ns::UInt64
    children::Vector{TimingNode}
    child_elapsed_seconds::Float64
end

"""
    _TimingState

Internal task-local timing state.

Each Julia task gets an independent timing stack and root-node collection so
timed scopes do not leak across unrelated task execution.
"""
mutable struct _TimingState
    stack::Vector{_TimingFrame}
    roots::Vector{TimingNode}
end

"""
    TimingReport

Immutable snapshot of the current task's collected timing tree.

`roots` contains top-level `TimingNode`s in execution order.
"""
struct TimingReport
    roots::Vector{TimingNode}
end

function Base.show(io::IO, report::TimingReport)
    print(io, "TimingReport(nroots=", length(report.roots), ")")
end

# Internal process-global timing configuration.
const _TIMING_CONFIG = Ref(TimingConfig(true, true, 10.0, 0.1))

function _timing_state()
    storage = task_local_storage()
    if !haskey(storage, :gaussletbases_timing_state)
        state = _TimingState(_TimingFrame[], TimingNode[])
        task_local_storage(:gaussletbases_timing_state, state)
        return state
    end
    return storage[:gaussletbases_timing_state]
end

"""
    timing_enabled() -> Bool

Return whether `@timeg` currently records timed scopes.
"""
timing_enabled() = _TIMING_CONFIG[].enabled

"""
    timing_live_enabled() -> Bool

Return whether completed timed scopes are printed immediately as they finish.
"""
timing_live_enabled() = _TIMING_CONFIG[].live_print

"""
    set_timing!(enabled::Bool) -> Bool

Enable or disable `@timeg` collection globally for the current process.
"""
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

"""
    set_timing_live!(enabled::Bool) -> Bool

Enable or disable live printing for completed timed scopes.
"""
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

"""
    set_timing_thresholds!(; expand = ..., drop = ...) -> TimingConfig

Set the elapsed-time thresholds used by `timing_report` and live printing.

- `expand`: minimum elapsed time required before child nodes are recursively
  shown in the rendered report
- `drop`: minimum elapsed time required for a node to appear at all
"""
function set_timing_thresholds!(;
    expand::Real = _TIMING_CONFIG[].expand_threshold_seconds,
    drop::Real = _TIMING_CONFIG[].drop_threshold_seconds,
)
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

"""
    reset_timing_report!()

Clear the current task's accumulated timing tree and active timing stack.
"""
function reset_timing_report!()
    state = _timing_state()
    empty!(state.stack)
    empty!(state.roots)
    return nothing
end

"""
    current_timing_report() -> TimingReport

Return a snapshot of the current task's recorded timing tree.

This requires all timed scopes to be closed; otherwise an `ArgumentError` is
thrown.
"""
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

"""
    @timeg label expr

Evaluate `expr` while recording a named timing scope.

When timing is disabled, `expr` is evaluated directly with no report entry.
When timing is enabled, the scope is added to the current task's timing tree
and may also be printed immediately if live timing is enabled.
"""
macro timeg(label, expr)
    time_scope = GlobalRef(@__MODULE__, :_time_scope)
    return esc(:($time_scope(() -> $(expr), $(label))))
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

"""
    timing_report(io::IO, report::TimingReport = current_timing_report())
    timing_report(report::TimingReport = current_timing_report()) -> String

Render a human-readable timing report for the current task or a supplied
`TimingReport`.
"""
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

end
