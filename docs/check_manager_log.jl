module ManagerLogPolicy

const LIVE_LOG_MAX_LINES = 2_000
const LIVE_LOG_PATH = joinpath(@__DIR__, "src", "developer", "pqs_manager_running_log.md")

function file_line_count(path::AbstractString)
    isfile(path) || throw(ArgumentError("manager log does not exist: $(path)"))
    return open(path, "r") do io
        count(_ -> true, eachline(io))
    end
end

function check_live_log(
    path::AbstractString = LIVE_LOG_PATH;
    max_lines::Integer = LIVE_LOG_MAX_LINES,
)
    max_lines > 0 || throw(ArgumentError("manager-log line limit must be positive"))
    lines = file_line_count(path)
    lines <= max_lines || error(
        "live manager log has $(lines) lines; rotate it before the $(max_lines)-line limit",
    )
    return lines
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) || ARGS == ["--check"] ||
        error("usage: julia docs/check_manager_log.jl [--check]")
    lines = ManagerLogPolicy.check_live_log()
    println("manager_log_lines=$(lines) limit=$(ManagerLogPolicy.LIVE_LOG_MAX_LINES)")
end
