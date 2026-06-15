include("cartesian_driver_ladder_lib.jl")

options = cartesian_ladder_options()

if options.list_lines
    print_cartesian_ladder_lines()
    exit(0)
end

isnothing(options.line) && begin
    println("usage: julia --project=. tools/run_cartesian_line_ladder.jl --line=<line>")
    println()
    print_cartesian_ladder_lines()
    exit(2)
end

cases = cartesian_ladder_cases_for_line(options.line)

run_cartesian_ladder(
    cases;
    title = "Cartesian temporary line ladder: $(options.line)",
    stop_after_stage = options.stop_after_stage,
    dry_run = options.dry_run,
)
