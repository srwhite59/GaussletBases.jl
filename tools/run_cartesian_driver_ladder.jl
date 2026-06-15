include("cartesian_driver_ladder_lib.jl")

options = cartesian_ladder_options()

if options.list_lines
    print_cartesian_ladder_lines()
    exit(0)
end

run_cartesian_ladder(
    CARTESIAN_MATRIX_CASES;
    title = "Cartesian driver 2x2x2 matrix",
    stop_after_stage = options.stop_after_stage,
    dry_run = options.dry_run,
)
