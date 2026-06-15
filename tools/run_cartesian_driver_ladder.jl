# Minimal glass-box Cartesian driver matrix for the demolition branch.
#
# This is intentionally not a Test.jl suite. It runs the 2x2x2 driver matrix
# in separate Julia processes, records the driver stage reached by each case,
# and reports the earliest failing driver stage across the matrix.

using Dates

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const DRIVER_STAGES = [
    :cartesian_system,
    :cartesian_recipe,
    :cartesian_parent,
    :cartesian_shells,
    :cartesian_units,
    :cartesian_transforms,
    :cartesian_pair_terms,
    :cartesian_assembly,
    :cartesian_report,
    :cartesian_materialization,
    Symbol("cartesian_print/save"),
]

const CASES = [
    (;
        name = "he_wl",
        input = "test/driver_inputs/he_wl_q5_pure_gausslet_h1.jl",
    ),
    (;
        name = "h2_wl",
        input = "test/driver_inputs/h2_wl_q5_pure_gausslet_h1.jl",
    ),
    (;
        name = "he_wl_gto",
        input = "test/driver_inputs/he_wl_q5_gto_h1.jl",
    ),
    (;
        name = "h2_wl_gto",
        input = "test/driver_inputs/h2_wl_q5_gto_h1.jl",
    ),
    (;
        name = "he_pqs",
        input = "test/driver_inputs/he_pqs_q5_wlmap.jl",
    ),
    (;
        name = "h2_pqs",
        input = "test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl",
    ),
    (;
        name = "he_pqs_gto",
        input = "test/driver_inputs/he_pqs_q5_gto.jl",
    ),
    (;
        name = "h2_pqs_gto",
        input = "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl",
    ),
]

const DRIVER_ARGS = [
    "save_artifact=false",
    "save_tsv=false",
]

function _argument_value(prefix)
    for argument in ARGS
        startswith(argument, prefix) || continue
        return argument[length(prefix) + 1:end]
    end
    return nothing
end

const STOP_AFTER_STAGE = _argument_value("--stop-after=")

function _driver_command(input)
    driver_args = copy(DRIVER_ARGS)
    if !isnothing(STOP_AFTER_STAGE)
        push!(driver_args, "stop_after_stage=$(repr(Symbol(STOP_AFTER_STAGE)))")
    end
    expression = """
        empty!(ARGS)
        append!(ARGS, $(repr(vcat([input], driver_args))))
        include("bin/cartesian_ham_builder.jl")
        """
    return Cmd(vcat(Base.julia_cmd().exec, ["--project=$(REPO_ROOT)", "-e", expression]))
end

function _print_command(cmd)
    show(stdout, cmd)
    println()
end

function _stage_index(stage)
    isnothing(stage) && return 0
    index = findfirst(==(stage), DRIVER_STAGES)
    return isnothing(index) ? 0 : index
end

function _stage_from_line(line)
    prefix = "[driver stage] "
    startswith(line, prefix) || return nothing
    return Symbol(strip(line[length(prefix) + 1:end]))
end

function _last_driver_stage(output)
    last_stage = nothing
    for line in split(output, '\n')
        stage = _stage_from_line(strip(line))
        isnothing(stage) || (last_stage = stage)
    end
    return last_stage
end

function _first_failure_message(output)
    lines = split(output, '\n')
    for line in lines
        stripped = strip(line)
        startswith(stripped, "ERROR:") && return stripped
    end
    for line in reverse(lines)
        stripped = strip(line)
        isempty(stripped) && continue
        return stripped
    end
    return "(driver failed without output)"
end

function _useful_final_line(output)
    for line in reverse(split(output, '\n'))
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "[driver stage]") && continue
        return stripped
    end
    return "(driver produced no output)"
end

function _run_case(case)
    cmd = _driver_command(case.input)
    output_path = tempname()
    passed = true
    elapsed = @elapsed try
        open(output_path, "w") do io
            run(pipeline(cmd; stdout = io, stderr = io))
        end
    catch
        passed = false
    end
    output = read(output_path, String)
    rm(output_path; force = true)
    stage = _last_driver_stage(output)
    return (;
        case.name,
        case.input,
        elapsed,
        passed,
        stage,
        stage_index = _stage_index(stage),
        note = passed ? _useful_final_line(output) : _first_failure_message(output),
        output,
    )
end

function _summarize(results)
    failures = [result for result in results if !result.passed]
    if isempty(failures)
        println("earliest failing driver stage: none")
        println("all cases passed")
        return nothing
    end

    earliest_index = minimum(result.stage_index for result in failures)
    earliest_stage =
        earliest_index == 0 ? :before_driver_stage_print :
        DRIVER_STAGES[earliest_index]
    failing_there = [result.name for result in failures if result.stage_index == earliest_index]
    passing_that_stage =
        [result.name for result in results if result.stage_index > earliest_index || result.passed]
    later_failures =
        [result.name for result in failures if result.stage_index > earliest_index]

    println("earliest failing driver stage: ", earliest_stage)
    println("cases failing there: ", isempty(failing_there) ? "(none)" : join(failing_there, ", "))
    println(
        "cases passing that stage: ",
        isempty(passing_that_stage) ? "(none)" : join(passing_that_stage, ", "),
    )
    println("later failures: ", isempty(later_failures) ? "(none)" : join(later_failures, ", "))
    return nothing
end

function run_ladder(; dry_run::Bool = false)
    println("Cartesian driver 2x2x2 matrix")
    println("started = ", Dates.now())
    println("repo = ", REPO_ROOT)
    println("stop_after_stage = ", isnothing(STOP_AFTER_STAGE) ? "none" : STOP_AFTER_STAGE)
    println()

    results = NamedTuple[]
    for (index, case) in enumerate(CASES)
        cmd = _driver_command(case.input)
        println("case $index / $(length(CASES)): ", case.name)
        println("input: ", case.input)
        print("command: ")
        _print_command(cmd)
        if dry_run
            println("elapsed_s: 0.0")
            println("result: dry-run")
            println("last_driver_stage: not_run")
            println("note: not_run")
            println()
            continue
        end

        result = _run_case(case)
        push!(results, result)
        println("driver_output_begin")
        print(result.output)
        isempty(result.output) || endswith(result.output, "\n") || println()
        println("driver_output_end")
        println("elapsed_s: ", result.elapsed)
        println("result: ", result.passed ? "pass" : "fail")
        println("last_driver_stage: ", isnothing(result.stage) ? "none" : result.stage)
        println(result.passed ? "useful_final_output: " : "first_error_line: ", result.note)
        println()
    end

    dry_run && return nothing
    _summarize(results)
    any(result -> !result.passed, results) && exit(1)
    return nothing
end

run_ladder(dry_run = "--dry-run" in ARGS)
