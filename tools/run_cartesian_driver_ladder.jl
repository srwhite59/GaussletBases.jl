# Minimal glass-box Cartesian driver ladder for the demolition branch.
#
# This is intentionally not a Test.jl suite. It runs selected driver inputs in
# order, one Julia process per case, and stops at the first failing driver run.

using Dates

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const CASES = [
    (;
        name = "he_wl_q5_pure_gausslet_h1",
        input = "test/driver_inputs/he_wl_q5_pure_gausslet_h1.jl",
    ),
    (;
        name = "h2_pqs_q5_independent_source_box_r4",
        input = "test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl",
    ),
    (;
        name = "h2_pqs_q5_independent_source_box_r4_final_basis",
        input = "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl",
    ),
    (;
        name = "h2_pqs_q5_independent_source_box_r4_h1",
        input = "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl",
    ),
    (;
        name = "h2_pqs_q5_independent_source_box_r4_h1_j",
        input = "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl",
    ),
    (;
        name = "h2_pqs_q5_independent_source_box_r4_private_rhf",
        input = "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl",
    ),
    (;
        name = "h2_pqs_q5_independent_source_box_r4_supplement_preflight",
        input = "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl",
    ),
]

const DRIVER_ARGS = [
    "save_artifact=false",
    "save_tsv=false",
]

function _driver_command(input)
    expression = """
        empty!(ARGS)
        append!(ARGS, $(repr(vcat([input], DRIVER_ARGS))))
        include("bin/cartesian_ham_builder.jl")
        """
    return Cmd(vcat(Base.julia_cmd().exec, ["--project=$(REPO_ROOT)", "-e", expression]))
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

function _print_command(cmd)
    show(stdout, cmd)
    println()
end

function run_ladder(; dry_run::Bool = false)
    println("Cartesian driver ladder")
    println("started = ", Dates.now())
    println("repo = ", REPO_ROOT)
    println()

    passed = String[]
    for (index, case) in enumerate(CASES)
        cmd = _driver_command(case.input)
        println("case $index / $(length(CASES)): ", case.name)
        println("input: ", case.input)
        print("command: ")
        _print_command(cmd)
        if dry_run
            println("result: dry-run")
            println()
            continue
        end

        output_path = tempname()
        elapsed = @elapsed process = open(output_path, "w") do io
            run(pipeline(ignorestatus(cmd); stdout = io, stderr = io))
        end
        output = read(output_path, String)
        rm(output_path; force = true)
        if process.exitcode == 0
            push!(passed, case.name)
            println("elapsed_s: ", elapsed)
            println("result: pass")
            println()
            continue
        end

        println("elapsed_s: ", elapsed)
        println("result: fail")
        println("first_failure_message: ", _first_failure_message(output))
        println()
        println("stopped_after_passed_cases: ", isempty(passed) ? "(none)" : join(passed, ", "))
        exit(process.exitcode)
    end

    dry_run && return nothing
    println("all_passed: ", join(passed, ", "))
    return nothing
end

run_ladder(dry_run = "--dry-run" in ARGS)
