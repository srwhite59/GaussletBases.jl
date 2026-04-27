using GaussletBases
using LinearAlgebra
using Printf

const SELECTED_TIMING_LABELS = [
    "diatomic.fixed_source.total",
    "diatomic.fixed_source.axis_bundles",
    "diatomic.fixed_source.source_assembly",
    "diatomic.source.total",
    "diatomic.source.shared_shell_construction",
    "diatomic.source.child_sequence_builds",
    "diatomic.source.final_sequence_merge",
    "diatomic.sequence.packet",
    "diatomic.packet.total",
    "diatomic.packet.setup",
    "diatomic.packet.base.overlap",
    "diatomic.packet.base.kinetic",
    "diatomic.packet.base.position_x",
    "diatomic.packet.base.position_y",
    "diatomic.packet.base.position_z",
    "diatomic.packet.base.x2_x",
    "diatomic.packet.base.x2_y",
    "diatomic.packet.base.x2_z",
    "diatomic.packet.pair_terms",
    "diatomic.packet.gaussian_terms",
    "diatomic.fixed_block.contraction",
]

mutable struct TimingStats
    elapsed_seconds::Float64
    self_seconds::Float64
    call_count::Int
    max_elapsed_seconds::Float64
end

TimingStats() = TimingStats(0.0, 0.0, 0, 0.0)

function _parse_bool(value::AbstractString)
    lowered = lowercase(value)
    lowered in ("1", "true", "yes", "on") && return true
    lowered in ("0", "false", "no", "off") && return false
    throw(ArgumentError("expected boolean value, got \"$value\""))
end

function _parse_args(args::Vector{String})
    values = Dict{String,String}()
    for arg in args
        occursin("=", arg) ||
            throw(ArgumentError("arguments must use key=value form, got \"$arg\""))
        key, value = split(arg, "="; limit = 2)
        values[key] = value
    end
    return values
end

function _add_timing!(stats::Dict{String,TimingStats}, node)
    current = get!(stats, node.label, TimingStats())
    current.elapsed_seconds += node.elapsed_seconds
    current.self_seconds += node.self_seconds
    current.call_count += node.call_count
    current.max_elapsed_seconds = max(current.max_elapsed_seconds, node.elapsed_seconds)
    for child in node.children
        _add_timing!(stats, child)
    end
    return stats
end

function _timing_stats(report::GaussletBases.TimeG.TimingReport)
    stats = Dict{String,TimingStats}()
    for root in report.roots
        _add_timing!(stats, root)
    end
    return stats
end

function _print_selected_timings(stats::Dict{String,TimingStats})
    println("selected_timings_begin")
    println("label,calls,elapsed_s,self_s,max_call_elapsed_s,status")
    for label in SELECTED_TIMING_LABELS
        if haskey(stats, label)
            timing = stats[label]
            @printf(
                "%s,%d,%.9f,%.9f,%.9f,ok\n",
                label,
                timing.call_count,
                timing.elapsed_seconds,
                timing.self_seconds,
                timing.max_elapsed_seconds,
            )
        else
            println(label, ",0,NaN,NaN,NaN,missing")
        end
    end
    println("selected_timings_end")
    return nothing
end

function _run_kernel(
    basis,
    expansion::CoulombGaussianExpansion,
    packet_kernel::Symbol;
    nside::Int,
)
    GC.gc()
    reset_timing_report!()
    set_timing!(true)
    set_timing_live!(false)
    set_timing_thresholds!(expand = 0.0, drop = 0.0)

    result = nothing
    wall_seconds = @elapsed begin
        result = bond_aligned_diatomic_nested_fixed_block(
            basis;
            expansion = expansion,
            nside = nside,
            packet_kernel = packet_kernel,
        )
    end
    report = current_timing_report()
    set_timing!(false)

    return (
        packet_kernel = packet_kernel,
        wall_seconds = wall_seconds,
        result = result,
        report = report,
        stats = _timing_stats(report),
    )
end

function _warmup_kernel(
    basis,
    expansion::CoulombGaussianExpansion,
    packet_kernel::Symbol;
    nside::Int,
)
    set_timing!(false)
    set_timing_live!(false)
    bond_aligned_diatomic_nested_fixed_block(
        basis;
        expansion = expansion,
        nside = nside,
        packet_kernel = packet_kernel,
    )
    GC.gc()
    return nothing
end

function _print_run_summary(run)
    source = run.result.source
    fixed_block = run.result.fixed_block
    println("run_begin")
    println("packet_kernel=", run.packet_kernel)
    @printf("wall_s=%.9f\n", run.wall_seconds)
    println("fixed_dimension=", size(fixed_block.coefficient_matrix, 2))
    println("source_support_count=", length(source.sequence.support_indices))
    println("source_working_box=", source.sequence.working_box)
    println("did_split=", source.geometry.did_split)
    println("shared_shell_count=", length(source.shared_shell_layers))
    println("child_sequence_count=", length(source.child_sequences))
    _print_selected_timings(run.stats)
    println("run_end")
    return nothing
end

function main(args::Vector{String})
    values = _parse_args(args)

    timing_context = get(values, "timing_context", "fresh_process")
    timing_context in ("fresh_process", "warm_persistent_session") || throw(
        ArgumentError(
            "timing_context must be fresh_process or warm_persistent_session, got \"$timing_context\"",
        ),
    )
    warmup = _parse_bool(get(values, "warmup", "true"))
    nside = parse(Int, get(values, "nside", "5"))
    bond_length = parse(Float64, get(values, "bond_length", "1.4"))
    core_spacing = parse(Float64, get(values, "core_spacing", "0.5"))
    xmax_parallel = parse(Float64, get(values, "xmax_parallel", "6.0"))
    xmax_transverse = parse(Float64, get(values, "xmax_transverse", "4.0"))
    bond_axis = Symbol(get(values, "bond_axis", "z"))

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = bond_length,
        core_spacing = core_spacing,
        xmax_parallel = xmax_parallel,
        xmax_transverse = xmax_transverse,
        bond_axis = bond_axis,
    )

    println("diatomic_packet_kernel_benchmark")
    println("timing_context=", timing_context)
    println("warmup=", warmup)
    println("case_begin")
    println("bond_length=", bond_length)
    println("core_spacing=", core_spacing)
    println("xmax_parallel=", xmax_parallel)
    println("xmax_transverse=", xmax_transverse)
    println("bond_axis=", bond_axis)
    println("nside=", nside)
    println("term_count=", length(expansion.coefficients))
    println("case_end")

    kernels = (:support_reference, :factorized_direct)
    if warmup
        for kernel in kernels
            _warmup_kernel(basis, expansion, kernel; nside = nside)
        end
    end

    runs = [_run_kernel(basis, expansion, kernel; nside = nside) for kernel in kernels]
    for run in runs
        _print_run_summary(run)
    end

    support = runs[1].result.fixed_block
    direct = runs[2].result.fixed_block
    println("comparison_begin")
    println("fixed_dimension_match=", size(support.coefficient_matrix, 2) == size(direct.coefficient_matrix, 2))
    @printf(
        "coefficient_matrix_inf_residual=%.9e\n",
        norm(support.coefficient_matrix - direct.coefficient_matrix, Inf),
    )
    @printf("overlap_inf_residual=%.9e\n", norm(support.overlap - direct.overlap, Inf))
    @printf("kinetic_inf_residual=%.9e\n", norm(support.kinetic - direct.kinetic, Inf))
    @printf("gaussian_sum_inf_residual=%.9e\n", norm(support.gaussian_sum - direct.gaussian_sum, Inf))
    @printf("pair_sum_inf_residual=%.9e\n", norm(support.pair_sum - direct.pair_sum, Inf))
    println("comparison_end")

    return nothing
end

main(ARGS)
