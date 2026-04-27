using GaussletBases
using LinearAlgebra
using Printf

const SELECTED_TIMING_LABELS = [
    "diatomic.fixed_source.total",
    "diatomic.fixed_source.axis_bundles",
    "diatomic.fixed_source.source_assembly",
    "diatomic.source.total",
    "diatomic.source.split_geometry.initial",
    "diatomic.source.shared_shell_construction",
    "diatomic.source.adaptive_shell_retention",
    "shell_layer.nonpacket",
    "diatomic.source.child_sequence_builds",
    "diatomic.source.child_sequence",
    "diatomic.source.child_sequence.shell_construction",
    "diatomic.source.child_sequence.core_block",
    "diatomic.source.child_sequence.sequence_merge",
    "diatomic.source.nonuniform_core_block",
    "diatomic.source.nonuniform_core.coefficient_merge",
    "diatomic.source.final_sequence_merge",
    "sequence_merge.nonpacket",
    "diatomic.sequence.coefficient_concat",
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
    println("source_cutset_timings_begin")
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
    println("source_cutset_timings_end")
    return nothing
end

_format_range(range::UnitRange{Int}) = string(first(range), ":", last(range), "(", length(range), ")")

function _format_box(box)
    return string("(", join((_format_range(axis_range) for axis_range in box), ", "), ")")
end

function _format_widths(widths)
    return @sprintf("(%.6g, %.6g, %.6g)", widths[1], widths[2], widths[3])
end

function _format_ranges(ranges)
    return string("[", join((_format_range(range) for range in ranges), ", "), "]")
end

function _run_kernel(
    basis,
    expansion::CoulombGaussianExpansion,
    route::Symbol;
    nside::Int,
)
    route in (:default, :support_reference, :factorized_direct) ||
        throw(ArgumentError("packet-kernel benchmark route must be :default, :support_reference, or :factorized_direct"))
    GC.gc()
    reset_timing_report!()
    set_timing!(true)
    set_timing_live!(false)
    set_timing_thresholds!(expand = 0.0, drop = 0.0)

    result = nothing
    wall_seconds = @elapsed begin
        if route == :default
            result = bond_aligned_diatomic_nested_fixed_block(
                basis;
                expansion = expansion,
                nside = nside,
            )
        else
            result = bond_aligned_diatomic_nested_fixed_block(
                basis;
                expansion = expansion,
                nside = nside,
                packet_kernel = route,
            )
        end
    end
    report = current_timing_report()
    set_timing!(false)

    return (
        route = route,
        wall_seconds = wall_seconds,
        result = result,
        report = report,
        stats = _timing_stats(report),
    )
end

function _warmup_kernel(
    basis,
    expansion::CoulombGaussianExpansion,
    route::Symbol;
    nside::Int,
)
    set_timing!(false)
    set_timing_live!(false)
    if route == :default
        bond_aligned_diatomic_nested_fixed_block(
            basis;
            expansion = expansion,
            nside = nside,
        )
    else
        bond_aligned_diatomic_nested_fixed_block(
            basis;
            expansion = expansion,
            nside = nside,
            packet_kernel = route,
        )
    end
    GC.gc()
    return nothing
end

function _print_run_summary(run)
    source = run.result.source
    fixed_block = run.result.fixed_block
    println("run_begin")
    println("route=", run.route)
    println("packet_kernel=", run.route == :default ? "default" : String(run.route))
    @printf("wall_s=%.9f\n", run.wall_seconds)
    println("fixed_dimension=", size(fixed_block.coefficient_matrix, 2))
    println("source_support_count=", length(source.sequence.support_indices))
    println("source_working_box=", source.sequence.working_box)
    println("did_split=", source.geometry.did_split)
    println("shared_shell_count=", length(source.shared_shell_layers))
    println("child_sequence_count=", length(source.child_sequences))
    println("parent_box=", _format_box(source.geometry.parent_box))
    println("split_working_box=", _format_box(source.geometry.working_box))
    println("split_index=", source.geometry.split_index)
    println(
        "child_boxes=",
        string("[", join((_format_box(box) for box in source.geometry.child_boxes), ", "), "]"),
    )
    println(
        "child_physical_widths=",
        string(
            "[",
            join((_format_widths(widths) for widths in source.geometry.child_physical_widths), ", "),
            "]",
        ),
    )
    println("child_column_ranges=", _format_ranges(source.child_column_ranges))
    println("midpoint_slab_column_range=", source.midpoint_slab_column_range)
    _print_selected_timings(run.stats)
    println("run_end")
    return nothing
end

function _print_residuals(left_label::AbstractString, left, right_label::AbstractString, right)
    println("comparison_begin")
    println("left=", left_label)
    println("right=", right_label)
    println("fixed_dimension_match=", size(left.coefficient_matrix, 2) == size(right.coefficient_matrix, 2))
    @printf(
        "coefficient_matrix_inf_residual=%.9e\n",
        norm(left.coefficient_matrix - right.coefficient_matrix, Inf),
    )
    @printf("overlap_inf_residual=%.9e\n", norm(left.overlap - right.overlap, Inf))
    @printf("kinetic_inf_residual=%.9e\n", norm(left.kinetic - right.kinetic, Inf))
    @printf("gaussian_sum_inf_residual=%.9e\n", norm(left.gaussian_sum - right.gaussian_sum, Inf))
    @printf("pair_sum_inf_residual=%.9e\n", norm(left.pair_sum - right.pair_sum, Inf))
    println("comparison_end")
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
    nside = parse(Int, get(values, "nside", "6"))
    bond_length = parse(Float64, get(values, "bond_length", "1.4"))
    core_spacing = parse(Float64, get(values, "core_spacing", "0.5"))
    xmax_parallel = parse(Float64, get(values, "xmax_parallel", "18.0"))
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

    routes = (:default, :support_reference, :factorized_direct)
    if warmup
        for route in routes
            _warmup_kernel(basis, expansion, route; nside = nside)
        end
    end

    runs = [_run_kernel(basis, expansion, route; nside = nside) for route in routes]
    for run in runs
        _print_run_summary(run)
    end

    default_block = runs[1].result.fixed_block
    support_block = runs[2].result.fixed_block
    direct_block = runs[3].result.fixed_block
    _print_residuals("default", default_block, "factorized_direct", direct_block)
    _print_residuals("support_reference", support_block, "factorized_direct", direct_block)

    return nothing
end

main(ARGS)
