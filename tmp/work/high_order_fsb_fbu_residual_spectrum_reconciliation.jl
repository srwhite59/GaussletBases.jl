using Dates
using LinearAlgebra
using Printf

using GaussletBases

struct ReconciliationCase
    label::String
    count::Int
    doside::Int
    sides::Vector{Int}
end

function timed_value(f::Function)
    timed = @timed f()
    return (
        value = timed.value,
        seconds = Float64(timed.time),
        bytes = Int(timed.bytes),
        gctime = Float64(timed.gctime),
    )
end

function mapping_for_scale(scale_factor::Union{Nothing,Float64})
    if isnothing(scale_factor)
        return IdentityMapping()
    end
    return AsinhMapping(c = 0.2, s = scale_factor * sqrt(0.4), tail_spacing = 10.0)
end

function mapping_label(scale_factor::Union{Nothing,Float64})
    if isnothing(scale_factor)
        return "identity"
    end
    return @sprintf("s_over_s0=%.1f", scale_factor)
end

function metric_gram_stats(
    coefficients::AbstractMatrix{<:Real},
    parent_overlap::AbstractMatrix{<:Real};
    rank_tol::Float64 = 1.0e-12,
)
    gram = Matrix{Float64}(transpose(coefficients) * parent_overlap * coefficients)
    eigenvalues = Float64[Float64(value) for value in eigvals(Symmetric(gram))]
    positive = sort(Float64[value for value in eigenvalues if value > 0.0])
    return (
        trace = isempty(positive) ? 0.0 : sum(positive),
        largest_eigenvalue = isempty(positive) ? 0.0 : positive[end],
        smallest_positive_eigenvalue = isempty(positive) ? 0.0 : positive[1],
        metric_rank = count(>(rank_tol), positive),
        count_above_1e8 = count(>(1.0e-8), positive),
        count_above_1e10 = count(>(1.0e-10), positive),
        count_above_1e12 = count(>(1.0e-12), positive),
    )
end

function centered_local_cube_coefficients(parent_side::Int, side::Int)
    interval = GaussletBases._experimental_high_order_centered_interval(parent_side, side)
    parent_dimension = parent_side^3
    column_count = side^3
    coefficients = zeros(Float64, parent_dimension, column_count)
    column = 0
    for ix in interval, iy in interval, iz in interval
        column += 1
        flat = GaussletBases._cartesian_flat_index(ix, iy, iz, (parent_side, parent_side, parent_side))
        coefficients[flat, column] = 1.0
    end
    return coefficients
end

function project_out_or_identity(
    target_coefficients::AbstractMatrix{<:Real},
    accumulated_coefficients::AbstractMatrix{<:Real},
    parent_overlap::AbstractMatrix{<:Real},
)
    if size(accumulated_coefficients, 2) == 0
        return Matrix{Float64}(target_coefficients)
    end
    return Matrix{Float64}(
        GaussletBases._experimental_high_order_metric_project_out(
            target_coefficients,
            accumulated_coefficients,
            parent_overlap,
        ),
    )
end

function route_side_targets(
    axis_data,
    parent_overlap::AbstractMatrix{<:Real},
    case::ReconciliationCase,
)
    parent_side = length(axis_data.centers)
    physical_targets = Dict{Int,NamedTuple}()
    debug_targets = Dict{Int,NamedTuple}()
    true_cube_targets = Dict{Int,NamedTuple}()
    for side in case.sides
        physical = GaussletBases._experimental_high_order_physical_full_block_3d(
            axis_data,
            side;
            doside = case.doside,
        )
        debug = GaussletBases._experimental_high_order_tensor_shell_3d(
            axis_data,
            side;
            doside = case.doside,
        )
        true_cube = centered_local_cube_coefficients(parent_side, side)
        physical_full = Matrix{Float64}(physical.shell.full_block_coefficients)
        physical_shell = Matrix{Float64}(physical.shell.shell_coefficients)
        debug_full = Matrix{Float64}(debug.full_block_coefficients)
        debug_shell = Matrix{Float64}(debug.shell_coefficients)
        physical_targets[side] = (
            full = physical_full,
            shell = physical_shell,
            raw_column_count = size(physical_full, 2),
            gram_stats = metric_gram_stats(physical_full, parent_overlap),
        )
        debug_targets[side] = (
            full = debug_full,
            shell = debug_shell,
            raw_column_count = size(debug_full, 2),
            gram_stats = metric_gram_stats(debug_full, parent_overlap),
        )
        true_cube_targets[side] = (
            full = true_cube,
            raw_column_count = size(true_cube, 2),
            gram_stats = metric_gram_stats(true_cube, parent_overlap),
        )
    end
    return (physical = physical_targets, debug = debug_targets, true_cube = true_cube_targets)
end

function build_route_fsb(
    route_label::String,
    route_targets::Dict{Int,NamedTuple},
    parent_overlap::AbstractMatrix{<:Real},
    parent_weights::AbstractVector{<:Real},
    all_targets,
    case::ReconciliationCase,
)
    accumulated = zeros(Float64, size(parent_overlap, 1), 0)
    rows = Any[]
    for (side_index, side) in enumerate(case.sides)
        route_target = route_targets[side]
        before_dim = size(accumulated, 2)
        if side_index == 1
            next_accumulated = Matrix{Float64}(route_target.full)
            GaussletBases._experimental_high_order_sign_fix_columns!(next_accumulated, parent_weights)
        else
            shell_increment = GaussletBases._experimental_high_order_metric_project_out(
                route_target.shell,
                accumulated,
                parent_overlap,
            )
            cleaned_increment = GaussletBases._experimental_high_order_lowdin_cleanup(
                shell_increment,
                parent_overlap;
                sign_vector = parent_weights,
            )
            next_accumulated = Matrix{Float64}(hcat(accumulated, cleaned_increment))
        end
        after_dim = size(next_accumulated, 2)
        for (target_kind, target_dict) in (
            ("reduced_physical_transformed_block", all_targets.physical),
            ("debug_u_transformed_block", all_targets.debug),
            ("true_local_distorted_cube", all_targets.true_cube),
        )
            target = target_dict[side]
            residual = project_out_or_identity(target.full, next_accumulated, parent_overlap)
            residual_stats = metric_gram_stats(residual, parent_overlap)
            push!(rows, (
                route_label = route_label,
                side = side,
                target_kind = target_kind,
                raw_target_column_count = target.raw_column_count,
                target_gram_stats = target.gram_stats,
                accumulated_before_dimension = before_dim,
                accumulated_after_dimension = after_dim,
                residual_stats = residual_stats,
            ))
        end
        accumulated = next_accumulated
    end
    union_raw = Matrix{Float64}(hcat([route_targets[side].full for side in case.sides]...))
    union_raw_stats = metric_gram_stats(union_raw, parent_overlap)
    cleaned_union = GaussletBases._experimental_high_order_lowdin_cleanup(
        union_raw,
        parent_overlap;
        sign_vector = parent_weights,
    )
    cross_overlap = Matrix{Float64}(transpose(accumulated) * parent_overlap * cleaned_union)
    cross_overlap_error = if size(cross_overlap, 1) == size(cross_overlap, 2)
        norm(transpose(cross_overlap) * cross_overlap - I, Inf)
    else
        NaN
    end
    return (
        final_fsb = accumulated,
        rows = rows,
        union_raw_column_count = size(union_raw, 2),
        union_raw_stats = union_raw_stats,
        cleaned_union_dimension = size(cleaned_union, 2),
        cleaned_fsb_dimension = size(accumulated, 2),
        cross_overlap_error = cross_overlap_error,
    )
end

function audit_case(
    case::ReconciliationCase,
    scale_factor::Union{Nothing,Float64};
    backend::Symbol,
    expansion,
    z_value::Float64,
)
    timing_total = timed_value() do
        mapping = mapping_for_scale(scale_factor)
        basis_data = timed_value() do
            build_basis(
                MappedUniformBasisSpec(
                    :G10;
                    count = case.count,
                    mapping = mapping,
                    reference_spacing = 1.0,
                ),
            )
        end
        axis_data = timed_value() do
            GaussletBases._experimental_high_order_axis_data_1d(
                basis_data.value;
                backend = backend,
                one_body_exponents = expansion.exponents,
                one_body_center = 0.0,
            )
        end
        parent_data = timed_value() do
            GaussletBases._experimental_high_order_parent_one_body_data(
                basis_data.value;
                axis_data = axis_data.value,
                backend = backend,
                expansion = expansion,
                Z = z_value,
                include_parent_projection_data = true,
                include_reference_energy = false,
            )
        end
        parent_overlap = parent_data.value.parent_overlap
        parent_weights = GaussletBases._experimental_high_order_parent_weights_3d(axis_data.value)
        targets_data = timed_value() do
            route_side_targets(axis_data.value, parent_overlap, case)
        end
        physical_route = timed_value() do
            build_route_fsb(
                "physical_x",
                targets_data.value.physical,
                parent_overlap,
                parent_weights,
                targets_data.value,
                case,
            )
        end
        debug_route = timed_value() do
            build_route_fsb(
                "debug_u",
                targets_data.value.debug,
                parent_overlap,
                parent_weights,
                targets_data.value,
                case,
            )
        end
        return (
            case = case,
            scale_factor = scale_factor,
            basis_size_parent = length(basis_data.value)^3,
            physical_route = physical_route.value,
            debug_route = debug_route.value,
            timing_basis = basis_data,
            timing_axis = axis_data,
            timing_parent = parent_data,
            timing_targets = targets_data,
            timing_physical_route = physical_route,
            timing_debug_route = debug_route,
        )
    end
    return merge(timing_total.value, (timing_total = timing_total,))
end

function write_tsv(path::AbstractString, case_rows)
    header = [
        "case_label",
        "count",
        "doside",
        "sides",
        "mapping_label",
        "route_label",
        "side",
        "target_kind",
        "raw_target_column_count",
        "target_metric_rank",
        "target_smallest_positive_eigenvalue",
        "target_largest_eigenvalue",
        "target_count_above_1e8",
        "target_count_above_1e10",
        "target_count_above_1e12",
        "fbu_union_raw_column_count",
        "fbu_union_metric_rank",
        "fbu_union_smallest_positive_eigenvalue",
        "fbu_union_largest_eigenvalue",
        "fbu_union_count_above_1e8",
        "fbu_union_count_above_1e10",
        "fbu_union_count_above_1e12",
        "fbu_cleaned_dimension",
        "fsb_cleaned_dimension",
        "accumulated_before_dimension",
        "accumulated_after_dimension",
        "final_cross_overlap_error",
        "residual_trace",
        "residual_largest_eigenvalue",
        "residual_count_above_1e8",
        "residual_count_above_1e10",
        "residual_count_above_1e12",
        "time_total_s",
        "bytes_total",
    ]
    open(path, "w") do io
        println(io, join(header, '\t'))
        for case_row in case_rows
            for route_data in (case_row.physical_route, case_row.debug_route)
                for residual_row in route_data.rows
                    tstats = residual_row.target_gram_stats
                    ustats = route_data.union_raw_stats
                    rstats = residual_row.residual_stats
                    values = [
                        case_row.case.label,
                        string(case_row.case.count),
                        string(case_row.case.doside),
                        join(case_row.case.sides, ","),
                        mapping_label(case_row.scale_factor),
                        residual_row.route_label,
                        string(residual_row.side),
                        residual_row.target_kind,
                        string(residual_row.raw_target_column_count),
                        string(tstats.metric_rank),
                        @sprintf("%.16e", tstats.smallest_positive_eigenvalue),
                        @sprintf("%.16e", tstats.largest_eigenvalue),
                        string(tstats.count_above_1e8),
                        string(tstats.count_above_1e10),
                        string(tstats.count_above_1e12),
                        string(route_data.union_raw_column_count),
                        string(ustats.metric_rank),
                        @sprintf("%.16e", ustats.smallest_positive_eigenvalue),
                        @sprintf("%.16e", ustats.largest_eigenvalue),
                        string(ustats.count_above_1e8),
                        string(ustats.count_above_1e10),
                        string(ustats.count_above_1e12),
                        string(route_data.cleaned_union_dimension),
                        string(route_data.cleaned_fsb_dimension),
                        string(residual_row.accumulated_before_dimension),
                        string(residual_row.accumulated_after_dimension),
                        @sprintf("%.16e", route_data.cross_overlap_error),
                        @sprintf("%.16e", rstats.trace),
                        @sprintf("%.16e", rstats.largest_eigenvalue),
                        string(rstats.count_above_1e8),
                        string(rstats.count_above_1e10),
                        string(rstats.count_above_1e12),
                        @sprintf("%.6f", case_row.timing_total.seconds),
                        string(case_row.timing_total.bytes),
                    ]
                    println(io, join(values, '\t'))
                end
            end
        end
    end
end

function write_summary(path::AbstractString, case_rows)
    open(path, "w") do io
        println(io, "High-order full-shell basis (FSB) / full-block union (FBU) target-rank audit")
        println(io, "backend = :pgdg_localized_experimental")
        println(io, "Routes: physical_x, debug_u")
        println(io)
        for case_row in case_rows
            println(io, string(case_row.case.label, " / ", mapping_label(case_row.scale_factor)))
            println(io, @sprintf("  parent basis size: %d", case_row.basis_size_parent))
            println(io, @sprintf("  wall time total: %.3fs", case_row.timing_total.seconds))
            for route_data in (case_row.physical_route, case_row.debug_route)
                println(io, string("  route ", route_data.rows[1].route_label))
                println(io, @sprintf("    FBU raw columns: %d", route_data.union_raw_column_count))
                println(io, @sprintf("    FBU raw metric rank (>1e-12): %d", route_data.union_raw_stats.metric_rank))
                println(io, @sprintf("    FBU cleaned dimension: %d", route_data.cleaned_union_dimension))
                println(io, @sprintf("    FSB cleaned dimension: %d", route_data.cleaned_fsb_dimension))
                println(io, @sprintf("    final cross-overlap error: %.3e", route_data.cross_overlap_error))
                for residual_row in route_data.rows
                    tstats = residual_row.target_gram_stats
                    rstats = residual_row.residual_stats
                    println(
                        io,
                        @sprintf(
                            "    side %-2d %-33s target raw=%4d rank=%4d residual trace=% .3e largest=% .3e counts=(%d,%d,%d)",
                            residual_row.side,
                            residual_row.target_kind,
                            residual_row.raw_target_column_count,
                            tstats.metric_rank,
                            rstats.trace,
                            rstats.largest_eigenvalue,
                            rstats.count_above_1e8,
                            rstats.count_above_1e10,
                            rstats.count_above_1e12,
                        ),
                    )
                end
            end
            println(io)
        end
    end
end

function main()
    backend = :pgdg_localized_experimental
    expansion = coulomb_gaussian_expansion(doacc = false)
    z_value = 2.0
    cases = Any[
        (ReconciliationCase("count11_doside5", 11, 5, [5, 7, 9, 11]), Union{Nothing,Float64}[nothing, 1.0]),
        (ReconciliationCase("count13_doside5", 13, 5, [5, 7, 9, 11, 13]), Union{Nothing,Float64}[1.0]),
    ]
    case_rows = Any[]
    for (case, scales) in cases, scale_factor in scales
        push!(case_rows, audit_case(case, scale_factor; backend = backend, expansion = expansion, z_value = z_value))
    end
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    tsv_path = joinpath(@__DIR__, "high_order_fsb_fbu_residual_spectrum_reconciliation_$timestamp.tsv")
    summary_path = joinpath(@__DIR__, "high_order_fsb_fbu_residual_spectrum_reconciliation_$timestamp.txt")
    write_tsv(tsv_path, case_rows)
    write_summary(summary_path, case_rows)
    println("Wrote:")
    println(tsv_path)
    println(summary_path)
end

main()
