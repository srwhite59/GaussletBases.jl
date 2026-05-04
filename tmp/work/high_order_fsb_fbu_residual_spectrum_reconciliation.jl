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

function gram_stats(
    residual_coefficients::AbstractMatrix{<:Real},
    parent_overlap::AbstractMatrix{<:Real};
    thresholds = (1.0e-8, 1.0e-10, 1.0e-12),
)
    gram = Matrix{Float64}(transpose(residual_coefficients) * parent_overlap * residual_coefficients)
    eigenvalues = eigvals(Symmetric(gram))
    positive = Float64[value for value in eigenvalues if value > 0.0]
    return (
        trace = isempty(positive) ? 0.0 : sum(positive),
        largest_eigenvalue = isempty(positive) ? 0.0 : maximum(positive),
        count_above_1e8 = count(>(thresholds[1]), positive),
        count_above_1e10 = count(>(thresholds[2]), positive),
        count_above_1e12 = count(>(thresholds[3]), positive),
        raw_rank = length(positive),
    )
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

function build_physical_full_shell_basis(
    side_data::AbstractVector,
    parent_overlap::AbstractMatrix{<:Real},
    parent_weights::AbstractVector{<:Real},
)
    accumulated = zeros(Float64, size(parent_overlap, 1), 0)
    rows = Any[]
    for (side_index, data) in enumerate(side_data)
        side = data.side
        full_block = data.full_block
        shell_block = data.shell_block
        before_dim = size(accumulated, 2)

        full_before_residual = project_out_or_identity(full_block, accumulated, parent_overlap)
        shell_before_residual = project_out_or_identity(shell_block, accumulated, parent_overlap)

        if side_index == 1
            next_accumulated = Matrix{Float64}(full_block)
            GaussletBases._experimental_high_order_sign_fix_columns!(next_accumulated, parent_weights)
        else
            shell_increment = GaussletBases._experimental_high_order_metric_project_out(
                shell_block,
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
        full_after_residual = project_out_or_identity(full_block, next_accumulated, parent_overlap)

        push!(rows, (
            side = side,
            target_kind = "full_block_before_add",
            target_dimension = size(full_block, 2),
            accumulated_before_dimension = before_dim,
            accumulated_after_dimension = after_dim,
            stats = gram_stats(full_before_residual, parent_overlap),
        ))
        push!(rows, (
            side = side,
            target_kind = "shell_only_before_add",
            target_dimension = size(shell_block, 2),
            accumulated_before_dimension = before_dim,
            accumulated_after_dimension = after_dim,
            stats = gram_stats(shell_before_residual, parent_overlap),
        ))
        push!(rows, (
            side = side,
            target_kind = "full_block_after_add",
            target_dimension = size(full_block, 2),
            accumulated_before_dimension = before_dim,
            accumulated_after_dimension = after_dim,
            stats = gram_stats(full_after_residual, parent_overlap),
        ))

        accumulated = next_accumulated
    end
    return accumulated, rows
end

function build_physical_full_block_union(
    side_data::AbstractVector,
    parent_overlap::AbstractMatrix{<:Real},
    parent_weights::AbstractVector{<:Real},
)
    union_coefficients = Matrix{Float64}(hcat([data.full_block for data in side_data]...))
    return GaussletBases._experimental_high_order_lowdin_cleanup(
        union_coefficients,
        parent_overlap;
        sign_vector = parent_weights,
    )
end

function fbu_ground_state_coefficients(fbu_heplus_data)
    decomposition = eigen(Symmetric(Matrix{Float64}(fbu_heplus_data.projected_hamiltonian)))
    return Float64[Float64(value) for value in decomposition.vectors[:, 1]]
end

function fsb_capture_deficiency(parent_overlap, fsb_coefficients, fbu_coefficients, fbu_ground_coefficients)
    cross_overlap = Matrix{Float64}(transpose(fsb_coefficients) * parent_overlap * fbu_coefficients)
    projected = cross_overlap * fbu_ground_coefficients
    deficiency = 1.0 - dot(projected, projected)
    return max(0.0, deficiency), cross_overlap
end

function reconcile_case(
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
        physical_blocks = timed_value() do
            map(case.sides) do side
                physical = GaussletBases._experimental_high_order_physical_full_block_3d(
                    axis_data.value,
                    side;
                    doside = case.doside,
                )
                (
                    side = side,
                    full_block = Matrix{Float64}(physical.shell.full_block_coefficients),
                    shell_block = Matrix{Float64}(physical.shell.shell_coefficients),
                )
            end
        end
        fsb_data = timed_value() do
            build_physical_full_shell_basis(
                physical_blocks.value,
                parent_overlap,
                parent_weights,
            )
        end
        fbu_data = timed_value() do
            build_physical_full_block_union(
                physical_blocks.value,
                parent_overlap,
                parent_weights,
            )
        end
        fsb_coefficients = fsb_data.value[1]
        residual_rows = fsb_data.value[2]
        fbu_coefficients = fbu_data.value
        fsb_heplus = timed_value() do
            GaussletBases._experimental_high_order_doside_heplus_data(parent_data.value, fsb_coefficients)
        end
        fbu_heplus = timed_value() do
            GaussletBases._experimental_high_order_doside_heplus_data(parent_data.value, fbu_coefficients)
        end
        fbu_ground_coefficients = fbu_ground_state_coefficients(fbu_heplus.value)
        capture_deficiency, cross_overlap = fsb_capture_deficiency(
            parent_overlap,
            fsb_coefficients,
            fbu_coefficients,
            fbu_ground_coefficients,
        )
        cross_overlap_error = if size(cross_overlap, 1) == size(cross_overlap, 2)
            norm(transpose(cross_overlap) * cross_overlap - I, Inf)
        else
            NaN
        end
        return (
            case = case,
            scale_factor = scale_factor,
            basis_size_parent = length(basis_data.value)^3,
            fsb_dimension = size(fsb_coefficients, 2),
            fbu_dimension = size(fbu_coefficients, 2),
            energy_gap = fsb_heplus.value.ground_energy - fbu_heplus.value.ground_energy,
            capture_deficiency = capture_deficiency,
            cross_overlap_error = cross_overlap_error,
            fsb_overlap_error = fsb_heplus.value.overlap_error,
            fbu_overlap_error = fbu_heplus.value.overlap_error,
            residual_rows = residual_rows,
            timing_basis = basis_data,
            timing_axis = axis_data,
            timing_parent = parent_data,
            timing_blocks = physical_blocks,
            timing_fsb = fsb_data,
            timing_fbu = fbu_data,
            timing_heplus_fsb = fsb_heplus,
            timing_heplus_fbu = fbu_heplus,
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
        "target_kind",
        "side",
        "target_dimension",
        "accumulated_before_dimension",
        "accumulated_after_dimension",
        "fbu_final_dimension",
        "fsb_final_dimension",
        "energy_gap",
        "capture_deficiency",
        "final_cross_overlap_error",
        "residual_trace",
        "largest_eigenvalue",
        "count_above_1e8",
        "count_above_1e10",
        "count_above_1e12",
        "raw_rank",
        "time_basis_s",
        "time_axis_s",
        "time_parent_s",
        "time_blocks_s",
        "time_fsb_s",
        "time_fbu_s",
        "time_heplus_fsb_s",
        "time_heplus_fbu_s",
        "time_total_s",
        "bytes_total",
    ]
    open(path, "w") do io
        println(io, join(header, '\t'))
        for case_row in case_rows
            for residual_row in case_row.residual_rows
                stats = residual_row.stats
                values = [
                    case_row.case.label,
                    string(case_row.case.count),
                    string(case_row.case.doside),
                    join(case_row.case.sides, ","),
                    mapping_label(case_row.scale_factor),
                    residual_row.target_kind,
                    string(residual_row.side),
                    string(residual_row.target_dimension),
                    string(residual_row.accumulated_before_dimension),
                    string(residual_row.accumulated_after_dimension),
                    string(case_row.fbu_dimension),
                    string(case_row.fsb_dimension),
                    @sprintf("%.16e", case_row.energy_gap),
                    @sprintf("%.16e", case_row.capture_deficiency),
                    @sprintf("%.16e", case_row.cross_overlap_error),
                    @sprintf("%.16e", stats.trace),
                    @sprintf("%.16e", stats.largest_eigenvalue),
                    string(stats.count_above_1e8),
                    string(stats.count_above_1e10),
                    string(stats.count_above_1e12),
                    string(stats.raw_rank),
                    @sprintf("%.6f", case_row.timing_basis.seconds),
                    @sprintf("%.6f", case_row.timing_axis.seconds),
                    @sprintf("%.6f", case_row.timing_parent.seconds),
                    @sprintf("%.6f", case_row.timing_blocks.seconds),
                    @sprintf("%.6f", case_row.timing_fsb.seconds),
                    @sprintf("%.6f", case_row.timing_fbu.seconds),
                    @sprintf("%.6f", case_row.timing_heplus_fsb.seconds),
                    @sprintf("%.6f", case_row.timing_heplus_fbu.seconds),
                    @sprintf("%.6f", case_row.timing_total.seconds),
                    string(case_row.timing_total.bytes),
                ]
                println(io, join(values, '\t'))
            end
        end
    end
end

function write_summary(path::AbstractString, case_rows)
    open(path, "w") do io
        println(io, "High-order full-shell basis (FSB) / full-block union (FBU) residual-spectrum reconciliation")
        println(io, "backend = :pgdg_localized_experimental")
        println(io, "Route = current physical-coordinate polynomial construction")
        println(io)
        for case_row in case_rows
            println(io, string(case_row.case.label, " / ", mapping_label(case_row.scale_factor)))
            println(io, @sprintf("  parent basis size: %d", case_row.basis_size_parent))
            println(io, @sprintf("  final FSB dim: %d", case_row.fsb_dimension))
            println(io, @sprintf("  final FBU dim: %d", case_row.fbu_dimension))
            println(io, @sprintf("  E_FSB - E_FBU: %.3e", case_row.energy_gap))
            println(io, @sprintf("  FBU-state capture deficiency: %.3e", case_row.capture_deficiency))
            println(io, @sprintf("  final cross-overlap error: %.3e", case_row.cross_overlap_error))
            println(io, @sprintf("  wall time total: %.3fs", case_row.timing_total.seconds))
            for residual_row in case_row.residual_rows
                stats = residual_row.stats
                println(
                    io,
                    @sprintf(
                        "  side %-2d %-22s trace=% .3e largest=% .3e counts=(%d,%d,%d) dims %d -> %d target=%d",
                        residual_row.side,
                        residual_row.target_kind,
                        stats.trace,
                        stats.largest_eigenvalue,
                        stats.count_above_1e8,
                        stats.count_above_1e10,
                        stats.count_above_1e12,
                        residual_row.accumulated_before_dimension,
                        residual_row.accumulated_after_dimension,
                        residual_row.target_dimension,
                    ),
                )
            end
            println(io)
        end
    end
end

function main()
    backend = :pgdg_localized_experimental
    expansion = coulomb_gaussian_expansion(doacc = false)
    z_value = 2.0
    scale_factors = Union{Nothing,Float64}[nothing, 1.0, 0.8, 0.6, 0.4, 0.2]
    cases = Any[
        (ReconciliationCase("count11_doside5", 11, 5, [5, 7, 9, 11]), scale_factors),
        (ReconciliationCase("count13_doside5_s1", 13, 5, [5, 7, 9, 11, 13]), Union{Nothing,Float64}[1.0]),
    ]
    case_rows = Any[]
    for (case, scales) in cases, scale_factor in scales
        push!(case_rows, reconcile_case(case, scale_factor; backend = backend, expansion = expansion, z_value = z_value))
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
