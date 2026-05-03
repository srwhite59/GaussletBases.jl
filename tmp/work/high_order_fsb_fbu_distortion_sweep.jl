using Dates
using LinearAlgebra
using Printf

using GaussletBases

struct SweepCase
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

function debug_full_block_union(axis_data, sides::AbstractVector{<:Integer}, doside::Integer)
    return GaussletBases._experimental_high_order_orthonormalized_full_block_union_coefficients(
        axis_data,
        sides;
        doside = doside,
    )
end

function debug_full_shell_basis(axis_data, sides::AbstractVector{<:Integer}, doside::Integer)
    side_values = Int[Int(side) for side in sides]
    parent_overlap = GaussletBases._experimental_high_order_parent_overlap_3d(axis_data)
    parent_weights = GaussletBases._experimental_high_order_parent_weights_3d(axis_data)
    first_shell = GaussletBases._experimental_high_order_tensor_shell_3d(
        axis_data,
        first(side_values);
        doside = doside,
    )
    coefficients = Matrix{Float64}(first_shell.full_block_coefficients)
    GaussletBases._experimental_high_order_sign_fix_columns!(coefficients, parent_weights)
    for side in side_values[2:end]
        shell = GaussletBases._experimental_high_order_tensor_shell_3d(
            axis_data,
            side;
            doside = doside,
        )
        projected_out = GaussletBases._experimental_high_order_metric_project_out(
            shell.shell_coefficients,
            coefficients,
            parent_overlap,
        )
        cleaned = GaussletBases._experimental_high_order_lowdin_cleanup(
            projected_out,
            parent_overlap;
            sign_vector = parent_weights,
        )
        coefficients = Matrix{Float64}(hcat(coefficients, cleaned))
    end
    return coefficients
end

function physical_full_block_union(axis_data, sides::AbstractVector{<:Integer}, doside::Integer)
    parent_overlap = GaussletBases._experimental_high_order_parent_overlap_3d(axis_data)
    parent_weights = GaussletBases._experimental_high_order_parent_weights_3d(axis_data)
    full_blocks = Matrix{Float64}[]
    for side in sides
        physical = GaussletBases._experimental_high_order_physical_full_block_3d(
            axis_data,
            Int(side);
            doside = doside,
        )
        push!(full_blocks, Matrix{Float64}(physical.shell.full_block_coefficients))
    end
    union_coefficients = Matrix{Float64}(hcat(full_blocks...))
    return GaussletBases._experimental_high_order_lowdin_cleanup(
        union_coefficients,
        parent_overlap;
        sign_vector = parent_weights,
    )
end

function physical_full_shell_basis(axis_data, sides::AbstractVector{<:Integer}, doside::Integer)
    side_values = Int[Int(side) for side in sides]
    parent_overlap = GaussletBases._experimental_high_order_parent_overlap_3d(axis_data)
    parent_weights = GaussletBases._experimental_high_order_parent_weights_3d(axis_data)
    first_full = GaussletBases._experimental_high_order_physical_full_block_3d(
        axis_data,
        first(side_values);
        doside = doside,
    )
    coefficients = Matrix{Float64}(first_full.shell.full_block_coefficients)
    GaussletBases._experimental_high_order_sign_fix_columns!(coefficients, parent_weights)
    for side in side_values[2:end]
        shell = GaussletBases._experimental_high_order_physical_shell_3d(
            axis_data,
            side;
            doside = doside,
        )
        projected_out = GaussletBases._experimental_high_order_metric_project_out(
            shell.shell.shell_coefficients,
            coefficients,
            parent_overlap,
        )
        cleaned = GaussletBases._experimental_high_order_lowdin_cleanup(
            projected_out,
            parent_overlap;
            sign_vector = parent_weights,
        )
        coefficients = Matrix{Float64}(hcat(coefficients, cleaned))
    end
    return coefficients
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

function run_route(
    route_label::String,
    build_fsb::Function,
    build_fbu::Function,
    case::SweepCase,
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
        fsb_data = timed_value() do
            build_fsb(axis_data.value, case.sides, case.doside)
        end
        fbu_data = timed_value() do
            build_fbu(axis_data.value, case.sides, case.doside)
        end
        fsb_heplus = timed_value() do
            GaussletBases._experimental_high_order_doside_heplus_data(parent_data.value, fsb_data.value)
        end
        fbu_heplus = timed_value() do
            GaussletBases._experimental_high_order_doside_heplus_data(parent_data.value, fbu_data.value)
        end
        fbu_ground_coefficients = fbu_ground_state_coefficients(fbu_heplus.value)
        capture_deficiency, cross_overlap = fsb_capture_deficiency(
            parent_data.value.parent_overlap,
            fsb_data.value,
            fbu_data.value,
            fbu_ground_coefficients,
        )
        cross_overlap_error = if size(cross_overlap, 1) == size(cross_overlap, 2)
            norm(transpose(cross_overlap) * cross_overlap - I, Inf)
        else
            NaN
        end
        return (
            route_label = route_label,
            case = case,
            scale_factor = scale_factor,
            basis_size_parent = length(basis_data.value)^3,
            fsb_dimension = size(fsb_data.value, 2),
            fbu_dimension = size(fbu_data.value, 2),
            fsb_overlap_error = fsb_heplus.value.overlap_error,
            fbu_overlap_error = fbu_heplus.value.overlap_error,
            fsb_ground_energy = fsb_heplus.value.ground_energy,
            fbu_ground_energy = fbu_heplus.value.ground_energy,
            energy_gap = fsb_heplus.value.ground_energy - fbu_heplus.value.ground_energy,
            capture_deficiency = capture_deficiency,
            cross_overlap_error = cross_overlap_error,
            timing_basis = basis_data,
            timing_axis = axis_data,
            timing_parent = parent_data,
            timing_fsb = fsb_data,
            timing_fbu = fbu_data,
            timing_heplus_fsb = fsb_heplus,
            timing_heplus_fbu = fbu_heplus,
        )
    end
    return merge(timing_total.value, (timing_total = timing_total,))
end

function write_tsv(path::AbstractString, rows)
    header = [
        "route_label",
        "case_label",
        "count",
        "doside",
        "sides",
        "mapping_label",
        "basis_size_parent",
        "fsb_dimension",
        "fbu_dimension",
        "fsb_ground_energy",
        "fbu_ground_energy",
        "energy_gap",
        "capture_deficiency",
        "cross_overlap_error",
        "fsb_overlap_error",
        "fbu_overlap_error",
        "time_basis_s",
        "time_axis_s",
        "time_parent_s",
        "time_fsb_s",
        "time_fbu_s",
        "time_heplus_fsb_s",
        "time_heplus_fbu_s",
        "time_total_s",
        "bytes_basis",
        "bytes_axis",
        "bytes_parent",
        "bytes_fsb",
        "bytes_fbu",
        "bytes_heplus_fsb",
        "bytes_heplus_fbu",
        "bytes_total",
    ]
    open(path, "w") do io
        println(io, join(header, '\t'))
        for row in rows
            values = [
                row.route_label,
                row.case.label,
                string(row.case.count),
                string(row.case.doside),
                join(row.case.sides, ","),
                mapping_label(row.scale_factor),
                string(row.basis_size_parent),
                string(row.fsb_dimension),
                string(row.fbu_dimension),
                @sprintf("%.16f", row.fsb_ground_energy),
                @sprintf("%.16f", row.fbu_ground_energy),
                @sprintf("%.16e", row.energy_gap),
                @sprintf("%.16e", row.capture_deficiency),
                @sprintf("%.16e", row.cross_overlap_error),
                @sprintf("%.16e", row.fsb_overlap_error),
                @sprintf("%.16e", row.fbu_overlap_error),
                @sprintf("%.6f", row.timing_basis.seconds),
                @sprintf("%.6f", row.timing_axis.seconds),
                @sprintf("%.6f", row.timing_parent.seconds),
                @sprintf("%.6f", row.timing_fsb.seconds),
                @sprintf("%.6f", row.timing_fbu.seconds),
                @sprintf("%.6f", row.timing_heplus_fsb.seconds),
                @sprintf("%.6f", row.timing_heplus_fbu.seconds),
                @sprintf("%.6f", row.timing_total.seconds),
                string(row.timing_basis.bytes),
                string(row.timing_axis.bytes),
                string(row.timing_parent.bytes),
                string(row.timing_fsb.bytes),
                string(row.timing_fbu.bytes),
                string(row.timing_heplus_fsb.bytes),
                string(row.timing_heplus_fbu.bytes),
                string(row.timing_total.bytes),
            ]
            println(io, join(values, '\t'))
        end
    end
end

function write_summary(path::AbstractString, rows)
    grouped = Dict{Tuple{String,String},Vector{Any}}()
    for row in rows
        push!(get!(grouped, (row.route_label, row.case.label), Any[]), row)
    end
    open(path, "w") do io
        println(io, "Bounded high-order full-shell basis (FSB) versus full-block union (FBU) distortion sweep")
        println(io, "backend = :pgdg_localized_experimental")
        println(io, "route_label = physical_x is the intended physical-coordinate polynomial route")
        println(io, "route_label = debug_u is the older compatibility/debug route retained only as a control")
        println(io)
        for key in sort(collect(keys(grouped)))
            route_label, case_label = key
            case_rows = grouped[key]
            max_gap = maximum(abs(row.energy_gap) for row in case_rows)
            max_capture = maximum(row.capture_deficiency for row in case_rows)
            max_cross = maximum(abs(row.cross_overlap_error) for row in case_rows if isfinite(row.cross_overlap_error))
            println(io, string(route_label, " / ", case_label))
            println(io, @sprintf("  max |E_FSB - E_FBU|: %.3e", max_gap))
            println(io, @sprintf("  max FBU-state capture deficiency: %.3e", max_capture))
            println(io, @sprintf("  max cross-overlap subspace error: %.3e", max_cross))
            for row in case_rows
                println(
                    io,
                    @sprintf(
                        "  %-12s fsb=%4d fbu=%4d gap=% .3e capture=% .3e total=%.3fs",
                        mapping_label(row.scale_factor),
                        row.fsb_dimension,
                        row.fbu_dimension,
                        row.energy_gap,
                        row.capture_deficiency,
                        row.timing_total.seconds,
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
    cases = SweepCase[
        SweepCase("count11_doside5", 11, 5, [5, 7, 9, 11]),
        SweepCase("count13_doside5", 13, 5, [5, 7, 9, 11, 13]),
    ]
    routes = [
        ("physical_x", physical_full_shell_basis, physical_full_block_union),
        ("debug_u", debug_full_shell_basis, debug_full_block_union),
    ]
    rows = Any[]
    for case in cases, scale_factor in scale_factors, (route_label, build_fsb, build_fbu) in routes
        push!(
            rows,
            run_route(
                route_label,
                build_fsb,
                build_fbu,
                case,
                scale_factor;
                backend = backend,
                expansion = expansion,
                z_value = z_value,
            ),
        )
    end
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    tsv_path = joinpath(@__DIR__, "high_order_fsb_fbu_distortion_sweep_$timestamp.tsv")
    summary_path = joinpath(@__DIR__, "high_order_fsb_fbu_distortion_sweep_$timestamp.txt")
    write_tsv(tsv_path, rows)
    write_summary(summary_path, rows)
    println("Wrote:")
    println(tsv_path)
    println(summary_path)
end

main()
